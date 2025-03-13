library("wesanderson")
library(grDevices)
library("corrplot")
corrplot.bis<-corrplot
library(vegan)
library(ade4)
library(Matrix)
#library(mlbench) ## Smooth plots
library(cluster)
library("RColorBrewer")
library("phyloseq")
library("ggplot2")
library("ggsignif")
library(psych)
library(dunn.test)
library(DESeq2)
library(CoDaSeq)
library(imputeTS)
library(ggordiplots)
library(GMPR)
set.seed(1234)
#################################
######### Read infile FUNCTIONS #
#################################

dada2genusTable <- function(abundance.file="", tax.file="", min.num.seqs = 10000, min.rar="min", plot_rar = TRUE){
	########################### Phyloseq object ####################
	#### 1: Read object 
	physeq.ent<-dada2_to_phyloseq(abundance.file, tax.file,sep=" ")
	join.physeq <- subset_taxa(physeq.ent,  Kingdom != "k_Eukaryota" &  Kingdom != "uc_" & Phylum != "NA" & Family  != "f_mitochondria" & Class   != "c_Chloroplast") 

	#### 2: Summarize to the Genus level ####
	join.physeq.genus <- tax_glom(join.physeq, taxrank = 'Genus') ### Could take a while

	#### 3: Rarefaction curve ####  Check is the rarefaction cut-off min.num.seqs  is good enough for sampling the bacterial diversity
	if(plot_rar == TRUE){
		samples2use<-colnames(otu_table( join.physeq.genus ))
		join.physeq.rarcuve = prune_samples(samples2use, join.physeq.genus)
		sample_sums(join.physeq.rarcuve)
		OTUtable<-data.frame(otu_table( join.physeq.rarcuve ))
		pdf("rarCurve.pdf", width=14, height=9)
			rarecurve(as.matrix(t(OTUtable)), step = 1000,  xlab = "Sample Size", ylab = "Genus", label = TRUE )
		dev.off()
	}
	#### 4: Create the Genus level table ####
	if(min.rar=="min"){min.num.seqs <- min(sample_sums(join.physeq.genus))}
	
	num_reads<-sample_sums(join.physeq.genus)[sample_sums(join.physeq.genus) < min.num.seqs]
	join.physeq.genus = prune_samples(names(which(sample_sums(join.physeq.genus) > min.num.seqs)), join.physeq.genus)
	Genus = data.frame(otu_table(join.physeq.genus)) 
	rownames(Genus) <-  tax_table(join.physeq.genus)[,6]

	#### 5: Dataset filtering and rarefaction ####
	join.physeq.genus.rar = rarefy_even_depth(join.physeq.genus, replace=FALSE,sample.size = min.num.seqs)
	Genus.rar = data.frame(otu_table(join.physeq.genus.rar)) 
	rownames(Genus.rar) <-  tax_table(join.physeq.genus.rar)[,6]
	ret.list<-list(Genus, Genus.rar)
	names(ret.list)<-c("Genus","Genus.rar")
	return(ret.list)
}


dada2_to_phyloseq = function(abundance.file,tax.file,sep=" "){
	library("phyloseq"); packageVersion("phyloseq")
	library("data.table")
	#system.time(data<-read.table(abundance.file, header=T, row.names=1, dec=".", sep=sep))
	system.time(data <- fread(abundance.file, header=T, fill = TRUE, dec=".", sep=sep) )
	data.colnames<-colnames(data)[1:length(colnames(data))-1]
	data <-data.frame(data, row.names=1)
	colnames(data) = data.colnames	
		
	data<-t(data)
	tax<-read.table(tax.file, header=T, row.names=1, dec=".", sep=sep)
	table(rownames(data) == rownames(tax))
	tax<-as.matrix(tax)

	for(i in rownames(tax)){
		temp.line<-tax[i,]
		temp.line[1]<-paste("k",temp.line[1],sep="_")
		temp.line[2]<-paste("p",temp.line[2],sep="_")
		temp.line[3]<-paste("c",temp.line[3],sep="_")
		temp.line[4]<-paste("o",temp.line[4],sep="_")
		temp.line[5]<-paste("f",temp.line[5],sep="_")
		temp.line[6]<-paste("g",temp.line[6],sep="_")
		#print(temp.line)
		#temp.line<- tax[1321,]
		#if( length(grep("T", is.na(temp.line))) > 0 ){ #### Looks for the NA (missing taxonomic labels)
		if(length(grep("_NA",  temp.line )) > 0 ){ #### Looks for the NA (missing taxonomic labels)
			na.ind<-grep("_NA",  temp.line )  ### Locates in the row where are such values
			na.lev<-min( na.ind )                   ### Locates the first NA taxonomic label
			last.tax<-temp.line[na.lev-1]           ### Takes the last-known taxonomic category
			###last.tax<-as.matrix(last.lev)[1]         
			temp.line[ na.ind ] <- rep( paste( "uc",last.tax, sep="_" ),length(temp.line[ na.ind ]) ) ### Add an uc_LAST_KNOWN_TAX_LEVEL to all the unkwon taxonomic levels, ie: 
			tax[i,]<-temp.line     ### Change the line for the one whith the new taxonomic info
		} else{
			tax[i,]<-temp.line     ### Change the line for the one whith the new taxonomic info
		}
			
	}

	OTU = otu_table(data, taxa_are_rows = TRUE)
	TAX = tax_table(tax[rownames(OTU),])
	physeq = phyloseq(OTU, TAX)
	return(physeq)

}

# CLR.transformation <- function(in.table = data.frame() ,min.reads=0, min.prop = 0.001, cutoff = 0){
# 	##### CLR transformation #######
# 	# It uses all the genus data to perfom the transformation
# 	min.reads <- 0  #all will be kept
# 	min.prop = 0.001 #OTUs with abundance of at least 0.001 (default) #this shouldn't have too much of an effect, and this way we avoid imputing too many 0s
# 	cutoff = 0  #i'll filter low abundant taxa afterwards (don't want to remove them from the matrix before doing the clr transformation)
# 
# 	matrix = t(in.table) 
# 	matrix.f = codaSeq.filter(matrix, min.reads=min.reads, min.occurrence=cutoff, min.prop=min.prop, samples.by.row=TRUE) 
# 
# 	matrix.f.n0 = estimate0.min(matrix.f) ### < - 
# 	matrix.f.n0.clr <- codaSeq.clr(matrix.f.n0, samples.by.row=F)
# 	pdf("CLR_hist.pdf")
# 		hist(rowSums(matrix.f.n0.clr))
# 	dev.off()
# 	in.genus.clr<-matrix.f.n0.clr
# 	return(in.genus.clr)
# }


CLR.transformation <- function(in.table = data.frame() ,min.reads=0, min.prop = 0.001, cutoff = 0, GMPR_aproximation =  F, GMPR_aproximation_log=F, GMPR_aproximation_square_root=T, min_ct = 2, intersect_no = 4){

	##### Two different ways to estimate the CLR transformation depending on the sparseness of the data, 
	##### high sparseness data should be normalized using the method of the 

	##### CLR transformation #######
	# It uses all the genus data to perfom the transformation
	min.reads <- min.reads  #all will be kept
	min.prop = min.prop #OTUs with abundance of at least 0.001 (default) #this shouldn't have too much of an effect, and this way we avoid imputing too many 0s
	cutoff = 0  #i'll filter low abundant taxa afterwards (don't want to remove them from the matrix before doing the clr transformation)

	if( GMPR_aproximation ==  TRUE ){

		#### Method for very sparse data sets 
		### Returns the median of all the different geometric mean samples
		# "We first estimate the size factors based on the OTU-level data and 
		# the genus-level counts are divided by the size factors to produce normalized genus-level abundances. "

		# "In contrast, the size factor-based approaches are capable of capturing the invariant part of the taxa 
		# counts and address the compositional problem efficiently through normalization by the size factors. 
		# The size factors could be naturally included as offsets in count-based parametric models to address 
		# uneven sequencing depth (Chen et al., 2018)."
		# Application 1: Counts are normalized by size factors to reduce the variation due to different library sizes
		# The normalized counts are subject to further downstream analysis such as ordination (PCA, PCoA), clustering,
		# and other multivariate methods. Note that further data transformation such as VST transformation (DESeq2)
		# may be needed in order to reveal patterns. Here shows an example of BC distance based ordination 
		#### paramenters:
		##### OTU table, colnames samples rownames species
		# OTUmatrix: An OTU table matrix, where OTUs arranged in columns and samples in rows.
		# min_ct: The minimal number of OTUs counts. Default 2.
		# intersect_no: The minimal number of shared OTUs between samples.  Default 4.

		# For DAA, one major challenge is to address the compositional problem. Rarefaction has a limited ability in this regard since the total sum constraint still exists after rarefaction. In addition, it suffers from a great power loss due to the discard of a large number of reads (McMurdie & Holmes, 2014). In contrast, the size factor-based approaches are capable of capturing the invariant part of the taxa counts and address the compositional problem efficiently through normalization by the size factors. The size factors could be naturally included as offsets in count-based parametric models to address uneven sequencing depth (Chen et al., 2018).



		cat("GMPR aproximation for sparse data","\n")
		cat("OTUmatrix: An OTU table matrix, where OTUs arranged in columns and samples in rows","\n")
		cat( "Colnames: ",colnames(t(in.table))[1:4], "\n")
		cat( "Rownames: ",rownames(t(in.table))[1:4], "\n")
		cat("Otherwise, transpose the matrix","\n")

		##  The first step is to calculate rij, which is the median count ratio of nonzero counts between sample i and j, 
		## The second step is to calculate the size factor si for a given sample i. It is the geometric mean of all the pairwise ratios
		gmpr.size.factor <- GMPR(t(in.table), min_ct = min_ct, intersect_no = intersect_no)
		in.genus.GMPR <- t(in.table) / gmpr.size.factor
		if(GMPR_aproximation_log==TRUE){
			cat("GMPR Log normalization ","\n")
			in.genus.GMPR<-log(in.genus.GMPR + 0.0000001)
		}
		if(GMPR_aproximation_square_root==TRUE){
			cat("GMPR Log normalization ","\n")
			in.genus.GMPR<-sqrt(in.genus.GMPR)
		}
		return(in.genus.GMPR)
		
	} else{
	

		matrix = t(in.table) 
		matrix.f = codaSeq.filter(matrix, min.reads=min.reads, min.occurrence=cutoff, min.prop=min.prop, max.prop=1.1, samples.by.row=TRUE) 
	
		matrix.f.n0 = estimate0.min(matrix.f) ###
		matrix.f.n0.clr <- codaSeq.clr(matrix.f.n0, samples.by.row=F)
		pdf("CLR_hist.pdf")
			hist(rowSums(matrix.f.n0.clr))
		dev.off()
		in.genus.clr<-matrix.f.n0.clr

		return(in.genus.clr)
	}
}


phyloseq2table<-function(phyloseq_file = ""){
	G.table<-load(phyloseq_file)
	if(G.table == "join.physeq.genus.rar"){
		Genus.rar<-data.frame( otu_table(join.physeq.genus.rar) )
		tax_vec_genus<-tax_table(join.physeq.genus.rar)[,6]
		rownames(Genus.rar) <- tax_vec_genus[rownames(Genus.rar)]
		return(Genus.rar)
	} else if(G.table == "join.physeq.genus"){
		Genus<-data.frame( otu_table(join.physeq.genus) )
		tax_vec_genus<-tax_table(join.physeq.genus)[,6]
		rownames(Genus) <- tax_vec_genus[rownames(Genus)]
		return(Genus)
	}
}


read.infile.data <- function(infileGenus = "", dada2taxfile="", dada2file="",format="",  min.rar = TRUE,min.num.seqs = 10000, Transpose_Matrix=F){

	min.rar.num<-min.num.seqs 
	if(format=="table"){
		# ########################### Genus table format #########################
		Genus<-read.table(infileGenus, header=T, row.names=1, dec=".", sep="\t")
		if(Transpose_Matrix==T){Genus<-t(Genus)}
		if(min.rar == TRUE){min.rar.num = min(colSums(Genus))}
		Genus.rar <- t(rrarefy(t(Genus), sample=min.rar.num))

	} else if(format=="phyloseq"){
		# ########################### phyloseq format #########################
		Genus<-phyloseq2table(phyloseq_file)
		if(min.rar == TRUE){min.rar.num = min(colSums(Genus))}
		Genus.rar <- t(rrarefy(t(Genus), sample=min.rar.num))
	} else if(format=="DADA2"){
		# ########################### DADA2 format #########################
		genus.abund.list<-dada2genusTable( abundance.file=dada2file, tax.file=dada2taxfile, min.num.seqs = min.num.seqs, min.rar = min.rar, plot_rar = TRUE )
		Genus<-genus.abund.list[["Genus"]]
		Genus.rar<-genus.abund.list[["Genus.rar"]]
	}
	ret.list<-list(Genus, Genus.rar)
	names(ret.list)<-c("Genus","Genus.rar")
	return(ret.list)
	
}


#################################
######### Filter FUNCTIONS ######
#################################

filter.table <-function(in.table = "", min.percentage = 0.01){
	in.table.prop <- (prop.table(as.matrix(in.table),2)) * 100
	median.vec<-apply(in.table.prop,1,function(x){mean(x)})
	median.vec<-median.vec[median.vec >= min.percentage]
	out.table<-t(in.table[names(median.vec),])
	return(out.table)
}


filter.prevalence<-function(in.table = "", max.prev =0.8){
	print(paste("You have",nrow(in.table),"samples. If this is not correct, transpose matrix!"))
	print(head(rownames(in.table)))
	esp.prev<-apply(in.table, 2, function(x){
			prev<-(length(grep("TRUE",x<=0))/length(x))
			return(prev)})

	in.table<-in.table[ , names(esp.prev[esp.prev <=max.prev])  ]
	return(in.table)
}

estimate0.min = function(matrix.test) { 
  print(paste("You have",ncol(matrix.test),"samples. If this is not correct, transpose matrix!"))
  matrix.test.p = t(t(matrix.test)/rowSums(t(matrix.test)))
  samplesums = colSums(matrix.test)
  
  matrix.f.n0 = matrix.test
  for (i in 1:nrow(matrix.f.n0)) {
    min = min(matrix.test.p[i,][matrix.test.p[i,] > 0])
    for (j in 1:ncol(matrix.f.n0 )) {
      if (matrix.f.n0 [i,j] == 0)
        matrix.f.n0 [i,j] = min*samplesums[j]
    }
  }
  return(matrix.f.n0)
}

median_per_condition<-function(in.mat="",condition=""){
	
	if( ncol(in.mat) != length(condition)){
		paste0("Number of columns ", ncol(in.mat)," do not match with the conditions ",length(condition))
	}
	names(condition)<-colnames(in.mat)
	cond2comp <- unique(sort(condition))
	out.median<-matrix(0,nrow(in.mat),length(cond2comp))
	rownames(out.median) <- rownames(in.mat)
	colnames(out.median) <- cond2comp
	for(i in 1:nrow(in.mat)){
		for( cond in cond2comp ){		
			#out.median[i,cond] <- median(in.mat[i,grep(cond,condition)])
			out.median[i,cond] <- median( in.mat[i,condition==cond] )
		}
	}
	return(out.median)
}

##############################################
######### ADONS and beta diversity FUNCTIONS #
##############################################


adonis.func<-function(sub.metadata = "" , temp.data.dist="", permutations = 1000){
	names.metadata<-colnames(sub.metadata)
	matrix.res <- matrix("",length(names.metadata), 4)
	rownames(matrix.res) <-names.metadata
	for(i in names.metadata ){
		#i<-"Succes"
		if(any(is.na(sub.metadata[,i])) == TRUE){
			no.naSamples<-rownames( sub.metadata[!is.na(sub.metadata[,i]) , ] )
			temp.matrix<-as.matrix(temp.data.dist)
			temp.data.dist.noNas<-as.dist(temp.matrix[no.naSamples,no.naSamples])
			temp.sub.metadta.noNAs<-sub.metadata[no.naSamples , i ]

			adononis.res<-adonis(temp.data.dist.noNas~temp.sub.metadta.noNAs,data=sub.metadata,permutations = permutations)
			p.value.temp<-adononis.res$aov.tab[ 6]$Pr[1]
			R2.temp<-adononis.res$aov.tab$R2[1]
			FModel.temp<-adononis.res$aov.tab$F.Model[1]
			matrix.res[i,] <-c(i,FModel.temp,R2.temp,p.value.temp)
			cat(i,"  ",FModel.temp,"  ",R2.temp,"  ",p.value.temp,"\n")
			next
			
		}
		
		adononis.res<-adonis(temp.data.dist~sub.metadata[,i],data=sub.metadata,permutations = permutations)
		p.value.temp<-adononis.res$aov.tab[ 6]$Pr[1]
		R2.temp<-adononis.res$aov.tab$R2[1]
		FModel.temp<-adononis.res$aov.tab$F.Model[1]
		matrix.res[i,] <-c(i,FModel.temp,R2.temp,p.value.temp)
		cat(i,"  ",FModel.temp,"  ",R2.temp,"  ",p.value.temp,"\n")
	
	}
	colnames(matrix.res) <- c("var","FModel","R2","adonis p-value")
	return(matrix.res)
}

ADONIS_func <- function(in.matrix="",Distance="",in.Metadata="",prefix="", permutations = 1000){

	in.matrix<-in.matrix[rownames(in.Metadata),]
	data.dist<-vegdist(in.matrix, method=Distance)

	#res<-adonis(data.dist~Day:Feed ,data= in.matrix ,permutations = 2000)
	#adonisDayFeed<-cbind(paste( "Join var", rownames(res$aov.tab)[1:1]), res$aov.tab$F.Model[1:1], res$aov.tab$R2[1:1], res$aov.tab[1:1, 6])
	#rownames(adonisDayFeed) <- adonisDayFeed[,1]

	#res<-adonis(data.dist~Day:PJ_cycles ,data= in.Metadata[rownames(as.matrix(data.dist)),] ,permutations = 2000)
	#adonisDayPJ_cycles<-cbind(paste( "Join var", rownames(res$aov.tab)[1:1]), res$aov.tab$F.Model[1:1], res$aov.tab$R2[1:1], res$aov.tab[1:1, 6])
	#rownames(adonisDayPJ_cycles) <- adonisDayPJ_cycles[,1]

	#res<-adonis(data.dist~Day:Butyrate ,data= in.Metadata[rownames(as.matrix(data.dist)),] ,permutations = 2000)
	#adonisDayButyrate<-cbind(paste( "Join var", rownames(res$aov.tab)[1:1]), res$aov.tab$F.Model[1:1], res$aov.tab$R2[1:1], res$aov.tab[1:1, 6])
	#rownames(adonisDayButyrate) <- adonisDayButyrate[,1]

	#res<-adonis(data.dist~Day:IL2 ,data= in.Metadata[rownames(as.matrix(data.dist)),] ,permutations = 2000)
	#adonisDayIL2<-cbind(paste( "Join var", rownames(res$aov.tab)[1:1]), res$aov.tab$F.Model[1:1], res$aov.tab$R2[1:1], res$aov.tab[1:1, 6])
	#rownames(adonisDayIL2) <- adonisDayIL2[,1]

	res.table.adonis<-adonis.func(sub.metadata = in.Metadata[rownames(as.matrix(data.dist)),] , temp.data.dist=data.dist)
	#res.table.adonis<-rbind(  res.table.adonis)
	BH.adj.p.value<-p.adjust(as.numeric(res.table.adonis[,4]),method="BH")
	res.table.adonis<-cbind(res.table.adonis,BH.adj.p.value)

	colnames(res.table.adonis) <- c("Variable","Fmodel","R2","p-value","BH adj p-value")

	adonis_table_tile<-paste0(prefix,"_adonis.tsv")
	write.table(res.table.adonis,adonis_table_tile,col.names=T,row.names = F,quote=FALSE,sep = "\t")
	return(data.frame(res.table.adonis))
}



median.table<-function( matrix.inf=matrix.inf, metadata=Metadata){ 
	median.matrix.inf<-matrix(0,0,dim(matrix.inf)[2])
	vec.names<-c()

	for( i in unique(sort(Metadata$Day)) ){
		#if(i==1){next}
		subdf<-subset(Metadata, Day==i)

		if( length(grep( "TRUE" , rownames(matrix.inf) %in% rownames(subdf))) == 0 ){next}

		vec<-apply( matrix.inf[ rownames(subdf) , ] , 2 , median )
		median.matrix.inf<-rbind(median.matrix.inf , vec)
		vec.names<-c(vec.names,i)
	}
	rownames(median.matrix.inf) <- vec.names

	### Remouve low abundace featires, It must have at least 4 values
	above.zero.features<-apply( median.matrix.inf, 2 , function(x) length(x[x!=0]) > 3 )
	median.matrix.inf<-median.matrix.inf[,above.zero.features]
	return(median.matrix.inf)
}

timePointdiff<-function( TimeSeries=""){   #### Input, vector data 
	ret.diff<-c()
	for( ind in 2:length(TimeSeries) ) {
		if(ind==2){
			abs(TimeSeries[1]-TimeSeries[ind])
			ret.diff<-c(ret.diff, abs(TimeSeries[1]-TimeSeries[ind]))
			next
		}
		ret.diff<-c(ret.diff, abs(TimeSeries[ind-1]-TimeSeries[ind]))
		#print(TimeSeries[ind])
	}
	return(ret.diff)
}


###### Format of the boxplot
format2boxplot<-function( data.table = "", OmicName=""){
	if(is.vector(data.table)){
		df.ret<-data.frame( Omic=rep(OmicName,length(data.table)), Varaible=names(data.table), adjr2=data.table )

	}
	if(is.data.frame(data.table)){
		df.ret<-cbind( rep(OmicName,nrow(data.table)), rownames(data.table), data.table )
	}

	colnames(df.ret) <- c("Omic","Varaible","adjr2")
	return(df.ret)

}




###### Format capscale_cum_variance
capscale_cum_variance <- function(in.Metadata ="", in.matrix="", Distance="", prefix = "", permutations = 1000, adj.pval.cutof = 0.1){
	Distance <<- Distance
	in.matrix<-in.matrix[rownames(in.Metadata),]

	all <- c()   #capscale(X ~ Y + Condition(Z))
	for (i in 1:ncol(in.Metadata)) { 
	  capsc <- capscale(in.matrix ~ in.Metadata[,i], distance = Distance, na.action=na.omit, permutations= 999999) #add=TRUE
	  an <- anova.cca(capsc, permutations = permutations) #permutations =999999 takes long
	  pval <- an["Pr(>F)"][[1]][[1]]
	  Fa <- an["F"][[1]][[1]]
	  r2 <- RsquareAdj(capsc)[[1]] #RsquaredAdj gives 2 ouputs [1] r.squared, [2] adj.r.squared
	  adjr2 <- RsquareAdj(capsc)[[2]]
	  all <- rbind(all,cbind(Fa,r2,adjr2,pval))
	  print(all)
	}

	colnames(all) <- c("F","r2","adjr2","p.value")
	row.names(all) <- colnames(in.Metadata)
	qval <- p.adjust(as.numeric(all[,"p.value"]),method="BH")
	all<-data.frame(cbind(all,qval))
	allRet<-all
	write.table(all,paste0(prefix,"_capscale_result.txt"),quote=F,col.names=TRUE,row.names=TRUE,sep="\t")
	in.Metadata.sub<-data.frame(in.Metadata[,rownames(all[all[,5]<adj.pval.cutof,])] )
	colnames(in.Metadata.sub) <- rownames(all[all[,5]<adj.pval.cutof,])
	rownames(in.Metadata.sub)<- rownames(in.Metadata)

	if(dim(all[all[,5]<adj.pval.cutof,])[1] == 0){
		print("Non significant variables")
		return("Non significant variables")
		

	}else if(dim(all[all[,5]<adj.pval.cutof,])[1] == 1){
		all<-all[all[,5]<adj.pval.cutof,]
	        #R2.adj Df      AIC        F    Pr..F.
		#cond 0.05528187  2 366.0042 3.165121 9.999e-05
		# F  r2  adjr2  p.value  qval
		vaf2plot<-c(as.numeric(all["adjr2"]),NA,NA,NA,as.numeric(all["qval"]))
		names(vaf2plot)<-c("R2.adj","Df","AIC","F","Pr..F.")
		vaf2plot<-t(data.frame(vaf2plot))
		rownames(vaf2plot) <- rownames(all)
		write.table(vaf2plot,paste0(prefix,"_non_redundant_variables.tsv"),col.names=T,row.names = T,quote=FALSE,sep = "\t")
		retList <- list(vaf2plot, allRet)
		names(retList) <- c("non_redundant","all")
		return( retList )
	}else{
	 
		#print("Pass the single capscale")
		# forward stepwise regression
		attach(in.Metadata.sub)
		str(in.Metadata.sub)
		dim(in.matrix)  ### numbers of rows should be the same!!!!!!
		sink(paste0(prefix,"_forward_stepwise_regression_results.txt"),append = FALSE, type = c("output", "message"), split = FALSE)
		# forward stepwise regression with all the variables without frozen cell count
		mod0=capscale(in.matrix ~ 1, data=in.Metadata.sub, distance=Distance) # Model with intercept only. 
		mod0
		#cat("Pass the single capscale mod0 ",Distance,"\n")

		mod1=capscale(in.matrix ~ ., data=in.Metadata.sub, distance=Distance,  na.action = na.exclude) # Model with all explanatory variables; add=T if u want to remove the (-) eigenvalues, u can add a constant
		mod1
		#cat("Pass the single capscale mod1 ",Distance,"\n")
		step.res<-ordiR2step(mod0, scope=formula(mod1), data=in.Metadata.sub, direction="forward", Pin = 0.1, R2scope = FALSE, permutations = permutations,
			trace = T, na.action=na.exclude ) #stepwise forward option (adding variables until best model)
		cat("Pass the single capscale ordiR2step ",Distance,"\n")
		anova(step.res,permutations = permutations) # Summary table for best model
		step.res$anova
		sink()
	
		vaf2plot<-data.frame(step.res$anova)
		vaf2plot<-subset(vaf2plot, Pr..F. <= adj.pval.cutof)
		rownames(vaf2plot)<-gsub("[+] ","",rownames(vaf2plot))
		write.table(vaf2plot,paste0(prefix,"_non_redundant_variables.tsv"),col.names=T,row.names = T,quote=FALSE,sep = "\t")

		retList <- list(vaf2plot, allRet)

		names(retList) <- c("non_redundant","all")
		return( retList )
	}
}

##### Non-redundance variance estimation
nonredundant_R2_matrix <- function(namesvars="",in.table=""){
	vec_var<-in.table[,1]
	All_vars<-matrix(0,length(namesvars)+1,2); rownames(All_vars) <- c(namesvars,"Unconstrained")
	colnames(All_vars) <- c("R2","cumR2")
	vec_var <- c(vec_var,1)
	names(vec_var) <- c( rownames(in.table) , "Unconstrained") 
	names_var<-names(vec_var)

	nonredundant_covariates <- c(vec_var[1])
	
	for(i in 2:length(vec_var)){
		#print(vec_var[i])
		diff_var<- vec_var[i] - vec_var[i-1]
		#cat(names_var[i-1], " ",  vec_var[i-1], " " ,  names_var[i], " ", vec_var[i] , "Non Redundant var ", diff_var ,"\n")
		nonredundant_covariates<-c(nonredundant_covariates, diff_var)
		
	}
	res.matrix<-cbind(nonredundant_covariates, vec_var)
	All_vars[rownames(res.matrix),] <- res.matrix
	return(All_vars)
}


nonredundant2plot <- function(non.redundant="",all.var=""){

	non_redundant<-data.frame(non.redundant)
	all_var<-all.var[order(all.var$adjr2,decreasing=T),]

	####################################################
	####### Non redundant effect all variables  ########
	####################################################

	ret_matrix <- matrix(0,0,3)
	colnames(ret_matrix) <- c("Variable","Variance","R2")
	#rownames(ret_matrix) <- rownames(all_var)
	sum<-0
	for(index in 1:length(rownames(all_var))){

		i<-rownames(all_var)[index]
		print(i)
		sum<-sum+1  		
		if(sum==1){
			temp_ret_matrix<-c(i, "R2", all_var[i,"adjr2"])
			ret_matrix<-rbind(ret_matrix,temp_ret_matrix)
			temp_ret_matrix<-c(i, "Cum_R2", all_var[i,"adjr2"])
			ret_matrix<-rbind(ret_matrix,temp_ret_matrix)

			temp_ret_matrix<-c(i, "Unconstrain",(1- all_var[i,"adjr2"]) )
			ret_matrix<-rbind(ret_matrix,temp_ret_matrix)
			 			

			next
		}
		temp_ret_matrix<-c(i, "R2", all_var[i,"adjr2"])
		ret_matrix<-rbind(ret_matrix,temp_ret_matrix)

		if( i %in% rownames(non_redundant) == FALSE ){ ##### The total varaince is estimated by substracting the previous time point varainec to the current time-point
			#cumVar=
			prev.var <- rownames(all_var)[index-1]
			ind.ex<-grep(paste0(prev.var,"$"), ret_matrix[,"Variable"] )
			tempDF<-data.frame( ret_matrix[ind.ex,] )
			oldR2<-as.character(subset( tempDF, Variance == "Cum_R2" )$R2)
			
			temp_ret_matrix<-c(i, "Cum_R2", oldR2 )
			ret_matrix<-rbind(ret_matrix,temp_ret_matrix)

			temp_ret_matrix<-c(i, "Unconstrain",(1- as.numeric(oldR2) ) )
			ret_matrix<-rbind(ret_matrix,temp_ret_matrix)
		
			next
		}

		temp_ret_matrix<-c(i, "Cum_R2", non_redundant[i,"R2.adj"])
		ret_matrix<-rbind(ret_matrix,temp_ret_matrix)	
		##### The non-redundant varaince is estimated by substracting 1 - the constrain varaince
		temp_ret_matrix<-c(i, "Unconstrain",(1- non_redundant[i,"R2.adj"] ) )
		ret_matrix<-rbind(ret_matrix,temp_ret_matrix)			

	}
	non_redundant_covar<-data.frame(Variable=ret_matrix[,1], Variance = ret_matrix[,2], R2=as.numeric(ret_matrix[,3]) ) 
	non_redundant_covar$Variance <-factor(non_redundant_covar$Variance, levels = c("R2", "Cum_R2","Unconstrain"))
	non_redundant_covar$Variable <-factor(non_redundant_covar$Variable, levels = rownames(all_var) )

	##############################################################
	####### Non redundant effect significative variables  ########
	##############################################################
	
	ret_matrix_sig <- matrix(0,0,3)
	colnames(ret_matrix_sig) <- c("Variable","Variance","R2")
	#rownames(ret_matrix) <- rownames(all_var)
	sum<-0
	for(index in 1:length(rownames(non_redundant))){
		i<-rownames(non_redundant)[index]
		print(i)
		sum<-sum+1  		
		if(sum==1){
			temp_ret_matrix<-c(i, "R2", non_redundant[i,"R2.adj"])
			ret_matrix_sig<-rbind(ret_matrix_sig,temp_ret_matrix)
			temp_ret_matrix<-c(i, "Cum_R2", non_redundant[i,"R2.adj"])
			ret_matrix_sig<-rbind(ret_matrix_sig,temp_ret_matrix)
			temp_ret_matrix<-c(i, "Unconstrain",(1- non_redundant[i,"R2.adj"]) )
			ret_matrix_sig<-rbind(ret_matrix_sig,temp_ret_matrix)
			 			

			next
		}
		temp_ret_matrix<-c(i, "R2", all_var[i,"adjr2"])
		ret_matrix_sig<-rbind(ret_matrix_sig,temp_ret_matrix)

		temp_ret_matrix<-c(i, "Cum_R2", non_redundant[i,"R2.adj"])
		ret_matrix_sig<-rbind(ret_matrix_sig,temp_ret_matrix)	
		##### The non-redundant varaince is estimated by substracting 1 - the constrain varaince
		temp_ret_matrix<-c(i, "Unconstrain",(1- non_redundant[i,"R2.adj"] ) )
		ret_matrix_sig<-rbind(ret_matrix_sig,temp_ret_matrix)			
	}
	non_redundant_covar_sig<-data.frame(Variable=ret_matrix_sig[,1], Variance = ret_matrix_sig[,2], R2=as.numeric(ret_matrix_sig[,3]) ) 
	non_redundant_covar_sig$Variance <-factor(non_redundant_covar_sig$Variance, levels = c("R2", "Cum_R2","Unconstrain"))
	non_redundant_covar_sig$Variable <-factor(non_redundant_covar_sig$Variable, levels = rownames(non_redundant) )

	ret_list<-list(non_redundant_covar,non_redundant_covar_sig)
	names(ret_list) <- c("nr","nr_sig")
	return(ret_list)	
}




