library("wesanderson")
library(grDevices)
library("corrplot")
corrplot.bis<-corrplot
library(vegan)
library(ade4)
#library(Matrix)
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
library(ggrepel)
library(microbiome)
set.seed(1234)
#################################
######### Read infile FUNCTIONS #
#################################

filter_otu_table <- function(in_phylo = phyloseq(), Samples = data.frame(), prev = 0.2, Tax_level = "Genus"){

	in_phylo <- prune_samples( Samples , in_phylo)
	Prevalence <- prevalence(in_phylo)
	in_phylo_filter <- prune_taxa( names(Prevalence[Prevalence >= prev]) , in_phylo)
	table <- tax_table(in_phylo_filter)
	otable <- as.matrix( otu_table(in_phylo_filter))
	otable <-  otable[ match( rownames(table)  ,  rownames(otable)  ) , ]
	if(Tax_level == "Genus"){
		rownames(otable) <- as.character(table[,6])	
	}
	if(Tax_level == "Species"){
		rownames(otable) <- as.character(table[,7])	
	}
	taxa_matrix <- t(otable)
	taxa_matrix <-taxa_matrix[match( Samples , rownames(taxa_matrix)),]
	return(taxa_matrix)
}

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
		gmpr.size.factor <- GMPR(data.frame(t(in.table)), min_ct = min_ct, intersect_no = intersect_no)
		in.genus.GMPR <- t(in.table) / gmpr.size.factor
		if(GMPR_aproximation_log==TRUE){
			cat("GMPR Log normalization ","\n")
			#in.genus.GMPR<-log(in.genus.GMPR + 0.0000001)
			in.genus.GMPR<-log(in.genus.GMPR + 1)
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


filter_low_prevalence <- function(in.table = matrix(), Categories_groups = data.frame(), max_percentage_0 = 80 ){

	cat("\n")	
	cat("Samples ",colnames(in.table)[1:3], "\n")
	cat("Taxas ",rownames(in.table)[1:3], "\n")

	if(length(Categories_groups)== 0){
		cat("Max percentage of zero values per taxa = ", max_percentage_0,"%","\n\n")
		percentage0<-apply(in.table == 0, 1, function(x){Num0<-length(grep("TRUE",x));total<-length( x );percentage <- (Num0/total) * 100})
		in.table <- in.table[names(percentage0[percentage0 < max_percentage_0]),]
		rownames(in.table)<-gsub("/",".",rownames(in.table))
		return(in.table)
	} else{

		cat("\n\nIf the filtering of prevalence is done by groups, those taxa will be taken in which at least one of the groups, for said taxa, has at least a lower prevalence of zeros than the cutoff. It is advisable to go down to avoid taking taxa with many zeros\n\n")


		Categories_groups$cond <- Categories_groups[,1]
		cat("Num categoreis =" ,length(levels(Categories_groups$cond)), "Categories =", levels(Categories_groups$cond),"\n\n")
		matrix0percentage <- matrix(0, dim(in.table)[1] , length(levels(Categories_groups$cond)) )
		colnames(matrix0percentage) <- levels(Categories_groups$cond)
		rownames(matrix0percentage) <- rownames(in.table)

		for(i in levels(Categories_groups$cond)){
			print(i)
			tdf <- rownames(subset(Categories_groups, cond == i))
			percentage0<-apply(in.table[,tdf] == 0, 1, function(x){Num0<-length(grep("TRUE",x));total<-length( x );percentage <- (Num0/total) * 100})
			matrix0percentage[names(percentage0),i] <- percentage0
		}
		#### Now subset filter all the taxas
		congruenceCat<- apply(matrix0percentage<max_percentage_0,1,function(x){AllPer<-any(x);return(AllPer)})
		in.table<-in.table[names(congruenceCat[congruenceCat]),]
		rownames(in.table)<-gsub("/",".",rownames(in.table))
		return(in.table)

	}
}


filter_low_prevalence <- function(in.table = matrix(), Categories_groups = data.frame(), max_percentage_0 = 80 ){

	cat("\n")	
	cat("Samples ",colnames(in.table)[1:3], "\n")
	cat("Taxas ",rownames(in.table)[1:3], "\n")

	if(length(Categories_groups)== 0){
		cat("Max percentage of zero values per taxa = ", max_percentage_0,"%","\n\n")
		percentage0<-apply(in.table == 0, 1, function(x){Num0<-length(grep("TRUE",x));total<-length( x );percentage <- (Num0/total) * 100})
		in.table <- in.table[names(percentage0[percentage0 < max_percentage_0]),]
		rownames(in.table)<-gsub("/",".",rownames(in.table))
		return(in.table)
	} else{

		cat("\n\nIf the filtering of prevalence is done by groups, those taxa will be taken in which at least one of the groups, for said taxa, has at least a lower prevalence of zeros than the cutoff. It is advisable to go down to avoid taking taxa with many zeros\n\n")


		Categories_groups$cond <- Categories_groups[,1]
		cat("Num categoreis =" ,length(levels(Categories_groups$cond)), "Categories =", levels(Categories_groups$cond),"\n\n")
		matrix0percentage <- matrix(0, dim(in.table)[1] , length(levels(Categories_groups$cond)) )
		colnames(matrix0percentage) <- levels(Categories_groups$cond)
		rownames(matrix0percentage) <- rownames(in.table)

		for(i in levels(Categories_groups$cond)){
			print(i)
			tdf <- rownames(subset(Categories_groups, cond == i))
			percentage0<-apply(in.table[,tdf] == 0, 1, function(x){Num0<-length(grep("TRUE",x));total<-length( x );percentage <- (Num0/total) * 100})
			matrix0percentage[names(percentage0),i] <- percentage0
		}
		#### Now subset filter all the taxas
		congruenceCat<- apply(matrix0percentage<max_percentage_0,1,function(x){AllPer<-any(x);return(AllPer)})
		in.table<-in.table[names(congruenceCat[congruenceCat]),]
		rownames(in.table)<-gsub("/",".",rownames(in.table))
		return(in.table)

	}
}



##############################################
######### ADONS and beta diversity FUNCTIONS #
##############################################


adonis.func<-function(sub.metadata = "" , temp.data.dist="", permutations = 1000){

	sub.metadata <- sub.metadata[match(rownames(as.matrix(temp.data.dist)), rownames(sub.metadata)),]
	names.metadata<-colnames(sub.metadata)
	matrix.res <- matrix("",length(names.metadata), 5)
	rownames(matrix.res) <-names.metadata
	for(i in names.metadata ){
		# i<-"center"
		if(any(is.na(sub.metadata[,i])) == TRUE){
			no.naSamples<-rownames( sub.metadata[!is.na(sub.metadata[, colnames(sub.metadata) == i]) , ] )

			#### Filter the rows and columns in the distance matrix 
			temp.matrix<-as.matrix(temp.data.dist)
			ind.rows <- match(  no.naSamples   ,  rownames(temp.matrix)   )
			ind.col <- match(  no.naSamples   ,  colnames(temp.matrix)   )			
			temp.data.dist.noNas<-as.dist(temp.matrix[ind.rows,ind.col])
			
			#### Filter the rows and columns in the metadata matrix 
			temp.sub.metadta.noNAs<-sub.metadata[match(no.naSamples,rownames(sub.metadata)), ]
			temp.sub.metadta.noNAs<- temp.sub.metadta.noNAs[ ,colnames(temp.sub.metadta.noNAs) == i ]

			#### Adonis test 
			adononis.res<-adonis2(temp.data.dist.noNas ~ temp.sub.metadta.noNAs,permutations = permutations)
			
			#### Filter the results
			p.value.temp <- adononis.res$`Pr(>F)`[1]
			R2.temp <- adononis.res$R2[1]
			FModel.temp <- adononis.res$F[1]
			N <- length(temp.sub.metadta.noNAs)
			
			#### Return the matrix
			matrix.res[i,] <-c(i,FModel.temp,R2.temp,p.value.temp,N)
			cat(i,"  ",FModel.temp,"  ",R2.temp,"  ",p.value.temp,N,"\n")
			
			rm(no.naSamples,temp.matrix,ind.rows,ind.col,temp.data.dist.noNas,
				temp.sub.metadta.noNAs,adononis.res,p.value.temp,R2.temp,FModel.temp,N )
			next
			
		}
		
		#### Filter the rows and columns in the metadata matrix 
		temp.sub.metadta.noNAs <- sub.metadata[ ,colnames(sub.metadata) == i ]
		
		#### Adonis test 					
		adononis.res<-adonis2(temp.data.dist~temp.sub.metadta.noNAs,permutations = permutations)

		#### Filter the results
		p.value.temp <- adononis.res$`Pr(>F)`[1]
		R2.temp <- adononis.res$R2[1]
		FModel.temp <- adononis.res$F[1]
		N <- length(temp.sub.metadta.noNAs)

		#### Return the matrix		
		matrix.res[i,] <-c(i,FModel.temp,R2.temp,p.value.temp,N)
		cat(i,"  ",FModel.temp,"  ",R2.temp,"  ",p.value.temp,N,"\n")
		rm(temp.sub.metadta.noNAs,adononis.res,p.value.temp,R2.temp,FModel.temp,N )
	
	}
	colnames(matrix.res) <- c("var","FModel","R2","adonis p-value","N")
	return(matrix.res)
}

ADONIS_func <- function(in.matrix="",Distance="",in.Metadata="",prefix="", permutations = 1000){

	in.matrix<- in.matrix[match(rownames(in.Metadata), rownames(in.matrix)),]
	data.dist<-vegdist(in.matrix, method=Distance)

	if(ncol(in.Metadata) == 1){
		sample_names <- rownames(as.matrix(data.dist))
		colnames <- colnames(in.Metadata)
		in.Metadata <- data.frame(in.Metadata[sample_names,] )
		colnames(in.Metadata) <- colnames; rownames(in.Metadata) <- sample_names
	}else{
		in.Metadata <- in.Metadata[match(rownames(as.matrix(data.dist)), rownames(in.matrix)),]
	}

	res.table.adonis <- adonis.func(sub.metadata = in.Metadata, temp.data.dist=data.dist)
	#res.table.adonis<-rbind(  res.table.adonis)
	BH.adj.p.value<-p.adjust(as.numeric(res.table.adonis[,4]),method="BH")
	res.table.adonis<-cbind(res.table.adonis,BH.adj.p.value)

	colnames(res.table.adonis) <- c("Variable","Fmodel","R2","p-value","N","BH adj p-value")
	res.table.adonis <- res.table.adonis[, c("Variable","Fmodel","R2","p-value","BH adj p-value","N")]

	adonis_table_tile<-paste0(prefix,"_adonis.tsv")
	write.table(res.table.adonis,adonis_table_tile,col.names=T,row.names = F,quote=FALSE,sep = "\t")

	res.table.adonis <- data.frame(res.table.adonis)	
	res.table.adonis$BH.adj.p.value <- as.numeric(as.character(res.table.adonis$BH.adj.p.value))
	res.table.adonis$p.value <- as.numeric(as.character(res.table.adonis$p.value))
	res.table.adonis$R2 <- as.numeric(as.character(res.table.adonis$R2))
	res.table.adonis$Fmodel <- as.numeric(as.character(res.table.adonis$Fmodel))
	res.table.adonis$N <- as.numeric(as.character(res.table.adonis$N))
	res.table.adonis <- res.table.adonis[order(res.table.adonis$BH.adj.p.value),]
	
	return(res.table.adonis)
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
capscale_cum_variance <- function(in.Metadata ="", in.matrix="", Distance="", prefix = "", permutations = 1000, adj.pval.cutof = 0.1, ordiR2step_Pin = 0.1  ){
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
	colnames(in.Metadata.sub) <- rownames(all[all[,5]< adj.pval.cutof,])
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
		step.res<-ordiR2step(mod0, scope=formula(mod1), data=in.Metadata.sub, direction="forward", Pin = ordiR2step_Pin, R2scope = FALSE, permutations = permutations,
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

##############################################
######### ADONIS plot functions FUNCTIONS ####
##############################################


ord_labels <-  function(ord){
    ev <- vegan::eigenvals(ord)
    tol <- -(1e-07)*ev[1]
    ord.labels <- rep("", length(ev))
    if ((any(is.na(ev))) | (any(ev < tol))) {
      for ( i in 1:length(ev)) {
        ord.labels[i] <- paste("DIM", i, sep = "")
      }
    }
    else {
      ev.pc <- round(100*(ev/sum(ev)), 2)
      axis.names <- names(ev)
      if (is.null(axis.names)) {
        for ( i in 1:length(ev.pc)) {
          ord.labels[i] <- paste("DIM", i, " ", sprintf(ev.pc[i], fmt = '%#.1f'), "%", sep="")
        }
      } else {
        for (i in 1:length(ev.pc)){
          ord.labels[i] <- paste(axis.names[i], " ", ev.pc[i],"%", sep="")
        }
      }
    }
    return(ord.labels)
  }

scale_arrow <- function(arrows, data, at = c(0, 0), fill = 0.75) {
  u <- c(range(data[,1], range(data[,2])))
  u <- u - rep(at, each = 2)
  r <- c(range(arrows[, 1], na.rm = TRUE), range(arrows[, 2], na.rm = TRUE))
  rev <- sign(diff(u))[-2]

  if (rev[1] < 0) {
    u[1:2] <- u[2:1]
  }
  if (rev[2] < 0) {
    u[3:4] <- u[4:3]
  }
  u <- u/r
  u <- u[is.finite(u) & u > 0]
  invisible(fill * min(u))
}



ggplot_envfit <- function(df_ord,ord,fit,Metadata,alpha_pval,choices=c(1,2)){


	df_ord <- as.data.frame(df_ord)
	axis.labels <- ord_labels(ord)[choices]
	#xlab <- axis.labels[1]
	#ylab <- axis.labels[2]

	##### Create the arrows
	df_arrows <- as.data.frame(scores(fit, "vectors"))

	 dim(df_arrows)[1]	
	colnames(df_ord) <- c("x", "y")

	if(dim(df_arrows)[1] > 0){  
		mult <- scale_arrow(df_arrows, df_ord[, c("x", "y")])
		df_arrows <- mult * df_arrows
		df_arrows$var <- rownames(df_arrows)
		df_arrows$p.val <- fit$vectors$pvals
		df_arrows$q.val <- p.adjust(fit$vectors$pvals,method="BH")
		colnames(df_arrows) <- c("x", "y", "var", "p.val")
		df_arrows <- df_arrows[df_arrows$p.val <= alpha_pval, ]
	} else{
		df_arrows<-data.frame(x=NULL, y=NULL, var=NULL, p.val=NULL)
		#colnames(df_arrows) <- c("x", "y", "var", "p.val")
	}
	
	##### Create the factors
	df_factors <- as.data.frame(scores(fit, "factors"))
	mult <- scale_arrow(df_factors, df_ord[, c("x", "y")])
	df_factors <- mult * df_factors
	#df_factors$var <- gsub("Fake_metadata[.]","",rownames(df_factors)) ### You can modify the names of the variable
	df_factors$var <- rownames(df_factors) ### You can modify the names of the variable
	fit$factors$qvals <- p.adjust(fit$factors$pvals,method="BH")
	q_values <- c(rep(1, length(df_factors$var) ))
	for(i in names(fit$factors$qvals)){
		#df_factors$var[grep(i,df_factors$var)]
		q_values[grep(i,df_factors$var)] <- fit$factors$qvals[i]
	}

	df_factors$q.val <- q_values
	colnames(df_factors) <- c("x", "y", "var", "q.val")
	df_factors <- df_factors[df_factors$q.val <= alpha_pval, ]

	if(dim(Metadata)[2] == 1){
		sample_names <- rownames(df_ord)
		colnames <- colnames(Metadata)
		Metadata <- data.frame(Metadata[sample_names,] )
		colnames(Metadata) <- colnames; rownames(Metadata) <- sample_names
		df_ord<-cbind(df_ord,Metadata)
	}else{
		df_ord<-cbind(df_ord,Metadata[rownames(df_ord),])
	}




	ret_list <- list(df_ord,df_arrows,df_factors)
	names(ret_list) <- c("df_ord","df_arrows","df_factors")
	return(ret_list)

}

vegan_PCoA_envfit <- function(in.matrix=NULL, Metadata2enfit=NULL, choices = c(1,2),  scaling = 1, perm = 999, distance = "euclidean"){

#	PCoA<-vegan::dbrda(in.matrix ~ 1, dist=distance, metaMDS = F) ### you can choose between the euclidean and the bray (bray curtis) dissimilarity
#	xlab=paste("Comp1", format(100*PCoA$CA$eig[1]/sum(abs(PCoA$CA$eig)), digits=4), "%", sep=" ")
#	ylab=paste("Comp2", format(100*PCoA$CA$eig[2]/sum(abs(PCoA$CA$eig)), digits=4), "%", sep=" ")
	PCoA<-vegan::capscale(in.matrix ~ 1, dist=distance, metaMDS = F,na.rm = T) ### you can choose between the euclidean and the bray (bray curtis) dissimilarity
	Sum_dbrda<-summary(PCoA)
	xlab=paste("PC1:", format( 100*Sum_dbrda$cont$importance["Proportion Explained",1], digits=4), "%", sep=" ")
	ylab=paste("PC2:", format( 100*Sum_dbrda$cont$importance["Proportion Explained",2], digits=4), "%", sep=" ")
	fit <- vegan::envfit(PCoA, Metadata2enfit, choices = choices, perm = perm , na.rm = TRUE)
	df_ord <- vegan::scores(PCoA, display = "sites", choices = choices,  scaling = scaling)

	list_return<-list(PCoA,xlab,ylab,fit,df_ord)
	names(list_return) <- c("PCoA","xlab","ylab","fit","df_ord")
	return(list_return)
}


filter_low_prevalence <- function(in.table = matrix(), Categories_groups = data.frame(), max_percentage_0 = 80 ){

	cat("\n")	
	cat("Samples ",colnames(in.table)[1:3], "\n")
	cat("Taxas ",rownames(in.table)[1:3], "\n")

	if(length(Categories_groups)== 0){
		cat("Max percentage of zero values per taxa = ", max_percentage_0,"%","\n\n")
		percentage0<-apply(in.table == 0, 1, function(x){Num0<-length(grep("TRUE",x));total<-length( x );percentage <- (Num0/total) * 100})
		in.table <- in.table[names(percentage0[percentage0 < max_percentage_0]),]
		rownames(in.table)<-gsub("/",".",rownames(in.table))
		return(in.table)
	} else{

		cat("\n\nIf the filtering of prevalence is done by groups, those taxa will be taken in which at least one of the groups, for said taxa, has at least a lower prevalence of zeros than the cutoff. It is advisable to go down to avoid taking taxa with many zeros\n\n")


		Categories_groups$cond <- Categories_groups[,1]
		cat("Num categoreis =" ,length(levels(Categories_groups$cond)), "Categories =", levels(Categories_groups$cond),"\n\n")
		matrix0percentage <- matrix(0, dim(in.table)[1] , length(levels(Categories_groups$cond)) )
		colnames(matrix0percentage) <- levels(Categories_groups$cond)
		rownames(matrix0percentage) <- rownames(in.table)

		for(i in levels(Categories_groups$cond)){
			print(i)
			tdf <- rownames(subset(Categories_groups, cond == i))
			percentage0<-apply(in.table[,tdf] == 0, 1, function(x){Num0<-length(grep("TRUE",x));total<-length( x );percentage <- (Num0/total) * 100})
			matrix0percentage[names(percentage0),i] <- percentage0
		}
		#### Now subset filter all the taxas
		congruenceCat<- apply(matrix0percentage<max_percentage_0,1,function(x){AllPer<-any(x);return(AllPer)})
		in.table<-in.table[names(congruenceCat[congruenceCat]),]
		rownames(in.table)<-gsub("/",".",rownames(in.table))
		return(in.table)

	}
}
