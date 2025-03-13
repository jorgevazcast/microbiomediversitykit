###################################################################################
######################		IMPORTANT         #######################################################################################
###### IMPORTANT!!!! The Wilcox-test has been replaced by the T-student given the lows sample SIZE, this in the pairwise boxplots ######
#########################################################################################################################################

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
set.seed(1234)

#######################################################################################################
########################################### Functions #################################################
#######################################################################################################

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

CLR.transformation <- function(in.table = data.frame() ,min.reads=0, min.prop = 0.001, cutoff = 0){
	##### CLR transformation #######
	# It uses all the genus data to perfom the transformation
	min.reads <- 0  #all will be kept
	min.prop = 0.001 #OTUs with abundance of at least 0.001 (default) #this shouldn't have too much of an effect, and this way we avoid imputing too many 0s
	cutoff = 0  #i'll filter low abundant taxa afterwards (don't want to remove them from the matrix before doing the clr transformation)

	matrix = t(in.table) 
	matrix.f = codaSeq.filter(matrix, min.reads=min.reads, min.occurrence=cutoff, min.prop=min.prop, samples.by.row=TRUE) 

	matrix.f.n0 = estimate0.min(matrix.f) ### < - 
	matrix.f.n0.clr <- codaSeq.clr(matrix.f.n0, samples.by.row=F)
	pdf("CLR_hist.pdf")
		hist(rowSums(matrix.f.n0.clr))
	dev.off()
	in.genus.clr<-matrix.f.n0.clr
	return(in.genus.clr)
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


read.infile.data <- function(infileGenus = "", dada2taxfile="", dada2file="",format="",  min.rar = TRUE,min.num.seqs = 10000){

	min.rar.num<-min.num.seqs 
	if(format=="table"){
		# ########################### Genus table format #########################
		Genus<-read.table(infileGenus, header=T, row.names=1, dec=".", sep="\t")
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

adonis.func<-function(sub.metadata = "" , temp.data.dist=""){
	names.metadata<-colnames(sub.metadata)
	matrix.res <- matrix("",length(names.metadata), 4)
	rownames(matrix.res) <-names.metadata
	for(i in names.metadata ){
		
		p.value.temp<-adonis(temp.data.dist~sub.metadata[,i],data=sub.metadata,permutations = 2000)$aov.tab[ 6]$Pr[1]
		R2.temp<-adonis(temp.data.dist~sub.metadata[,i],data=sub.metadata,permutations = 2000)$aov.tab$R2[1]
		FModel.temp<-adonis(temp.data.dist~sub.metadata[,i],data=sub.metadata,permutations = 2000)$aov.tab$F.Model[1]
		matrix.res[i,] <-c(i,FModel.temp,R2.temp,p.value.temp)
		cat(i,"  ",FModel.temp,"  ",R2.temp,"  ",p.value.temp,"\n")
	
	}
	colnames(matrix.res) <- c("var","FModel","R2","adonis p-value")
	return(matrix.res)
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

ADONIS_func <- function(in.matrix="",Distance="",in.Metadata="",prefix=""){
	data.dist<-vegdist(in.matrix, method=Distance)
	in.matrix=in.Metadata[rownames(as.matrix(data.dist)),]

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

capscale_cum_variance <- function(in.Metadata ="", in.matrix="", Distance="", prefix = ""){
	Distance <<- Distance
	all <- c()   #capscale(X ~ Y + Condition(Z))
	for (i in 1:ncol(in.Metadata)) { 
	  capsc <- capscale(in.matrix ~ in.Metadata[,i], distance = Distance, na.action=na.omit, permutations= 999999) #add=TRUE
	  an <- anova.cca(capsc, permutations = 1000) #permutations =999999 takes long
	  pval <- an["Pr(>F)"][[1]][[1]]
	  Fa <- an["F"][[1]][[1]]
	  r2 <- RsquareAdj(capsc)[[1]] #RsquaredAdj gives 2 ouputs [1] r.squared, [2] adj.r.squared
	  adjr2 <- RsquareAdj(capsc)[[2]]
	  all <- rbind(all,cbind(Fa,r2,adjr2,pval))
	  print(all)
	}

	colnames(all) <- c("F","r2","adjr2","p-value")
	row.names(all) <- colnames(in.Metadata)
	qval <- p.adjust(all[,"p-value"],method="BH")
	all<-data.frame(cbind(all,qval))
	allRet<-all
	write.table(all,paste0(prefix,"_capscale_result.txt"),quote=F,col.names=TRUE,row.names=TRUE,sep="\t")
	in.Metadata.sub<-data.frame(in.Metadata[,rownames(all[all[,5]<0.1,])] )
	colnames(in.Metadata.sub) <- rownames(all[all[,5]<0.1,])
	rownames(in.Metadata.sub)<- rownames(in.Metadata)

	if(dim(all[all[,5]<0.1,])[1] == 0){
		print("Non significant variables")
		return("Non significant variables")
		

	}else if(dim(all[all[,5]<0.1,])[1] == 1){
		all<-all[all[,5]<0.1,]
	        #R2.adj Df      AIC        F    Pr..F.
		#cond 0.05528187  2 366.0042 3.165121 9.999e-05
		# F  r2  adjr2  p.value  qval
		vaf2plot<-c(as.numeric(all["adjr2"]),NA,NA,NA,as.numeric(all["qval"]))
		names(vaf2plot)<-c("R2.adj","Df","AIC","F","Pr..F.")
		vaf2plot<-t(data.frame(vaf2plot))
		rownames(vaf2plot) <- rownames(all)
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
		mod1=capscale(in.matrix ~ ., data=in.Metadata.sub, distance=Distance) # Model with all explanatory variables; add=T if u want to remove the (-) eigenvalues, u can add a constant
		mod1
		#cat("Pass the single capscale mod1 ",Distance,"\n")
		step.res<-ordiR2step(mod0, scope=formula(mod1), data=in.Metadata.sub, direction="forward", Pin = 0.1, R2scope = FALSE, permutations = 10000,trace = T) #stepwise forward option (adding variables until best model)
		cat("Pass the single capscale ordiR2step ",Distance,"\n")
		anova(step.res,permutations = 1000) # Summary table for best model
		step.res$anova
		sink()
	
		vaf2plot<-data.frame(step.res$anova)
		vaf2plot<-subset(vaf2plot, Pr..F. <= 0.1)
		rownames(vaf2plot)<-gsub("[+] ","",rownames(vaf2plot))
		vaf2plot[,1]
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

		if( i %in% rownames(non_redundant) == FALSE ){
			#cumVar=
			prev.var <- rownames(all_var)[index-1]
			ind.ex<-grep(prev.var, ret_matrix[,"Variable"] )
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

		temp_ret_matrix<-c(i, "Unconstrain",(1- non_redundant[i,"R2.adj"] ) )
		ret_matrix<-rbind(ret_matrix,temp_ret_matrix)			

	}
	resDF<-data.frame(Variable=ret_matrix[,1], Variance = ret_matrix[,2], R2=as.numeric(ret_matrix[,3]) ) 

	resDF$Variance <-factor(resDF$Variance, levels = c("R2", "Cum_R2","Unconstrain"))
	resDF$Variable <-factor(resDF$Variable, levels = rownames(all_var) )
	return(resDF)
}


PCoA_plot <- function(in.table = "",  in.Metadata="", Distance="", metaMDS=T, PDF_title="",plot_title="",colors="",pchPJ="",PC1_PDF_title="",treat="", remove_var_envfit="none",pch2PCoA="", legend_lab = c(), legend_col=c(), tretColors = c(), PCHlab=c(), PCH_title="" ){

	PCoA<-dbrda(in.table ~ 1, dist=Distance, metaMDS = metaMDS)
	colo2PCoA<-colors
	#pch2PCoA <- pchPJ[as.character(in.Metadata[rownames(PCoA$CA$u),]$PJ_cycles)]
	
	#### Envfit ####
	if( remove_var_envfit=="none"){
		enfit.pcoa<-envfit(PCoA,in.Metadata[rownames(PCoA$CA$u),])
	} else{
		enfit.Metadata<-in.Metadata[,colnames(in.Metadata)!=remove_var_envfit]
		enfit.pcoa<-envfit(PCoA,enfit.Metadata[rownames(PCoA$CA$u),])
	}

	xlab=paste("Comp1", format(100*PCoA$CA$eig[1]/sum(abs(PCoA$CA$eig)), digits=4), "%", sep=" ")
	ylab=paste("Comp2", format(100*PCoA$CA$eig[2]/sum(abs(PCoA$CA$eig)), digits=4), "%", sep=" ")
	#xlab<-paste("Comp1", format(summary(PCoA)$cont$importance[2,1] * 100 , digits=4), "%", sep=" ")
	#ylab<-paste("Comp2", format(summary(PCoA)$cont$importance[2,2] * 100 , digits=4), "%", sep=" ")
	### Define the limits
	xlim<-sort(c(enfit.pcoa$factors$centroids[,1],PCoA$CA$u[,1]))
	ylim<-sort(c(enfit.pcoa$factors$centroids[,2],PCoA$CA$u[,2]))	

#plot(PCoA$CA$u, cex = 4.5, xlab=xlab, ylab=ylab, main=plot_title, col = colors[as.character(in.Metadata[rownames(PCoA$CA$u),]$Day)] , pch = pch2PCoA, cex.main=2, cex.axis=1, cex.lab = 2 )	

	pdf(PDF_title,width=14,height=12)
		par(mar=c(5.7, 6.3, 4.1, 12.9), xpd=TRUE) # bottom, left, top, and right. 
		plot(PCoA$CA$u, cex = 2.5, xlab=xlab, ylab=ylab, main=plot_title, col = colors[rownames(PCoA$CA$u)] , 
			pch=pch2PCoA[rownames(PCoA$CA$u)], cex.main=2, cex.axis=1, cex.lab = 2, xlim=c( min(xlim)-0.1, max(xlim)+0.1 ),
			ylim=c( min(ylim)-0.1, max(ylim)+0.1 ))
		plot(enfit.pcoa , p.max = 0.05, col = "black", cex=1)

	treat<-treat[rownames(PCoA$CA$u)]

	for(i in levels(treat)) {
		print(i)
		color=tretColors[i]		
		ordiellipse(PCoA$CA$u[grep("TRUE",i==treat),1:2], groups=as.character(treat[treat==i]),col=color,draw="polygon",alpha=20, conf=0.70, ### Its confidence is set to the 80%
		show.groups=i) 
		ordispider(PCoA$CA$u[grep("TRUE",i==treat),], groups=treat[treat==i],col=color,lty=2,lwd=0.7) 
	} 



		#col2plot<-colo2PCoA[unique(names(colo2PCoA))]; pch2plot<-pch2PCoA[unique(sort(names(pch2PCoA)))]
		#legend(x=0.58,y=0.3, inset=c(-0.2,0), pch=19,legend=names(col2plot),col=col2plot, title="Day",bty = "n")
		legend("topright",  pch=19,inset=c(-0.23,0.3),legend=names(tretColors),col=tretColors, title="Condition",bty = "n", cex=2)
		legend("topright", inset=c(-0.23,0.55), pch=PCHlab,legend=names(PCHlab), title=PCH_title,bty = "n", cex=2)
		#legend(x=0.58,y=-0.03, inset=c(-0.2,0), pch=pch2plot,pt.cex = 1.5,legend=names(pch2plot), title="Feed",bty = "n")
		#if(length(pch2PCoA) > 0){
		#	legend("topright", inset=c(-0.12,0.4), pch=pch2plot,pt.cex = 1.5,legend=names(col2plot), title="Feed",bty = "n")
		#}
	dev.off()

	#pdf(PC1_PDF_title)
	#	plot(cbind(PCoA$CA$u[,1],rep(1,length(PCoA$CA$u[,1]))),col=coloresday[as.character(in.Metadata[rownames(PCoA$CA$u[,1:2]) , ]$Day)], pch = 20, cex=4, 
	#		xlab=paste("PC1", format(100*PCoA$CA$eig[1]/sum(abs(PCoA$CA$eig)), digits=4), "%", sep=" ")  )
			#text(c(pca$x[,1],rep(1,21)), as.character(Metadata[rownames(pca$x[,1]) , ]$Day) , cex=1)
			#legend(x=0.36,y=0.3, inset=c(-0.2,0), pch=19,legend=names(coloresday),col=coloresday, title="Day",bty = "n")
	#dev.off()

}






###############################################################################################################################################################
###############################################################################################################################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###############################################################################################################################################################
###############################################################################################################################################################

#################################################################################
########################	READ DATA		#########################
#################################################################################



#setworkingDir = "/home/luna.kuleuven.be/u0120466/Postdoc_Raes/Projects/Prodigest/Analysis/Dysbiosis_and_deseases/Arthritis/SPoA/Di_Paola_M/Biomarker_discovery"
#setwd(setworkingDir)

### Infiles
Metadata.File <- "/home/luna.kuleuven.be/u0120466/Scripts/R/biomarker_discovery/Prefilter/HIV/Metadata_HIV.tsv"
Genus.File <- "/home/luna.kuleuven.be/u0120466/Scripts/R/biomarker_discovery/Prefilter/HIV/otu_table_no_chimerasL6_raw.in_enteroscript"
### Variable to compare
#metadata_variable = "cond"  ## Description Histological_Inflammation Histology

#### Rarefaction and transformation
min.num.seqs = 2000
max.prev = 0.3 #### Max percentage of "0" counts (0-1)

#### Samples to exlude
#location2use <- "colon"
#Metadata<-read.table(Metadata.File, header=T, row.names=1, dec=".", sep="\t")
#samples2exlude<-rownames(subset(Metadata, Location != location2use))
samples2exlude<-c()

### Output directory
#dir.create(file.path("./", location2use), showWarnings = FALSE)
#setwd(paste0("./", location2use))

####################################################################################################################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~		Infiles		~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
####################################################################################################################################################

### Metadata file
Metadata<-read.table(Metadata.File, header=T, row.names=1, dec=".", sep="\t")
#Metadata<-subset(Metadata, Location == location2use)
#excSamples<-rownames(subset(Metadata,Location=="neg"))
#IncludeSamples<-setdiff(   rownames(Metadata), excSamples )
#Metadata.2<-Metadata[,setdiff( colnames(Metadata), c("id.old", "oldBC", "PlateName", "ID", "id2",  "Paper", "id1")   ) ]
#rownames(Metadata.2)<-rownames(Metadata)
Metadata.2<-Metadata

#Metadata.2$Description <- factor(Metadata.2$Description)
#Metadata.2$Histological_Inflammation <- factor(Metadata.2$Histological_Inflammation)
#Metadata.2$Histology <- factor(Metadata.2$Histology)
#Metadata.2$LocaHis <- factor(Metadata.2$LocaHis)


#infileGenus <- "simulated_viral_data.tsv"
#dada2taxfile="/home/luna.kuleuven.be/u0120466/Postdoc_Raes/Projects/Prodigest/Analysis/Dysbiosis_and_deseases/Arthritis/SpA_paper_biopsies/dada2_170_120/taxa_SV.tsv"
#dada2file="/home/luna.kuleuven.be/u0120466/Postdoc_Raes/Projects/Prodigest/Analysis/Dysbiosis_and_deseases/Arthritis/SpA_paper_biopsies/dada2_170_120/sequence_table_SV.tsv"
#phyloseq_file = "/home/luna.kuleuven.be/u0120466/Postdoc_Raes/Projects/Prodigest/Analysis/Dysbiosis_and_deseases/Arthritis/SpA_paper_biopsies/dada2_170_120/Genus_above.25000.physeq.RData"
#format="table" ## possible values "table", "DADA2", "phyloseq"
#min.num.seqs = 97105

## Infile table
ret.list<-read.infile.data(
	infileGenus = Genus.File,
	format="table", ## phyloseq, DADA2, table
	min.num.seqs = min.num.seqs,
	min.rar = TRUE ## TRUE (use the minimun number of samples for the rarefactions) FALSE (use the value in the variable min.num.seqs )
)
Genus<-ret.list[["Genus"]]
Genus.rar<-ret.list[["Genus.rar"]]

colSums(Genus)
colSums(Genus.rar)

### Exclude samples
if( length(samples2exlude) > 0 ){
	subsetSamples<-setdiff( colnames(Genus.rar), samples2exlude )
	Genus.rar<-Genus.rar[,subsetSamples]
	Genus<-Genus[,subsetSamples]
	Genus<-Genus[rowSums(Genus) != 0,]
}

samples2use<- colnames(Genus)[colnames(Genus) %in% rownames(Metadata)]
Genus = Genus[,samples2use] 
Genus.rar = Genus.rar[,samples2use] 

### Set the working directory
#setwd(file.path("./", metadata_variable))

####################################################################################################################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~		Data format	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
####################################################################################################################################################

#### Putative contaminants list
# https://bmcbiol.biomedcentral.com/articles/10.1186/s12915-014-0087-z
contamination<-c("g_Afipia", "g_Aquabacterium", "g_Asticcacaulis", "g_Aurantimonas", "g_Beijerinckia", "g_Bosea", "g_Bradyrhizobium", "g_Brevundimonas", "g_Caulobacter", "g_Craurococcus", "g_Devosia", "g_Hoeflea", "g_Mesorhizobium", "g_Methylobacterium", "g_Novosphingobium", "g_Ochrobactrum", "g_Paracoccus", "g_Pedomicrobium", "g_Phyllobacterium", "g_Rhizobium", "g_Roseomonas", "g_Sphingobium", "g_Sphingomonas", "g_Sphingopyxis", "g_Acidovorax", "g_Azoarcus", "g_Azospira", "g_Burkholderiaf", "g_Comamonas", "g_Cupriavidus", "g_Curvibacter", "g_Delftia", "g_Duganella", "g_Herbaspirillum", "g_Janthinobacterium", "g_Kingella", "g_Leptothrix", "g_Limnobacter", "g_Massilia", "g_Methylophilus", "g_Methyloversatilis", "g_Oxalobacter", "g_Pelomonas", "g_Polaromonas", "g_Ralstonia", "g_Schlegelella", "g_Sulfuritalea", "g_Undibacterium", "g_Variovorax", "g_Acinetobacter", "g_Enhydrobacter", "g_Nevskia", "g_Pseudomonas", "g_Pseudoxanthomonas", "g_Psychrobacter", "g_Stenotrophomonas", "g_Xanthomonas", "g_Aeromicrobium", "g_Arthrobacter", "g_Beutenbergia", "g_Brevibacterium", "g_Corynebacterium", "g_Curtobacterium", "g_Dietzia", "g_Geodermatophilus", "g_Janibacter", "g_Kocuria", "g_Microbacterium", "g_Micrococcus", "g_Microlunatus", "g_Patulibacter", "g_Propionibacteriume", "g_Rhodococcus", "g_Tsukamurella", "g_Abiotrophia", "g_Bacillus", "g_Brevibacillus", "g_Brochothrix", "g_Facklamia", "g_Paenibacillus", "g_Streptococcus", "g_Chryseobacterium", "g_Dyadobacter", "g_Flavobacterium", "g_Hydrotalea", "g_Niastella", "g_Olivibacter", "g_Pedobacter", "g_Wautersiella", "g_Enterobacter", "g_Escherichia","g_Kineococcus","g_Ornithinimicrobium","g_Phenylobacterium","g_Escherichia/Shigella","g_Escherichia.Shigella","g_Tumebacillus","g_Escherichia_Shigella")

### In case of having specific negative controls
CNsigcontamination<-read.table("/home/luna.kuleuven.be/u0120466/Postdoc_Raes/Projects/Prodigest/Analysis/Dysbiosis_and_deseases/Arthritis/SPoA/SpA_paper_biopsies/Contamination/Genus_cont_sig.tsv", header=T, row.names=1, dec=".", sep="\t")
more.cont <- rownames(subset(CNsigcontamination,Condition=="neg" & adj.pval< 0.1 ))
contamination<-unique(sort(c(contamination,more.cont)))

################################################################
########################### DATA TRANSFORMATION ################
#### The input is a data.frame: rows species and columns samples
##### CLR transformation #######
in.genus.clr<-CLR.transformation(
	in.table = Genus, 
	min.reads <- 0, # min.num.seqs #all will be kept
	min.prop = 0.001, #OTUs with abundance of at least 0.01 (default) #this shouldn't have too much of an effect, and this way we avoid imputing too many 0s
	cutoff = 0  #i'll filter low abundant taxa afterwards (don't want to remove them from the matrix before doing the clr transformation)

)

################################################################
########################### Filter data ########################

#### Filter prevalence
### Combine the genus from the CLR and the ones from the total samples
in.genus <- filter.prevalence(in.table = t(Genus), max.prev = max.prev) ### Samples are rows
Genus2use<-colnames(in.genus)[colnames(in.genus) %in% colnames(in.genus.clr)]
Genus2use <- setdiff(Genus2use ,contamination )

################################################################
########################### Final input data ###################
##### All filtering #######
#IncludeSamples<-intersect(rownames(in.genus.clr),IncludeSamples)
IncludeSamples<-rownames(in.genus.clr)

in.genus<-Genus[Genus2use, IncludeSamples]

##### Rar filtering #######
in.genus.rar<-Genus.rar[Genus2use, IncludeSamples]

##### CLR filtering #######
in.genus.clr<-t(in.genus.clr[IncludeSamples,Genus2use])
# matrix.f.n0["g_Clostridium_XI",]
# matrix.f.n0.clr["g_Clostridium_XI",]
#Genus.rar[ "g_Holdemania" ,]
# setdiff( IncludeSamples , colnames(Genus) )


##### Metadata filtering #######
#Metadata.2 <- Metadata.2[IncludeSamples,]
#Metadata_all_samp <-  Metadata.2[,c("Cons", "Description","LocaHis", "Histological_Inflammation", "Histology")]
#Metadata_cases <-  subset(Metadata.2 , is.na(ASDAS) == FALSE )
#Metadata_cases <-  Metadata_cases[,c("Cons","LocaHis", "Histological_Inflammation", "Histology", "NSAID","Age", "BMI", "Sy.Du", "CRP","BASDAI","ASDAS" )]

#cases.genus.rar<-in.genus.rar[,rownames(Metadata_cases)]
#cases.genus.clr<-in.genus.clr[,rownames(Metadata_cases)]


#################################################################################
########################	ADONIS		#################################
#################################################################################
out_dir  = "./"
#####  bray.table
bray.table.ADONIS <- ADONIS_func( in.matrix = t(in.genus.rar), Distance = "bray", in.Metadata = Metadata , prefix = paste0(out_dir,"/Common_metadata_bray") )

#cases_bray.table.ADONIS <- ADONIS_func( in.matrix = t(cases.genus.rar), Distance = "bray", in.Metadata = Metadata_cases, prefix = paste0(out_dir,"/Cases_metadata_bray") )


#####  clr.table
clr.table.ADONIS <- ADONIS_func( in.matrix = t(in.genus.clr), Distance = "euclidean", in.Metadata = Metadata, prefix = paste0(out_dir,"/Common_metadata_euclidean"))

#cases_clr.table.ADONIS <- ADONIS_func( in.matrix = t(cases.genus.clr), Distance = "euclidean", in.Metadata = Metadata_cases, prefix = paste0(out_dir,"/Cases_metadata_euclidean"))

#################################################################################
########################	ordiR2step		#########################
#################################################################################
#attach(Metadata.2)
#####  bray.table
capscale_bray<-capscale_cum_variance( in.matrix = t(in.genus.rar), Distance = "bray", in.Metadata = Metadata.2, prefix = "Common_bray")
capscale_bray[["non_redundant"]]
capscale_bray[["all"]]

#cases_capscale_bray<-capscale_cum_variance( in.matrix = t(cases.genus.rar), Distance = "bray", in.Metadata = Metadata_cases, prefix = "Cases_bray")
#cases_capscale_bray[["non_redundant"]]
#cases_capscale_bray[["all"]]


#####  clr.table
capscale_clr<-capscale_cum_variance(in.matrix = t(in.genus.clr),Distance = "euclidean",in.Metadata = Metadata.2,prefix = "Common_clr")
capscale_clr[["non_redundant"]]
capscale_clr[["all"]]

#cases_capscale_clr<-capscale_cum_variance(in.matrix = t(cases.genus.clr),Distance = "euclidean",in.Metadata = Metadata_cases,prefix = "Cases_clr")
#cases_capscale_clr[["non_redundant"]]
#cases_capscale_clr[["all"]]


#################################################################################
########################	Plot the nonredundant covariates	#########
#################################################################################

### Estimate the difference in covariance
#Bray_nrR2 <- nonredundant_R2_matrix(namesvars=colnames(Metadata.2), in.table= capscale_bray)
#clr_nrR2 <- nonredundant_R2_matrix(namesvars=colnames(Metadata.2), in.table= capscale_clr)

capscale_bray2plot<-nonredundant2plot(non.redundant = capscale_bray[["non_redundant"]], all.var = capscale_bray[["all"]] )

P<- ggplot(data=capscale_bray2plot, aes(x=Variable, y=R2, fill=Variance)) +
	geom_bar(stat="identity", position=position_dodge()) +
	 ggtitle("B-C explained variance") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("BC_ordiR2step.pdf",plot=P,width = 12)	


P<- ggplot(data=capscale_bray2plot[capscale_bray2plot$Variable %in%rownames(capscale_bray[["non_redundant"]]),], aes(x=Variable, y=R2, fill=Variance)) +
	geom_bar(stat="identity", position=position_dodge()) +
	 ggtitle("B-C") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("BC_ordiR2step_sig.pdf",plot=P,width = 7)	


#cases_capscale_bray2plot<-nonredundant2plot(non.redundant = cases_capscale_bray[["non_redundant"]], all.var = cases_capscale_bray[["all"]] )

#P<- ggplot(data=cases_capscale_bray2plot, aes(x=Variable, y=R2, fill=Variance)) +
#	geom_bar(stat="identity", position=position_dodge()) +
#	 ggtitle("Cases B-C SpoA")
#ggsave("Cases_BC_ordiR2step.pdf",plot=P)	


capscale_clr2plot<-nonredundant2plot(non.redundant = capscale_clr[["non_redundant"]], all.var = capscale_clr[["all"]] )

P<- ggplot(data=capscale_clr2plot, aes(x=Variable, y=R2, fill=Variance)) +
	geom_bar(stat="identity", position=position_dodge()) +
	 ggtitle("CLR") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("CLR_ordiR2step.pdf",plot=P,width = 12)	

P<- ggplot(data=capscale_clr2plot[capscale_clr2plot$Variable %in%rownames(capscale_clr[["non_redundant"]]),], aes(x=Variable, y=R2, fill=Variance)) +
	geom_bar(stat="identity", position=position_dodge()) +
	 ggtitle("CLR explained variance") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("CLR_ordiR2step_sig.pdf",plot=P,width = 7)	

#cases_capscale_clr2plot<-nonredundant2plot(non.redundant = cases_capscale_clr[["non_redundant"]], all.var = cases_capscale_clr[["all"]] )

#P<- ggplot(data=cases_capscale_clr2plot, aes(x=Variable, y=R2, fill=Variance)) +
#	geom_bar(stat="identity", position=position_dodge()) +
#	 ggtitle("Cases CLR SpoA")
#ggsave("Cases_CLR_ordiR2step.pdf",plot=P)	


#################################################################################
########################  	PCoAs 	        #################################
#################################################################################

PChs<-rep(0,nrow(Metadata.2))
names(PChs) <- rownames(Metadata.2)
#PChs[Metadata.2$Group == "A" ] = 15
PChs[Metadata.2$transmision == "hsh" ] = 17
#PChs[Metadata.2$Group == "C" ] = 18
PChs[Metadata.2$transmision == "htx" ] = 19	

PCHlab = c(17,19)
names(PCHlab) <- c("hsh","htx")
PCH_title = "Transmision"
colors <- rep("",nrow(Metadata.2))
names(colors) <- rownames(Metadata.2)
colors[Metadata.2$Group == "A" ] = "red"
colors[Metadata.2$Group == "B" ] = "darkgreen"
colors[Metadata.2$Group == "C" ] = "orange"
colors[Metadata.2$Group == "D" ] = "blue"

tretColors<-c("red","darkgreen","orange","blue")
names(tretColors) <- c("A","B","C","D")

treat<-Metadata.2$Group
names(treat) <- rownames(Metadata.2)


PCoA_plot(in.table = t(in.genus.clr),
	in.Metadata=Metadata.2,
	Distance="euclidean",
	plot_title = "CLR PCoA",
	PDF_title  = "CLR_PCoA.pdf",
	colors=colors,
	treat=treat,
	tretColors = tretColors,
	pch2PCoA = PChs,
	PCHlab=PCHlab,
	PCH_title = "Transmision"	
)


PCoA_plot(in.table = t(in.genus.rar),
	in.Metadata=Metadata.2,
	Distance="bray",
	plot_title = "BC PCoA",
	PDF_title  = "BC_PCoA.pdf",
	colors=colors,
	treat=treat,
	tretColors = tretColors,
	pch2PCoA = PChs,
	PCHlab=PCHlab,
	PCH_title = "Transmision"	
	
	
)
















