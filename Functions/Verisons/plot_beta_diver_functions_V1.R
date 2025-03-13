set.seed(12345)
library(phyloseq)
library("ggplot2")
library(ggpubr)
source("/home/luna.kuleuven.be/u0141268/Postdoc_Raes/Projects/DISCOVERIE/Functions/beta_diver_functions_V9.R")
library(gridExtra)
library(RColorBrewer)
library(microbiome)
library("viridis") 

PCoA_grapper <- function(Var = "", tempMetadata = data.frame(), table.ADONIS=data.frame(),
		taxa_matrix=matrix(), distance="bray", Ellipses = F, Label = F,colors = c() ){

	tempMetadata$VarTemp <- tempMetadata[,Var]
	tempMetadata <- tempMetadata[complete.cases(tempMetadata$VarTemp),]
	taxa_matrix <- taxa_matrix[match( rownames(tempMetadata), rownames(taxa_matrix)),]
	#print(table(is.na(taxa_matrix)))

	N_samples <- nrow(tempMetadata)
	pval <-  round(  as.numeric(as.character(table.ADONIS[Var,"p.value"])) ,digits=3)
	FDR <- round(  as.numeric(as.character(table.ADONIS[Var,"BH.adj.p.value"])) ,digits=3)
	R2 <- round(  as.numeric(as.character(table.ADONIS[Var,"R2"])) ,digits=3)
	PCoA_title <- Var
	PCoA_subtitle <- paste0( "p-val = ", pval,"; FDR = ", FDR,"; R2 = ", R2,"; N = ", N_samples)
		
	moreColors <- brewer.pal(n = 8, name = "Dark2")

	if(length(colors) == 0){

#	if( length(sort(unique(tempMetadata$VarTemp))) == 2){colors<- c("red","blue");names(colors) <- sort(unique(tempMetadata$VarTemp))}
#	if( length(sort(unique(tempMetadata$VarTemp))) == 3){colors<- c("red","blue","darkgreen");names(colors) <- sort(unique(tempMetadata$VarTemp))}
#	if( length(sort(unique(tempMetadata$VarTemp))) == 4){colors<- moreColors[1:4];names(colors) <- sort(unique(tempMetadata$VarTemp))}
		if( length(sort(unique(tempMetadata$VarTemp))) == 2){colors<- moreColors[1:2];names(colors) <- sort(unique(tempMetadata$VarTemp))}
		if( length(sort(unique(tempMetadata$VarTemp))) == 3){colors<- moreColors[1:3];names(colors) <- sort(unique(tempMetadata$VarTemp))}
		if( length(sort(unique(tempMetadata$VarTemp))) == 4){colors<- moreColors[1:4];names(colors) <- sort(unique(tempMetadata$VarTemp))}
	}
		
	list_ordination <- vegan_PCoA_envfit(in.matrix= taxa_matrix, Metadata2enfit=data.frame(tempMetadata), distance = distance)
	PCoA<-list_ordination[["PCoA"]]; xlab<-list_ordination[["xlab"]]; ylab<-list_ordination[["ylab"]]; fit<-list_ordination[["fit"]]
	df_ord<-list_ordination[["df_ord"]]

	ret_list <- ggplot_envfit(df_ord, ord=PCoA,fit,Metadata=data.frame(tempMetadata),alpha_pval=0.99 )
		df_ord <- ret_list[["df_ord"]]
		df_arrows <- ret_list[["df_arrows"]]
		df_factors <- ret_list[["df_factors"]]

	if(Var == "Enterotype"){
	
		colors_ent = c("#f28118", "#ba1015", "#059571","#4b59a1")
		names(colors_ent) <- c( "Bacteroides 1" , "Bacteroides 2", "Prevotella", "Ruminococcus")
		PCoAplot <- ggplot(data = df_ord, aes(x = x, y = y, color = VarTemp )) + theme_bw() + 
			geom_point( aes(color = VarTemp), size =3) + 
			xlab(xlab) + ylab(ylab)  + scale_color_manual(values=colors_ent) + 
			stat_ellipse(aes(x = x, y = y,  group=VarTemp, fill=VarTemp ), 
			linetype = 2 ,type = "norm" , geom="polygon",level=0.8,alpha=0.1, show.legend=F) +	
			labs(title = PCoA_title, subtitle = PCoA_subtitle,color = Var) 				
			
	}else if( (class(tempMetadata$VarTemp) == "character" | class(tempMetadata$VarTemp) == "factor") & length(sort(unique(tempMetadata$VarTemp))) < 4 ){
		PCoAplot <- ggplot(data = df_ord, aes(x = x, y = y, color = VarTemp )) + theme_bw() + 
			geom_point( aes(color = VarTemp), size =3) + 
			xlab(xlab) + ylab(ylab)  + 
			#scale_color_manual(values=colors) + 
			labs(title = PCoA_title, subtitle = PCoA_subtitle,color = Var) 
			
			if(Ellipses == T){
				PCoAplot <- PCoAplot + 
					stat_ellipse(aes(x = x, y = y,  group=VarTemp, fill=VarTemp ), 
					linetype = 2 ,type = "norm" , geom="polygon",level=0.8,alpha=0.1, show.legend=F) 
			}						
			
			if(Label == T){
				subdf_factors <- df_factors[grep( Var,  df_factors$var ),]
				subdf_factors$var <- gsub(Var,"",subdf_factors$var)
							
				PCoAplot <- PCoAplot + 
					geom_text(data = subdf_factors,  aes(x = x, y = y, label = var), 
					#inherit.aes = FALSE, color =colors,  hjust = "outward", size =4 )
					inherit.aes = FALSE,  hjust = "outward", size =5 )					
			}									
			if(length(colors) != 0){
				
				PCoAplot <- PCoAplot + scale_color_manual(values=colors[levels(factor(tempMetadata$VarTemp))]) 
			}			
			
	}else if( (class(tempMetadata$VarTemp) == "character" | class(tempMetadata$VarTemp) == "factor") & length(sort(unique(tempMetadata$VarTemp))) >= 4 ){
		PCoAplot <- ggplot(data = df_ord, aes(x = x, y = y, color = VarTemp )) + theme_bw() + 
			geom_point( aes(color = VarTemp), size =3) + 
			xlab(xlab) + ylab(ylab)   + 
			labs(title = PCoA_title, subtitle = PCoA_subtitle,color = Var) 
			
			if(Ellipses == T){
				PCoAplot <- PCoAplot + 
					stat_ellipse(aes(x = x, y = y,  group=VarTemp, fill=VarTemp ), 
					linetype = 2 ,type = "norm" , geom="polygon",level=0.8,alpha=0.1, show.legend=F) 
			}						
			
			if(Label == T){
			
				subdf_factors <- df_factors[grep( Var,  df_factors$var ),]
				subdf_factors$var <- gsub(Var,"",subdf_factors$var)
				PCoAplot <- PCoAplot + 
					geom_text(data = subdf_factors,  aes(x = x, y = y, label = var), 
					inherit.aes = FALSE,  hjust = "outward", size =4 )
			}
			
			if(length(colors) != 0){
				
				PCoAplot <- PCoAplot + scale_color_manual(values=colors[levels(factor(tempMetadata$VarTemp))]) 
			}			
						
			
	
	}else{
		yourname_color_palette <- c("#74869c", "#6daddd", "#83adbb", "#b3d17c", "#ddaa7b","#ab5548")
		colors2use <- colorRampPalette(yourname_color_palette)(length(sort(unique(tempMetadata$VarTemp))))
		PCoAplot <- ggplot(data = df_ord, aes(x = x, y = y, color = VarTemp )) + theme_bw() + 
			geom_point( aes(color = VarTemp), size =3) + 
			xlab(xlab) + ylab(ylab)  +	
			labs(title = PCoA_title, subtitle = PCoA_subtitle,color = Var)  +   scale_colour_gradientn(colors=colors2use)	
			
		angle = 20;len = 0.5;unit = "cm"
		PCoAplot <- PCoAplot + geom_segment(data = df_arrows[grep( Var,  df_arrows$var ),], aes(x = 0,  xend = x, y = 0, yend = y), 
			arrow = arrow(angle = angle, length = unit(len, unit)), color = "black", inherit.aes = FALSE) + 
			geom_text(data = df_arrows[grep( Var,  df_arrows$var ),],  
			aes(x = x, y = y, label = var), inherit.aes = FALSE, color = "black",  hjust = "outward", size =4)

			
	}
	return(PCoAplot)	
}	


