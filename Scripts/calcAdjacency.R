### this script will load the correlation matrixes for networks then identify a correlation threshold based upon a structural metric
###this script will also output the network node and edge attributes

library(igraph)
threshold_cor_sparCC_cor <- function(matrix,cor_threshold, connex){

	matrix[lower.tri(matrix, diag = TRUE)] <- 0
	rowz <- which(abs(matrix) >= cor_threshold , arr.ind = TRUE)
	rowz.nr <- as.data.frame(rowz[rowz[,1] != rowz[,2],])

	linksDat <- data.frame(fromNode = rownames(rowz.nr), toNode = NA,	 cor = NA, stringsAsFactors = FALSE)

	linksDat[,2] <- rownames(matrix[rowz.nr[1:dim(rowz.nr)[1],2],])
	linksDat[,3] <- diag(as.matrix(matrix[rowz.nr[1:dim(rowz.nr)[1],1],rowz.nr[1:dim(rowz.nr)[1],2]]))
	toColz <- strsplit(linksDat$toNode, "\\.(?=\\d)", perl = TRUE)
	linksDat[,2] <- sapply(toColz, "[[",1)

return(linksDat)
} #http://www.biomedcentral.com/1756-0500/2/240


calcAll <- function(corFile){
	gc()
	#betRange <- seq(from = 0.19, to = 0.22, by= 0.01)
	betRange <- seq(from = 0.06, to = 0.5, by= 0.01)
	storeCoef <- data.frame(matrix(ncol = 4, nrow = length(betRange), data = 0 ))
	storeCoef[,1] <- betRange
	colnames(storeCoef) <- c("Threshold","Metric", "noEdges", "noNodes")

		for (i in 1:length(betRange)){
			cat(paste0("Running : ", c, " cor : ",betRange[i], "\n"))
			linksDat.nr <- threshold_cor_sparCC_cor(corFile, betRange[i],"SparCC" )
			ig <-graph.data.frame(linksDat.nr, directed=F)
			lcc = induced.subgraph(ig, V(ig)[clusters(ig)$membership == 1])


		component  <- get.data.frame(lcc,what=c("edges"))

			#storeCoef[i,2] <- round(transitivity(lcc), digits = 3) #maximal.cliques.count(ig)
			z <- multilevel.community(lcc, weights = NULL)
			storeCoef[i,2] <- round(modularity(z),digits=3)
			storeCoef[i,3] <- dim(component)[1]
			storeCoef[i,4] <- length(union(unique(component$to),unique(component$from)))
			cat(paste0( "edges = ", storeCoef[i,3], " nodes =",storeCoef[i,4], "\n"))

		}
		cat(paste0("completed: ", c))
	return(storeCoef)
	}




#c <- "CTRL"
#corFile <- read.table(paste0("/home/gaz/MAARS_p2/Scripts/compositionalCor_v2/SparCC/","CTRL","/sparCC_out_", "CTRL", "_union_i20.txt"),header=1,sep="\t",row.names = 1,check.names= FALSE)

####using the pvals



cohorts <- c("AD_LES","CTRL","PSO_LES")
cnt <- 1
clusterings <-list()
for (c in cohorts){

	#sparCC.cor <- read.table(paste0("/home/gaz/MAARS_p2/Scripts/compositionalCor_v2/SparCC/AD_LES/sparCC_out_CTRL.txt"),header=1,sep="\t",row.names = 1,check.names= FALSE)
	sparCC.cor <- read.table(paste0("/home/gaz/MAARS_p2/Scripts/compositionalCor_v2/SparCC/",c,"/sparCC_out_", c, "_union_i20.txt"),header=1,sep="\t",row.names = 1,check.names= FALSE)
	sparCC.sig <- read.table(paste0("/home/gaz/MAARS_p2/Scripts/compositionalCor_v2/SparCC/",c, "/", c, "_pvals_two_sided.txt"),header=1,sep="\t",row.names = 1,check.names= FALSE)

	sparCC.sig2 <- apply(sparCC.sig, 2,as.numeric)
	sparCC.cor2 <- sparCC.cor

	#correct for multiple testing
	#AD_sparCC.sig.adj <- p.adjust(AD_sparCC.sig,method = "BH")



	sparCC.cor2[sparCC.sig2 >= 0.05] <- 0

	sparCC.cor <- sparCC.cor2

	#sparCC.cor2 <- sparCC.cor
	#sparCC.cor[sparCC.sig <0.1] <- 0
	# linksDat.nr <- threshold_cor_sparCC_cor(sparCC.cor, 0.21,"SparCC" )


	# ig <-graph.data.frame(linksDat.nr, directed=F)
	# maximal.cliques.count(ig)s

	clusterings[[cnt]] <- calcAll(sparCC.cor)
	plot(clusterings[[cnt]]$Threshold,clusterings[[cnt]]$Metric, main = cohorts[[cnt]],xlab = "Threshold", ylab = "Metric")
	abline(v = clusterings[[cnt]][clusterings[[cnt]]$Metric ==max(clusterings[[cnt]]$Metric),1])
	cnt <- cnt + 1
}

res.clust <- data.frame(cbind(clusterings[[1]][c(2,3,4)],clusterings[[2]][c(2,3,4)],clusterings[[3]][c(2,3,4)]))

row.names(res.clust) <- clusterings[[1]]$Threshold



##write out best
c <- "CTRL"
	sparCC.cor <- read.table(paste0("/home/gaz/MAARS_p2/Scripts/compositionalCor_v2/SparCC/",c,"/sparCC_out_", c, "_union_i20.txt"),header=1,sep="\t",row.names = 1,check.names= FALSE)
	sparCC.sig <- read.table(paste0("/home/gaz/MAARS_p2/Scripts/compositionalCor_v2/SparCC/",c, "/", c, "_pvals_two_sided.txt"),header=1,sep="\t",row.names = 1,check.names= FALSE)

	sparCC.sig2 <- apply(sparCC.sig, 2,as.numeric)
	sparCC.cor2 <- sparCC.cor

	#correct for multiple testing
	sparCC.cor2[sparCC.sig2 >= 0.05] <- 0 


linksdat.out <- threshold_cor_sparCC_cor(sparCC.cor2, 0.22,"SparCC" )
	igg <-graph.data.frame(linksdat.out, directed=F)
			lccc = induced.subgraph(igg, V(igg)[clusters(igg)$membership == 1])


component1.n  <- get.data.frame(lccc,what=c("edges"))
length(union(unique(component1.n$to),unique(component1.n$from)))



cat(paste0( "edges = ", dim(component1.n)[1], " nodes =",length(union(unique(component1.n$to),unique(component1.n$from))), "\n"))

write.table(component1.n, file = paste0("/home/gaz/MAARS_p2/Scripts/compositionalCor_v2/SparCC/mainNetworks_calc/",c, "_links.txt"), sep = "\t", quote = FALSE,row.names = FALSE)


zee <- multilevel.community(lccc, weights = NULL)
round(modularity(zee),digits=3)