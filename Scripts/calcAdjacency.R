### this script will load the correlation matrixes for networks then identify a correlation threshold based upon a structural metric
###this script will also output the network node and edge attributes
source("/home/gaz/MAARS_p2/Scripts/compositionalCor_v3/Scripts/inference_functions.R")



cohorts <- c("AD_LES","CTRL","PSO_LES")
cnt <- 1
clusterings <-list()
for (c in cohorts){

	#sparCC.cor <- read.table(paste0("/home/gaz/MAARS_p2/Scripts/compositionalCor_v2/SparCC/AD_LES/sparCC_out_CTRL.txt"),header=1,sep="\t",row.names = 1,check.names= FALSE)
	sparCC.cor <- read.table(paste0("/home/gaz/MAARS_p2/Scripts/compositionalCor_v3/SparCC/",c,"/sparCC_out_", c, "_union_i20.txt"),header=1,sep="\t",row.names = 1,check.names= FALSE)
	sparCC.sig <- read.table(paste0("/home/gaz/MAARS_p2/Scripts/compositionalCor_v3/SparCC/",c, "/", c, "_pvals_two_sided.txt"),header=1,sep="\t",row.names = 1,check.names= FALSE)

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
	sparCC.cor <- read.table(paste0("/home/gaz/MAARS_p2/Scripts/compositionalCor_v3/SparCC/",c,"/sparCC_out_", c, "_union_i20.txt"),header=1,sep="\t",row.names = 1,check.names= FALSE)
	sparCC.sig <- read.table(paste0("/home/gaz/MAARS_p2/Scripts/compositionalCor_v3/SparCC/",c, "/", c, "_pvals_two_sided.txt"),header=1,sep="\t",row.names = 1,check.names= FALSE)

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

#write.table(component1.n, file = paste0("/home/gaz/MAARS_p2/Scripts/compositionalCor_v2/SparCC/mainNetworks_calc/",c, "_links.txt"), sep = "\t", quote = FALSE,row.names = FALSE)


zee <- multilevel.community(lccc, weights = NULL)
round(modularity(zee),digits=3)
