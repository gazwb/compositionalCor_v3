# function calls for comparison of networks
setwd("~/MAARS_p2/MB")

#########################################################################################
# set of functions for fixing the module names and realigning them to the control modules
#########################################################################################



# create overlap matrix i.e. the otus which overlap over the modules over two networks
defOverlapMatrix <- function(net1, net2){

	ovlpMatrix <- matrix(ncol = max(unique(net1$Module)), nrow = max(unique(net2$Module)), data = 0)
	for (i in 1:nrow(ovlpMatrix)){
		for (j in 1:ncol(ovlpMatrix)){
			ovlpMatrix[i,j] <- sum(net2[net2$Module==i,1] %in% net1[net1$Module==j,1]) 
		}
		
	}
	return(ovlpMatrix)
}

# define the jaccard matrix 

defJaccardMatrix <- function(net1, net2,ovlpMatrix){

	jacMatrix <- matrix(ncol = max(unique(net1$Module)), nrow = max(unique(net2$Module)), data = 0)
	for (i in 1:nrow(ovlpMatrix)){
		for (j in 1:ncol(ovlpMatrix)){
			jacMatrix[i,j] <- sum(net2[net2$Module==i,1] %in% net1[net1$Module==j,1]) / length(union(net2[net2$Module==i,1],net1[net1$Module==j,1])) 
		}
	}
	return(jacMatrix)
}

# realign the modules
reAlign <- function(fx, jacMatrix){
	passIndx <- as.data.frame(which(jacMatrix >= 0.2,arr.ind = TRUE)) #row is thesecondary network, col is the primary
	passIndx[duplicated(passIndx$col),2] <- paste0(passIndx[duplicated(passIndx$col),2],"99")

	for (y in 1:dim(passIndx)[1]){
		fx[fx$Module==passIndx[y,1],10] <- passIndx[y,2]
	}
	colnames(fx)[10] <- "module.fixed"

	changeTable <- data.frame(matrix(ncol = length(unique(fx[is.na(fx$module.fixed),3])), nrow =length(unique(fx[is.na(fx$module.fixed),3])),data = 0 ))
	changeTable[,1] <- unique(fx[is.na(fx$module.fixed),3])
	changeTable[,2] <- order(unique(fx[is.na(fx$module.fixed),3]),decreasing =  FALSE)
	fx[is.na(fx$module.fixed),10] <- changeTable[match(fx[is.na(fx$module.fixed),3], changeTable[,1]),2] + dim(jacMatrix)[2]
return(fx)
}

reAlign.ctrl <- function(fx){
	fx$module.fixed <- fx$Module
	return(fx)
}

write_realignment <- function(x,y,z){
	write.table(x,file = "/home/gaz/MAARS_p2/Scripts/compositionalCor_v3/SparCC/mainNetworks_calc/AD_LES_nodes_realign.txt", quote = FALSE, sep ="\t", row.names = FALSE)
	write.table(y,file = "/home/gaz/MAARS_p2/Scripts/compositionalCor_v3/SparCC/mainNetworks_calc/CTRL_nodes_realign.txt", quote = FALSE, sep ="\t", row.names = FALSE)
	write.table(z,file = "/home/gaz/MAARS_p2/Scripts/compositionalCor_v3/SparCC/mainNetworks_calc/PSO_LES_nodes_realign.txt", quote = FALSE, sep ="\t", row.names = FALSE)

}


# basic network statistics 
# calcStats is the main call and calls several functions to get networks statistics
calcStats <- function(nodeDataL,edgeDataL,netOrder) {
	stat.df <- data.frame(matrix(nrow = length(netOrder), ncol = 5,dat = 0))
	rownames(stat.df) <- netOrder
	colnames(stat.df) <- c("OTU_in_Network", "Edges_in_Network", "Clustering_Coefficient", "No_Modules", "Mean_sparCC")

	######populate
	stat.df$OTU_in_Network <-  unlist(lapply(nodeDataL, function(x) dim(x)[1]))
	stat.df$Edges_in_Network <-  unlist(lapply(edgeDataL, function(x) dim(x)[1]))
	stat.df$Clustering_Coefficient <- calcClusteringCoef(edgeDataL)[,1]
	stat.df$No_Modules <- unlist(lapply(nodeDataL, function(x) length(unique(x$Module))))
	stat.df$Mean_sparCC <- unlist(lapply(edgeDataL, function(x) mean(abs(x$cor))))
	stat.df$Modularity <- calcModularity(edgeDataL)[,1]

	o2Proportions <- calcO2Proportion(nodeDataL)
	stat.dfo2 <- cbind(stat.df,o2Proportions)
return(stat.dfo2)
}


createIgraphObject <- function(edgeData) {
	#ig <-createIgraphObject(edgeDataL[[3]])
	ig <-graph.data.frame(edgeData, directed=F)
	return(ig)
}


calcClusteringCoef <- function(edgeDataL) {
	edgeData_list <- edgeDataL
	coefTabs <- data.frame(matrix(ncol= 1, nrow= length(edgeData_list), dat = 0))
	colnames(coefTabs) <-  c("ClusteringCoef")
	rownames(coefTabs) <- netOrder

	for (j in 1:length(edgeData_list)) { 
			ig <-createIgraphObject(edgeData_list[[j]])
			coefTabs[j,1] <- round(transitivity(ig), digits = 3) 

	} 
	return(coefTabs)
}


calcO2Proportion <- function(nodeDataL){
	nodeData_list <- nodeDataL
	countTabs <- data.frame(matrix(ncol= 3, nrow= length(nodeData_list), dat = 0))
	colnames(countTabs) <-  c("Aer","AnAer", "FacAnAer")
	rownames(countTabs) <- netOrder
	propTabs <- list()
		for (i in 1:length(nodeData_list)){
				propTabs[[i]] <- table(nodeData_list[[i]]$o2_tolerance)
				countTabs[i,] <- as.numeric(propTabs[[i]][names(propTabs[[i]]) %in% colnames(countTabs)])
		}
	countTabs$Aerobic_Proportion <- round(countTabs$Aer / (countTabs$Aer + countTabs$AnAer),digits = 3) 

return(countTabs)
}

calcModularity <- function(edgeDataL) {
	edgeData_list <- edgeDataL
	coefTabs <- data.frame(matrix(ncol= 1, nrow= length(edgeData_list), dat = 0))
	colnames(coefTabs) <-  c("Modularity")
	rownames(coefTabs) <- netOrder

	for (j in 1:length(edgeData_list)) { 
			ig <-createIgraphObject(edgeData_list[[j]])
			z <- multilevel.community(ig, weights = NULL)
			coefTabs[j,1] <- round(modularity(z), digits = 3) 

	} 
	return(coefTabs)
}

#need this
getUnionNodes <- function(){

	node <- readAttributes()
	n1 <- node[[1]]
	n2 <- node[[2]]
	n3 <- node[[3]]

	unionNodes <- union(union(n1$nodeName,n2$nodeName),n3$nodeName)
	return(unionNodes)
}
# pairwise differential connectivity
# a funtion which calculates differential connectivity between two networks
# pass this two networks from a function that chooses which ones to pass

diffConnectivity <- function(edgeDataL){
	
net.ADL <- edgeDataL[[1]]
net.CTRL <- edgeDataL[[2]]
net.PSOL <- edgeDataL[[3]]

	ig.net.ADL <- createIgraphObject(net.ADL)
	ig.net.CTRL <- createIgraphObject(net.CTRL)
	ig.net.PSOL <- createIgraphObject(net.PSOL)

	ADL.degree <- degree(ig.net.ADL)
	CTRL.degree <- degree(ig.net.CTRL)
	PSOL.degree <- degree(ig.net.PSOL)

	unionNodes <- getUnionNodes()
	degreeTable <- data.frame(matrix(nrow = length(unionNodes), ncol = 9, data = 0))
	colnames(degreeTable) <- c("ADL.degree", "CTRL.degree", "PSOL.degree", "Diffc-ad","Diffc-pso", "Diffad-pso", "p1","p2","p3")
	row.names(degreeTable) <- unionNodes
# populate
	degreeTable[,1] <- ADL.degree[match(unionNodes, names(ADL.degree))]
	degreeTable[,2] <- CTRL.degree[match(unionNodes, names(CTRL.degree))]
	degreeTable[,3]<- PSOL.degree[match(unionNodes, names(PSOL.degree))]
	degreeTable[,4] <- (degreeTable[,1] - degreeTable[,2])
	degreeTable[,5] <- (degreeTable[,3] - degreeTable[,2])
	degreeTable[,6] <- (degreeTable[,3] - degreeTable[,1])
	degreeTable[is.na(degreeTable)] <- 0


return(degreeTable)
}


#######################################################################
# set of functions for calculating p value of differential connectivity
#######################################################################

## calc the pairwise differential connectivitiy between two networks ## 
pairwisediffConnectivity <- function(network1, network2, unionNodes) {
n1 <- createIgraphObject(network1)
n2 <- createIgraphObject(network2)

	n1.degree <- degree(n1)
	n2.degree <- degree(n2)

	degreeTable <- data.frame(matrix(nrow = length(unionNodes), ncol = 3, data = 0))
	colnames(degreeTable) <- c("n1.degree", "n2.degree", "diff")
	degreeTable[,1] <- n1.degree[match(unionNodes, names(n1.degree))]
	degreeTable[,2] <- n2.degree[match(unionNodes, names(n2.degree))]
	#degreeTable[is.na(degreeTable)] <- 0
	degreeTable[,3] <- (degreeTable[,2] - degreeTable[,1])
	return(degreeTable[,3])

}

# read in permuted networks 
#############################################################################################
# I will keep this as i expect i will need to do this again for the differential connectivity
#############################################################################################
readBootstraps <- function() {

	bootstrapStore <- list()
	bs.ADL.dir <- paste0("/home/gaz/MAARS_p2/Scripts/compositionalCor_v3/SparCC/AD_LES/Bootstraps/", list.files("/home/gaz/MAARS_p2/Scripts/compositionalCor_v3/SparCC/AD_LES/Bootstraps/"))
	bs.ADL <- lapply(bs.ADL.dir, function(x) read.table(x, header = TRUE, check.names = FALSE, row.names =1))
	cat("ADL bootstraps read in \n")

	bs.CTRL.dir <- paste0("/home/gaz/MAARS_p2/Scripts/compositionalCor_v3/SparCC/CTRL/Bootstraps/", list.files("/home/gaz/MAARS_p2/Scripts/compositionalCor_v3/SparCC/CTRL/Bootstraps/"))
	bs.CTRL <- lapply(bs.CTRL.dir, function(x) read.table(x, header = TRUE, check.names = FALSE,row.names =1)) 
	cat("CTRL bootstraps read in \n")

	bs.PSOL.dir <- paste0("/home/gaz/MAARS_p2/Scripts/compositionalCor_v3/SparCC/PSO_LES/Bootstraps/", list.files("/home/gaz/MAARS_p2/Scripts/compositionalCor_v3/SparCC/PSO_LES/Bootstraps/"))
	bs.PSOL <- lapply(bs.PSOL.dir, function(x) read.table(x, header = TRUE, check.names = FALSE,row.names =1)) 
	cat("PSOL bootstraps read in \n")

	bootstrapStore[[1]] <- bs.ADL
	bootstrapStore[[2]] <- bs.CTRL
	bootstrapStore[[3]] <- bs.PSOL
return(bootstrapStore)
}

# index of union set in bootstrap
# pass this ONLY ONE bootstrap to get the index
getUnionIndex <- function(bootstrap) { 
	unionNodes <- getUnionNodes()
	rcIndex <- which(colnames(bootstrap) %in% unionNodes, arr.ind = TRUE)

return(rcIndex)
}


calcModuleStats <- function(fx){
modStats <- list()
mNodes <- list()
mdNames <- as.numeric(unique(fx$module.fixed))[order(as.numeric(unique(fx$module.fixed)))]
mdNames2 <- mdNames
mdNames2[mdNames > 100] <- max(mdNames[mdNames < 100]) + 1
md2 <- mdNames
	for (i in 1:length(mdNames)){
		mNodes[[i]] <- fx[fx$module.fixed== mdNames[i] ,]
	}

modSizes <- data.frame(module.fixed = mdNames, noGenes = 0, stringsAsFactors = FALSE )
	for (r in 1:dim(modSizes)[1]){
		modSizes[r,2] <- dim(mNodes[[r]])[1]
	}

	modStats[[1]] <- mdNames ## module names
	modStats[[2]] <- mNodes ## module bugs
	modStats[[3]] <- modSizes ## module sizes
return(modStats)
}

#exec one sparcc
