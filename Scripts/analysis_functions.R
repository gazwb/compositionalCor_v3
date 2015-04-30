# function calls for comparison of networks
setwd("~/MAARS_p2/MB")
library(reshape)
library(ggplot2)
library("RColorBrewer")
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

defOverlapMatrix.fx <- function(net1, net2){

	ovlpMatrix <- matrix(ncol = max(unique(net1$module.fixed)), nrow = max(unique(net2$module.fixed)), data = 0)
	for (i in 1:nrow(ovlpMatrix)){
		for (j in 1:ncol(ovlpMatrix)){
			ovlpMatrix[i,j] <- sum(net2[net2$module.fixed==i,1] %in% net1[net1$module.fixed==j,1]) 
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
	passIndx <- which(jacMatrix >= 0.1) #row is thesecondary network, col is the primary
	


	passIndx[duplicated(passIndx$col),2] <- max(as.numeric(passIndx$col))

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

# realign the modules
reAlign.2 <- function(fx, jacMatrix){
	arr.max <- NULL
	for (i in 1:dim(jacMatrix)[2]){
		arr.max[i] <- which(jacMatrix[,i] == apply(jacMatrix,2,function(x) max(x))[i],arr.ind = TRUE)
}

ctrlMods <- 1:dim(jacMatrix)[2]
identityMatrix <- jacMatrix

if (any(duplicated(arr.max) == TRUE)) {  

	conflictCtrlModules <- c(which(arr.max == arr.max[duplicated(arr.max)]))
	conflictDisModule <- arr.max[duplicated(arr.max)]
	trueIdentCtrl <- ctrlMods[-conflictCtrlModules]
	trueIdentDis <- arr.max[-conflictCtrlModules]
	} else { 
	conflictCtrlModules = 0
	conflictDisModule = 0 
	trueIdentCtrl <- ctrlMods
	trueIdentDis <- arr.max
} 

diag(identityMatrix[trueIdentDis,trueIdentCtrl]) <- 1
identityMatrix[conflictDisModule[which.max(jacMatrix[conflictDisModule,conflictCtrlModules])],conflictCtrlModules[which.max(jacMatrix[conflictDisModule,conflictCtrlModules])]] <- 1
identityMatrix[identityMatrix != 1] <- 0
identityLink <- which(identityMatrix == 1, arr.ind = TRUE)

for (y in 1:dim(identityLink)[1]){
		fx[fx$Module==identityLink[y,1],10] <- identityLink[y,2]
	}

remainders <- cbind(unique(fx[is.na(fx$V10),3])[order(unique(fx[is.na(fx$V10),3]))],seq(1:length(unique(fx[is.na(fx$V10),3]))) + dim(jacMatrix)[2])

for (y in 1:dim(remainders)[1]){
		fx[fx$Module==remainders[y,1],10] <- remainders[y,2]
	}
colnames(fx)[10] <- "module.fixed"

return(fx)
}

reAlign.ctrl <- function(fx){
	fx$module.fixed <- fx$Module
	return(fx)
}

plotNetworkHubs <- function(nets){
	ad.deg <- getHubNodes(nets[[1]])
	ctrl.deg <- getHubNodes(nets[[2]])
	pso.deg <- getHubNodes(nets[[3]])

adTop <- ad.deg[1:5]
psoTop <- pso.deg[1:5]
ctrlTop <- ctrl.deg[1:5]

unionHubs <- union(union(names(adTop), names(psoTop)),names(ctrlTop))

hubTable <- data.frame(matrix(nrow = length(unionHubs), ncol = 3, data = 0))
rownames(hubTable) <- unionHubs#
colnames(hubTable) <- c("ADL", "CTRL", "PSOL")

hubTable[,1] <- as.numeric(ad.deg[match(rownames(hubTable), names(ad.deg))])
hubTable[,2] <- as.numeric(ctrl.deg[match(rownames(hubTable), names(ctrl.deg))])
hubTable[,3] <- as.numeric(pso.deg[match(rownames(hubTable), names(pso.deg))])

}
getHubNodes <- function(network){
	ig <- createIgraphObject(network)
	deg.net  <- degree(ig)
	deg.net <- deg.net[order(deg.net, decreasing = TRUE)]
	return(deg.net)
}

write_realignment <- function(x,y,z){
	write.table(x,file = "/home/gaz/MAARS_p2/Scripts/compositionalCor_v3/SparCC/mainNetworks_calc/AD_LES_nodes_realign.txt", quote = FALSE, sep ="\t", row.names = FALSE)
	write.table(y,file = "/home/gaz/MAARS_p2/Scripts/compositionalCor_v3/SparCC/mainNetworks_calc/CTRL_nodes_realign.txt", quote = FALSE, sep ="\t", row.names = FALSE)
	write.table(z,file = "/home/gaz/MAARS_p2/Scripts/compositionalCor_v3/SparCC/mainNetworks_calc/PSO_LES_nodes_realign.txt", quote = FALSE, sep ="\t", row.names = FALSE)

}

getDegree <- function(network){
	ig <- graph.data.frame(network, directed=F)
	topDegree <- degree(ig)[order(degree(ig),decreasing = TRUE)]
	return(topDegree)
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



	ADL.degree <- (degree(ig.net.ADL))#/length(degree(ig.net.ADL))
	CTRL.degree <- (degree(ig.net.CTRL))#/length(degree(ig.net.CTRL))
	PSOL.degree <- (degree(ig.net.PSOL))#/length(degree(ig.net.PSOL))

	unionNodes <- getUnionNodes()
	degreeTable <- data.frame(matrix(nrow = length(unionNodes), ncol = 6, data = 0))
	colnames(degreeTable) <- c("ADL.degree", "CTRL.degree", "PSOL.degree", "Diffcvad","Diffcvpso", "Diffadvpso")
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
pairwisediffConnectivity <- function(network1, network2) {
unionNodes <- getUnionNodes()
n1 <- createIgraphObject(network1)
n2 <- createIgraphObject(network2)

	n1.degree <- (degree(n1))#/length(degree(n1))
	n2.degree <- (degree(n2))#/length(degree(n2))

	degreeTable <- data.frame(matrix(nrow = length(unionNodes), ncol = 3, data = 0))
	colnames(degreeTable) <- c("n1.degree", "n2.degree", "diff")
	row.names(degreeTable) <- unionNodes
	degreeTable[,1] <- n1.degree[match(unionNodes, names(n1.degree))]
	degreeTable[,2] <- n2.degree[match(unionNodes, names(n2.degree))]
	#degreeTable[is.na(degreeTable)] <- 0
	degreeTable[,3] <- (degreeTable[,1] - degreeTable[,2])
	return(degreeTable[,3])

}

difConexPermutation <- function(permutations.netdis,permutations.netCtrl){
unionNodes <- getUnionNodes()

permDf <- data.frame(matrix(nrow = length(unionNodes),ncol = length(permutations.netdis)))
rownames(permDf) <- unionNodes
for (i in 1:length(permutations.netdis)){
	dis.net <- threshold_cor_sparCC_cor(permutations.netdis[[i]],0.21,"Sparcc")
	ctrl.net <- threshold_cor_sparCC_cor(permutations.netCtrl[[i]],0.22,"Sparcc")
	cat(paste0("perm " ,i, " done \n"))
	permDf[,i] <- pairwisediffConnectivity(dis.net,ctrl.net)
}

return(permDf)
}

calculatePermutationPval <- function(diffConect, diffPerm.ADL,diffPerm.PSOL) {
diffPerm.ADL[is.na(diffPerm.ADL)] <- 0
diffPerm.PSOL[is.na(diffPerm.PSOL)] <- 0

ad.p <- 0

for (i in 1:dim(diffConect)[1]) {
	if (diffConect[i,4] == 0){
		ad.p[i] <- 1
		}else if (diffConect[i,4] < 0) {
		ad.p[i] <-  sum(diffConect[i,4] > diffPerm.ADL[i,]) / ncol(diffPerm.ADL)
			}else {
			ad.p[i] <-  sum(diffConect[i,4] < diffPerm.ADL[i,]) / ncol(diffPerm.ADL)
			}
		}
	

	pso.p <- 0


for (i in 1:dim(diffConect)[1]) {
	if (diffConect[i,5] == 0){
		pso.p[i] <- 1
		}else if (diffConect[i,4] < 0) {
		pso.p[i] <-  sum(diffConect[i,5] > diffPerm.PSOL[i,]) / ncol(diffPerm.PSOL)
			}else {
			pso.p[i] <-  sum(diffConect[i,5] < diffPerm.PSOL[i,]) / ncol(diffPerm.PSOL)
			}
		}

	diffConect$ad.p <- ad.p
	diffConect$pso.p <- pso.p

	return(diffConect)
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


readPermutations <- function() {

	permutationStore <- list()
	perm.ADL.dir <- paste0("/home/gaz/MAARS_p2/Scripts/compositionalCor_v2/SparCC/permutations/ADL/solved/", list.files("/home/gaz/MAARS_p2/Scripts/compositionalCor_v2/SparCC/permutations/ADL/solved/"))
	perm.ADL <- lapply(perm.ADL.dir, function(x) read.table(x, header = TRUE, check.names = FALSE, row.names =1))
	cat("ADL permutations read in \n")

	perm.PSOL.dir <- paste0("/home/gaz/MAARS_p2/Scripts/compositionalCor_v2/SparCC/permutations/PSOL/solved/", list.files("/home/gaz/MAARS_p2/Scripts/compositionalCor_v2/SparCC/permutations/PSOL/solved/"))
	perm.PSOL <- lapply(perm.PSOL.dir, function(x) read.table(x, header = TRUE, check.names = FALSE,row.names =1)) 
	cat("PSOL  permutations read in \n")

# we have two sets of permutations for CTRL 
	perm.CTRL.PSOL.dir <- paste0("/home/gaz/MAARS_p2/Scripts/compositionalCor_v2/SparCC/permutations/CTRL/PSOL/solved/", list.files("/home/gaz/MAARS_p2/Scripts/compositionalCor_v2/SparCC/permutations/CTRL/PSOL/solved/"))
	perm.CTRL.PSOL <- lapply(perm.CTRL.PSOL.dir, function(x) read.table(x, header = TRUE, check.names = FALSE,row.names =1)) 
	cat("CTRL PSO permutations read in\n")

	perm.CTRL.ADL.dir <- paste0("/home/gaz/MAARS_p2/Scripts/compositionalCor_v2/SparCC/permutations/CTRL/ADL/solved/", list.files("/home/gaz/MAARS_p2/Scripts/compositionalCor_v2/SparCC/permutations/CTRL/ADL/solved/"))
	perm.CTRL.ADL <- lapply(perm.CTRL.ADL.dir, function(x) read.table(x, header = TRUE, check.names = FALSE,row.names =1)) 
	cat("CTRL AD permutations read in\n")

# NEEDS TO BE FIXED WHEN ADDING IN AD 
permutationStore[[1]] <- perm.ADL
permutationStore[[2]] <- perm.CTRL.ADL
permutationStore[[3]] <- perm.PSOL
permutationStore[[4]] <- perm.CTRL.PSOL
return(permutationStore)
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
mdNames <- 1:max(fx$module.fixed)
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


###################################
# This is where the plot scripts go
###################################

# plot of o2 proportions by module
o2PropPlot <- function(nodeDataL){ # pass a list of node atributes 
p <- list()
for (m in 1:length(nodeDataL)) {
	nodeData <- nodeDataL[[m]]
		moduleO2Counts <- data.frame(matrix(ncol =max(nodeData$module.fixed), nrow = 3, dat =0),check.names = 0)
	rownames(moduleO2Counts) <- c("Aer","AnAer", "FacAnAer")
	#colnames(moduleO2Counts) <- paste0("m",as.numeric(unique(nodeData$module.fixed)[order(as.numeric(unique(nodeData$module.fixed)))]))
	colnames(moduleO2Counts) <- paste0("m",1:max(nodeData$module.fixed))
	count <- 1
	for (l in 1:max(nodeData$module.fixed)){
		nodO2Tab <- table(nodeData[nodeData$module.fixed== l,8])
		moduleO2Counts[,count] <- data.frame(nodO2Tab[names(nodO2Tab) %in% rownames(moduleO2Counts)])
		count = count + 1
	}

	s.moduleO2Counts<- moduleO2Counts
	s.moduleO2Counts$rn <- rownames(s.moduleO2Counts)
	melt.moduleO2Counts <- melt(s.moduleO2Counts, id = "rn")

	p[[m]] <- ggplot(melt.moduleO2Counts ,aes(x=factor(variable),y=value,fill=factor(rn))) + 
    geom_bar(stat="identity",width=0.8) +labs(title=netOrder[m]) + scale_y_continuous("Frequency") + scale_x_discrete("")

	}
	mp <- multiplot(p[[1]],p[[2]],p[[3]])
	return(mp)
}

o2PropNetworkPlot <- function(nodeDataL){
o2netDf <- data.frame(matrix(nrow =length(nodeDataL), ncol = 4,data = 0))
colnames(o2netDf) <- c("Aer","AnAer", "FacAnAer","Other")
	for (m in 1:length(nodeDataL)) {
		nodeData <- nodeDataL[[m]]
		o2netDf[m,1] <- sum(nodeData$o2_tolerance == "Aer",na.rm = TRUE) / dim(nodeData)[1]
		o2netDf[m,2]<- sum(nodeData$o2_tolerance == "AnAer",na.rm = TRUE) / dim(nodeData)[1] 
		o2netDf[m,3]<- sum(nodeData$o2_tolerance == "FacAnAer",na.rm = TRUE) / dim(nodeData)[1]
		o2netDf[m,4]<- 1-sum(o2netDf[m,1],o2netDf[m,2],o2netDf[m,3])

	}
 o2netDf$rn <- c("ADL","CTRL","PSOL")
	melt.o2netdf <- melt(o2netDf)
	(p <- ggplot(melt.o2netdf ,aes(x=factor(rn),y=value,fill=factor(variable))) + 
    geom_bar(stat="identity",width=0.8) +labs(title="Proportion of O2 tolerating bacteria") + scale_y_continuous("Frequency") + scale_x_discrete(""))

return(p)
}

orderNetworkPlot <- function(nodeDatatL){
	


	orderTables <- list()

		for (m in 1:length(nodeDataL)) {
			nodeData <- nodeDataL[[m]]
							orderTables[[m]] <- table(nodeData$Order)
		}
	allNames.Order <- unique(unlist(lapply(orderTables, function(x) names(x))))
	allNames.Order <- allNames.Order[allNames.Order != "", drop = TRUE]

orderNetDf <- data.frame(matrix(nrow =length(allNames.Order), ncol = 3,data = 0))
rownames(orderNetDf) <- allNames.Order

		for (n in 1:length(nodeDataL)) {
			nodeData <- nodeDataL[[n]]
			orderNetDf[,n] <- table(nodeData$Order)[match(allNames.Order, names(table(nodeData$Order)))]
		}
orderNetDf[is.na(orderNetDf)] <- 0
# calculate proportion ### I HAD TO HACK THIS FOR SOME REASON?
orderNetDf[,1] <- orderNetDf[,1] / colSums(orderNetDf)[1]
orderNetDf[,2] <- orderNetDf[,2] / colSums(orderNetDf)[2]
orderNetDf[,3] <- orderNetDf[,3] / colSums(orderNetDf)[3]


colnames(orderNetDf) <- c("ADL","CTRL","PSOL")
orderNetDf$rn <- rownames(orderNetDf)

#melt and ggplot
melt.orderNetDf <- melt(orderNetDf)
(p <- ggplot(melt.orderNetDf,aes(x=factor(variable),y=value,fill=factor(rn))) + 
    geom_bar(stat="identity",width=0.8) +labs(title="Order proportion by network") + scale_y_continuous("Proportion") + scale_x_discrete(""))

return(p)
}

taxNetworkPlot <- function(nodeDatatL){
	
	orderTables <- list()

		for (m in 1:length(nodeDataL)) {
			nodeData <- nodeDataL[[m]]
						fullTax <- grab_taxonomy(nodeData$nodeName)
		nodeData$Phylum <- fullTax[match(nodeData$nodeName,fullTax$OTU_id),3]
		nodeData$Class <- fullTax[match(nodeData$nodeName,fullTax$OTU_id),4]
			orderTables[[m]] <- table(nodeData$Class)
		}
	allNames.Order <- unique(unlist(lapply(orderTables, function(x) names(x))))
	allNames.Order <- allNames.Order[allNames.Order != "", drop = TRUE]

orderNetDf <- data.frame(matrix(nrow =length(allNames.Order), ncol = 3,data = 0))
rownames(orderNetDf) <- allNames.Order

		for (n in 1:length(nodeDataL)) {
			nodeData <- nodeDataL[[n]]
			nodeDataL[[m]]
			fullTax <- grab_taxonomy(nodeData$nodeName)
		nodeData$Phylum <- fullTax[match(nodeData$nodeName,fullTax$OTU_id),3]
		nodeData$Class <- fullTax[match(nodeData$nodeName,fullTax$OTU_id),4]
			orderNetDf[,n] <- table(nodeData$Class)[match(allNames.Order, names(table(nodeData$Class)))]
		}
orderNetDf[is.na(orderNetDf)] <- 0
# calculate proportion ### I HAD TO HACK THIS FOR SOME REASON?
orderNetDf[,1] <- orderNetDf[,1] / colSums(orderNetDf)[1]
orderNetDf[,2] <- orderNetDf[,2] / colSums(orderNetDf)[2]
orderNetDf[,3] <- orderNetDf[,3] / colSums(orderNetDf)[3]


colnames(orderNetDf) <- c("ADL","CTRL","PSOL")
orderNetDf$rn <- rownames(orderNetDf)

#melt and ggplot
melt.orderNetDf <- melt(orderNetDf)
(p <- ggplot(melt.orderNetDf,aes(x=factor(variable),y=value,fill=factor(rn))) + 
    geom_bar(stat="identity",width=0.8) +labs(title="Class proportion by network") + scale_y_continuous("Proportion") + scale_x_discrete(""))

return(p)
}


# plot of module similarity

 moduleDegreePlot <- function(){
 	 modDegrees <- extractModuleDegrees(nets,node.fx)
	p <- list()
 	for (q in 1:length(netOrder)) {

	 	networkDegree <-  lapply(modDegrees[[q]],function(x) mean(x))
	
		netDegree.ul <- as.data.frame(unlist(networkDegree))
		netDegree.ul$rn <- rownames(netDegree.ul)

		melt.moduleDegree <- melt(netDegree.ul, id = "rn")

		p[[q]] <- ggplot(melt.moduleDegree ,aes(x=factor(rn),y=value,fill=factor(rn))) + 
	    geom_bar(stat="identity",width=0.8) +labs(title=netOrder[q]) + scale_y_continuous("Degree") + scale_x_discrete("")
	}
	mp <- multiplot(p[[1]],p[[2]],p[[3]])
	return(mp)
 }

moduleOverlapPlot <- function(node.fx){

	net1  <- node.fx[[1]]
	net2  <- node.fx[[2]]
	net3  <- node.fx[[3]]
	##change the 199 to lower numbers
	net1$module.fixed[as.numeric(net1$module.fixed) > 100] <- as.numeric(max(net1$module.fixed[as.numeric(net1$module.fixed) < 100] )) + 1
	net2$module.fixed[as.numeric(net2$module.fixed) > 100] <- as.numeric(max(net2$module.fixed[as.numeric(net2$module.fixed) < 100] )) + 2
	net3$module.fixed[as.numeric(net3$module.fixed) > 100] <- as.numeric(max(net3$module.fixed[as.numeric(net3$module.fixed) < 100] )) + 3

	net1$Module <- as.numeric(net1$module.fixed)
	net2$Module <- as.numeric(net2$module.fixed)
	net3$Module <- as.numeric(net3$module.fixed)


	olm.1 <- defOverlapMatrix(net1,net2)
	jm1 <- defJaccardMatrix(net1,net2,olm.1)

	olm.1.melt <-  melt(jm1)


	a <- (ggplot(data =  olm.1.melt, aes(x = X1, y = X2)) + 
	    geom_tile(aes(fill = value), colour = "white") + scale_fill_gradient(low='white', high='darkred',limits=c(0,1)) + 
	    scale_x_continuous("CTRL",breaks = 1:dim(olm.1)[1]) + scale_y_continuous("AD",breaks = 1:dim(olm.1)[2]))

	olm.2 <- defOverlapMatrix(net3,net2)
	jm2 <- defJaccardMatrix(net3,net2,olm.1)
	olm.2.melt <-  melt(jm2)

	b <- (ggplot(data =  olm.2.melt, aes(x = X1, y = X2)) + 
	    geom_tile(aes(fill = value), colour = "white") + scale_fill_gradient(low='white', high='darkred',limits=c(0,1)) + 
	    scale_x_continuous("CTRL", breaks = 1:dim(olm.2)[1]) + scale_y_continuous("PSO", breaks = 1:dim(olm.2)[2]))


mp <- multiplot(a,b)
return(mp)
}

moduleParallelPlot <- function(node.fx){
	net1  <- node.fx[[1]]
	net2  <- node.fx[[2]]
	net3  <- node.fx[[3]]



olp.ad_ctrl <- defOverlapMatrix.fx(net1,net2)
olp.pso_ctrl <-  defOverlapMatrix.fx(net3,net2)


# module sizes
modNodeNames <- c(paste0("A",modStats.ADL[[3]]$module.fixed),paste0("C",modStats.CTRL[[3]]$module.fixed),paste0("P",modStats.PSOL[[3]]$module.fixed))
modNodeSizes <- c(modStats.ADL[[3]]$noGenes,modStats.CTRL[[3]]$noGenes,modStats.PSOL[[3]]$noGenes)

nodeDataDF <- data.frame(matrix(ncol = 2, nrow = length(modNodeNames), data = 0))
nodeDataDF[,1] <- modNodeNames
nodeDataDF[,2] <- modNodeSizes
colnames(nodeDataDF) <- c("node", "size")


adCtrl.olp <- as.data.frame(olp.ad_ctrl)
adCtrl.olp1 <- adCtrl.olp[,-3]
colnames(adCtrl.olp1) <- paste0("A",modStats.ADL[[3]]$module.fixed)
rownames(adCtrl.olp1) <- paste0("C",modStats.CTRL[[3]]$module.fixed)
adCtrl.olp1$rn <- rownames(adCtrl.olp1)
melt.ad.df <- melt(adCtrl.olp1,id = "rn")

psoCtrl.olp <- as.data.frame(olp.pso_ctrl)
colnames(psoCtrl.olp) <- paste0("P",modStats.PSOL[[3]]$module.fixed)
rownames(psoCtrl.olp) <- paste0("C",modStats.CTRL[[3]]$module.fixed)
psoCtrl.olp$rn <- rownames(psoCtrl.olp)
melt.pso.df <- melt(psoCtrl.olp,id = "rn")

boundEdges <- rbind(melt.pso.df,melt.ad.df)

# scaling
boundEdges$value <- boundEdges$value * 10
nodeDataDF$size <- nodeDataDF$size * 10


write.table(boundEdges , file = "/home/gaz/MAARS_p2/Scripts/compositionalCor_v3/SparCC/randOutput/edges.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(nodeDataDF, file = "/home/gaz/MAARS_p2/Scripts/compositionalCor_v3/SparCC/randOutput/nodes.txt", sep = "\t", quote = FALSE, row.names = FALSE)


melt.df <- melt(olp.pso_ctrl)
colnames(melt.df) <- c("ctrlnet","psonet", "freq")


 alluvial( melt.df[,1:2], freq=melt.df$freq, border=NA, hide = melt.df$freq < quantile(melt.df$freq, .5))


ovlpDF <- data.frame(matrix(ncol = 4, nrow = 5, data = 0))
colnames(ovlpDF) <- c("module", "adnet","ctrlnet", "ad")
 ovlpDF[,1] <- 1:dim(ovlpDF)[1]
 ovlpDF[,2] <- diag(olp.ad_ctrl)
 ovlpDF[,3] <- modStats.CTRL[[3]]
 ovlpDF[,4] <- diag(olp.pso_ctrl)
 melt.ovlpDF <- melt(ovlpDF, id = "module")

modStats.CTRL[[3]]
modStats.ADL[[3]]


 alluvial(  melt.ovlpDF[,1:2], freq= melt.ovlpDF$value, border=NA, hide =  melt.ovlpDF$value < quantile( melt.ovlpDF$value, .5))


dat <- read.table("/home/gaz/Documents/rand/alluvial.csv", header = TRUE)
 #melt.ovlpDF <- melt(dat,id = "network")


 alluvial(  dat[,1:4], freq= dat$prop, border=NA, hide =  dat$prop < quantile( dat$prop, .5))



}

# module similarity 
# use correlation between eigenvectors to show correlation between modules NOT COMPLETE

moduleSimilarity <- function(node.fx){


MBdat <- readRawProportionalCountDat()
ADL.pcd <- MBdat[[1]]
CTRL.pcd <- MBdat[[2]]
PSOL.pcd <- MBdat[[3]]

	net1  <- node.fx[[1]]
	net2  <- node.fx[[2]]
	net3  <- node.fx[[3]]


corList <- list()
network.moduleList <- list()
module.pca.list <- list()
for (n in 1:length(node.fx)){ 

noc <- node.fx[[n]]

mods <- unique(noc$module.fixed)[order(unique(noc$module.fixed))]
for (i in mods) {
	noc$module.fixed == i
	network.moduleList[[i]] <- noc[noc$module.fixed == i,]
}



module.pca.table <- data.frame(matrix(ncol = max(mods), nrow =dim(MBdat[[n]])[2], data = 0 ))

	for (j in mods){		
		module.subTable <- t(MBdat[[n]][match(network.moduleList[[j]]$nodeName,rownames(MBdat[[n]])),])
		pca.module <- princomp(module.subTable)
		pca1.module <- pca.module$scores[,1]
		module.pca.table[,j] <- pca1.module
	}
module.pca.list[[n]] <- module.pca.table
corList[[n]] <- cor(module.pca.table,method = "spearman")
	}


}

# module intercorrelation
moduleSimilarity <- function(node.fx,nets){
	ADL.net <- nets[[1]]
	CTRL.net <- nets[[2]]
	PSOL.net <- nets[[3]]

	ADL.node  <- node.fx[[1]]
	CTRL.node  <- node.fx[[2]]
	PSOL.node  <- node.fx[[3]]

	
	network.moduleList.all <- list()
	for (n in 1:length(node.fx)){ 
		network.moduleList <- list()
		noc <- node.fx[[n]]
		mods <- unique(noc$module.fixed)[order(unique(noc$module.fixed))]
		for (i in mods) {
			network.moduleList[[i]] <- noc[noc$module.fixed == i,]
		}
		network.moduleList.all[[n]] <- network.moduleList
	}
		moduleCorList <- list()
		for (k in 1:length(nets)){ 
			moduleCorDf <- data.frame(matrix(nrow = length(network.moduleList.all[[k]]),ncol = 1, dat = 0))
			for (j in 1:length(network.moduleList.all[[k]])){

				module.tmp <- extractModuleLinks(nets[[k]], network.moduleList.all[[k]][[j]])
				moduleCorDf[j,1] <- mean(abs(module.tmp$cor))
			}
	moduleCorList[[k]] <- moduleCorDf
	}
}

moduleAttributes <- network.moduleList.all[[1]][[1]] # tmp 
net <- net.ADL

extractModuleLinks <- function(net,moduleAttributes){

	tmpfrom <- net[which(net$from %in% moduleAttributes$nodeName),]
	tmpto <- net[which(net$to %in% moduleAttributes$nodeName),]
	merged <- merge(tmpfrom,tmpto)
return(merged)
}

net.1 <- CTRL.net
net.2 <- CTRL.net

net.1.nodes <- CTRL.node
net.2.nodes <- CTRL.node


extractModuleExtraCorrelation <- function(net.1, net.1.nodes){
modules <- unique(net.1.nodes$module.fixed)[order(unique(net.1.nodes$module.fixed))]
posMatrix <- data.frame(matrix(nrow = max(modules), ncol = max(modules), data = 0))
negMatrix <- data.frame(matrix(nrow = max(modules), ncol = max(modules), data = 0))
ovrMatrix <- data.frame(matrix(nrow = max(modules), ncol = max(modules), data = 0))

m <- 1
n <- 1

for (m in 1:max(modules)){ 
	moduleNodes <- net.1.nodes[net.1.nodes$module.fixed == m,1]
	s1 <- net.1[which(net.1$from %in% moduleNodes),]
	s2 <- net.1[which(net.1$to %in% moduleNodes),]

	bound <- rbind(s1,s2) 
	bound$fromModule <- net.1.nodes[match(bound$from,net.1.nodes$nodeName),10]
	bound$toModule <- net.1.nodes[match(bound$to,net.1.nodes$nodeName),10]
	for (n in 1:max(modules)){
		if (m == n){
			tmpNeg <- 0
			tmpPos <- 0
			tmpOvr <- 0
			}else {
			interCor <- bound[paste0(bound$fromModule, bound$toModule) != paste0(m,m),]
			tmpCors <- c(interCor[interCor$fromModule == n,3],interCor[interCor$toModule == n,3])

			tmpPos <- mean(tmpCors[tmpCors >= 0])
			tmpNeg <- mean(tmpCors[tmpCors <= 0])
			tmpOvr <- mean(abs(tmpCors))
			tmpNeg[is.na(tmpNeg)] <- 0
			tmpPos[is.na(tmpPos)] <- 0
			tmpOvr[is.na(tmpOvr)] <- 0
			} 
		posMatrix[n,m] <- tmpPos
		negMatrix[n,m] <- tmpNeg
		ovrMatrix[n,m] <- tmpOvr

	}
	
}

}
# this didnt really work archived
moduleEigenvectorCorrelation <- function(){

# for a module in a network get expression and calclate eigenvector, put in data frame for correlation

modulePC_lists <- list()


for (i in 1:length(node.fx)) {

	nodeOC <- node.fx[[i]]
	rawOC <- aa[[i]]
	netRaw.tmp <- rawOC[match(rawOC$nodeName, rownames(rawOC)),]
	noModules.tmp <- max(nodeOC$module.fixed)
	networkPcList.tmp <- data.frame(matrix(ncol = noModules.tmp,nrow = dim(rawOC)[2], data =0))
	modStat.tmp <- calcModuleStats(node.fx[[i]])
	for (j in 1:noModules.tmp){
			
			
			moduleOfChoice <- modStat.tmp[[2]][[j]]$nodeName
			if (length(moduleOfChoice) <1) {
				networkPcList.tmp[,j] <- 0
			} else {
			networkPcList.tmp[,j] <- prcomp(t(rawOC[match(moduleOfChoice, rownames(rawOC)),]))$x[,1]
		}
	rownames(networkPcList.tmp) <- colnames(rawOC)
	modulePC_lists[[i]] <- networkPcList.tmp
	}

}
 


}


extractModuleDegrees <- function(nets,node.fx){

network.moduleDegrees <- list()
for (n in 1:length(nets)){
	net <- createIgraphObject(nets[[n]])
	net.degree <- (degree(net))
	moduleDegrees <- list()
		for (i in 1:max(node.fx[[n]]$module.fixed)){
			moduleDegrees[[i]] <- net.degree[which(node.fx[[n]][match(names(net.degree), node.fx[[n]]$nodeName),10] == i)]
		}
	network.moduleDegrees[[n]] <- moduleDegrees

	}
return(network.moduleDegrees)
}

keystoneCtrlSpecies <- function(nets){
topBetweenTable <- data.frame(matrix(nrow = 10,ncol =3, data = 0))
colnames(topBetweenTable) <- netOrder
ig.PSOL <- createIgraphObject(net.PSOL)
	ig.CTRL <- createIgraphObject(net.CTRL)
		ig.ADL <- createIgraphObject(net.ADL)

	between.ADL<-betweenness(ig.ADL,normalized = TRUE,directed = FALSE)#evcent(ig)$vector# 
	between.CTRL<-betweenness(ig.CTRL,normalized = TRUE,directed = FALSE)#evcent(ig)$vector# 
	between.PSOL<-betweenness(ig.PSOL,normalized = TRUE,directed = FALSE)#evcent(ig)$vector# 

between.CTRL <- between.CTRL[order(between.CTRL,decreasing = TRUE)]
keyStone <- between.CTRL[1:10]

return(keyStone)


}


plotBetweenessCentrality <- function(nets){
topBetweenTable <- data.frame(matrix(nrow = 10,ncol =3, data = 0))
colnames(topBetweenTable) <- netOrder
ig.PSOL <- createIgraphObject(net.PSOL)
	ig.CTRL <- createIgraphObject(net.CTRL)
		ig.ADL <- createIgraphObject(net.ADL)

	between.ADL<-betweenness(ig.ADL,normalized = TRUE,directed = FALSE)#evcent(ig)$vector# 
	between.CTRL<-betweenness(ig.CTRL,normalized = TRUE,directed = FALSE)#evcent(ig)$vector# 
	between.PSOL<-betweenness(ig.PSOL,normalized = TRUE,directed = FALSE)#evcent(ig)$vector# 

between.CTRL <- between.CTRL[order(between.CTRL,decreasing = TRUE)]
between.CTRL[1:10]
rownames(topBetweenTable) <- names(between.CTRL[1:10])
topBetweenTable[,2] <- between.CTRL[1:10]
topBetweenTable[,1] <- between.ADL[match(rownames(topBetweenTable), names(between.ADL))]
topBetweenTable[,3] <- between.PSOL[match(rownames(topBetweenTable), names(between.PSOL))]
topBetweenTable[is.na(topBetweenTable)] <- 0

topBetweenTable$rn <- rownames(topBetweenTable)
melt.topBetweenTable <- melt(topBetweenTable)
melt.topBetweenTable$rn <- factor(as.character(melt.topBetweenTable$rn), levels = unique(melt.topBetweenTable$rn))
melt.topBetweenTable$join <- grab_taxonomy(melt.topBetweenTable$rn)$join
melt.topBetweenTable$join <- factor(as.character(melt.topBetweenTable$join), levels = unique(melt.topBetweenTable$join))

plot <- ggplot(data=melt.topBetweenTable,aes(x=as.factor(join), y=value, fill=variable)) +  geom_bar(stat="identity", position = "dodge") + labs(title="Betweeness Centrality") + scale_y_continuous("Betweeness Centrality") + scale_x_discrete("")
return(plot)
}


keystoneOccurence <- function(nets){

	keyStone <- keystoneCtrlSpecies(nets)
b <- 2

	fromKeystoneLinks <- nets[[b]][which(nets[[b]]$from %in% names(keyStone)),]
	toKeystoneLinks <- nets[[b]][which(nets[[b]]$to %in% names(keyStone)),]
	boundLinks <- rbind(fromKeystoneLinks,toKeystoneLinks)

	boundLinks.r1 <- boundLinks[which(boundLinks$from %in% names(keyStone)),]
	boundLinks.r2 <- boundLinks.r1[which(boundLinks.r1$to %in% names(keyStone)),]
	boundLinks.r2 <- boundLinks.r2[!duplicated(boundLinks.r2),]

corMatrix <- data.frame(matrix(nrow = length(keyStone), ncol = length(keyStone), data = 0))
rownames(corMatrix) <- names(keyStone)
colnames(corMatrix) <- names(keyStone)

write.table(boundLinks.r2, file = "/home/gaz/MAARS_p2/Scripts/compositionalCor_v3/Paper/Figures/presentation/CTRLKeystone.txt", row.names = FALSE, quote = FALSE, sep = "\t")
}