###### read in networks and create node attributes #####
source("/home/gaz/MAARS_p2/Scripts/compositionalCor_v3/Scripts/inference_functions.R")

setwd("/home/gaz/MAARS_p2/MB")

datList <- load_MB_data()
otuTable <- datList[[1]]
otuLabel <- datList[[2]]
otudescriptionFile <-datList[[3]]
descriptionFile <- datList[[4]]
OTU_O2_tolerance <- datList[[5]]

############sort out O2 tolerance #########
OTU_O2_tolerance[OTU_O2_tolerance$O2_tolerance == "strict-Aer",2] <- "Aer"
OTU_O2_tolerance[OTU_O2_tolerance$O2_tolerance == "strict-AnAer",2] <- "AnAer"
OTU_O2_tolerance[OTU_O2_tolerance$O2_tolerance == "Aer?",2] <- "Aer"

nets <- readNetworks()
modules <- assignCommunities(nets) # pass the networks to the community assignment function
# read in these networks 
net.ADL <- nets[[1]]
# net.ADNL <- read.table(file="/home/gaz/MAARS_p2/Scripts/compositionalCor/mainNetworks/edgeData_AD_NON_LES", header = TRUE)
net.CTRL <- nets[[2]] 
# net.PSONL <- read.table(file="/home/gaz/MAARS_p2/Scripts/compositionalCor/mainNetworks/edgeData_PSO_NON_LES", header = TRUE)
net.PSOL <- nets[[3]]

# get every node in all networks to calculate the intersectlistTaonomy
allNodes <- union(union(union(net.ADL$to, net.ADL$from) , union(net.PSOL$to, net.PSOL$from)) , union(net.CTRL$to, net.CTRL$from))

intersectListTaxonomy <- grab_taxonomy(allNodes)
# how many networks are being used? 

cohorts <- list(net.ADL,net.CTRL,net.PSOL)

# get attributes for all networks and store them in list
nodeDataLists <- list()
for (i in 1:length(cohorts)){

	assign("network",cohorts[[i]])
	nodes.network <- union(network$from, network$to)
	nodeData.network <- assign_node_attributes(nodes.network, modules[[i]])
	nodeDataLists[[i]] <- nodeData.network
}

# write these node data tables out 

write.table(nodeDataLists[[1]], file = "/home/gaz/MAARS_p2/Scripts/compositionalCor_v3/SparCC/mainNetworks_calc/AD_LES_nodes.txt", quote = FALSE, sep = "\t", row.names = FALSE)
write.table(nodeDataLists[[2]], file = "/home/gaz/MAARS_p2/Scripts/compositionalCor_v3/SparCC/mainNetworks_calc/CTRL_nodes.txt", quote = FALSE, sep = "\t", row.names = FALSE)
write.table(nodeDataLists[[3]], file = "/home/gaz/MAARS_p2/Scripts/compositionalCor_v3/SparCC/mainNetworks_calc/PSOL_nodes.txt", quote = FALSE, sep = "\t", row.names = FALSE)