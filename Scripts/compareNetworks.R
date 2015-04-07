setwd("~/MAARS_p2/MB")
source("/home/gaz/MAARS_p2/Scripts/compositionalCor_v3/Scripts/inference_functions.R")
source("/home/gaz/MAARS_p2/Scripts/compositionalCor_v3/Scripts/analysis_functions.R")

nets <- readNetworks() 
# read in these networks 
net.ADL <- nets[[1]]
# net.ADNL <- read.table(file="/home/gaz/MAARS_p2/Scripts/compositionalCor/mainNetworks/edgeData_AD_NON_LES", header = TRUE)
net.CTRL <- nets[[2]] 
# net.PSONL <- read.table(file="/home/gaz/MAARS_p2/Scripts/compositionalCor/mainNetworks/edgeData_PSO_NON_LES", header = TRUE)
net.PSOL <- nets[[3]]


node <- readAttributes() 
# read in these networks 
node.ADL <- node[[1]]
# net.ADNL <- read.table(file="/home/gaz/MAARS_p2/Scripts/compositionalCor/mainNetworks/edgeData_AD_NON_LES", header = TRUE)
node.CTRL <- node[[2]] 
# net.PSONL <- read.table(file="/home/gaz/MAARS_p2/Scripts/compositionalCor/mainNetworks/edgeData_PSO_NON_LES", header = TRUE)
node.PSOL <- node[[3]]


# call functions for module realigning
olm.CTRL_PSOL <- defOverlapMatrix(node.CTRL ,node.PSOL)
olm.CTRL_ADL <- defOverlapMatrix(node.CTRL ,node.ADL)
# pass module memberships for networks and overlap
jac.CTRL_PSOL <- defJaccardMatrix(node.CTRL ,node.PSOL,olm.CTRL_PSOL) 
jac.CTRL_ADL <- defJaccardMatrix(node.CTRL ,node.ADL,olm.CTRL_ADL) 

# realign
node.PSOL.fx <- reAlign(node.PSOL,jac.CTRL_PSOL)
node.ADL.fx <- reAlign(node.ADL,jac.CTRL_ADL)
node.CTRL.fx <- reAlign.ctrl(node.CTRL)

# write out to file
write_realignment(node.ADL.fx,node.CTRL.fx,node.PSOL.fx)


#define the order of the lists 
netOrder <- c("AD L","CTRL","PSO L")
# call network stats
netStats <- calcStats(node,nets, netOrder)
# call differential connectivity
diffConect <- diffConnectivity(nets)


###get module statistics
modStats.ADL <- calcModuleStats(node.ADL.fx)
modStats.CTRL <- calcModuleStats(node.CTRL.fx)
modStats.PSOL <- calcModuleStats(node.PSOL.fx)












rawUnionDat <- readRawProportionalCountDat()
AD.raw <- t(rawUnionDat[[1]])
CTRL.raw <- t(rawUnionDat[[2]])
PSO.raw <- t(rawUnionDat[[3]])

# exec one SparCC

#sparCC.PSO <- sparcc(PSO.raw) 