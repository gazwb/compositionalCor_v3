setwd("~/MAARS_p2/MB")
source("/home/gaz/MAARS_p2/Scripts/compositionalCor_v3/Scripts/inference_functions.R")
source("/home/gaz/MAARS_p2/Scripts/compositionalCor_v3/Scripts/analysis_functions.R")
source("/home/gaz/rand/multiplot/multiplot.R")

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

mbdat <- load_MB_data()
getHubNodes(net.ADL)


# call functions for module realigning
olm.CTRL_PSOL <- defOverlapMatrix(node.CTRL ,node.PSOL)
olm.CTRL_ADL <- defOverlapMatrix(node.CTRL ,node.ADL)
# pass module memberships for networks and overlap
jac.CTRL_PSOL <- defJaccardMatrix(node.CTRL ,node.PSOL,olm.CTRL_PSOL) 
jac.CTRL_ADL <- defJaccardMatrix(node.CTRL ,node.ADL,olm.CTRL_ADL) 

# realign
node.PSOL.fx <- reAlign.2(node.PSOL,jac.CTRL_PSOL)
node.ADL.fx <- reAlign.2(node.ADL,jac.CTRL_ADL)
node.CTRL.fx <- reAlign.ctrl(node.CTRL)

# create a new list of the network attributes with the fixed nodes
node.fx <- list(node.ADL.fx,node.CTRL.fx,node.PSOL.fx)

# write out to file
#write_realignment(node.ADL.fx,node.CTRL.fx,node.PSOL.fx)


#define the order of the lists 
netOrder <- c("AD L","CTRL","PSO L")
# call network stats
netStats <- calcStats(node,nets, netOrder)
# call differential connectivity
diffConect <- diffConnectivity(nets)
diffConectivityPerms <- readPermutations() 

perm.ADL <- diffConectivityPerms[[1]]
perm.CTRL.ADL <- diffConectivityPerms[[2]]
perm.PSOL <- diffConectivityPerms[[3]]
perm.CTRL.PSOL <- diffConectivityPerms[[4]]

diffPerm.ADL <- difConexPermutation(perm.ADL,perm.CTRL.ADL)
diffPerm.PSOL <- difConexPermutation(perm.PSOL,perm.CTRL.PSOL)

# function to calculate p values
diffConnect.P <- calculatePermutationPval(diffConect, diffPerm.ADL,diffPerm.PSOL)
tax <- mbdat[[2]][match(rownames(diffConnect.P),mbdat[[2]]$OTU_id),]

#fix tax
tax[tax$Family == "",6] <- "F."
tax[tax$Genus == "",7] <- "G."
tax[tax$Species == "",8] <- "sp."
tax$joinName <- as.factor(paste0(tax$Genus, " ", tax$Species))

## what modules are these nodes in?
cms <- node.CTRL.fx[match(rownames(diffConnect.P),node.CTRL.fx$nodeName),10]
pms <- node.PSOL.fx[match(rownames(diffConnect.P),node.PSOL.fx$nodeName),10]
ams <- node.ADL.fx[match(rownames(diffConnect.P),node.ADL.fx$nodeName),10]

diffConnect.tax.P <- cbind(tax[,c(6,7,8,10)],ams,cms,pms, diffConnect.P)
rownames(diffConnect.tax.P) <- rownames(diffConnect.P)

# diffConnect.tax.P[diffConnect.tax.P$ad.p < 0.05,]
# diffConnect.tax.P[diffConnect.tax.P$pso.p < 0.05,]


fold.changes <- 	calcDifferentialAbundance(diffConnect.tax.P)
diffConnect.tax.P.fc <- cbind(diffConnect.tax.P ,fold.changes)

# plot the differential connectivity by module
AD.diffconnect <- diffConnect.tax.P.fc[!is.na(diffConnect.tax.P.fc$ams),]
PSO.diffconnect <- diffConnect.tax.P.fc[!is.na(diffConnect.tax.P.fc$pms),]

AD.diffconnect$ad.BH <- p.adjust(AD.diffconnect$ad.p,method = "BH")
PSO.diffconnect$pso.BH <- p.adjust(PSO.diffconnect$pso.p,method = "BH")
# set color scale
PSO.diffconnect$pms <- as.factor(PSO.diffconnect$pms)
myColors <- brewer.pal(8,"Dark2")
names(myColors) <- levels(PSO.diffconnect$pms)
colScale <- scale_colour_manual(name = "pms",values = myColors)

# plot
qplot(Diffcvpso, PSOvCTRL.FC, data=PSO.diffconnect, colour=as.factor(pms),xlab = "Differential Connectivity", ylab = "log10(Fold Change)", main = "PSO Differential connectivity vs Fold change") + geom_point(aes(size = 1)) + colScale
qplot(Diffcvad, ADvCTRL.FC, data=AD.diffconnect, colour=as.factor(ams),xlab = "Differential Connectivity", ylab = "log10(Fold Change)",main = "AD Differential connectivity vs Fold change")+ geom_point(aes(size = 1)) + colScale

# get significant differentially connected bacterias

sig.dc.PSO <-  PSO.diffconnect[PSO.diffconnect$pso.BH < 0.1,]
sig.dc.AD <-  AD.diffconnect[AD.diffconnect$ad.BH < 0.1,]


# factor for order 
sig.dc.PSO <- sig.dc.PSO[order(sig.dc.PSO$Diffcvpso,decreasing = FALSE),]
sig.dc.PSO$joinName<- as.character(sig.dc.PSO$joinName)
sig.dc.PSO$joinName <- factor(sig.dc.PSO$joinName, levels=unique(sig.dc.PSO$joinName))

# factor for order 
sig.dc.AD<- sig.dc.AD[order(sig.dc.AD$Diffcvad,decreasing = FALSE),]
sig.dc.AD$joinName<- as.character(sig.dc.AD$joinName)
sig.dc.AD$joinName <- factor(sig.dc.AD$joinName, levels=unique(sig.dc.AD$joinName))
sig.dc.AD$rn <- rownames(sig.dc.AD)
sig.dc.AD$rn <- factor(sig.dc.AD$rn, levels=unique(sig.dc.AD$rn))

ggplot(data=sig.dc.PSO, aes(x=joinName, y=Diffcvpso, fill=joinName)) + geom_bar(colour="black", stat="identity") + guides(fill=FALSE) + labs(y = "Differential Connectivity", x = "") + ggtitle("PSO v CTRL") + coord_flip()
ggplot(data=sig.dc.AD, aes(x=rn, y=Diffcvad, fill=joinName)) + geom_bar(colour="black", stat="identity") + guides(fill=FALSE) + scale_x_discrete(labels = sig.dc.AD$joinName)+ labs(y = "Differential Connectivity", x = "") + ggtitle("AD v CTRL") + coord_flip() 



# more plots
moduleDegreePlot()
plotBetweenessCentrality(nets)

# extract positive and negative subgraphs
# 1 is positive links, 2 is negative links, 3 iis positive nodes, 4 is negative nodes
directionalSubgraphs <- extractDirectionalSubgraphs(nets,node.fx)



###get module statistics
modStats.ADL <- calcModuleStats(node.ADL.fx)
modStats.CTRL <- calcModuleStats(node.CTRL.fx)
modStats.PSOL <- calcModuleStats(node.PSOL.fx)


o <- o2PropNetworkPlot(node.fx)
# plot of o2 proportions per modules
p <- o2PropPlot(node.fx) # works now

# plot of module similarity
j <- moduleOverlapPlot(node.fx)

k <- orderNetworkPlot(node.fx)
l <- taxNetworkPlot(node.fx)




readPermutations




rawUnionDat <- readRawProportionalCountDat()
AD.raw <- t(rawUnionDat[[1]])
CTRL.raw <- t(rawUnionDat[[2]])
PSO.raw <- t(rawUnionDat[[3]])

# exec one SparCC

#sparCC.PSO <- sparcc(PSO.raw) 