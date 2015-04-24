
setwd("~/MAARS_p2/MB")
library("limma")
library("igraph")
library("minet")
library(SpiecEasi)
home <- "~/MAARS_p2/Scripts/"
folder <- "compositionalCor_v3"
#############################################################################################################
###################################### PARAMETERS ###########################################################
#############################################################################################################
#########LOAD THE AFFY DATA
descriptionFile <- read.table(file="~/MAARS_p2/Affy/descriptionFile140414.csv", header = TRUE, fill = TRUE)
#########LOAD THE AFFY DATA
affydescriptionFile <- read.table(file="~/MAARS_p2/Affy/descriptionFile140414.csv", header = TRUE, fill = TRUE)
affyTable <- read.table(file="~/MAARS_p2/Affy/MAARS_normTranscriptome_618samples_16042014.txt", header = TRUE)


#########LOAD THE OTU DATA ## THIS OTU DATA IS UN NORMALISED
otuTable <- read.table("/home/gaz/MAARS_p2/MB/OTU_tableFri_Apr_04_14_23_44_CEST_2014", header = TRUE,row.names = 1, stringsAsFactors = FALSE)
#otuLabel <- read.table(file="OTU_labelFri_Apr_04_14:23:44_CEST_2014", header = TRUE, fill = TRUE,stringsAsFactors = FALSE)
#otuTable <- read.table("/home/gaz/MAARS_p2/MB/otu_table_proportions.csv", header = TRUE,row.names = 1, stringsAsFactors = FALSE)
otuLabel <- read.table(file="otuLabel_exid.txt", header = TRUE, fill = TRUE,stringsAsFactors = FALSE)
otudescriptionFile <- read.table(file="descriptionFile040414v3.csv", header = TRUE)
descriptionFile <- read.table(file="descriptionFile040414v3.csv", header = TRUE)
OTU_O2_tolerance <- read.table(file="/home/gaz/MAARS_p2/MB/OTU_O2_tolerance2.csv", header = TRUE, fill = TRUE)


##parameters for ARGS

cohort <- "PSO_NON_LES"

inference <- 2 #1 for pearson, 2 for spearman, 3 for mutual information

############sort out O2 tolerance #########
OTU_O2_tolerance[OTU_O2_tolerance$O2_tolerance == "strict-Aer",2] <- "Aer"
OTU_O2_tolerance[OTU_O2_tolerance$O2_tolerance == "strict-AnAer",2] <- "AnAer"
OTU_O2_tolerance[OTU_O2_tolerance$O2_tolerance == "Aer?",2] <- "Aer"

#############################################################
#
### Select OTU cohorts
#
############################################################
if (cohort == "CTRL"){
v_otuDescription <- otudescriptionFile[otudescriptionFile$clinical_group == cohort,] ## create a variable of the otu description file
} else if (cohort == "PSO_LES" ){
	v1_otuDescription <- otudescriptionFile[otudescriptionFile$clinical_group == "PSO",] 
	v_otuDescription <- v1_otuDescription[v1_otuDescription$lesional == "LES",] 
}else if (cohort == "AD_LES" ){
	v1_otuDescription <- otudescriptionFile[otudescriptionFile$clinical_group == "AD",] 
	v_otuDescription <- v1_otuDescription[v1_otuDescription$lesional == "LES",] 
}else if (cohort == "AD_NON_LES" ){
	v1_otuDescription <- otudescriptionFile[otudescriptionFile$clinical_group == "AD",] 
	v_otuDescription <- v1_otuDescription[v1_otuDescription$lesional == "NON_LES",] 
}else if (cohort == "PSO_NON_LES" ){
	v1_otuDescription <- otudescriptionFile[otudescriptionFile$clinical_group == "PSO",] 
	v_otuDescription <- v1_otuDescription[v1_otuDescription$lesional == "NON_LES",] 
}

bu <- otuTablez




calculateOtuProportions <- function(otuTablez){
otuTablez <- apply(otuTablez,2, function(x) round((x/sum(x,na.rm = TRUE) * 1000), digits = 2))

return(otuTablez)
}
 

###create a cohort based upon clinical group --- pass a description file
OTU_createClinicalGroupCohort <- function(desc, print = TRUE) {
	cutTab <- otuTable[,match(desc$sample_id,colnames(otuTable))]
return(cutTab)
}


v_otuTable <- OTU_createClinicalGroupCohort(v_otuDescription) 

assign_node_attributes <- function(otuVector){
	nodeData <- data.frame(nodeName = otuVector)
	nodeData$otuID <- intersectListTaxonomy[match(otuVector, intersectListTaxonomy$OTU_id),1]
	nodeData$Module <- moduleAssignments.fixed[match(otuVector, moduleAssignments.fixed$node),3]
	nodeData$Symbol <- paste0(substr(intersectListTaxonomy$Genus,0,5), ".", intersectListTaxonomy$Species)
	nodeData$Order <- intersectListTaxonomy[match(otuVector, intersectListTaxonomy$OTU_id),5]
	nodeData$Family <- intersectListTaxonomy[match(otuVector, intersectListTaxonomy$OTU_id),6]
	nodeData$Genus <- intersectListTaxonomy[match(otuVector, intersectListTaxonomy$OTU_id),7]
	nodeData$Species <- intersectListTaxonomy[match(otuVector, intersectListTaxonomy$OTU_id),8]
	nodeData$o2_tolerance <-OTU_O2_tolerance[match(intersectListTaxonomy[match(uniqueNodes, intersectListTaxonomy$OTU_id),1],OTU_O2_tolerance$OTU_ID),2]
	nodeData$pubName <- paste0(substr(nodeData$Genus,1,1),".", " ", nodeData$Species)

	#o2 tollerance
	#degree
return(nodeData)
}

filter_low_abundant_otu <- function(otuTable, minNonZeroOtu){

#######filter redundant OTU
	nonRedundant<- rowSums(abs(contrast[1:dim(contrast)[1],1:dim(contrast)[2]-1])>0)>minNonZeroOtu 
	contrastRem<- contrast[nonRedundant,]

	data <- apply(contrastRem,2, as.numeric)
	rownames(data) <- rownames(contrastRem)
return(data)
}


grab_taxonomy <- function(otuVector){
	intersectListTaxonomy <- data.frame(stringsAsFactors = FALSE)
	intersectListTaxonomy <- otuLabel[match(otuVector, otuLabel$OTU_id),]
	intersectListTaxonomy$Family[intersectListTaxonomy$Family == ""] = "F."
	intersectListTaxonomy$Genus[intersectListTaxonomy$Genus == ""] = "G."
	intersectListTaxonomy$Species[intersectListTaxonomy$Species == ""] = "sp."

return(intersectListTaxonomy)
}


threshold_cor_sparCC_cor <- function(matrix,cor_threshold, connex){

	matrix[lower.tri(matrix, diag = TRUE)] <- 0
	rowz <- which(abs(matrix) >= cor_threshold , arr.ind = TRUE)
	rowz.nr <- as.data.frame(rowz[rowz[,1] != rowz[,2],])

	linksDat <- data.frame(fromNode = rownames(rowz.nr), toNode = NA,	 cor = NA, stringsAsFactors = FALSE)
		for (g in 1:dim(rowz.nr)[1]){
			linksDat[g,2] <- rownames(matrix[rowz.nr[g,2],])
			linksDat[g,3] <- matrix[rowz.nr[g,1],rowz.nr[g,2]]
			#linksDat[g,4] <- sparCC.sig[rowz.nr[g,1],rowz.nr[g,2]]
		}

return(linksDat)
} #http://www.biomedcentral.com/1756-0500/2/240


set.seed(1)
##### create a network for "otuTable1"
contrast<- v_otuTable

dat0 <- filter_low_abundant_otu(contrast,round(dim(contrast)[2]/20))

dat0 <- calculateOtuProportions(dat0)

#res is just for exporting to sparcc Format
res <- cbind(Row.Names = rownames(dat0), dat0)
#res <- dat0
write.table(res, file = paste0(home,folder,"/SparCC/",cohort,"/raw_dat_",cohort), sep ="\t", quote = FALSE,row.names = FALSE)

system(paste0("sh /home/gaz/MAARS_p2/Scripts/compositionalCor_v3/SparCC/",cohort,"/execSparC.sh")) ## joblib for parralelisation
#system(paste0("python /home/gaz/Documents/SparCC/SparCC.py ", paste0(home,folder,"/sparCC/",cohort)," -i 1 --cor_file=", paste0("/home/gaz/MAARS_p2/Scripts/compositionalCor/sparCC/",cohort,"_sparcc.txt")))

#system(paste0("python /home/gaz/Documents/SparCC/MakeBootstraps.py ", paste0(home,folder,"/sparCC/",cohort, "_sparcc.txt")," -n 10 -o", paste0("/home/gaz/MAARS_p2/Scripts/compositionalCor/sparCC/Resamplings/boot")))


# library(igraph)
# path="/home/gaz/MAARS_p2/Scripts/compositionalCor/sparCC/pvals_two_sided.txt"
# pvals=read.table(path,header=TRUE,sep="\t")
# pvals.mat=pvals[,2:ncol(pvals)]
# # set p-values of 0 to a non-zero, small p-value so we can take the logarithm
# pvals.mat[pvals.mat==0]=0.000000001
# # convert into significance
# sig.mat=-1*log10(pvals.mat) 
# # remove all edges with significance below 1
# sig.mat[sig.mat<1]=0
# sig.mat=as.matrix(sig.mat)
# # convert adjacency matrix into a graph
# sparcc.graph=graph.adjacency(sig.mat,mode="undirected")
# # display the graph
# layout=layout.spring
# plot(sparcc.graph, layout=layout)


# colnames(sig.mat) <- substr(colnames(sig.mat),2,nchar(colnames(sig.mat)))
# rownames(sig.mat) <- colnames(sig.mat)


# linksDat.nr <- threshold_cor_sparCC(sig.mat, 1,"SparCC")

###### what about the correlation coefficient
sparCC.cor <- read.table(paste0("/home/gaz/MAARS_p2/Scripts/compositionalCor/",cohort, "/SparCC/",cohort,"_sparcc.txt"),header=1,sep="\t",row.names = 1,check.names= FALSE)

#sparCC.sig <- read.table(paste0("/home/gaz/MAARS_p2/Scripts/compositionalCor/",cohort, "/sparCC/",cohort,"_pvals_two_sided.txt"),header =1,sep="\t",row.names = 1,check.names= FALSE)

sparCC.cor2 <- sparCC.cor
#sparCC.cor[sparCC.sig <0.1] <- 0
linksDat.nr <- threshold_cor_sparCC_cor(sparCC.cor, 0.22,"SparCC" )




mapper <- as.data.frame(union(linksDat.nr$fromNode,linksDat.nr$toNode))

colnames(dat0) <- 0:(dim(dat0)[2]-1)
mapper[,2] <- 0:(dim(mapper)[1]-1)
colnames(mapper) <- c("OTU_ID","OTU_new")

edgeData <- linksDat.nr
edgeData[,3] <- round(abs(edgeData[,3]*100))

edgeData$fromNode <- mapper[match(edgeData$fromNode,mapper$OTU_ID),2]
edgeData$toNode <- mapper[match(edgeData$toNode,mapper$OTU_ID),2]

edgeData <- edgeData[order(edgeData$fromNode),]

####new cliques####

	ig <-graph.data.frame(edgeData, directed=F)
	 maximal.cliques.count(ig)






write.table(edgeData, file = paste0(home,folder,"/Output/edgeData_map_",cohort), sep ="\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(mapper, file =  paste0(home,folder,"/mappers/mapper_", cohort, ".txt"), sep = "\t", quote = FALSE)


###exec louvain
system(paste0("/home/gaz/Community_latest/./convert -i ", paste0(home,folder,"/Output/edgeData_map_",cohort), " -o ", paste0(home,folder,"/Output/edgeData_map_",cohort,".bin -w " ), paste0(home,folder,"/Output/edgeData_map_",cohort,".weights")))

system(paste0("/home/gaz/Community_latest/./community ", paste0(home,folder,"/Output/edgeData_map_",cohort,".bin"), " -l -1 -w " , paste0(home,folder,"/Output/edgeData_map_",cohort,".weights "),"> ", paste0(home,folder,"/Output/edgeData_map_",cohort,".tree")))


system(paste0("/home/gaz/Community_latest/./hierarchy ", paste0(home,folder,"/Output/edgeData_map_",cohort,".tree"), " -l 2 > ", paste0(home,folder,"/Output/edgeData_map_",cohort,".l3")))

modules <- read.table(paste0(home,folder,"/Output/edgeData_map_",cohort,".l3"), sep = " ", header= FALSE)



########get a list of all nodes in the network
uniqueNodes <- unique(union(linksDat.nr$fromNode, linksDat.nr$toNode))

###### this is a list of taxonomy for the unique node set
intersectListTaxonomy <- grab_taxonomy(uniqueNodes)



#load modules and fix
moduleAssignments <- modules
moduleAssignments$OTU_id <- mapper[match(moduleAssignments$V1, mapper$OTU_new),1]
moduleAssignments$v3 <- moduleAssignments$v3 + 1

moduleAssignments.fixed <- data.frame(node = moduleAssignments$OTU_id, OTU_new =moduleAssignments[,1],  module = moduleAssignments$v3)

######asign node attribues
nodeData <- assign_node_attributes(uniqueNodes) 

#write nodes and edges.
#write.table(linksDat.nr, file = paste0(home,folder,"/mainNetworks/edgeData_",cohort), quote = FALSE, sep = "\t", row.names = FALSE)
#write.table(nodeData, file = paste0(home,folder,"/mainNetworks/nodeData_",cohort), quote = FALSE, sep = "\t", row.names = FALSE)



