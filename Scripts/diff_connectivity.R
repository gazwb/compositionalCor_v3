# differential connectivity



setwd("~/MAARS_p2/MB")
source("/home/gaz/MAARS_p2/Scripts/compositionalCor_v3/Scripts/inference_functions.R")
source("/home/gaz/MAARS_p2/Scripts/compositionalCor_v3/Scripts/analysis_functions.R")

#permutation functions

rawList <- readRawProportionalCountDat()

raw.ADL <- rawList[[1]]
raw.CTRL <- rawList[[2]]
raw.PSOL <- rawList[[3]]


set.seed(1)

# do ADL and CTRL first
# join ADL and CTRL together


createRandomNetwork <- function(raw.n1, raw.n2){ 
	join.n1.n2 <- cbind(raw.n1, raw.n2)
	# randomly select the same numbers of samples used to infer true network
	randOrder <- sample(dim(join.n1.n2)[2])

	n1.sampIndex <- randOrder[1:dim(raw.n1)[2]]
	n2.sampIndex <- randOrder[(dim(raw.n1)[2]+1):length(randOrder)]

	tmp.n1samp <- join.n1.n2[,n1.sampIndex]
	tmp.n2samp <- join.n1.n2[,n2.sampIndex]

	permTemp <- list(tmp.n1samp, tmp.n2samp)
return(permTemp)
}

createPermutations <- function(raw.n1,raw.n2, no.perm, c1, c2){

if ((c1 == "ADL") && (c2 == "CTRL")){
	folder <- "ADL"
} else if ((c1 == "PSOL") && (c2 == "CTRL")){
	folder <- "PSOL"
} else if ((c1 == "AD") && (c2 == "PSOL")){
	folder <- "PSOL"

cat("folder")
	for (i in 1:no.perm){ 
		tmpperm <- createRandomNetwork(raw.n1,raw.n2)
		write.table(tmpperm[[1]], paste0("/home/gaz/MAARS_p2/Scripts/compositionalCor_v2/SparCC/permutations/",folder, "/perm_",i), sep = "\t", quote = FALSE)

	} 

}