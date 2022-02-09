matrixPath <- "HCC_Lu/matrix.mtx"
latentSpacePath <- "HCC_Lu/1541_15Patterns.rds"
bcs_mergedPath <- "HCC_Lu/1541_bcs_merged.rdata"

fullMat <- Matrix::readMM(matrixPath)
features <- read.table("~/FertigLab/Visium/HCC_Lu/features.tsv")
visiumCoGAPS <- readRDS(CoGAPSpath)
patternList <- colnames(visiumCoGAPS@sampleFactors)
coGAPSMat <- visiumCoGAPS@featureLoadings %*% t(visiumCoGAPS@sampleFactors)

load(bcs_mergedPath)
spotCoords <- bcs_merge[,c("barcode","sample","imagerow","imagecol")]  
colnames(spotCoords)[colnames(spotCoords)=="imagerow"]="y"
colnames(spotCoords)[colnames(spotCoords)=="imagecol"]="x"

spatialData<- assignCoGAPSPatternstoSpots(spotCoords, visiumCoGAPS)

assignPatternstoSpots <- function(coords, cgObj){
  spData <- cbind(coords, as.data.frame(cgObj@sampleFactors))  
  return(spData)
}








fullMat <- as.matrix(fullMat)
rownames(fullMat) <- features[,2]
features<- features[features[,2] %in% rownames(visiumCoGAPS@featureLoadings),2]
fullMat <- fullMat[rownames(visiumCoGAPS@featureLoadings),]
fullMat <- log2(1+fullMat)
residualMat <-  fullMat - coGAPSMat

  