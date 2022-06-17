## author: Atul Deshpande
## email: adeshpande@jhu.edu
rm(list = ls())
setwd('.')
source("./R/preprocessing.R")
source("./R/getSpatialParameters.R")
source("./R/getInteractingGenes.R")
source('./R/find_genes_of_interest_nonparametric_fast.R')
## Specify data folder paths here: 
## Expected structure: parent_folder
##                              |____ VisiumDir (10x output format)
##                              |           |____ patient_id1
##                              |           |____ patient_id2
##                              |           .   .   .   .   .
##                              |____ CoGAPS_Analysis
##                                          |____ patient_id1
##                                          |           |____ cogapsFilePattern1
##                                          |           |____ cogapsFilePattern2
##                                          |             .   .   .   .   .
##                                          |____ patient_id2
##                                                      |____ cogapsFilePattern3
##                                                      |____ cogapsFilePattern4
##                                                        .   .   .   .   .
patient_id1 <- '1541'
visiumDir <- "./VisiumData/"
cogapsDir <- "./CoGAPS_Analysis/"
cogapsFilePattern <- ".rds"

## Set these parameters
# SpInMarkersMode: defaut mode is "residual". You can also set "DE" mode for Differential Expression mode.
SpInMarkersMode = "residual"  
# SpinMarkersRefPattern is the pattern whose "interaction" with every other pattern we want to study. If refPattern is not explicitly assigned, the code assumes Pattern_1 to be refPattern.
SpinMarkersRefPattern = "Pattern_8" 

## Loading data
pVisiumPath <- paste0(visiumDir,patient_id1)
fullMat <- load10XExpr(pVisiumPath)
good_gene_threshold <- 3
goodGenes <- rownames(fullMat)[apply(fullMat,1,function(x) sum(x>0)>=good_gene_threshold)]
fullMat <- fullMat[goodGenes,]
spCoords <- load10XCoords(pVisiumPath)
rownames(spCoords) <- spCoords$barcode
pCoGAPSPath <- paste0(cogapsDir,patient_id1)

## Matching formats
cogapsFilePath <- dir(pCoGAPSPath,cogapsFilePattern,full.names = T)
CoGAPS_Result <- readRDS(cogapsFilePath)
features <- intersect(rownames(fullMat),rownames(CoGAPS_Result@featureLoadings))
barcodes <- intersect(colnames(fullMat),rownames(CoGAPS_Result@sampleFactors))
fullMat <- fullMat[features,barcodes]
cgMat <- CoGAPS_Result@featureLoadings[features,] %*% t(CoGAPS_Result@sampleFactors[barcodes,])
spCoords <- spCoords[barcodes,]
spPatterns <- cbind(spCoords,CoGAPS_Result@sampleFactors[barcodes,])

## Running scripts
optParams <- getSpatialParameters(spPatterns)
SpInMarkers <- getInteractingGenes(data = fullMat, reconstruction = cgMat, spatialPatterns = spPatterns, refPattern = SpinMarkersRefPattern, mode = SpInMarkersMode)
SpInMarkers$optParams <- optParams
for (i in seq(1,length(SpInMarkers$interacting_genes)))
{
    filename <- paste0(gsub(x = cogapsFilePath,pattern = '.rds',replacement = '_SpInMarkers_'),names(SpInMarkers$interacting_genes[[i]])[1],".txt")
    write.table(rownames(SpInMarkers$interacting_genes[[i]]),file = filename,row.names = F,col.names = F, quote = F)
}  

filename <- gsub(x = cogapsFilePath,pattern = '.rds',replacement = '_patternHotspots.txt')
write.table((!is.na(SpInMarkers$hotspot.regions))*1,file = filename,row.names = T,col.names = T, quote = F)

saveRDS(SpInMarkers, file = gsub(x = cogapsFilePath,pattern = '.rds', replacement = '_SpInMarkers.rds'))
