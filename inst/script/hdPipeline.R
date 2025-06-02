
data_dir <- "~/0_Projects/02_BTC/DPT/HDsample/binned_outputs/square_016um/"
coords <- load10XCoords(data_dir)
rownames(coords) <- coords$barcode

spPatterns <- getSpatialFeatures("~/Downloads/rctd_cell_types-2.csv")
barcodes <- intersect(rownames(coords), rownames(spPatterns))
spPatterns <- cbind(coords[barcodes,],spPatterns[barcodes,])
optParams <- getSpatialParameters(spPatterns,visiumDir=data_dir)
sigmaPair <- optParams[,1]
patnames <- setdiff(colnames(spPatterns),c("x","y","barcode"))


patthresholds <- calcAllThresholds(spPatterns)
patHotspots <- findAllHotspots.value(spPatterns, threshold=patthresholds)

spInfluence <- calcInfluence(spPatterns,optParams)
infthresholds <- calcAllThresholds(spInfluence)
infHotspots <- findAllHotspots.value(spInfluence, threshold=infthresholds)

#create a table of pattern pairs
patternPairs <- t(combn(patnames,2))

data <- load10XExpr(data_dir)
data <- data[,barcodes]


lrlist <- data(lrlist)
lrgenes <- Reduce(union,lrlist)
lrgenes <- intersect(lrgenes,rownames(data))
data <- data[lrgenes,]

# Calculate interaction scores for all pattern pairs
IMscores <- calcAllIMscores.HD(data, patHotspots, infHotspots, patternpairs)

LR_df <- readRDS("SpaceMarkers_Test_Inputs/LR_df.rds")
LR_df <- LR_df %>% dplyr::select(ligand.symbol,receptor.symbol)
LR_df$LR <- paste0(LR_df$ligand.symbol,", ",LR_df$receptor.symbol)
LR_list <- strsplit(LR_df$LR,", ")
LRscores <- getGeneSetScore(IMscores,genes = LR_list)
