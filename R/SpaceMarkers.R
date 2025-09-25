#' @title Main dispathcher for SpaceMarkers
#' @description This function is the main dispatcher for SpaceMarkers. It takes a features object/path and a cell by gene matrix and returns IMScores.
#' @param features A path to a features file.
#' @param data A path to a Visium directory.
#' @param supervised Type of analysis to run, supervised (TRUE) or unsupervised (FALSE).
#' 
SpaceMarkers <- function(features, counts, supervised = FALSE, ...) {
  if (supervised)  {
    return(.supervised_SpaceMarkers(features, data, ...))
  } else if (!supervised) {
    return(.unsupervised_SpaceMarkers(features, data, ...))
  } else {
    stop("Check input: supervised must be TRUE or FALSE.")
  }
}

.unsupervised_SpaceMarkers <- function(features, data, cpus,...) {
    #load spatial coords from tissue positions, deconvolved patterns, and expression
    message("Loading data...")
    coords <- load10XCoords(data)
    dataMatrix <- load10XExpr(data)
    features <- getSpatialFeatures(features)

    #add spatial coordinates to deconvolved data, only use barcodes present in data
    message("Preparing data...")
    spPatterns <- merge(coords, features, by.x = "barcode", by.y = "row.names")
    spPatterns <- spPatterns[which(spPatterns[,"barcode"] %in% colnames(dataMatrix)),]

    #remove genes with low expression, only barcodes present in spatial data
    message("Filtering data...")
    keepGenes <- which(apply(dataMatrix, 1, sum) > 10)
    keepBarcodes <- which(colnames(dataMatrix) %in% spPatterns[,"barcode"])
    dataMatrix <- dataMatrix[keepGenes, keepBarcodes]

    #compute optimal parameters for spatial patterns
    message("Computing optimal parameters...")
    optParams <- getSpatialParameters(spPatterns, visiumDir=data)

    #find hotspots in spatial patterns
    message("Finding hotspots...")
    hotspots <- findAllHotspots(spPatterns)

    #find regions of overlapping spatial patterns
    message("Finding overlapping regions...")
    overlaps <- getOverlapScores(hotspots)

    #find genes that are differentially expressed in spatial patterns
    message("Finding interacting genes...")
    spaceMarkers <- getPairwiseInteractingGenes(data = dataMatrix,
                                                optParams = optParams,
                                                spPatterns = spPatterns,
                                                hotspots = hotspots,
                                                mode = "DE",
                                                analysis = "enrichment",
                                                minOverlap = 10,
					                            workers=cpus)

    #save SpaceMarkers Interaction Scores
    message("Calculating Interaction Scores...")
    IMScores <- getIMScores(spaceMarkers)

    return(IMScores)
}

.supervised_SpaceMarkers <- function(features, data, ...) {
    # limit the analysis to ligand-receptor genes
    use.ligand.receptor.genes <- TRUE

    coords <- load10XCoords(data)
    rownames(coords) <- coords[["barcode"]]

    spPatterns <- getSpatialFeatures(features)
    barcodes <- intersect(rownames(coords), rownames(spPatterns))
    spPatterns <- cbind(coords[barcodes,], spPatterns[barcodes, ])
    optParams <- getSpatialParameters(spPatterns, visiumDir = data)
    sigmaPair <- optParams[,1]
    patnames <- setdiff(colnames(spPatterns), c("x", "y", "barcode"))

    # Calculate hotspots and influence for each pattern
    print("Calculating hotspots and influence for each pattern...")
    patthresholds <- calcAllThresholds(spPatterns, minvals = 0.1, maxvals = 0.8)
    patHotspots <- findAllHotspots.value(spPatterns, threshold = patthresholds)

    spInfluence <- calcInfluence(spPatterns, optParams)
    infthresholds <- calcAllThresholds(spInfluence, minvals = 0.01, maxvals = 0.5)
    infHotspots <- findAllHotspots.value(spInfluence, threshold = infthresholds)

    #create a table of pattern pairs
    patternPairs <- t(combn(patnames, 2))

    data <- load10XExpr(data)
    data <- data[, barcodes]

    data(lrdf)
    lrpairs <- lrdf[["interaction"]][,c("ligand.symbol","receptor.symbol")]
    ligands <- sapply(lrpairs[["ligand.symbol"]],function(i) strsplit(i,split=", "))
    receptors <- sapply(lrpairs[["receptor.symbol"]],function(i) strsplit(i,split=", "))

    lrgenes <- union(unlist(ligands), unlist(receptors))

    if (use.ligand.receptor.genes) {
    print("Limiting data to ligand-receptor genes...")
        # Filter data to include only ligand-receptor genes
        data <- data[rownames(data) %in% lrgenes, ]
    } else{
        # Use "good" genes
        goodgenes <- which(apply(data, 1, sum) > goodgeneThreshold )
        data <- data[goodgenes,]
    }

    # Calculate interaction scores for all pattern pairs
    IMscores <- calcAllIMscores.HD(data = data,
                                   patHotspots = patHotspots,
                                   infHotspots = infHotspots,
                                   patternPairs = patternPairs,
                                   avoid_confounders = TRUE)



    return(IMscores)
}
