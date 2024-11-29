## author: Atul Deshpande
## email: adeshpande@jhu.edu

#===================
#' @title Identify hotspots of spatial pattern influence
#' @description This function calculates 'hotspots' which are regions of high
#'  spatial influence based on an outlier threshold from a null distribution. 
#' @export
#' @family getIntGenes
#' @param    spPatterns    A data frame that contains the spatial coordinates 
#' and metrics for spatial features (cell types/cell processes). The column 
#' names must include 'x' and 'y' as well as the spatially varying features.
#' @param    params    a named vector of the optimal sigma and threshold for a 
#' given spatial pattern. The names are should be 'sigmaOpt' and 'threshOpt'.
#' The default value is NULL.
#' @param patternName     a character string that specifies the pattern of 
#' interest
#' @param    outlier    a character string specifying whether to apply the 
#' outlier threshold to the kernel density distribution in a one-sided manner
#' (specify 'positive' the default) or in a two sided manner (specify 
#' 'two.sided').
#' @param nullSamples a numeric values specifying the number of spatial patterns
#' to randomly sample for a null distribution.
#' @param    includeSelf    a logic value specifying whether to consider the
#' spatial influence the pattern has on surrounding regions only (set to FALSE),
#' or whether to also consider the influence of the pattern itself (set to TRUE
#' , the default).
#' @param ... Arguments passed to methods
#' @return a character vector with the spatial feature name if the spatial 
#' influence exceeded the threshold for that spot/cell, and NA otherwise
#' @examples 
#' library(SpaceMarkers)
#' #Visium data links
#' urls <- read.csv(system.file("extdata","visium_data.txt",
#'                            package="SpaceMarkers",mustWork = TRUE))
#' sp_url <- urls[["visium_url"]][2]
#' #Remove present Directories if any
#' unlink(basename(sp_url))
#' unlink("spatial", recursive = TRUE)
#' #Obtaining CoGAPS Patterns i.e Spatial Features
#' cogaps_result <- readRDS(system.file("extdata","CoGAPS_result.rds",
#'                                    package="SpaceMarkers",mustWork = TRUE))
#' spFeatures <- slot(cogaps_result,"sampleFactors")
#' #Obtaining Spatial Coordinates
#' download.file(sp_url, basename(sp_url), mode = "wb")
#' untar(basename(sp_url))
#' spCoords <- load10XCoords(visiumDir = ".", version = "1.0")
#' rownames(spCoords) <- spCoords$barcode
#' #Match Dimensions
#' barcodes <- intersect(rownames(spFeatures),spCoords$barcode)
#' spCoords <- spCoords[barcodes,]
#' spFeatures <- spFeatures[barcodes,]
#' spPatterns <- cbind(spCoords,spFeatures[barcodes,])
#' spPatterns<-spPatterns[c("barcode","y","x","Pattern_1","Pattern_5")]
#' data("optParams")
#' hotspots <- findPatternHotspots(
#' spPatterns = spPatterns,
#' patternName = "Pattern_1",
#' params = optParams[,"Pattern_1"],
#' outlier = "positive",nullSamples = 100,includeSelf = TRUE)
#' #Remove present Directories if any
#' unlink(basename(sp_url))
#' unlink("spatial", recursive = TRUE)
#' 

findPatternHotspots <- function(
        spPatterns, params = NULL, patternName = "Pattern_1",
        outlier = "positive",
    nullSamples = 100, includeSelf = TRUE,...){
    if (is.null(params)){
        sigmaPair <- 10
        kernelthreshold <- 4
    } else {
        sigmaPair <- params["sigmaOpt"]
        kernelthreshold <- params["threshOpt"]
    }
    
    allwin <- spatstat.geom::owin(
        xrange = c(min(spPatterns$x),max(spPatterns$x)),yrange =c(
            min(spPatterns$y),max(spPatterns$y)))
    patternVector <- as.matrix(spPatterns[,patternName])
    X <-spatstat.geom::ppp(
        x=spPatterns$x,y = spPatterns$y, window = allwin,marks = patternVector)
    Kact1 <- spatstat.explore::Smooth(
        X, at = "points", sigma = sigmaPair[1], ...)
    if (includeSelf == TRUE){
      Kact1 <- spPatterns[,patternName] + Kact1
    } else {
      Kact1 <- Kact1
    }
    Karr1 <- vapply(seq(1,nullSamples),function(i){
        Xr<-X;
        spatstat.geom::marks(Xr) <- sample(spatstat.geom::marks(X));
        temp <- spatstat.explore::Smooth(
                    Xr, at="points", sigma = sigmaPair[1], ...); 
                return(temp)}, numeric(length(Kact1)))
    Karr1 <- unlist(Karr1)
    mKvec <- mean(Karr1)
    sKvec <- sd(Karr1)
    upthresh <- mKvec+kernelthreshold*sKvec
    lothresh <- mKvec-kernelthreshold*sKvec
    if (outlier == "positive"){
        ind1 <- which(Kact1 > upthresh[1])
    }
    else if (outlier == "two.sided")
    {
        ind1 <- which((Kact1 > upthresh)|(Kact1 < lothresh))
    }
    region <- array(NA, length(Kact1))
    region[ind1] <- patternName
    return(region)
}

gettestMat <- function(data,reconstruction,mode){
    if ("DE" %in% mode)
        testMat <- data
    else if ("residual" %in% mode)
        if (is.null(reconstruction))
            stop("Reconstruction matrix not provided for residual mode.")
        else if (all(dim(data)!=dim(reconstruction)))
            stop("Original and reconstruction have mismatched dimensions.")
        else if (any(dimnames(data)[[1]]!=dimnames(reconstruction)[[1]]))
            stop("Original and reconstruction have mismatched row names.")
        else if (any(dimnames(data)[[2]]!=dimnames(reconstruction)[[2]]))
            stop("Original and reconstruction have mismatched column names.")
        else 
            testMat <- data - reconstruction
    else
        stop("Invalid mode. mode must be either DE or residual.")
    return(testMat)
}

getSpaceMarkersMetric <- function(interacting.genes){
    interacting_genes <- lapply(interacting.genes, as.data.frame)
    if (length(interacting_genes) == 0)
    {
        message("No interacting genes found. Returning 
            result with only hotspots.")
        return(interacting_genes)
    }
        
    for (i in seq(1,length(interacting_genes))){
        interacting_genes[[i]]$KW.p.adj <- as.numeric(
            interacting_genes[[i]]$KW.p.adj)
        interacting_genes[[i]]$Dunn.pval_1_Int.adj <- as.numeric(
            interacting_genes[[i]]$Dunn.pval_1_Int.adj)
        interacting_genes[[i]]$Dunn.pval_2_Int.adj <- as.numeric(
            interacting_genes[[i]]$Dunn.pval_2_Int.adj)
    }
    for (i in seq(1,length(interacting_genes)))
    {
        if (all(dim(interacting_genes[[i]])>1))   {
            Zsign <- (2*(-1+((interacting_genes[[i]]$Dunn.zP1_Int<0)&
                        (interacting_genes[[i]]$Dunn.zP2_Int<0))*1)+1)
            Zmag <- abs((interacting_genes[[i]]$Dunn.zP1_Int)*
                    (interacting_genes[[i]]$Dunn.zP2_Int)/
                    (pmax(abs(interacting_genes[[i]]$Dunn.zP2_P1),1)))
                interacting_genes[[i]]$SpaceMarkersMetric <- Zsign*log2(Zmag+1)
                od <- order(
                    interacting_genes[[i]]$SpaceMarkersMetric,decreasing=TRUE)
                interacting_genes[[i]] <- interacting_genes[[i]][od,]
                interacting_genes[[i]] <- interacting_genes[[i]][!is.na(
                    interacting_genes[[i]]$SpaceMarkersMetric),]
        }
    }
    return(interacting_genes)
}


#===================
#' @title Calculate Interaction Regions and Associated Genes
#' @description This function calculates statistically significant genes using a 
#' non-parametric Kruskal-Wallis test for genes in any one region 
#' of influence and a post hoc Dunn's test is used for analysis of 
#' genes between regions.
#' @export
#' @family getIntGenes
#' @param    data    original spatial data matrix.
#' @param    reconstruction    reconstruction of the data matrix from latent 
#' spaces. Required for "residual" mode.
#' @param    spPatterns    A data frame that contains the spatial coordinates 
#' and metrics for spatial features (cell types/cell processes). The column 
#' names must include 'x' and 'y' as well as the spatially varying features.
#' @param    optParams    a matrix with dimensions 2 X N, where N is the number
#' of spatial patterns with optimal parameters. The first row contains the 
#' kernel width 'sigmaOpt' for each pattern, and the second row is the
#' threshOpt (outlier threshold) for each pattern. Users can also input their
#' preferred param values. The default value is NULL.
#' @param refPattern     a character string that specifies the pattern whose
#'  "interaction" with every other pattern we want
#' to study. The default value is "Pattern_1".
#' @param    mode    SpaceMarkers mode of operation. Possible values are
#'  "DE" (the default) or "residual".
#' @param    minOverlap    a number that specifies the minimum overlap between
#'  genes in two patterns to be considered for the statistical tests.
#'  The default is 50.
#' @param    hotspots    a vector that specifies the patterns to compare
#'  to the 'refPattern'. The default is NULL which indicates that all patterns
#'   would be compared to the 'refPattern'.
#' @param analysis a character string that specifies the type of downstream 
#' analysis to be performed. Possible values are "enrichment" (default) 
#' and "overlap". In enrichment mode, all genes are returned, ranked by 
#' the SpaceMarkers metric. In overlap mode, only the genes which are 
#' significantly overexpressed in the interaction region are returned.
#' @param ... Arguments passed to methods
#' @return a list of data frames with information about the interacting genes
#'  of the refPattern and each latent feature pattern matrix 
#'  (interacting_genes object). There is also a data frame with all of the
#'   regions of influence for any two of patterns (the hotspots object).
#' @examples 
#' library(SpaceMarkers)
#' #Visium data links
#' urls <- read.csv(system.file("extdata","visium_data.txt",
#' package="SpaceMarkers",mustWork = TRUE))
#' counts_url <- urls[["visium_url"]][1]
#' sp_url <- urls[["visium_url"]][2]
#' #Remove present Directories if any
#' unlink(basename(sp_url))
#' unlink("spatial", recursive = TRUE)
#' files <- list.files(".")[grepl(basename(counts_url),list.files("."))]
#' unlink(files)
#' download.file(counts_url,basename(counts_url), mode = "wb")
#' counts_matrix<-load10XExpr(visiumDir=".",h5filename = basename(counts_url))
#' #Obtaining CoGAPS Patterns
#' cogaps_result <- readRDS(system.file("extdata","CoGAPS_result.rds",
#' package="SpaceMarkers",mustWork = TRUE))
#' features <- intersect(rownames(counts_matrix),rownames(
#'     slot(cogaps_result,"featureLoadings")))
#' barcodes <- intersect(colnames(counts_matrix),rownames(
#'     slot(cogaps_result,"sampleFactors")))
#' counts_matrix <- counts_matrix[features,barcodes]
#' cogaps_matrix <- slot(cogaps_result,"featureLoadings")[features,]%*%
#'     t(slot(cogaps_result,"sampleFactors")[barcodes,])
#' #Obtaining Spatial Coordinates
#' download.file(sp_url, basename(sp_url), mode = "wb")
#' untar(basename(sp_url))
#' spCoords <- load10XCoords(visiumDir = ".", version = "1.0")
#' rownames(spCoords) <- spCoords$barcode
#' spCoords <- spCoords[barcodes,]
#' spPatterns <- cbind(spCoords,slot(cogaps_result,
#' "sampleFactors")[barcodes,])
#' data("curated_genes")
#' spPatterns<-spPatterns[c("barcode","y","x","Pattern_1","Pattern_5")]
#' counts_matrix <- counts_matrix[curated_genes,]
#' cogaps_matrix <- cogaps_matrix[curated_genes, ]
#' data("optParams")
#' SpaceMarkersMode <- "DE"
#' ref_Pattern <- "Pattern_1"
#' SpaceMarkers_test <- getInteractingGenes(
#'     data=counts_matrix,reconstruction=NULL,
#'     optParams = optParams,
#'     spPatterns = spPatterns,
#'     refPattern = "Pattern_1",
#'     mode="DE",analysis="overlap")
#' #Remove present Directories if any
#' unlink(basename(sp_url))
#' unlink("spatial", recursive = TRUE)
#' files <- list.files(".")[grepl(basename(counts_url),list.files("."))]
#' unlink(files)
#' 

getInteractingGenes <- function(data, spPatterns, refPattern="Pattern_1",
                                mode=c("DE","residual"), optParams=NULL,
                                reconstruction=NULL, hotspots=NULL,
                                analysis=c("enrichment","overlap"),
                                minOverlap=50,...) {

    testMat <- gettestMat(data,reconstruction,mode)
    pattList<- setdiff(colnames(spPatterns),c("barcode","x","y"))
    if (is.null(optParams)){
        message("optParams not provided. Calculating optParams.")
        optParams <- getSpatialParameters(spPatterns)
    } else {
        message("Using user provided optParams.")
        if (any(colnames(optParams)!=pattList)) stop(
            "Colnames of optParams must match Pattern names.")
        if (any(rownames(optParams)!=c("sigmaOpt","threshOpt"))) stop(
            "Rownames of optParams must match c(\"sigmaOpt\",\"threshOpt\")")
        if(any(!is.numeric(optParams))) stop("optParams must be numeric.")
    }
    if (is.null(hotspots)) {
        hotspots <- c()
        for (patternName in pattList) {
            hotspots <- cbind(
                hotspots,findPatternHotspots(
                    spPatterns=spPatterns,patternName=patternName,
                    params=optParams[,patternName],outlier = "positive") )
        }
        colnames(hotspots) <- pattList
    }
    else if (!((refPattern %in% colnames(hotspots)) && (nrow(
        hotspots)==ncol(data))))
        message("hotspots does not have refPattern column or dimension 
    does not match with data.")
    else 
        message("Using user provided hotspot regions.")
    interacting_genes <- list()
    for (pattern in setdiff(pattList,refPattern)){
        region <- hotspots[,refPattern]
        region <- ifelse(!is.na(region) & !is.na(
            hotspots[,pattern]),"Interacting",ifelse(
                !is.na(region),region,hotspots[,pattern]))
        region <- factor(region)
        if (length(levels(region))<3||any(table(region)<minOverlap))#def. 50
            message(refPattern," and ",pattern," do not sufficiently interact.
                Skipping statistical test for genes.")
        else {
            message("Calculating genes of interest for ",refPattern," and ",
                pattern)
            interacting_genes<-c(interacting_genes,find_genes_of_interest(
                    testMat = testMat, region=region,...))
        }
    }
    interacting_genes <- getSpaceMarkersMetric(interacting_genes)
    return(list(interacting_genes=interacting_genes,hotspots=hotspots))
}

#' ================
#' @title getPairwiseInteractingGenes
#' @description Performs pairwise analysis to find genes associated with 
#' spatial interaction between pairs of spatially varying patterns. 
#' @export
#' @family getIntGenes
#' @inheritParams getInteractingGenes
#' @param patternPairs A matrix of pattern pairs to be analyzed. Default is
#' @param workers (optional) Number of workers to 
#' be used for parallel processing.
#' @return a list of data frames for each pattern with 1) names of the patterns 
#' (patterns object) 2) data frame with the hotspots of influence for the
#'  two patterns (the hotspots object). 3) data frame with the genes 
#'  associated with the interaction between the two patterns
#'   (interacting genes object, empty if insufficient interaction). 
#' @examples
#' library(SpaceMarkers)
#' #Visium data links
#' urls <- read.csv(system.file("extdata","visium_data.txt",
#' package="SpaceMarkers",mustWork = TRUE))
#' counts_url <- urls[["visium_url"]][1]
#' sp_url <- urls[["visium_url"]][2]
#' #Remove present Directories if any
#' unlink(basename(sp_url))
#' unlink("spatial", recursive = TRUE)
#' files <- list.files(".")[grepl(basename(counts_url),list.files("."))]
#' unlink(files)
#' download.file(counts_url,basename(counts_url), mode = "wb")
#' counts_matrix<-load10XExpr(visiumDir=".",
#' h5filename = basename(counts_url))
#' #Obtaining CoGAPS Patterns
#' cogaps_result <- readRDS(system.file("extdata","CoGAPS_result.rds",
#' package="SpaceMarkers",mustWork = TRUE))
#' features <- intersect(rownames(counts_matrix),rownames(
#'     slot(cogaps_result,"featureLoadings")))
#' barcodes <- intersect(colnames(counts_matrix),rownames(
#'     slot(cogaps_result,"sampleFactors")))
#' counts_matrix <- counts_matrix[features,barcodes]
#' cogaps_matrix <- slot(cogaps_result,"featureLoadings")[features,]%*%
#'     t(slot(cogaps_result,"sampleFactors")[barcodes,])
#' #Obtaining Spatial Coordinates
#' download.file(sp_url, basename(sp_url), mode = "wb")
#' untar(basename(sp_url))
#' spCoords <- load10XCoords(visiumDir = ".", version = "1.0")
#' rownames(spCoords) <- spCoords$barcode
#' spCoords <- spCoords[barcodes,]
#' spPatterns <- cbind(spCoords,
#' slot(cogaps_result,"sampleFactors")[barcodes,])
#' data("curated_genes")
#' spPatterns<-spPatterns[c("barcode","y","x","Pattern_1",
#' "Pattern_3","Pattern_5")]
#' counts_matrix <- counts_matrix[curated_genes,]
#' cogaps_matrix <- cogaps_matrix[curated_genes, ]
#' optParams <- matrix(c(6, 2, 6, 2, 6, 2), nrow = 2)
#' rownames(optParams) <- c("sigmaOpt","threshOpt")
#' colnames(optParams) <- c("Pattern_1","Pattern_3","Pattern_5")
#' SpaceMarkersMode <- "DE"
#' patternPairs <- matrix(c("Pattern_1", "Pattern_1", 
#'            "Pattern_3", "Pattern_5"), nrow=2)
#' SpaceMarkers_test <- getPairwiseInteractingGenes(
#'     data=counts_matrix,reconstruction=NULL,
#'     optParams = optParams,
#'     spPatterns = spPatterns,
#'     mode="DE",analysis="enrichment", patternPairs=patternPairs)
#' #Remove present Directories if any
#' unlink(basename(sp_url))
#' unlink("spatial", recursive = TRUE)
#' files <- list.files(".")[grepl(basename(counts_url),list.files("."))]
#' unlink(files)
#' 

getPairwiseInteractingGenes <- function(data, spPatterns, 
                                    mode=c("DE","residual"), 
                                    optParams=NULL, reconstruction=NULL, 
                                    hotspots=NULL, minOverlap=50,
                                    analysis=c("enrichment","overlap"),
                                    patternPairs=NULL ,..., workers=NULL) {
    # save params in a list for easy passing
    argsParams <- list(spPatterns=spPatterns, mode=mode, 
                        optParams=optParams, hotspots=hotspots, 
                        analysis=analysis, minOverlap=minOverlap)
    if (!exists("patternList"))
        patternList <- setdiff(colnames(spPatterns), c("x","y","barcode"))
    else {
        if (!all(patternList %in% colnames(spPatterns)))
            stop("patternList contains patterns not present in spPatterns.")
    }
    
    # check patternPairs if provided, if not generate them
    patternPairs <- getPatternPairs(patternList=patternList, 
                                        patternPairs=patternPairs)

    input_list <- apply(patternPairs,1, pairArgsList, argsParams)
    
    # Check if BiocParallel is installed
    use_biocparallel <- requireNamespace("BiocParallel", quietly = TRUE) &&
                            (length(input_list)>1)
    
    if (use_biocparallel) {
        
        # Determine the number of workers to use
        if (is.null(workers)) {
            workers <- BiocParallel::multicoreWorkers()  # Use default 
        } else if (workers <= 0) {
            stop("Invalid number of workers. #
                Please provide a positive integer.")
        }
        
        # Register parallel backend
        bpparam <- BiocParallel::MulticoreParam(workers = workers)
        
        # Run the loop with BiocParallel
        results <- BiocParallel::bplapply(input_list, function(args) {
            # call getInteractingGenes for a pair 
            do.call(getInteractingGenes, c(list(data = data), 
                list(reconstruction=reconstruction), args, ...))
        }, BPPARAM = bpparam)  # Use the specified parallel backend
    } else {
        # Use a regular for loop with input_list
        results <- vector("list", length = length(input_list))
        
        for (ii in seq(1,length(input_list))){
            results[[ii]] <- do.call(getInteractingGenes, c(list(data = data), 
                list(reconstruction=reconstruction), input_list[[ii]],...))
        }
    }
    for (ii in seq(1,length(results)))
        results[[ii]]$patterns <- patternPairs[ii,]
    names(results) <- apply(patternPairs,1,paste,collapse="_")
    return(results)
}

# Create a list of arguments for passing two patterns to getInteractingGenes
pairArgsList <- function(patternPair, argsParams) {
    patcols <- c("x","y", "barcode",patternPair)
    argList <- list(
        spPatterns = argsParams$spPatterns[, patcols],
        hotspots = argsParams$hotspots[, patcols],
        optParams = argsParams$optParams[, patternPair],
        refPattern = patternPair[1],
        mode = argsParams$mode,
        analysis = argsParams$analysis,
        minOverlap = argsParams$minOverlap)
    return(argList)
}

getPatternPairs <- function(patternList, patternPairs) {
    if (is.null(patternPairs)){
        message("patternPairs not provided. Calculating all 
            possible pairs.")
        patternPairs <- t(utils::combn(patternList,2,simplify = TRUE))
    } else {
        if (is.list(patternPairs)) {
            if (!all(vapply(patternPairs,length, FUN.VALUE = 1)==2)) 
                stop("Each item in the patternPairs list 
                    must be a vector of 2 patterns.")
            patternPairs <- do.call(rbind, patternPairs)
        } else if (is.matrix(patternPairs)) {
            if (ncol(patternPairs) != 2) 
                stop("patternPairs matrix must have 2 columns.")
        } else if(is.character(patternPairs)) {
            if (length(patternPairs)!=2){
                stop("A single pair be a character vector of length 2.")
            } else
                patternPairs <- matrix(patternPairs, nrow=1)
        }
        else
            stop("patternPairs must either be a 2-column matrix 
                or a list of vectors.")
        mess <- paste0(setdiff(patternPairs, patternList)," ")
        if (!all(patternPairs %in% patternList)) 
            stop("Following are not pattern names: ",mess )
    }
    return(patternPairs)
}
