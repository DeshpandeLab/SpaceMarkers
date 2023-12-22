## author: Atul Deshpande
## email: adeshpande@jhu.edu

find_pattern_hotspots <- function(
        spPatterns, params = NULL,
        patternName = "Pattern_1",
        outlier = "positive",...){
    if (is.null(params)){
        sigmaPair <- 10
        kernelthreshold <- 2
    } else {
        sigmaPair <- params["sigmaOpt"]
        kernelthreshold <- params["threshOpt"]
    }
    
    allwin <- spatstat.geom::owin(
        xrange = c(min(spPatterns$x),max(spPatterns$x)),yrange =c(
            min(spPatterns$y),max(spPatterns$y)))
    patternVector <- as.matrix(spPatterns[,patternName])
    X <-spatstat.geom::ppp(
        x=spPatterns$x,y = spPatterns$y,window = allwin,marks = patternVector)
    Kact1 <- spatstat.explore::Smooth(
        X, at = "points", sigma = sigmaPair[1],leaveoneout = TRUE)
    Karr1 <- vapply(seq(1,100),function(i){
        Xr<-X;
        spatstat.geom::marks(Xr) <- spatstat.geom::marks(X)[pracma::randperm(
            seq_len(length(spatstat.geom::marks(
                X))))];temp <- spatstat.explore::Smooth(
                    Xr, at="points", sigma = sigmaPair[1],leaveoneout=TRUE); 
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
    if (mode=="residual"){
        if (is.null(reconstruction))
            stop("Reconstruction matrix not provided for residual mode.")
        else if (all(dim(data)!=dim(reconstruction)))
            stop("Original and reconstruction have mismatched dimensions.")
        else 
            testMat <- data - reconstruction
    }
    else if (mode=="DE") testMat <- data
    else
        stop("Invalid mode.")
    return(testMat)
}

getSpaceMarkersMetric <- function(interacting.genes){
    interacting_genes <- lapply(interacting.genes, as.data.frame)
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
            interacting_genes[[i]]$SpaceMarkersMetric <-  abs(
                interacting_genes[[i]]$Dunn.zP1_Int) * abs(
                    interacting_genes[[i]]$Dunn.zP2_Int) * (2*(-1)^(pmin(
                        interacting_genes[[i]]$Dunn.zP1_Int,
                        interacting_genes[[i]]$Dunn.zP2_Int)>=0)-1)
            od<-order(interacting_genes[[i]]$SpaceMarkersMetric,decreasing=TRUE)
            interacting_genes[[i]] <- interacting_genes[[i]][od,]
            interacting_genes[[i]] <- interacting_genes[[i]][!is.na(
                interacting_genes[[i]]$SpaceMarkersMetric),]
        }
    }
        return(interacting_genes)
}


#===================
#' getInteractingGenes
#' Calculate Interaction Regions and Associated Genes
#' This function calculates statistically significant genes using a 
#' non-parametric Kruskal-Wallis test for genes in any one region of influence
#'  and a post hoc Dunn's test is used for analysis of genes between regions.
#' @export
#'
#' @param    data    original spatial data matrix.
#' @param    reconstruction    reconstruction of the data matrix from latent 
#' spaces. Required for "residual" mode.
#' @param    spPatterns    A data frame that contains the spatial coordinates 
#' for each cell type. The column names must include 'x' and 'y' as well as a
#'  set of numbered columns named 'Pattern_1.....N'.
#' @param    optParams    a matrix with dimensions 2 X N, where N is the number
#'  of patterns with optimal parameters for outlier
#' detection calculated from function getSpatialParameters(). The first row
#'  contains the kernel width sigmaOpt for each
#' pattern, and the second row is the threshOpt (outlier threshold) for each
#'  pattern. Users can also input their
#' preferred param values.
#' The default value is NULL.
#' @param refPattern     a character string that specifies the pattern whose
#'  "interaction" with every other pattern we want
#' to study. The default value is "Pattern_1".
#' @param    mode    SpaceMarkers mode of operation. Possible values are
#'  "residual" (the default) or "DE".
#' @param    minOverlap    a number that specifies the minimum overlap between
#'  genes in two patterns to be considered for the statistical tests.
#'  The default is 50.
#' @param    hotspots    a vector that specifies the patterns to compare
#'  to the 'refPattern'. The default is NULL which indicates that all patterns
#'   would be compared to the 'refPattern'.
#' @param ... Arguments passed to methods
#' @return a list of data frames with information about the interacting genes
#'  of the refPattern and each latent feature pattern matrix 
#'  (interacting_genes object). There is also a data frame with all of the
#'   regions of influence for any two of patterns (the hotspots object).
#' @examples 
#' library(SpaceMarkers)
#' #Visium data links
#' urls <- read.csv("inst/extdata/visium_data.txt")
#' counts_url <- urls[1,1]
#' sp_url <- urls[2,1]
#' #Remove present Directories if any
#' unlink(basename(sp_url))
#' unlink("spatial", recursive = TRUE)
#' files <- list.files(".")[grepl(basename(counts_url),list.files("."))]
#' unlink(files)
#' download.file(counts_url,basename(counts_url))
#' counts_matrix<-load10XExpr(visiumDir=".",h5filename = basename(counts_url))
#' good_gene_threshold <- 1500
#' goodGenes <- rownames(counts_matrix)[apply(counts_matrix,1,function(x)
#' sum(x>0)>=good_gene_threshold)]
#' print(length(goodGenes))
#' counts_matrix <- counts_matrix[goodGenes,]
#' #Latent Feature Space
#' cogaps_factors <- read.table("inst/extdata/cogaps_factors.csv")
#' cogaps_features <- read.table("inst/extdata/cogaps_features.csv")
#' features <- intersect(rownames(counts_matrix),rownames(cogaps_features))
#' barcodes <- intersect(colnames(counts_matrix),rownames(cogaps_factors))
#' counts_matrix <- counts_matrix[features,barcodes]
#' cogaps_matrix<-data.matrix(cogaps_features[features,])%*%
#' t(data.matrix(cogaps_factors[barcodes,]))
#' #Get Breast Cancer data
#' optParams <- data.matrix(data.frame("Pattern_1" = c(7,2.1),
#' "Pattern_5"= c(6.4,1.2)))
#' rownames(optParams) <- c("sigmaOpt","threshOpt")
#' #Spatial Coordinates
#' download.file(sp_url, basename(sp_url))
#' untar(basename(sp_url))
#' spCoords <- load10XCoords(visiumDir = ".")
#' rownames(spCoords) <- spCoords$barcode
#' spCoords <- spCoords[barcodes,]
#' spPatterns <- cbind(spCoords,cogaps_factors[barcodes,])
#' spPatterns<-spPatterns[c("barcode","y","x","Pattern_1","Pattern_5")]
#' 
#' SpaceMarkersMode <- "DE"
#' ref_Pattern <- "Pattern_1"
#' SpaceMarkers_test <- getInteractingGenes(data = counts_matrix,
#' reconstruction = NULL,
#' optParams = optParams,
#' spPatterns = spPatterns,
#' refPattern = ref_Pattern,
#' mode = SpaceMarkersMode,
#' analysis = "enrichment")
#' 

getInteractingGenes <- function(data,spPatterns,refPattern="Pattern_1",
                                mode=c("residual","DE"),optParams=NULL,
                                reconstruction=NULL,hotspots=NULL,
                                minOverlap=50,...) {
    testMat <- gettestMat(data,reconstruction,mode)
    pattList<-colnames(spPatterns)[startsWith(colnames(spPatterns),"Pattern_")]
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
                hotspots,find_pattern_hotspots(
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
    interacting_genes <- list(); data <- as.matrix(data)
    for (pattern in setdiff(pattList,refPattern)){
        region <- hotspots[,refPattern]; 
        region <- ifelse(!is.na(region) & !is.na(
            hotspots[,pattern]),"Interacting",ifelse(
                !is.na(region),region,hotspots[,pattern]))
        region <- factor(region)
        if (length(levels(region))<3||any(table(region)<minOverlap))#default 50
            message(refPattern,"and",pattern,"do not sufficiently interact.
                Skipping statistical test for genes.")
        else
            interacting_genes<-c(interacting_genes,find_genes_of_interest(
                    testMat = testMat, region=region,...))
    }
    interacting_genes <- getSpaceMarkersMetric(interacting_genes)
    return(list(interacting_genes=interacting_genes,hotspots=hotspots))
}