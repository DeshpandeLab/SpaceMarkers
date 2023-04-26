## author: Atul Deshpande
## email: adeshpande@jhu.edu

find_pattern_hotspots <- function(spatialPatterns, params = NULL, patternName = "Pattern_1", outlier = "positive")
{
  if (is.null(params)){
    sigmaPair = 10
    kernelthreshold = 2
  }
  else{
    sigmaPair = params["sigmaOpt"]
    kernelthreshold = params["threshOpt"]
  }
  
  allwin = spatstat.geom::owin(xrange = c(min(spatialPatterns$x),max(spatialPatterns$x)), yrange = c(min(spatialPatterns$y),max(spatialPatterns$y)))
  patternVector = as.matrix(spatialPatterns[,patternName])
  X <-spatstat.geom::ppp(x = spatialPatterns$x, y = spatialPatterns$y, window = allwin, marks = patternVector)
  Kact1 = spatstat.explore::Smooth(X, at = "points", sigma = sigmaPair[1], leaveoneout = T)
  Karr1 = sapply(seq(1,100), function(i) {Xr = X; spatstat.geom::marks(Xr) = spatstat.geom::marks(X)[pracma::randperm(1:length(spatstat.geom::marks(X)))]; temp = spatstat.explore::Smooth(Xr, at = "points", sigma = sigmaPair[1], leaveoneout = T); return(temp)})
  Karr1 = unlist(Karr1)
  mKvec = mean(Karr1)
  sKvec = sd(Karr1)
  upthresh = mKvec+kernelthreshold*sKvec
  lothresh = mKvec-kernelthreshold*sKvec
  if (outlier == "positive"){
    ind1 = which(Kact1 > upthresh[1])
  }
  else if (outlier == "two.sided")
  {
    ind1 = which((Kact1 > upthresh)|(Kact1 < lothresh))
  }
  region = array(NA, length(Kact1))
  region[ind1] = patternName
  return(region)
}
#===================
#' getInteractingGenes
#' Calculate Interaction Regions and Associated Genes
#'
#' This function calculates statistically significant genes using a non-parametric Kruskal-Wallis test for genes in any one region of influence and a post hoc Dunn's test is used for analysis of genes between regions.
#'

#' @export
#'
#' @param data  original spatial data matrix.
#' @param reconstruction  reconstruction of the data matrix from latent spaces.
#' Required for "residual" mode.
#' @param spatialPatterns A data frame that contains the spatial coordinates for each cell type. The column names must include 'x' and 'y' as well as a set of numbered columns named 'Pattern_1.....N'.
#' @param optParams  a matrix with dimensions 2 X N, where N is the number of patterns with optimal parameters for outlier
#' detection calculated from function getSpatialParameters(). The first row contains the kernel width sigmaOpt for each
#' pattern, and the second row is the threshOpt (outlier threshold) for each pattern. Users can also input their
#' preferred param values.
#' The default value is NULL.
#' @param refPattern	 a character string that specifies the pattern whose "interaction" with every other pattern we want
#' to study. The default value is "Pattern_1".
#' @param mode  SpaceMarkers mode of operation. Possible values are "residual" (the default) or "DE".
#' @param minOverlap a number that specifies the minimum overlap between genes in two patterns to be considered for the statistical tests. The default is 50.
#' @param hotspotRegions a vector that specifies the patterns to compare to the 'refPattern'. The default is NULL which indicates that all patterns would be compared to the 'refPattern'.
#'
#'
#'
#' @return a list of data frames with information about the interacting genes of the refPattern and each latent feature pattern matrix (interacting_genes object). There is also a data frame with all of the regions of influence for any two of patterns (the hotspotRegions object).


getInteractingGenes <- function(data, reconstruction=NULL, spatialPatterns, optParams=NULL,
                                refPattern="Pattern_1", mode=c("residual","DE"), minOverlap = 50, hotspotRegions = NULL){
  
  if (mode=="residual"&&is.null(reconstruction)) stop("Reconstruction matrix not provided for residual mode.")
  if (mode=="residual"&&all(dim(data)!=dim(reconstruction))) stop("Original and reconstructed matrix do not have the same dimensions.")
  patternList <- colnames(spatialPatterns)[startsWith(colnames(spatialPatterns),"Pattern_")]
  if (is.null(optParams)){
    print("optParams not provided. Calculating optParams.")
    optParams <- getSpatialParameters(spatialPatterns)
  }
  else{
    print("Using user provided optParams.")
    if (any(colnames(optParams)!=patternList)) stop("Error: colnames of optParams must match Pattern names.")
    if (any(rownames(optParams)!=c("sigmaOpt","threshOpt"))) stop("Error: rownames of optParams must match c(\"sigmaOpt\",\"threshOpt\")")
    if(any(!is.numeric(optParams))) stop("Error: optParams must be numeric.")
  }
  if (is.null(hotspotRegions))
  {
    hotspotRegions = c()
    for (patternName in patternList)
    {
      hotspotRegions <- cbind(hotspotRegions,find_pattern_hotspots(spatialPatterns = spatialPatterns,
                                                                   patternName = patternName,
                                                                   params = optParams[,patternName],
                                                                   outlier = "positive") )
    }
    colnames(hotspotRegions) <- patternList
  }
  else if (!((refPattern %in% colnames(hotspotRegions)) && (nrow(hotspotRegions)== ncol(data))))
  {
    stop("Error: hotspotRegions does not have refPattern column or dimension does not match with data.")
  }
  else
    print("Using user provided hotspot regions.")
  
  
  interacting.genes <- list();
  data <- as.matrix(data)
  
  for (pattern in setdiff(patternList,refPattern)){
    region <- hotspotRegions[,refPattern];
    region <- ifelse(!is.na(region) & !is.na(hotspotRegions[,pattern]),"Interacting",ifelse(!is.na(region),region,hotspotRegions[,pattern]))
    region <- factor(region)
    if (length(levels(region))<3||any(table(region)<minOverlap)) #default 50
    {
      print(paste0(refPattern, " and ", pattern, " do not sufficiently interact. Skipping statistical test for genes."))
    } else {
      if (mode=="residual"){
        residualMat <- data - reconstruction
        interacting.genes <- c(interacting.genes,find_genes_of_interest_nonparametric_fast(testMat = residualMat, goodGenes = NULL, region=region))
      }
      else if (mode=="DE")
        interacting.genes <- c(interacting.genes,find_genes_of_interest_nonparametric_fast(testMat = data, goodGenes = NULL, region=region))
      else
        stop("Invalid mode.")
    }
  }
  
  interacting_genes <- lapply(interacting.genes, as.data.frame)
  for (i in seq(1,length(interacting_genes)))
    interacting_genes[[i]]$KW.p.adj <- as.numeric(interacting_genes[[i]]$KW.p.adj)
  
  for (i in seq(1,length(interacting_genes)))
  {
    if (all(dim(interacting_genes[[i]])>1))   {
      od <- order(interacting_genes[[i]]$KW.p.adj)
      interacting_genes[[i]] <- interacting_genes[[i]][od,]
    }
  }
  return(list(interacting_genes=interacting_genes,hotspotRegions=hotspotRegions))
}







