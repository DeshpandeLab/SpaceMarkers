#' @importFrom spatstat.geom owin ppp marks
#' @importFrom spatstat.explore Smooth
#' @importFrom ape Moran.I
#import description end
0



## author: Atul Deshpande
## email: adeshpande@jhu.edu

find_kernel_outliers_for_sensitivity <- function(pattern,locs,
                                                    pattern_threshold=0.15,
                                                sigma = 10,kernelthreshold = 2,
                                                method = "Pattern_Threshold",
                                                outlier = "positive",...)
{
    allwin<-spatstat.geom::owin(xrange = c(min(locs$x),max(locs$x)),
                                    yrange=c(min(locs$y),max(locs$y)))
    X<-spatstat.geom::ppp(x=locs$x,y=locs$y,window=allwin,marks=pattern)
    Kact<-spatstat.explore::Smooth(X,at ="points",sigma=sigma,...)
    Karr<-vapply(seq(1,100), function(i){Xr<-X;
    spatstat.geom::marks(Xr)<-sample(spatstat.geom::marks(X));
    temp<-spatstat.explore::Smooth(Xr,at="points",sigma=sigma,...);
        return(temp)}, numeric(length(Kact)))
    Kvec <- unlist(Karr)
    mKvec <- mean(Kvec)
    sKvec <- sd(Kvec)
    Kact<-(Kact-mKvec)/sKvec
    return(Kact)
}

getOptimalSigmaThresh <- function(pattern, locs, sigVec, threshVec,...){
    visium.dist <- as.matrix(dist(locs))
    visium.dist.inv <-1/visium.dist
    diag(visium.dist.inv) <- 0
    allwin<-spatstat.geom::owin(xrange=c(min(locs$x),max(locs$x)),
                                    yrange=c(min(locs$y),max(locs$y)))
    X<-spatstat.geom::ppp(x=locs$x,y=locs$y,window=allwin,marks=pattern)
    Ks<-vapply(sigVec,function(i) spatstat.explore::Smooth(X,at="points",
                                                            sigma=i,...),
                                                            numeric(X$n))
    mor_1<-vapply(seq(1,length(sigVec)),function(i) unlist(ape::Moran.I(
        spatstat.geom::marks(X)-Ks[,i], visium.dist.inv)),numeric(4))
    sigOpt1_ind <- which.min(abs(unlist(mor_1[1,])-unlist(mor_1[2,])))
    if (sigOpt1_ind>1&&sigOpt1_ind<length(sigVec)){
        smallsigVec<-seq(sigVec[sigOpt1_ind-1],sigVec[sigOpt1_ind+1],
                            (sigVec[sigOpt1_ind+1]-sigVec[sigOpt1_ind-1])/10)
    }else if (sigOpt1_ind==1){
        smallsigVec <- seq(sigVec[sigOpt1_ind],sigVec[sigOpt1_ind+1],
                            (sigVec[sigOpt1_ind+1] - sigVec[sigOpt1_ind])/10)
    }else{
        smallsigVec <- seq(sigVec[sigOpt1_ind-1],sigVec[sigOpt1_ind],
                            (sigVec[sigOpt1_ind] - sigVec[sigOpt1_ind-1])/10)
    }
    smallKs<-vapply(smallsigVec,function(i) spatstat.explore::Smooth(
            X,at="points", sigma = i, ...), numeric(X$n))
    smallmor_2<-vapply(seq(1,length(smallsigVec)), function(i) unlist(
            ape::Moran.I(spatstat.geom::marks(X)-smallKs[,i],visium.dist.inv)),
            numeric(4))
    sigOpt1_ind <- which.min(abs(unlist(smallmor_2[1,])-
                                unlist(smallmor_2[2,])))
    Kact2<-find_kernel_outliers_for_sensitivity(pattern=pattern,locs=locs,
                                                sigma=smallsigVec[sigOpt1_ind],
                                                method="Kernel2",
                                                kernelthreshold=0,
                                                outlier="positive")
    inds2<-(kronecker(matrix(1,1,length(threshVec)),Kact2)>threshVec)*1
    mor_2_ind<-vapply(seq(1,length(threshVec)),function(j){
        visium.dist.inv_1<- visium.dist.inv;visium.dist.inv_1[inds2[,j]==0]<-0;
        unlist(ape::Moran.I(spatstat.geom::marks(X)-smallKs[,sigOpt1_ind],
                            visium.dist.inv_1))}, numeric(4))
    threshOpt1_ind<-which.min(abs(unlist(mor_2_ind[1,])-unlist(mor_2_ind[2,])))
    return(data.frame(
        sigmaOpt=smallsigVec[sigOpt1_ind],threshOpt=threshVec[threshOpt1_ind]))
}

#===================
#' Calculate the optimal parameters from spatial kernel density for 
#' cell-cell interactions
#' @description This function uses Morans.I to calculate the optimal width of 
#' the kernel density (sigmaOpt) as well as the outlier threshold around the set
#' of spots (threshOpt) for a null distribution.
#' @export
#' @param spatialPatterns  A data frame that contains the spatial coordinates 
#' for each cell type. The column names must include 'x' and 'y' as well as a 
#' set of numbered columns named  'Pattern_1.....N'.
#' @param ... Arguments passed to methods
#' @return a numeric matrix of sigmaOpts - the optimal width of the gaussian 
#' distribution, and the threshOpt - outlier threshold around the set of spots 
#' for each pattern
#' @examples
#' library(SpaceMarkers)
#' # Create test data
#' cells <- c()
#' test_num <- 500
#' for(i in 1:test_num){
#'     cells[length(cells)+1] <- paste0("cell_",i)
#' }
#' spPatterns <- data.frame(barcode = cells,
#' y = runif(test_num, min=0, max=test_num),
#' x = runif(test_num, min=0, max=test_num),
#' Pattern_1 = runif(test_num, min=0, max=1),
#' Pattern_2 = runif(test_num, min=0, max=1) )
#' # Call the getSpatialParamsMoransI function with the test data
#' optParams <- getSpatialParamsMoransI(spPatterns)
#'
getSpatialParamsMoransI <- function(spatialPatterns,...){
    .Deprecated(new = "getSpatialParameters")
    sigmaRes <- max(floor(min(diff(range(spatialPatterns$x)),
                                diff(range(spatialPatterns$y)))/250),1)
    sigVec <- seq(2,40*sigmaRes,sigmaRes)
    threshVec <- seq(1,3,0.1)
    patternList <- setdiff(colnames(spatialPatterns),c("barcode","x","y"))
    optParams<-vapply(patternList,function(i) unlist(getOptimalSigmaThresh(
        pattern=spatialPatterns[,i],locs=data.frame(
        x=spatialPatterns$x,y=spatialPatterns$y),sigVec,threshVec,...)),numeric(2))
    return(optParams)
}

#===================
#' Read optimal parameters for spatial kernel density from user input or .json
#' file
#'
#' @description This function obtains the width of a spatial kernel density
#' (sigma) from either the user input or from a scale factors .json file. 
#' The outlier threshold around the set of spots (threshold) for each pattern is
#' specified  by the user (default is 4). 
#' @export
#' @param spatialPatterns  A data frame that contains the spatial coordinates 
#' for each cell type. The column names must include 'x' and 'y' as well as a 
#' set of numbered columns named  'Pattern_1.....N'.
#' @param visiumDir A string path specifying the location of the 10xVisium
#' directory
#' @param spatialDir A string path specifying the location of the spatial folder
#' containing the .json file of the spot characteristics
#' @param pattern A string specifying the name of the .json file
#' @param sigma A numeric value specifying the width of the kernel density 
#' estimate to be used for smoothing
#' @param threshold A numeric value specifying how many standard deviations 
#' above the mean of a null distribution to use an outlier threshold for
#' identifying 'hotspots'
#' @param resolution A string specifying image resolution to be used for spot 
#' diameter. Can take values of "fullres" (default), "lowres" or "hires".
#' @param ... Arguments passed to methods
#' @return a numeric matrix of sigmaOpts - the optimal width of the gaussian 
#' distribution, and the threshOpt - outlier threshold around the set of spots 
#' for each pattern
#' @examples
#' library(SpaceMarkers)
#' # Create test data
#' cells <- c()
#' test_num <- 500
#' for(i in 1:test_num){
#'     cells[length(cells)+1] <- paste0("cell_",i)
#' }
#' spPatterns <- data.frame(barcode = cells,
#' y = runif(test_num, min=0, max=test_num),
#' x = runif(test_num, min=0, max=test_num),
#' Pattern_1 = runif(test_num, min=0, max=1),
#' Pattern_2 = runif(test_num, min=0, max=1) )
#' # Call the getSpatialParameters function with the test data
#' optParams <- getSpatialParameters(spPatterns, sigma = 10)
#'

getSpatialParameters <- function(spatialPatterns,visiumDir = ".",
                                         spatialDir ="spatial",
                                         pattern = "scalefactors_json.json",
                                         sigma = NULL,threshold = 4,
                                         resolution = 
                                           c("fullres","lowres","hires"), ...) {
    #resolve resolution parameter
    resolution <- match.arg(resolution, several.ok = FALSE)
    message("resolution: ", resolution)
    patternList <- setdiff(colnames(spatialPatterns),c("barcode","x","y"))
    if (!is.null(sigma)) {
      sigmaOpt <- sigma
      threshOpt <- threshold
    } else if (is.null(sigma) & 
               file.exists(file.path(visiumDir,spatialDir,pattern))) {
        message("Reading spot diameter from specified .json file")
        scale_values <- jsonlite::read_json(file.path(visiumDir,spatialDir,
                                                    pattern))
      if ("lowres" %in% resolution){
        resolution <- "lowres"
      } else if ("hires" %in% resolution){
        resolution <- "hires"
      } else if ("fullres" %in% resolution){
        resolution <- "fullres"
      } else {
        stop("Resolution argument not recognized.
             Please supply either lowres, hires or fullres.")
      }
      
      threshOpt <- threshold
      sigmaOpt <- as.numeric(scale_values$spot_diameter_fullres)
      if (resolution != "fullres") {
        scale_factor <- scale_values[grepl(resolution,
                                           names(scale_values))][[1]]
        sigmaOpt <- sigmaOpt * as.numeric(scale_factor)
      } 
      
      
    } else  {
      stop("Please specify the sigma or correct path to .json")
    }
    optParams <-matrix(c(sigmaOpt,threshOpt),nrow = 2,
                       ncol = length(patternList))
    colnames(optParams) <- patternList
    rownames(optParams) <- c("sigmaOpt","threshOpt")
  return(optParams)
}
