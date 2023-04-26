#' @import spatstat
#' @import ape
#' @import pracma
#import description end
0



## author: Atul Deshpande
## email: adeshpande@jhu.edu

find_kernel_outliers_for_sensitivity <- function(pattern, locs, pattern_threshold = 0.15, sigma = 10, kernelthreshold = 2, method = "Pattern_Threshold", outlier = "positive")
{
  allwin <- spatstat.geom::owin(xrange = c(min(locs$x),max(locs$x)), yrange = c(min(locs$y),max(locs$y)))
  X <-spatstat.geom::ppp(x = locs$x, y = locs$y, window = allwin, marks = pattern)
  
  Kact = spatstat.explore::Smooth(X, at = "points", sigma = sigma, leaveoneout = T)
  Karr = sapply(seq(1,100), function(i) {Xr = X; spatstat.geom::marks(Xr) = spatstat.geom::marks(X)[pracma::randperm(1:length(spatstat.geom::marks(X)))]; temp = spatstat.explore::Smooth(Xr, at = "points", sigma = sigma, leaveoneout = T); return(temp)})
  Kvec = unlist(Karr)
  mKvec = mean(Kvec)
  sKvec = sd(Kvec)
  Kact<-(Kact-mKvec)/sKvec
  return(Kact)
}

getOptimalSigmaThresh <- function(pattern, locs, sigmaVec, threshVec){
  visium.dist <- as.matrix(dist(locs))
  visium.dist.inv <-1/visium.dist
  diag(visium.dist.inv) <- 0
  allwin <- spatstat.geom::owin(xrange = c(min(locs$x),max(locs$x)), yrange = c(min(locs$y),max(locs$y)))
  X <-spatstat.geom::ppp(x = locs$x, y = locs$y, window = allwin, marks = pattern)
  Ks <- sapply(sigmaVec,function(i) spatstat.explore::Smooth(X, at = "points", sigma = i, leaveoneout = T))
  mor_1 <- sapply(1:length(sigmaVec), function(i) ape::Moran.I(spatstat.geom::marks(X)-Ks[,i], visium.dist.inv))
  sigmaOpt1_ind <- which.min(abs(unlist(mor_1[1,])-unlist(mor_1[2,])))
  if (sigmaOpt1_ind>1&&sigmaOpt1_ind<length(sigmaVec)){
    smallSigmaVec <- seq(sigmaVec[sigmaOpt1_ind-1],sigmaVec[sigmaOpt1_ind+1],(sigmaVec[sigmaOpt1_ind+1] - sigmaVec[sigmaOpt1_ind-1])/10)
  }else if (sigmaOpt1_ind==1){
    smallSigmaVec <- seq(sigmaVec[sigmaOpt1_ind],sigmaVec[sigmaOpt1_ind+1],(sigmaVec[sigmaOpt1_ind+1] - sigmaVec[sigmaOpt1_ind])/10)
  }else{
    smallSigmaVec <- seq(sigmaVec[sigmaOpt1_ind-1],sigmaVec[sigmaOpt1_ind],(sigmaVec[sigmaOpt1_ind] - sigmaVec[sigmaOpt1_ind-1])/10)
  }
  
  smallKs <- sapply(smallSigmaVec,function(i) spatstat.explore::Smooth(X, at = "points", sigma = i, leaveoneout = T))
  smallmor_2 <- sapply(1:length(smallSigmaVec), function(i) ape::Moran.I(spatstat.geom::marks(X)-smallKs[,i], visium.dist.inv))
  sigmaOpt1_ind <- which.min(abs(unlist(smallmor_2[1,])-unlist(smallmor_2[2,])))
  Kact2 = find_kernel_outliers_for_sensitivity(pattern = pattern, locs = locs, sigma = smallSigmaVec[sigmaOpt1_ind], method = "Kernel2", kernelthreshold = 0, outlier = "positive")
  inds2<- (kronecker(matrix(1,1,length(threshVec)),Kact2)>threshVec)*1
  mor_2_ind <- sapply(1:length(threshVec), function(j){visium.dist.inv_1 <- visium.dist.inv; visium.dist.inv_1[inds2[,j]==0] = 0; ape::Moran.I(spatstat.geom::marks(X)-smallKs[,sigmaOpt1_ind],visium.dist.inv_1)})
  threshOpt1_ind <- which.min(abs(unlist(mor_2_ind[1,])-unlist(mor_2_ind[2,])))
  return(data.frame(sigmaOpt = smallSigmaVec[sigmaOpt1_ind], threshOpt = threshVec[threshOpt1_ind]))
}
#===================
#' getSpatialParameters
#' Calculate the Optimal Parameters for Interacting Cells
#'
#' This function calculates the optimal width of the gaussian distribution (sigmaOpt) as well as the outlier threshold around the set of spots (thresOpt) for each pattern from a latent feature space.
#'
#' @export
#'
#' @param spatialPatterns  A data frame that contains the spatial coordinates for each cell type. The column names must include 'x' and 'y' as well as a set of numbered columns named  'Pattern_1.....N'.
#'
#'

getSpatialParameters <- function(spatialPatterns){
  good_gene_threshold <- 3;
  sigmaRes <- max(floor(min(diff(range(spatialPatterns$x)),diff(range(spatialPatterns$y)))/250),1)
  sigmaVec <- seq(2,40*sigmaRes,sigmaRes)
  threshVec <- seq(1,3,0.1)
  
  patternList <- colnames(spatialPatterns)[startsWith(colnames(spatialPatterns),"Pattern_")]
  optParams<-sapply(patternList, function(i) unlist(getOptimalSigmaThresh(pattern = spatialPatterns[,i], locs = data.frame(x = spatialPatterns$x, y = spatialPatterns$y), sigmaVec, threshVec)))
  return(optParams)
}
