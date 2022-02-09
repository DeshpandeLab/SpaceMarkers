find_kernel_outliers_for_sensitivity <- function(pattern, locs, pattern_threshold = 0.15, sigma = 10, kernelthreshold = 2, method = "Pattern_Threshold", outlier = "positive")
{
  require("spatstat","pracma")
  allwin <- owin(xrange = c(min(locs$x),max(locs$x)), yrange = c(min(locs$y),max(locs$y)))
  X <-ppp(x = locs$x, y = locs$y, window = allwin, marks = pattern)
    
    Kact = Smooth(X, at = "points", sigma = sigma, leaveoneout = T)
    Karr = sapply(seq(1,100), function(i) {Xr = X; marks(Xr) = marks(X)[pracma::randperm(1:length(marks(X)))]; temp = Smooth(Xr, at = "points", sigma = sigma, leaveoneout = T); return(temp)})
    Kvec = unlist(Karr)
    mKvec = mean(Kvec)
    sKvec = sd(Kvec)
    Kact<-(Kact-mKvec)/sKvec
    return(Kact)
}
