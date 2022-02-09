find_pattern_hotspots <- function(sample2Exprs, params = NULL, patternName = "Pattern_1", outlier = "positive")
{
  require("spatstat","pracma")
  if (is.null(params)){
    sigmaPair = 10
    kernelthreshold = 2
  }
  else{
    sigmaPair = params["sigmaOpt"]
    kernelthreshold = params["threshOpt"]
  }

    allwin = owin(xrange = c(min(sample2Exprs$x),max(sample2Exprs$x)), yrange = c(min(sample2Exprs$y),max(sample2Exprs$y)))
    patternVector = as.matrix(sample2Exprs[,patternName])
    X <-ppp(x = sample2Exprs$x, y = sample2Exprs$y, window = allwin, marks = patternVector)
    Kact1 = Smooth(X, at = "points", sigma = sigmaPair[1], leaveoneout = T)
    Karr1 = sapply(seq(1,100), function(i) {Xr = X; marks(Xr) = marks(X)[pracma::randperm(1:length(marks(X)))]; temp = Smooth(Xr, at = "points", sigma = sigmaPair[1], leaveoneout = T); return(temp)})
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