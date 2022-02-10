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
  Kact1 = spatstat.core::Smooth(X, at = "points", sigma = sigmaPair[1], leaveoneout = T)
  Karr1 = sapply(seq(1,100), function(i) {Xr = X; spatstat.geom::marks(Xr) = spatstat.geom::marks(X)[pracma::randperm(1:length(spatstat.geom::marks(X)))]; temp = spatstat.core::Smooth(Xr, at = "points", sigma = sigmaPair[1], leaveoneout = T); return(temp)})
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


getInteractingGenes <- function(data, reconstruction=NULL, spatialPatterns, refPattern="Pattern_1", mode=c("residual","DE")){
    if (mode=="residual"&&is.null(reconstruction)) stop("Reconstruction matrix not provided for residual mode.")
    if (all(dim(data)!=dim(reconstruction))) stop("Original and reconstructed matrix do not have the same dimensions.")
    patternList <- colnames(spatialPatterns)[startsWith(colnames(spatialPatterns),"Pattern_")]
    hotspot.regions = c()
    for (patternName in patternList)
      {
        hotspot.regions <- cbind(hotspot.regions,find_pattern_hotspots(spatialPatterns = spatialPatterns, patternName = patternName, params = optParams[,patternName], outlier = "positive") )
    }
    colnames(hotspot.regions) <- patternList
    
    interacting.genes <- list();
    data <- as.matrix(data)
    
    for (pattern in setdiff(patternList,refPattern)){
      region <- hotspot.regions[,refPattern];
      region <- ifelse(!is.na(region) & !is.na(hotspot.regions[,pattern]),"Interacting",ifelse(!is.na(region),region,hotspot.regions[,pattern]))
      region <- factor(region)
      if (length(levels(region))<3||any(table(region)<50)) #default 50
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
      interacting_genes[[i]]$p.adj <- as.numeric(interacting_genes[[i]]$p.adj)
    
    for (i in seq(1,length(interacting_genes)))
    {
      if (all(dim(interacting_genes[[i]])>1))   {
        od <- order(interacting_genes[[i]]$p.adj)
        interacting_genes[[i]] <- interacting_genes[[i]][od,]
      }
    }
    return(list(interacting_genes=interacting_genes,hotspot.regions=hotspot.regions))
}
    
  
  
  
  
  

  