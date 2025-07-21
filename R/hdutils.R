#' @title Compute the spatial influence of a spatial feature
#' @description This function computes the spatial influence of a specified pattern
#' @param spPatterns A data frame containing x, y coordinates and pattern name
#' @param optParams A data frame with optimal parameters for the pattern
#' @param ... Additional parameters for the Smooth function
#' @return A data frame with the spatial influence of the specified pattern
#' @export
calcInfluence <- function(spPatterns, optParams,...) {
 
    patnames <- setdiff(colnames(spPatterns),
                       c("x", "y", "barcode"))

    allwin <- spatstat.geom::owin(
    range(spPatterns$x),
    range(spPatterns$y))
    X <- spatstat.geom::ppp(x = spPatterns$x, y = spPatterns$y,
                            window = allwin,
                            marks = spPatterns[,patnames[1]])

    spInfluence <- sapply(patnames, function(pat) {
    # Create a point pattern object for each pattern
    X <- spatstat.geom::ppp(x = spPatterns$x, y = spPatterns$y,
                            window = allwin,
                            marks = spPatterns[,pat])
    
    # Calculate the kernel for the specified pattern
    Kact1 <- spatstat.explore::Smooth(
      X, at = "points", sigma = optParams[1,pat],...)
    
    # Plot the K-function
    return(Kact1)
    })
    spInfluence <- as.data.frame(spInfluence)
    colnames(spInfluence) <- patnames
    spInfluence <- cbind(spPatterns[,c("barcode","x", "y")], spInfluence)

    return(spInfluence)
}

#' @title Compute the threshold for identifying outlier values or hotspots
#' @description This function computes the threshold for identifying outlier 
#' values or hotspots by fitting a normal mixture model to the data.
#' @param df A vector containing pattern values
#' @param minval Minimum value for quantile threshold
#' @param maxval Maximum value for quantile threshold
#' @return A list containing the computed thresholds
#' @importFrom mixtools normalmixEM
#' @export
calcThresholds <- function(df, minval = 0.01, maxval = 0.99) {

#Calculate the two compmonents of the normal mixture model
comps <- mixtools::normalmixEM(df,k=2)

# Identify smaller component
small <- which.min(comps$mu)
minthresh <- quantile(df,minval)
maxthresh <- quantile(df,maxval)
thresh <- min(max(comps$mu[small] + (comps$sigma[small] * 4), minthresh),maxthresh)

return(thresh)
}
#' @title Compute the thresholds for all columns in a data frame
#' @description This function computes the thresholds for all columns in a 
#' data frame. The data frame could be an spPatterns object or an spInfluence
#' object.
#' @param df A data frame with pattern values (optionally 
#' with x, y, barcode columns)
#' @param minvals Minimum value for quantile threshold
#' @param maxvals Maximum value for quantile threshold
#' @param ... Additional parameters to pass to lower level functions
#' @return A list containing the computed thresholds for each pattern
#' @export 
calcAllThresholds <- function(df, minvals = 0.01, maxvals = 0.99,...) {
  
  patnames <- setdiff(colnames(df), c("x", "y", "barcode"))
  #Check if minvals and maxvals are vectors
  if (length(minvals) == 1) minvals <- rep(minvals, length(patnames))
  names(minvals) <- patnames
  if (length(maxvals) == 1) maxvals <- rep(maxvals, length(patnames))
  names(maxvals) <- patnames
  # Check if minvals and maxvals are of the same length and match 
  # the number of columns in df
  if ((length(minvals) == length(maxvals)) && (ncol(df)==length(minvals)))
        stop("minvals and maxvals must be scalar or vectors of 
              the length ncol(df)")
    # Calculate thresholds for all patterns
    thresholds <- sapply(patnames, function(pat) {
        calcThresholds(df[,pat], minval=minvals[pat], maxval=maxvals[pat],...)
    })
    return(thresholds)
}


#' @title Find hotspots for all patterns or influences based on values
#' @description Convenience function to find hotspots for all spatial patterns 
#' or influence dataframes based on provided thresholds
#' @inheritParams calcAllThresholds
#' @param threshold a scalar or vector of thresholds for each column in the data frame.
#'  Either user provided or the output of @calcAllThresholds
#' @return a data frame with the same dimensions as the input data frame.
#' @export 
findAllHotspots.value <- function(df, threshold = 0.1,...){
    patnames <- setdiff(colnames(df),c("x","y","barcode"))
    if (length(threshold)==1){
        threshold <- rep(threshold,length(patnames))
    }
    if (length(threshold)!=length(patnames)){
        stop("Length of threshold must be 1 or equal to number of patterns.")
    }
    names(threshold) <- patnames
    
    hotspots <- matrix(NA, nrow=nrow(df), ncol=length(patnames))
    colnames(hotspots) <- patnames
    for (pat in patnames){
        hotspots[,pat] <- ifelse(
            df[,pat]>threshold[pat],pat,NA)
    }
    hotspots <- cbind(df[c("barcode","y","x")],hotspots)
    row.names(hotspots) <- hotspots$barcode
    hotspots <- as.data.frame(hotspots)
    return(hotspots)
}

classifySpots <- function(patHotspots, infHotspots, patternpair = NULL) {
    patnames <- setdiff(colnames(patHotspots), c("x", "y", "barcode"))
    infnames <- setdiff(colnames(infHotspots), c("x", "y", "barcode"))  
    #check if patHotspots and infHotspots have the same dimensions
    if (!all(dim(patHotspots) == dim(infHotspots))) {
        stop("patHotspots and infHotspots must have the same dimensions.")
    } else if (!all(patnames == infnames)) {
        stop("patHotspots and infHotspots must have the same column names.")
    }
    # Check if patternpair is NULL or is contained in patnames
    if (!is.null(patternpair) && !all(patternpair %in% patnames)) {
        stop("patternpair must be NULL or contained in patnames.")
    } else if (is.null(patternpair) && length(patnames) > 2) {
        stop("More than 2 patterns found. Please provide patternpair.")
    }
    
    region <- patHotspots[,patternpair[1]]
    pat1.inf2 <- ifelse(!is.na(region) & !is.na(infHotspots[,patternpair[2]]),"Interacting",
    ifelse(!is.na(region),region,NA))

    region <- patHotspots[,patternpair[2]]
    pat2.inf1 <- ifelse(!is.na(region) & !is.na(infHotspots[,patternpair[1]]),"Interacting",
    ifelse(!is.na(region),region,NA))
    df <- data.frame(pat1.inf2=pat1.inf2, pat2.inf1=pat2.inf1)
return(df)
}

#' @title Perform row-wise t-tests from scratch
#' @description This function iterates over the rows of a matrix and performs a
#' t-test comparing two groups of columns. It calculates the t-statistic, p-value,
#' and sample sizes without relying on `stats::t.test()` for the core logic.
#' @param in.data A numeric matrix. Rows represent features, columns represent samples.
#' @param region A factor or vector indicating the group membership for each column of `in.data`.
#'                Must have exactly two levels/unique values. Its length must equal `ncol(in.data)`.
#' @param min_bins Minimum number of non-missing observations required in each group to perform the t-test.
#' @param ... Additional arguments passed to the t-test function.
#' @return A matrix with rows corresponding to the features and columns:
#'          - `statistic`: The calculated t-statistic.
#'          - `p.value`: The calculated two-sided p-value.
#'         - `n1`: Number of non-missing observations in group 1 for that row.
#'         - `n2`: Number of non-missing observations in group 2 for that row.
#' @importFrom stats t.test
row.t.test <- function(in.data,region,min_bins=50,...){
    if (!is.factor(region)) {
        region <- factor(region)
    }
    if (nlevels(region) != 2) {
        stop("'region' must have exactly two levels (groups).")
    }
    group_levels <- levels(region)
    int <- which(as.character(group_levels)=="Interacting")
    interacting <- group_levels[int]
    patname <- group_levels[-int]
    idx_pat <- which(region == patname)
    idx_int <- which(region == interacting)
    if (length(idx_pat) >= min_bins || length(idx_int) >= min_bins) {
    temp <- sapply(rownames(in.data), function(r) {
        pat <- in.data[r, idx_pat]
        int <- in.data[r, idx_int]
    if (length(x) < min_bins || length(y) < min_bins) {
            return(c(statistic=NA, p.value=NA, n1=0, n2=0))
        }
        tmp <- t.test(x=int, y=pat, 
                      alternative = "two.sided", var.equal = FALSE, 
                      na.action = na.omit, ...)
        return(c(statistic=tmp$statistic, p.value=tmp$p.value, n1=length(x), n2=length(y)))
    })
    } else {
        temp <- matrix(NA, nrow = nrow(in.data), ncol = 4)
        return(c(statistic=array(NA, dim=nrow(in.data)),p.value=array(NA, dim=nrow(in.data)), n1=0, n2=0))
    }
    temp <- t(temp)
    return(temp)
}

calcIMscores.HD <- function(data, patHotspots, infHotspots, patternpair,...) {
        spotClass <- classifySpots(patHotspots, infHotspots, patternpair = patternpair)
        pat1 <- patternpair[1]
        pat2 <- patternpair[2]
        t1table <- row.t.test(data, spotClass[,1])
        colnames(t1table)[1] <- paste0("t_", pat1, "_near_", pat2)
        t2table <- row.t.test(data, spotClass[,2],...)
        colnames(t2table)[1] <- paste0("t_", pat2, "_near_", pat1)
        tscores <- cbind(t1table[,1], t2table[,1])
        colnames(tscores) <- c(paste0("t_", pat1, "_near_", pat2), paste0("t_", pat2, "_near_", pat1))
        return(tscores)
}

#' @title Calculate interaction scores for all pattern pairs
#' @description This function calculates interaction scores for all pattern pairs
#' using the `calcIMscores.HD` function. It can run in parallel if BiocParallel is available.
#' @param data A numeric matrix with genes as rows and barcodes as columns.
#' @param patHotspots A data frame with pattern hotspots, containing columns for x, y, and barcode.
#' @param infHotspots A data frame with influence hotspots, containing columns for x, y, and barcode.
#' @param patternPairs A data frame with pattern pairs to calculate interaction scores for. If NULL, 
#' all combinations of patterns in `patHotspots` will be used.
#' If provided, it should have two columns with pattern names. 
#' Each row should represent a pair of patterns for which interaction scores will be calculated.
#' @param ... Additional parameters to pass to lower level functions.
#' @return A data frame with interaction scores for all pattern pairs.
#' @export
calcAllIMscores.HD <- function(data, patHotspots, infHotspots, patternPairs=NULL,...) {
    if (is.null(patternPairs)) {
        patternPairs <- utils::combn(setdiff(colnames(patHotspots), c("x", "y", "barcode")), 2, simplify = FALSE)
    }
    if (requireNamespace("BiocParallel", quietly = TRUE) && BiocParallel::bpparam()$workers > 1) {
        bpp <- BiocParallel::bpparam()
        IMscores_list <- BiocParallel::bplapply(
            seq_len(nrow(patternPairs)),
            function(i) {
                patternpair <- patternPairs[i,]
                IMscores.pair <- calcIMscores.HD(data, patHotspots, infHotspots, patternpair)
                colnames(IMscores.pair) <- c(
                    paste0("t_", patternpair[1], "_near_", patternpair[2]),
                    paste0("t_", patternpair[2], "_near_", patternpair[1])
                )
                IMscores.pair
            },
            BPPARAM = bpp
        )
        IMscores <- do.call(cbind, IMscores_list)
        message("Processed all pattern pairs in parallel.\n")
    } else {
        IMscores <- c()
        for (i in 1:nrow(patternPairs)) {
            patternpair <- patternPairs[i,]
            IMscores.pair <- calcIMscores.HD(data, patHotspots, infHotspots, patternpair)
            IMscores <- cbind(IMscores, IMscores.pair)
            colnames(IMscores)[(ncol(IMscores)-1):ncol(IMscores)] <- c(
                paste0("t_", patternpair[1], "_near_", patternpair[2]),
                paste0("t_", patternpair[2], "_near_", patternpair[1])
            )
            message("Processed pattern pair:", patternpair[1], " and ", patternpair[2], "\n")
            if (i %% 10 == 0) {
                message("Processed", i, "pattern pairs so far.\n")
            }
        }
    }
    return(IMscores)
}

#' @title getOverlapScores.HD
#' @description Calculate the overlap scores between patterns in hotspots
#' @param patHotspots A data frame with columns x, y, barcode and pattern names
#' @param infHotspots A data frame with columns x, y, barcode and pattern names
#' @param method The method to calculate overlapping abundance scores. Options are
#' "relative-abundance", "differential-abundance" and "absolute"
#' @param patternList A character vector of pattern names to calculate overlap 
#' scores for. If NULL, all patterns in patHotspots and infHotspots will be used.
#' @details The function calculates the overlap scores between patterns hotspots
#' using the specified method. The default method is "relative-abundance"
#' @return A data frame with columns pattern, influence and overlapping abundance
#' @export
#' @examples
#' hotspots <- data.frame(x = c(1,2,3,4,5),
#'                         y = c(1,2,3,4,5),
#'                         barcode = c("A","B","C","D","E"),
#'                         pattern1 = c(1,0,1,0,1),
#'                         pattern2 = c(1,1,0,0,1))
#' getOverlapScores(hotspots)   
#' getOverlapScores(hotspots, c("pattern1","pattern2"))
#' @importFrom ggplot2 ggplot geom_tile geom_text theme_minimal
#' @importFrom reshape2 melt
#' @importFrom stats complete.cases
getOverlapScores.HD <- function(patHotspots, infHotspots,
                             patternList = NULL, method = c("relative-abundance",
                                                            "differential-abundance",
                                                             "absolute") ) {
    
    #warn if more than one method is supplied, do not warn by default
    if(length(method) > 1){
        method <- method[1]
        message("Only one method can be used at a time. Using ", method)}

    if (is.null(patternList)) {
        patternList <- setdiff(colnames(patHotspots),c("x","y","barcode"))
    } else if (!all(patternList %in% colnames(patHotspots))) {
        stop("Pattern names not found in hotspots")
    }

    patBinarized <- (!is.na(patHotspots[,patternList]))*1
    infBinarized <- (!is.na(infHotspots[,patternList]))*1
    intersects <- t(patBinarized) %*% infBinarized
    #nHotspots <- colSums(binarized)
    nHotsP1 <- t(t(colSums(patBinarized))) %*% array(1, length(patternList))
    nHotsP2 <- matrix(1, nrow=length(patternList),ncol=1) %*% colSums(infBinarized)
    colnames(nHotsP1) <- patternList
    colnames(nHotsP2) <- patternList
    rownames(nHotsP2) <- patternList
    overlapScore <- switch(method,
        "relative-abundance" = intersects/nHotsP2/(nHotsP1[,1]/nrow(patBinarized)),
        "differential-abundance" = intersects/nHotsP2 - (nHotsP1 - intersects)/(nrow(patBinarized) - nHotsP2),
        "absolute" = intersects,
        stop("Method not supported")
    )

    diag(overlapScore) <- NA
    colnames(overlapScore) <- paste0("near.",colnames(overlapScore))
    #overlapScore[upper.tri(overlapScore,diag = TRUE)] <- NA
  
    # Melt normalized Jaccard for output
    dfOverlap <- reshape2::melt(overlapScore)
    dfOverlap <- dfOverlap[stats::complete.cases(dfOverlap),]
    # Due to melting in lower triangular orientation, the column names are flipped
    colnames(dfOverlap) <- c("pattern", "influence", "relAbundance")
    return(dfOverlap)
}
