
#' @title getOverlapScores
#' @description Calculate the overlap scores between patterns in hotspots
#' @param hotspots A data frame with columns x, y, barcode and pattern names
#' @param patternList A character vector of pattern names to calculate overlap scores for
#' @return A data frame with columns pattern1, pattern2 and overlapScore
#' @export
#' @examples
#' hotspots <- data.frame(x = c(1,2,3,4,5), y = c(1,2,3,4,5), barcode = c("A","B","C","D","E"), pattern1 = c(1,0,1,0,1), pattern2 = c(1,1,0,0,1))
#' getOverlapScores(hotspots)   
#' getOverlapScores(hotspots, c("pattern1","pattern2"))
#' @import ggplot2
#' @import reshape2
getOverlapScores <- function(hotspots,
                             patternList = NULL) {
    if (is.null(patternList))
        patternList <- setdiff(colnames(hotspots),c("x","y","barcode"))
    else if (!all(patternList %in% colnames(hotspots)))
        stop("Pattern names not found in hotspots")
  binarized <- (!is.na(hotspots[,patternList]))*1
  intersects <- t(binarized) %*% binarized
  unions <- sapply(colnames(binarized),function(c) colSums((binarized[,c] + binarized)>0))
  Jaccard <- intersects/unions
  
  nHotspots <- colSums(binarized)
  nHotsP1 <- t(t(nHotspots)) %*% array(1, length(patternList))
  nHotsP2 <- t(nHotsP1)
  bestJaccard <- pmin(nHotsP1,nHotsP2)/pmax(nHotsP1,nHotsP2)
  
  Jaccard[upper.tri(Jaccard,diag = TRUE)] <- NA
  normJaccard <- Jaccard/bestJaccard
  normJaccard[upper.tri(normJaccard,diag = TRUE)] <- NA
  
  # Melt normalized Jaccard for output
  dfJacc <- reshape2::melt(normJaccard)
  dfJacc <- dfJacc[complete.cases(dfJacc),]
  # Due to melting in lower triangular orientation, the column names are flipped
  colnames(dfJacc) <- c("pattern2", "pattern1", "overlapScore")
  dfJacc <- dfJacc[,c(2,1,3)]

  return(dfJacc)
}

#' @title plotOverlapScores
#' @description Plot the overlap scores between patterns in hotspots
#' @param df A data frame with columns pattern1, pattern2 and overlapScore
#' @param title The title of the plot
#' @param output_path The path to save the plot
#' @param fontsize The font size of the plot
#' @return A ggplot object
#' @export
#' @examples
#' df <- data.frame(pattern1 = c("pattern1","pattern1","pattern2","pattern2"), pattern2 = c("pattern1","pattern2","pattern1","pattern2"), overlapScore = c(0.5,0.7,0.3,0.9))
#' plotOverlapScores(df)
#' plotOverlapScores(df, "Overlap Scores", "overlapScores.png", 15)
#' @import ggplot2
#'
plotOverlapScores <- function(df, title = "Spatial Overlap Scores", output_path = "overlapScores.png",fontsize = 15) {
    p <- ggplot2::ggplot(data = df, aes(pattern1, pattern2, fill = overlapScore)) +
        geom_tile(color = "black", size = 0.8) +
        geom_text(aes(label = round(overlapScore, 2)), size = 6) +  # Display values on the plot
        scale_fill_gradient2(low = "#FFF7EC", mid = "#FDBB84", high = "#D7301F", midpoint = 0.5) +
        scale_y_discrete(limits = rev, guide = guide_axis(angle = 45)) +
        theme_minimal() +
        theme(
            axis.text.x = element_text(angle = 45, vjust = 1, size = fontsize, hjust = 1),
            axis.text.y = element_text(size = fontsize)
        ) +
        coord_fixed() +
        labs(x = NULL, y = NULL, fill = "overlapScore", title = title) +
        theme(
            panel.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.border = element_blank()
        )

    return(p)
}

#' @title getIMScores
#' @description Get the interaction scores for SpaceMarkers
#' @param SpaceMarkers A list of SpaceMarkers objects
#' @return A data frame with columns Gene and SpaceMarkersMetric
#' @export
#' @examples
#' example(getPairwiseInteractingGenes)
#' getIMScores(SpaceMarkers)
#'
getIMScores <- function(SpaceMarkers){
    smi <- SpaceMarkers[which(sapply(SpaceMarkers, function(x) length(x[['interacting_genes']]))>0)]
    fields <- c('Gene', 'SpaceMarkersMetric')

    imscores <- lapply(seq_along(smi), function(x) {
        df <- smi[[x]][['interacting_genes']][[1]][,fields]
        #rename to metric to its parent item name
        setNames(df, c('Gene', names(smi)[x]))
    })

    imscores <- Reduce(function(x, y) {
                merge(x, y, by="Gene", all=TRUE)
            }, x=imscores, right=FALSE)
    
    if(is.null(imscores)) {
        imscores <- data.frame(Gene=character(0))
    }
    return(imscores)
}