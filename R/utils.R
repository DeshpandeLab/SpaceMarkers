
#' @title getOverlapScores
#' @description Calculate the overlap scores between patterns in hotspots
#' @param hotspots A data frame with columns x, y, barcode and pattern names
#' @param patternList A character vector of pattern names to calculate overlap 
#' scores for
#' @param method The method to calculate overlap scores. Options are
#' "Szymkiewicz–Simpson", "Jaccard", "Sørensen–Dice", "Ochiai" and "absolute"
#' @details The function calculates the overlap scores between patterns hotspots
#' using the specified method. The default method is "Szymkiewicz–Simpson" 
#' overlap coefficient.
#' @return A data frame with columns pattern1, pattern2 and overlapScore
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
getOverlapScores <- function(hotspots,
                             patternList = NULL, method = c("Szymkiewicz–Simpson",
                                                            "Jaccard", "Sørensen–Dice", 
                                                            "Ochiai", "absolute") ) {
    if (is.null(patternList))
        patternList <- setdiff(colnames(hotspots),c("x","y","barcode"))
    else if (!all(patternList %in% colnames(hotspots)))
        stop("Pattern names not found in hotspots")
    binarized <- (!is.na(hotspots[,patternList]))*1
    intersects <- t(binarized) %*% binarized
    nHotspots <- colSums(binarized)
    nHotsP1 <- t(t(nHotspots)) %*% array(1, length(patternList))
    nHotsP2 <- t(nHotsP1)
    if ("Szymkiewicz–Simpson" %in% method) {
        # Szymkiewicz–Simpson index
        overlapScore <- intersects/pmin(nHotsP1,nHotsP2)
    } else if ("Jaccard" %in% method) {
        # Jaccard index
        overlapScore <- intersects/(nHotsP1 + nHotsP2 - intersects)
    } else if ("Sørensen–Dice" %in% method) {
        # Sørensen–Dice index
        overlapScore <- 2*intersects/(nHotsP1 + nHotsP2)
    } else if ("Ochiai" %in% method) {
        # Ochiai index
        overlapScore <- intersects/sqrt(nHotsP1*nHotsP2)
    } else if ("absolute" %in% method) {
        # Absolute overlap
        overlapScore <- intersects
    } else
        stop("Method not supported")
  
    overlapScore[upper.tri(overlapScore,diag = TRUE)] <- NA
  
    # Melt normalized Jaccard for output
    dfOverlap <- reshape2::melt(overlapScore)
    dfOverlap <- dfOverlap[stats::complete.cases(dfOverlap),]
    # Due to melting in lower triangular orientation, the column names are flipped
    colnames(dfOverlap) <- c("pattern2", "pattern1", "overlapScore")
    dfOverlap <- dfOverlap[,c(2,1,3)]
    return(dfOverlap)
}

#' @title plotOverlapScores
#' @description Plot the overlap scores between patterns in hotspots
#' @param df A data frame with columns pattern1, pattern2 and overlapScore
#' @param title The title of the plot
#' @param fontsize The font size of the plot
#' @param out The output path for the plot
#' @return A ggplot object
#' @export
#' @examples
#' df <- data.frame(pattern1 = c("pattern1","pattern1","pattern2","pattern2"), pattern2 = c("pattern1","pattern2","pattern1","pattern2"), overlapScore = c(0.5,0.7,0.3,0.9))
#' plotOverlapScores(df)
#' plotOverlapScores(df, "Overlap Scores", "overlapScores.png", 15)
#' @import ggplot2
#'
plotOverlapScores <- function(df, title = "Spatial Overlap Scores", out = NULL,fontsize = 15) {
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
    if(!is.null(out)){
        ggplot2::ggsave(filename = out, plot = p)
    }
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
#' @importFrom stats setNames
#' 
getIMScores <- function(SpaceMarkers){
    smi <- SpaceMarkers[which(sapply(SpaceMarkers, function(x) length(x[['interacting_genes']]))>0)]
    fields <- c('Gene', 'SpaceMarkersMetric')

    imscores <- lapply(seq_along(smi), function(x) {
        df <- smi[[x]][['interacting_genes']][[1]][,fields]
        #rename to metric to its parent item name
        stats::setNames(df, c('Gene', names(smi)[x]))
    })

    imscores <- Reduce(function(x, y) {
                merge(x, y, by="Gene", all=TRUE)
            }, x=imscores, right=FALSE)
    
    if(is.null(imscores)) {
        imscores <- data.frame(Gene=character(0))
    }
    return(imscores)
}


#' @title plotIMScores
#' @description Plot the top SpaceMarkers IMScores
#' @param df A data frame with columns Gene and SpaceMarkersMetric
#' @param interaction The interaction to plot
#' @param cutOff The cut off value for the plot
#' @param nGenes The number of genes to plot
#' @param geneText The font size for the gene text
#' @param metricText The font size for the metric text
#' @param increments The increments for the y-axis
#' @param out The output path for the plot
#' @export
#' @examples 
#' example(getPairwiseInteractingGenes)
#' plotIMScores(getIMScores(SpaceMarkers), "Pattern_1_Pattern_3")
#' @import ggplot2
#' @importFrom stats reorder
#' @importFrom utils head
plotIMScores <- function(df, interaction, cutOff = 0, nGenes = 20,
    geneText = 12, metricText = 12, increments = 1, out = NULL) {
    df$genes <- df$Gene
    df <- df[order(df[[interaction]], decreasing = TRUE),]
    df <- utils::head(df,nGenes)
    df[[interaction]] <- round(df[[interaction]],2)
    
    p1 <- ggplot2::ggplot(df, ggplot2::aes(x = stats::reorder(Gene,!!sym(interaction)), y = !!sym(interaction))) +
                geom_bar(stat = "identity",
                show.legend = FALSE,
                fill = "blue",      # Background color
                color = "blue") + 
        theme(axis.text.x = element_text(size = metricText, angle=0), axis.text.y = element_text(angle = 0, size = geneText)) +
        xlab("Genes") +
        ylab("SpaceMarkers Metric") +
        ggtitle(interaction) + 
        scale_y_continuous(breaks= seq(0, max(df[[interaction]]),
            by = increments),
            limits = c(0, max(df[[interaction]]) + 0.2)) +
        coord_flip()
    if (!is.null(out)){
        ggplot2::ggsave(filename = out,plot = p1)
    }
    return(p1)
}