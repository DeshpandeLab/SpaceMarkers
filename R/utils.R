
#' @title getOverlapScores
#' @description Calculate the overlap scores between patterns in hotspots
#' @param hotspots A data frame with columns x, y, barcode and pattern names
#' @param patternList A character vector of pattern names to calculate overlap 
#' scores for
#' @param method The method to calculate overlap scores. Options are
#' "Szymkiewicz-Simpson", "Jaccard", "Sorensen-Dice", "Ochiai" and "absolute"
#' @details The function calculates the overlap scores between patterns hotspots
#' using the specified method. The default method is "Szymkiewicz-Simpson"
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
                             patternList = NULL, method = c("Szymkiewicz-Simpson",
                                                            "Jaccard", "Sorensen-Dice",
                                                            "Ochiai", "absolute") ) {
    
    #stop if more than one method is supplied, do not warn by default
    if(length(method) > 1){
        method <- method[1]
        message("Only one method can be used at a time. Using ", method)}

    if (is.null(patternList)) {
        patternList <- setdiff(colnames(hotspots),c("x","y","barcode"))
    } else if (!all(patternList %in% colnames(hotspots))) {
        stop("Pattern names not found in hotspots")
    }

    binarized <- (!is.na(hotspots[,patternList]))*1
    intersects <- t(binarized) %*% binarized
    nHotspots <- colSums(binarized)
    nHotsP1 <- t(t(nHotspots)) %*% array(1, length(patternList))
    nHotsP2 <- t(nHotsP1)
    overlapScore <- switch(method,
        "Szymkiewicz-Simpson" = intersects/pmin(nHotsP1,nHotsP2),
        "Jaccard" = intersects/(nHotsP1 + nHotsP2 - intersects),
        "Sorensen-Dice" = 2*intersects/(nHotsP1 + nHotsP2),
        "Ochiai" = intersects/sqrt(nHotsP1*nHotsP2),
        "absolute" = intersects,
        stop("Method not supported")
    )

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
#' df <- data.frame(pattern1 = c("pattern1","pattern1","pattern2","pattern2"), 
#'                  pattern2 = c("pattern1","pattern2","pattern1","pattern2"),
#'                  overlapScore = c(0.5,0.7,0.3,0.9))
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

#' @title getGeneSetScore
#' @description Calculate the mean interaction score for a set of genes
#' @param IMscores A matrix of interaction scores
#' @param genes A list of gene sets, where each set is a vector of gene names
#' @return A vector of mean interaction scores for each gene set
#' @export
getGeneSetScore <- function(IMscores, genes) {
    scores <- sapply(genes, function(geneSet) {
        if (length(geneSet) == 0) return(NA)
        geneSet <- intersect(geneSet, rownames(IMscores))
        if (length(geneSet) == 0) return(NA)
        if (length(geneSet) == 1) {
            return(IMscores[geneSet, , drop = FALSE])
        } else
            if (length(geneSet) > 1) {
                return(colMeans(IMscores[geneSet, , drop = FALSE], na.rm = TRUE))
            }
    })
    scores <- Reduce(rbind, scores)
    return(scores)
}


#' @title pickImage
#' @description The function picks the appropriate histology image file from the
#' spatial directory based on the specified resolution.
#' @param sp_dir path to the spatial directory
#' @param res a character string specifying the resolution of the image
#' @return a character string of the image file name
pickImage <- function(sp_dir, res) {
  files <- list.files(sp_dir, full.names = TRUE)
  low   <- grep("^tissue_.*lowres.*\\.(png|jpg|jpeg|tif|tiff)$",
                basename(files), ignore.case = TRUE, value = TRUE)
  high  <- grep("^tissue_.*hires.*\\.(png|jpg|jpeg|tif|tiff)$",
                basename(files), ignore.case = TRUE, value = TRUE)
  anyi  <- grep("\\.(png|jpg|jpeg|tif|tiff)$",
                basename(files), ignore.case = TRUE, value = TRUE)
  want <- switch(res,lowres = low,hires  = high,fullres = high)
  if (!length(want)) want <- if (res == "lowres") high else low
  if (!length(want)) want <- anyi
  if (!length(want)) stop("No histology image found in ", sp_dir)
  return(file.path(sp_dir, want[[1]]))
}

#' @title plotSpatialDataOverImage
#' @description This function plots spatial data over the complementary 
#' histology image of varing resolutions
#' @param visiumDir directory with a spatial folder containing 
#' scalefactors_json.json, images (lowres,or hires), and 
#' coordinates (tissue_positons_(list).csv or probe.csv)
#' @param df a dataframe with the features of interest. For example,
#' can behotspots (character), and/or influence (numeric))
#' @param feature_col feature to plot over spots, Default: NULL
#' @param barcode_col barcode column name to match with coordinates,
#' Default: 'barcode'
#' @param resolution Image resoultion to scale coordinates too, 
#' Default: c("lowres", "hires", "fullres")
#' @param version Visium version. Automatically infers from load10XCoords if 
#' NULL,Default: NULL
#' @param colors colors to be displayed over spots. If set to NULL, it
#' automatically colors the spots red for character values and uses viridis for 
#' numeric values. Default: NULL
#' @param point_size size of spots displayed on the plot, Default: 2.5
#' @param stroke thickness of spot outline, Default: 0.05
#' @param alpha Transparency of the spots, Default: 0.5
#' @param title Title displayed on the plot, Default: 'Spatial Heatmap'
#' @param bg_color background color of ggplot box, Default: NULL
#' @param crop crop spatial plot to a zoomed in window, Default: TRUE
#' @param text_size size of text on the plot, Default: 15
#' @return a ggplot object
#
#' @export 
#' @importFrom rlang sym
#' @importFrom SpaceMarkers load10XCoords
#' @importFrom readbitmap read.bitmap
#' @importFrom viridis scale_fill_viridis
#' @importFrom stats na.omit setNames
#' @importFrom RColorBrewer brewer.pal

plotSpatialDataOverImage <- function(
    visiumDir, df, feature_col, barcode_col="barcode",
    resolution=c("lowres","hires","fullres"), version=NULL, colors=NULL,
    point_size=2.5, stroke=0.05, alpha=0.5, title="Spatial Heatmap",
    bg_color=NULL, crop=TRUE, text_size = 15) {
  gg <- getNamespace("ggplot2"); dp <- getNamespace("dplyr")
  resolution <- match.arg(resolution)
  if (is.null(barcode_col)) {
    if (!is.null(rownames(df))) df$barcode <- rownames(df) else
      stop("barcode_col is NULL and df has no rownames.")
  } else if (barcode_col != "barcode") df <- dp$rename(
    df, barcode = !!rlang::sym(barcode_col))
  pos <- SpaceMarkers::load10XCoords(visiumDir, resolution, version)
  names(pos) <- c("barcode","y","x"); df <- merge(
    df[, c("barcode", feature_col)], pos, by="barcode")
  img <- readbitmap::read.bitmap(pickImage(
    file.path(visiumDir, "spatial"), resolution))
  if (crop) {xr <- range(df$x, na.rm=TRUE); yr <- range(df$y, na.rm=TRUE)
  xmin <- max(1L, floor(xr[1])); xmax <- min(ncol(img), ceiling(xr[2]))
  ymin <- max(1L, floor(yr[1])); ymax <- min(nrow(img), ceiling(yr[2]))
  img <- img[ymin:ymax, xmin:xmax,, drop=FALSE]
  df <- dp$mutate(df, x_c=x-xmin+1L, y_c=y-ymin+1L)
  xl <- c(0, xmax-xmin+1L); yl <- c(0, ymax-ymin+1L)
  } else {
    df <- dp$mutate(df, x_c=x, y_c=y);xl<-c(0, ncol(img));yl<-c(0, nrow(img))}
    p <- gg$ggplot() +
      gg$annotation_raster(as.raster(img), 0, diff(xl), 0, diff(yl)) +
      gg$geom_point(
        data = df,
        gg$aes(x_c, yl[2] - y_c, fill = .data[[feature_col]]),
        shape = 21, color = "black", size = point_size,
        stroke = stroke, alpha = alpha
      ) +
      gg$coord_fixed(xlim = xl, ylim = yl, expand = FALSE) +
      gg$labs(fill = feature_col, x = NULL, y = NULL, title = title) +
      gg$theme_bw(text_size) +
      gg$theme(
        plot.background = if (!is.null(bg_color))
          gg$element_rect(fill = bg_color, color = NA)
        else gg$element_blank()
      )
  if (is.numeric(df[[feature_col]])) {
    p <- p + (if (!is.null(colors)) gg$scale_fill_gradientn(colours=colors)
              else viridis::scale_fill_viridis())
  } else {
    vals <- unique(stats::na.omit(df[[feature_col]]))
    p <- p + if (!is.null(colors)) gg$scale_fill_manual(values=colors) else
      if (length(vals) > 1) gg$scale_fill_manual(values = stats::setNames(
        RColorBrewer::brewer.pal(max(3, min(length(vals), 9)), "Set1")
        [seq_along(vals)], vals)) else gg$scale_fill_manual(values="red")
  }
  return(p)
}

