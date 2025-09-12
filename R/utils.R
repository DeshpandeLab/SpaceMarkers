#' @title calcOverlapScores
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
#' calcOverlapScores(hotspots = hotspots)   
#' @importFrom ggplot2 ggplot geom_tile geom_text theme_minimal
#' @importFrom reshape2 melt
#' @importFrom stats complete.cases
calcOverlapScores <- function(hotspots, patternList = NULL,
                             method = c("Szymkiewicz-Simpson",
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
  return(overlapScore)
}


#' @title getOverlapScores
#' @description Calculate the overlap scores between patterns in hotspots
#' @param hotspots A data frame with columns x, y, barcode and pattern names
#' @param in_hotspots Hotspots from the influence of patterns, Default: NULL
#' @param pat_hotspots Hotspots from the pattern itself, Default: NULL
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
#' in_hotspots <- data.frame(x = c(1,2,3,4,5),
#'                           y = c(1,2,3,4,5),
#'                           barcode = c("A","B","C","D","E"),
#'                           pattern1 = c(1,0,1,0,1),
#'                           pattern2 = c(1,1,0,0,1))
#' pat_hotspots <- data.frame(x = c(1,2,3,4,5),y = c(1,2,3,4,5),
#'                           barcode = c("A","B","C","D","E"),
#'                            pattern1 = c(0.3,1,0,0,1),
#'                            pattern3 = c(0,NA,1,0,NA))
#' getOverlapScores(hotspots = NULL, in_hotspots = in_hotspots, 
#' pat_hotspots = pat_hotspots, method = "absolute")
#'
#' @importFrom ggplot2 ggplot geom_tile geom_text theme_minimal
#' @importFrom reshape2 melt
#' @importFrom stats complete.cases

getOverlapScores <- function(hotspots,in_hotspots = NULL,pat_hotspots = NULL,
                             patternList = NULL, method = c("Szymkiewicz-Simpson",
                                                            "Jaccard", "Sorensen-Dice",
                                                            "Ochiai", "absolute") ) {

    if (!is.null(hotspots) & (is.null(in_hotspots) | is.null(pat_hotspots))) {
        # Set upper triangle and diagonal to NA for better visualization
      message("Please add pass in/pat_hotspots to hotspots. Assuming symmetric
      overlap scores. Setting upper triangle and diagonal to NA.")
      overlapScore <- calcOverlapScores(hotspots = hotspots,
                                        patternList =patternList,method=method)
      overlapScore[upper.tri(overlapScore,diag = TRUE)] <- NA
    }  else {
      patternList_in <- setdiff(colnames(in_hotspots),c("x","y","barcode"))
      rownames(in_hotspots) <- in_hotspots$barcode
      in_hotspots <- in_hotspots[,patternList_in]
      colnames(in_hotspots) <- paste0(patternList_in,"_in")
      
      patternList_pat <- setdiff(colnames(pat_hotspots),c("x","y","barcode"))
      rownames(pat_hotspots) <- pat_hotspots$barcode
      hotspots <- cbind(pat_hotspots[,patternList_pat],in_hotspots[rownames(pat_hotspots),])
      overlapScore <- calcOverlapScores(hotspots = hotspots,method=method)
    }
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
    rownames(imscores) <- imscores$Gene
    imscores$Gene <- NULL
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
    df$Gene <- rownames(df)
    df <- df[order(df[[interaction]], decreasing = TRUE),]
    df <- utils::head(df,nGenes)
    df[[interaction]] <- round(df[[interaction]],2)
    
    p1 <- ggplot2::ggplot(df, ggplot2::aes(x = stats::reorder(Gene,!!sym(interaction)), y = !!sym(interaction))) +
                geom_bar(stat = "identity",
                show.legend = FALSE,
                fill = "blue",      # Background color
                color = "blue") + 
        theme(axis.text.x = element_text(size = metricText, angle=0),
              axis.text.y = element_text(angle = 0, size = geneText)) +
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
#' @title plotScoresHeatmap
#' @description Generates top IMscores heatmap based on the maximum score per 
#' row where the number of rows is defined by geneCutOff
#' @param df A data frame with genes as rownames and IMscores as columns
#' @param interactions a vector of column names to subset the heatmap
#' @param geneCol a character specifying the column with the gene names.
#' Default: NULL and assumes rownames
#' @param geneCutOff Number of genes to plot based on max score, Default: 20
#' @param out Filename of heatmap, Default: 'ScoresHeatmap.png'
#' @param fontsize_row Size of row text, Default: 7
#' @param fontsize_col Size of column text N, Default: 7
#' @param cluster_cols Whether or not to cluster columns on the heatmap,
#' Default: TRUE
#' @param cluster_rows Whether or not to clutser rows on the heatmap,
#'  Default: TRUE
#' @param ... Additional parameters to pass to pheatmap::pheatmap
#' @return A Pheatmap object
#' @examples 
#' example(getPairwiseInteractingGenes)
#' imScores <- getIMScores(SpaceMarkers)
#' plotScoresHeatmap(imScores, geneCutOff = 15, out = "IMscoresHeatmap.png")
#' @export 
#' @importFrom ggplotify as.ggplot
#' @importFrom pheatmap pheatmap

plotScoresHeatmap <- function(df,
                              interactions = NULL,
                              geneCol = NULL,
                              geneCutOff = 20,
                              out         = "ScoresHeatmap.png",
                              fontsize_row = 7,
                              fontsize_col = 7,
                              cluster_cols = TRUE,
                              cluster_rows = TRUE,
                              ...) {
  # If the gene column is null, then assume it is already the rownames and 
  # move forward
  if (!is.null(geneCol)){
    rownames(df) <- df[[geneCol]]
    df[[geneCol]] <- NULL
  }
  mtx <- as.matrix(df)
  # optionally subset columns
  if (!is.null(interactions)) {
    mtx <- mtx[, interactions, drop = FALSE]
  }
  
  # order by top genes
  ## Get Max value by row
  mtx <- mtx[order(apply(mtx,1,max), decreasing = TRUE),]
  mtx <- head(mtx,geneCutOff)
  
  # draw heatmap and save
  ph <- ggplotify::as.ggplot(pheatmap::pheatmap(
    mtx,
    filename = out,
    cluster_cols = cluster_cols,
    cluster_rows = cluster_rows,
    fontsize_row = fontsize_row,
    fontsize_col = fontsize_col,
    ...
  ))
  
  return(ph)
}



#' @title getGseaPathways
#' @description An fgsea wrapper to get enriched pathways for IMscores
#' @param imscores A data frame with genes as rownames and IMscores as columns
#' @param gene.sets a list where the names are pathway names and the values are
#' vectors of gene names
#' @param ... params to pass to fgsea::fgsea
#' @return A list of data frames with enriched pathways for each IMscore column
#' @export 
#' @importFrom fgsea fgsea
#' @importFrom dplyr rename arrange
getGseaPathways<- function(imscores, gene.sets,...){
  patternNames <- setdiff(colnames(imscores),"Gene")
  features <- rownames(imscores)
  if (any(grepl("_",features))) {
    stop("Gene names cannot contain underscores (_).
         Please remove them and try again.")
  }
  d <- lapply(
    patternNames, FUN = function(p) {
      mscores <- imscores[,p]
      names(mscores) <- features
      mscores <- sort(mscores,decreasing = TRUE)
      result <- fgsea::fgsea(pathways = gene.sets, stats = mscores,...)
      if (nrow(result) == 0) {
        message(paste0("No pathways enriched for ",p))
        return(result)
      } else {
        result <- dplyr::rename(result, gene.set = "pathway")
        result <- dplyr::arrange(result,padj, -NES, na.rm = TRUE)
        result$leadingEdge <- vapply(result$leadingEdge, FUN = toString, FUN.VALUE = character(1))
        result$log10P_adj <- -log10(result$padj)
        result$interaction <- p
        return(result)
      }
    })
  return(d)
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
#' @param resolution Image resoultion to scale coordinates too, Default: c("lowres", "hires", "fullres")
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
#' @return a ggplot object
#
#' @export 
#' @importFrom dplyr mutate rename
#' @importFrom SpaceMarkers load10XCoords
#' @importFrom readbitmap read.bitmap
#' @importFrom ggplot2 ggplot annotation_raster geom_point aes coord_fixed scale_x_continuous scale_y_continuous labs theme_bw theme element_blank element_rect scale_fill_gradientn scale_fill_manual
#' @importFrom viridis scale_fill_viridis
#' @importFrom stats na.omit setNames
#' @importFrom RColorBrewer brewer.pal

plotSpatialDataOverImage <- function(
    visiumDir,
    df,
    feature_col,
    barcode_col = "barcode",
    resolution  = c("lowres", "hires", "fullres"),
    version     = NULL,
    colors      = NULL,
    point_size  = 2.5,
    stroke      = 0.05,
    alpha       = 0.5,
    title       = "Spatial Heatmap",
    bg_color    = NULL,
    crop        = TRUE
) {
  # --- deps ---
  requireNamespace("ggplot2", quietly = TRUE)
  requireNamespace("grid", quietly = TRUE)
  requireNamespace("readbitmap", quietly = TRUE)
  requireNamespace("viridis", quietly = TRUE)
  requireNamespace("RColorBrewer", quietly = TRUE)
  requireNamespace("dplyr", quietly = TRUE)
  
  resolution <- match.arg(resolution)
  spatial_dir <- file.path(visiumDir, "spatial")
  if (!dir.exists(spatial_dir)) stop("Could not find 'spatial' at: ", spatial_dir)
  if (!is.data.frame(df))       stop("`df` must be a data.frame.")
  if (!feature_col %in% names(df)) stop("`feature_col` not found in `df`.")
  
  # normalize barcode column name
  if (is.null(barcode_col)) {
    if (!is.null(rownames(df))) {
      df <- dplyr::mutate(df, barcode = rownames(df))
    } else stop("barcode_col is NULL and df has no rownames.")
  } else if (barcode_col != "barcode") {
    df <- dplyr::rename(df, barcode = !!barcode_col)
  }
  
  # ---- coordinates (ALREADY SCALED to chosen resolution) ----
  if (!requireNamespace("SpaceMarkers", quietly = TRUE)) {
    stop("Need SpaceMarkers::load10XCoords()")
  }
  pos <- SpaceMarkers::load10XCoords(visiumDir = visiumDir,
                                     resolution = resolution,
                                     version    = version)
  names(pos) <- c("barcode","y","x")
  pos$x <- as.numeric(pos$x); pos$y <- as.numeric(pos$y)
  
  # keep only barcodes present in coords
  df  <- merge(df[, c("barcode", feature_col)], pos, by = "barcode")
  
  # ---- pick the matching image file for the same resolution ----
  pick_image <- function(sp_dir, res) {
    files <- list.files(sp_dir, full.names = TRUE)
    low   <- grep("^tissue_.*lowres.*\\.(png|jpg|jpeg|tif|tiff)$",
                  basename(files), ignore.case = TRUE, value = TRUE)
    high  <- grep("^tissue_.*hires.*\\.(png|jpg|jpeg|tif|tiff)$",
                  basename(files), ignore.case = TRUE, value = TRUE)
    anyi  <- grep("\\.(png|jpg|jpeg|tif|tiff)$",
                  basename(files), ignore.case = TRUE, value = TRUE)
    
    want <- switch(res,
                   lowres = low,
                   hires  = high,
                   fullres = high  # best available; coords will be unscaled (scale=1)
    )
    if (!length(want)) want <- if (res == "lowres") high else low
    if (!length(want)) want <- anyi
    if (!length(want)) stop("No histology image found in ", sp_dir)
    file.path(sp_dir, want[[1]])
  }
  
  image_path <- pick_image(spatial_dir, resolution)
  img        <- readbitmap::read.bitmap(image_path)
  H <- nrow(img); W <- ncol(img)
  
  # ---- compute crop window in the SAME pixel system as df/coords ----
  if (crop) {
    # Keep everything in the SAME top-origin pixel space as load10XCoords()
    xmin_i <- max(1L, floor(min(df$x)))
    xmax_i <- min(W,  ceiling(max(df$x)))
    ymin_i <- max(1L, floor(min(df$y)))
    ymax_i <- min(H,  ceiling(max(df$y)))
    
    # Crop directly using top-origin row/col indices (NO H - y here)
    img_crop <- img[ymin_i:ymax_i, xmin_i:xmax_i, , drop = FALSE]
    
    # Cropped dimensions
    Wc <- xmax_i - xmin_i + 1L
    Hc <- ymax_i - ymin_i + 1L
    
    # Re-express points in the cropped frame (still top-origin)
    df2 <- df |>
      dplyr::mutate(
        x_c = x - (xmin_i - 1L),
        y_c = y - (ymin_i - 1L)
      )
    
    # Plot ranges and raster for ggplot (we'll flip y at draw time only)
    xlim_use <- c(0, Wc)
    ylim_use <- c(0, Hc)
    
    raster_xmin <- 0; raster_xmax <- Wc
    raster_ymin <- 0; raster_ymax <- Hc
    bg_grob <- as.raster(img_crop)
    
  } else {
    # full image (unchanged)
    df2 <- dplyr::mutate(df, x_c = x, y_c = y)
    xlim_use <- c(0, W)
    ylim_use <- c(0, H)
    raster_xmin <- 0; raster_xmax <- W
    raster_ymin <- 0; raster_ymax <- H
    bg_grob <- as.raster(img)
  }
  
  
  # ---- plot (annotation_raster + consistent axis) ----
  p <- ggplot2::ggplot() +
    ggplot2::annotation_raster(
      bg_grob,
      xmin = raster_xmin, xmax = raster_xmax,
      ymin = raster_ymin, ymax = raster_ymax
    ) +
    ggplot2::geom_point(
      data = df2,
      ggplot2::aes(x = x_c, y = (ylim_use[2] - y_c), fill = .data[[feature_col]]),
      shape  = 21, color = "black", size = point_size, stroke = stroke, alpha = alpha
    ) +
    ggplot2::coord_fixed(xlim = xlim_use, ylim = ylim_use, expand = FALSE) +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::labs(fill = feature_col, x = NULL, y = NULL, title = title) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank(),
      axis.text  = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      plot.background = if (!is.null(bg_color)) ggplot2::element_rect(fill = bg_color, color = NA) else ggplot2::element_blank()
    )
  
  # color scale
  if (is.numeric(df[[feature_col]])) {
    p <- if (!is.null(colors)) p + ggplot2::scale_fill_gradientn(colours = colors)
    else                   p + viridis::scale_fill_viridis()
  } else {
    vals <- unique(stats::na.omit(df[[feature_col]]))
    if (!is.null(colors)) {
      p <- p + ggplot2::scale_fill_manual(values = colors)
    } else if (length(vals) > 1) {
      pal <- stats::setNames(
        RColorBrewer::brewer.pal(max(3, min(length(vals), 9)), "Set1")[seq_along(vals)],
        vals
      )
      p <- p + ggplot2::scale_fill_manual(values = pal)
    } else {
      p <- p + ggplot2::scale_fill_manual(values = "red")
    }
  }
  return(p)
}

