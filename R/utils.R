#' @title calc_overlap_scores
#' @description Calculate the overlap scores between patterns in hotspots
#' @param hotspots A data frame with columns x, y, barcode and pattern names
#' @param patternList A character vector of pattern names to calculate overlap 
#' scores for
#' @param method The method to calculate overlap scores. Options are
#' "Szymkiewicz-Simpson", "Jaccard", "Sorensen-Dice", "Ochiai" and "absolute"
#' @details The function calculates the overlap scores between patterns hotspots
#' using the specified method. The default method is "Szymkiewicz-Simpson"
#' overlap coefficient.
#' @return A matrix of overlap scores (patterns × patterns)
#' @export
#' @examples
#' hotspots <- data.frame(x = c(1,2,3,4,5),
#'                         y = c(1,2,3,4,5),
#'                         barcode = c("A","B","C","D","E"),
#'                         pattern1 = c(1,0,1,0,1),
#'                         pattern2 = c(1,1,0,0,1))
#' calc_overlap_scores(hotspots = hotspots)   
#' @importFrom ggplot2 ggplot geom_tile geom_text theme_minimal
#' @importFrom reshape2 melt
#' @importFrom stats complete.cases
calc_overlap_scores <- function(hotspots, patternList = NULL,
                                method = c("Szymkiewicz-Simpson",
                                           "Jaccard", "Sorensen-Dice",
                                           "Ochiai", "absolute") ) {
  # stop if more than one method is supplied, do not warn by default
  if(length(method) > 1){
    method <- method[1]
    message("Only one method can be used at a time. Using ", method)
  }
  
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
  overlapScore <- switch(
    method,
    "Szymkiewicz-Simpson" = intersects/pmin(nHotsP1,nHotsP2),
    "Jaccard"             = intersects/(nHotsP1 + nHotsP2 - intersects),
    "Sorensen-Dice"       = 2*intersects/(nHotsP1 + nHotsP2),
    "Ochiai"              = intersects/sqrt(nHotsP1*nHotsP2),
    "absolute"            = intersects,
    stop("Method not supported")
  )
  return(overlapScore)
}


#' @title get_overlap_scores
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
#' get_overlap_scores(hotspots = NULL, in_hotspots = in_hotspots, 
#'                    pat_hotspots = pat_hotspots, method = "absolute")
#'
#' @importFrom ggplot2 ggplot geom_tile geom_text theme_minimal
#' @importFrom reshape2 melt
#' @importFrom stats complete.cases
get_overlap_scores <- function(hotspots = NULL, in_hotspots = NULL, pat_hotspots = NULL,
                               patternList = NULL,
                               method = c("Szymkiewicz-Simpson", "Jaccard", "Sorensen-Dice",
                                          "Ochiai", "absolute")) {
  
  if (!is.null(hotspots) & (is.null(in_hotspots) | is.null(pat_hotspots))) {
    message("Assuming symmetric overlap scores. Setting upper triangle and diagonal to NA.")
    overlapScore <- calc_overlap_scores(hotspots = hotspots, patternList = patternList, method = method)
    overlapScore[upper.tri(overlapScore, diag = TRUE)] <- NA
    
    dfOverlap <- reshape2::melt(overlapScore)
    dfOverlap <- dfOverlap[stats::complete.cases(dfOverlap), ]
    colnames(dfOverlap) <- c("pattern2", "pattern1", "overlapScore")
    dfOverlap <- dfOverlap[, c(2, 1, 3)]
    
  } else {
    patternList_in <- setdiff(colnames(in_hotspots), c("x", "y", "barcode"))
    rownames(in_hotspots) <- in_hotspots$barcode
    in_hotspots <- in_hotspots[, patternList_in, drop = FALSE]
    # prefix influence columns with "near_"
    colnames(in_hotspots) <- paste0("near_", patternList_in)
    
    patternList_pat <- setdiff(colnames(pat_hotspots), c("x", "y", "barcode"))
    rownames(pat_hotspots) <- pat_hotspots$barcode
    
    hotspots <- cbind(
      pat_hotspots[, patternList_pat, drop = FALSE],
      in_hotspots[rownames(pat_hotspots), , drop = FALSE]
    )
    
    overlapScore <- calc_overlap_scores(hotspots = hotspots, method = method)
    
    dfOverlap <- reshape2::melt(overlapScore)
    dfOverlap <- dfOverlap[stats::complete.cases(dfOverlap), ]
    colnames(dfOverlap) <- c("pattern2", "pattern1", "overlapScore")
    dfOverlap <- dfOverlap[, c(2, 1, 3)]
    
    # filter using "near_" instead of "_in"
    dfOverlap <- dfOverlap[!grepl("^near_", dfOverlap$pattern1), ]
    dfOverlap <- dfOverlap[grepl("^near_",  dfOverlap$pattern2), ]
  }
  
  return(dfOverlap)
}

#' @title plot_overlap_scores
#' @description Plot the overlap scores between patterns in hotspots
#' @param df A data frame with columns pattern1, pattern2 and overlapScore
#' @param title The title of the plot
#' @param fontsize The font size of the plot
#' @param tilefontsize The font size of the tile labels
#' @param out The output path for the plot
#' @param influence If TRUE, hides self "pattern vs near_pattern" cells and
#'   strips the "near_" prefix from y-axis labels for readability.
#' @return A ggplot object
#' @export
#' @examples
#' df <- data.frame(pattern1 = c("pattern1","pattern1","pattern2","pattern2"),
#'                  pattern2 = c("near_pattern1","near_pattern2","near_pattern1","near_pattern2"),
#'                  overlapScore = c(0.5,0.7,0.3,0.9))
#' plot_overlap_scores(df, influence = TRUE)
#' @import ggplot2
plot_overlap_scores <- function(df, title = "Spatial Overlap Scores",
                                out = NULL, fontsize = 15, tilefontsize = 6,
                                influence = FALSE) {
  if (isTRUE(influence)) {
    p1 <- as.character(df$pattern1)
    p2 <- as.character(df$pattern2)
    # mask "pattern1 vs near_pattern1"
    mask <- p2 == paste0("near_", p1)
    df$overlapScore[mask] <- NA_real_
    df <- df[!is.na(df$overlapScore), ]
  }
  
  p <- ggplot2::ggplot(data = df, ggplot2::aes(pattern1, pattern2, fill = overlapScore)) +
    ggplot2::geom_tile(color = "black", linewidth = 0.8) +
    ggplot2::geom_text(ggplot2::aes(label = round(overlapScore, 2)), size = tilefontsize) +
    ggplot2::scale_fill_gradient2(low = "#FFF7EC", mid = "#FDBB84", high = "#D7301F", midpoint = 0.5) +
    ggplot2::scale_y_discrete(
      limits = rev,
      guide = ggplot2::guide_axis(angle = 45),
      labels = function(x) if (isTRUE(influence)) sub("^near_", "", x) else x
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, size = fontsize, hjust = 1),
      axis.text.y = ggplot2::element_text(size = fontsize),
      axis.title.x = if (isTRUE(influence)) ggplot2::element_text(size = fontsize) else ggplot2::element_blank(),
      axis.title.y = if (isTRUE(influence)) ggplot2::element_text(size = fontsize) else ggplot2::element_blank()
    ) +
    ggplot2::coord_fixed() +
    ggplot2::labs(
      x = if (isTRUE(influence)) "Pattern_Hotspots" else NULL,
      y = if (isTRUE(influence)) "Influence_Hotspots" else NULL,
      fill = "overlapScore",
      title = title
    ) +
    ggplot2::theme(
      panel.background = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank()
    )
  
  if (!is.null(out)) {
    ggplot2::ggsave(filename = out, plot = p)
  }
  return(p)
}


#' @title get_im_scores
#' @description Get the interaction scores for SpaceMarkers
#' @param SpaceMarkers A list of SpaceMarkers objects
#' @return A data frame with columns Gene and SpaceMarkersMetric
#' @export
#' @examples
#' example(get_pairwise_interacting_genes)
#' get_im_scores(SpaceMarkers)
#' @importFrom stats setNames
#' 
get_im_scores <- function(SpaceMarkers){
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


#' @title plot_im_scores_heatmap
#' @description A pheatmap wrapper to plot the top IMscores based on the maximum score
#'   per row, where the number of rows is defined by \code{geneCutOff}.
#'
#' @param df A data.frame with genes as rownames and IMscores as columns. If
#'   \code{geneCol} is not \code{NULL}, that column is used as gene names.
#' @param interactions A character vector of column names to subset the
#'   heatmap to specific interactions. Default: \code{NULL} (use all columns).
#' @param geneCol Character specifying the column with the gene names.
#'   Default: \code{"Gene"}. If \code{NULL}, assumes rownames already hold
#'   gene names.
#' @param geneCutOff Integer. Number of genes (rows) to plot, based on
#'   descending maximum score per gene. Default: \code{20}.
#' @param out Optional filename for saving the heatmap directly via
#'   \code{pheatmap}. Default: \code{"ScoresHeatmap.png"}. Set to \code{NULL}
#'   to avoid automatic saving.
#' @param fontsize_row Numeric. Size of row text. Default: \code{7}.
#' @param fontsize_col Numeric. Size of column text. Default: \code{7}.
#' @param cluster_cols Logical. Whether to cluster columns. Default: \code{TRUE}.
#' @param cluster_rows Logical. Whether to cluster rows. Default: \code{TRUE}.
#' @param ... Additional parameters passed to \code{pheatmap::pheatmap()}.
#'
#' @return A \code{pheatmap} object (invisibly). To draw it, use \code{print(ph)}.
#'
#' @examples
#' \dontrun{
#' # Basic usage with automatic saving (PNG)
#' ph <- plot_scores_heatmap(imScores, geneCutOff = 15)
#'
#' # No automatic file, then save a custom TIFF via tiff()/dev.off():
#' ph <- plot_scores_heatmap(imScores, geneCutOff = 15, out = NULL)
#' tiff("IMscoresHeatmap.tiff", width = 7, height = 7, units = "in", res = 300)
#' print(ph)
#' dev.off()
#' }
#'
#' @importFrom pheatmap pheatmap
#' @export
plot_im_scores_heatmap <- function(df,
                                interactions  = NULL,
                                geneCol       = "Gene",
                                geneCutOff    = 20,
                                out           = "ScoresHeatmap.png",
                                fontsize_row  = 7,
                                fontsize_col  = 7,
                                cluster_cols  = TRUE,
                                cluster_rows  = TRUE,
                                ...) {
  # If the gene column is not NULL, use it as rownames and drop the column
  if (!is.null(geneCol)) {
    if (!geneCol %in% colnames(df)) {
      stop("geneCol '", geneCol, "' not found in df.")
    }
    rownames(df) <- df[[geneCol]]
    df[[geneCol]] <- NULL
  }
  mtx <- as.matrix(df)
  # optionally subset columns
  if (!is.null(interactions)) {
    missing_cols <- setdiff(interactions, colnames(mtx))
    if (length(missing_cols)) {
      stop("The following interactions are not columns in df: ",
           paste(missing_cols, collapse = ", "))
    }
    mtx <- mtx[, interactions, drop = FALSE]
  }
  # order by top genes (max score per row)
  if (nrow(mtx) == 0L) stop("No rows to plot after subsetting.")
  mtx <- mtx[order(apply(mtx, 1, max, na.rm = TRUE), decreasing = TRUE), , drop = FALSE]
  mtx <- head(mtx, geneCutOff)
  # build argument list for pheatmap
  args <- list(
    mat           = mtx,
    cluster_cols  = cluster_cols,
    cluster_rows  = cluster_rows,
    fontsize_row  = fontsize_row,
    fontsize_col  = fontsize_col
  )
  if (!is.null(out)) {
    args$filename <- out
  } 
  ph <- do.call(pheatmap::pheatmap, c(args, list(...)))
  invisible(ph)
}
                                     


#' @title plot_im_scores
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
#' example(get_pairwise_interacting_genes)
#' plot_im_scores(get_im_scores(SpaceMarkers), "Pattern_1_Pattern_3")
#' @import ggplot2
#' @importFrom stats reorder
#' @importFrom utils head
plot_im_scores <- function(df, interaction, cutOff = 0, nGenes = 20,
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

#' @title get_enriched_pathways
#' @description An fgsea wrapper to get enriched pathways for IMscores
#' @param imscores A data frame or matrix with genes as rownames and IMscores as columns
#' @param gene.sets A list where the names are pathway names and the values are
#' vectors of gene names
#' @param method Character, one of \code{c("undirected","directed")}.
#' @param threshold Numeric FDR threshold (currently unused in filtering).
#' @param interactionNames Optional character vector of interaction names to test.
#'   If NULL, uses all IMscore columns except "Gene".
#' @param ... Additional parameters to pass to \code{fgsea::fgsea}
#' @return A list of data frames with enriched pathways for each IMscore column
#' @export 
#' @importFrom fgsea fgsea ora
#' @importFrom dplyr rename arrange
get_enriched_pathways <- function(imscores, gene.sets,method = c(
  "undirected","directed"),threshold = 0.05,interactionNames = NULL,
  ...){
  
  if (method == "undirected"){
    if (!is.null(interactionNames)){
      interactionNames <- interactionNames
    } else {
      interactionNames <- setdiff(colnames(imscores),"Gene")
    }
    features <- rownames(imscores)
    if (any(grepl("_",features))) {
      stop("Gene names cannot contain underscores (_).
         Please remove them and try again.")
    }
    d <- lapply(
      interactionNames, FUN = function(p) {
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
    
  } else if (method == "directed"){
    imscores <- unique(imscores[, c("gene","cell_interaction","effect_size")])
    imscores$effect_size[is.na(imscores$effect_size)] <- 0
    imscores <- reshape2::acast(
      imscores,
      gene ~ cell_interaction,
      value.var = "effect_size"
    )
    
    features <- rownames(imscores)
    if (any(grepl("_",features))) {
      stop("Gene names cannot contain underscores (_).
         Please remove them and try again.")
    }
    if (!is.null(interactionNames)){
      interactionNames <- interactionNames
    } else {
      interactionNames <- setdiff(colnames(imscores),"Gene")
    }
    
    d <- lapply(
      interactionNames, FUN = function(p) {
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
  }
  
  return(d)
}

                                     
#' Dotplot of enriched pathways across interactions
#'
#' @param enr_list A list of data.frames from \code{get_enriched_pathways()}.
#'   Must include: gene.set, interaction, leadingEdge, log10P_adj.
#' @param top_n Integer. Keep the top_n pathways *within each interaction*
#'   by log10P_adj; the plot shows the union. Default: 5.
#' @return A ggplot object.
plot_enriched_results <- function(enr_list, top_n = 5) {
  df <- Filter(NROW, enr_list)
  if (length(df) == 0L) stop("No enrichment results to plot.")
  df <- dplyr::bind_rows(df)
  req <- c("gene.set","interaction","leadingEdge","log10P_adj")
  miss <- setdiff(req, names(df))
  if (length(miss)) stop("Missing columns: ", paste(miss, collapse = ", "))
  
  n_lead <- function(x) ifelse(is.na(x) | trimws(x) == "", 0L, lengths(strsplit(x, "\\s*,\\s*")))
  df <- df |>
    dplyr::mutate(
      n_leading = n_lead(.data$leadingEdge),
      log10P_adj = ifelse(is.infinite(.data$log10P_adj),
                          max(.data$log10P_adj[is.finite(.data$log10P_adj)], 0, na.rm = TRUE),
                          .data$log10P_adj)
    )
  
  keep_sets <- df |>
    dplyr::group_by(.data$interaction) |>
    dplyr::filter(!is.na(.data$log10P_adj)) |>
    dplyr::slice_max(order_by = .data$log10P_adj, n = top_n, with_ties = TRUE) |>
    dplyr::ungroup() |>
    dplyr::pull(.data$gene.set) |>
    unique()
  df <- df[df$gene.set %in% keep_sets, , drop = FALSE]
  
  y_levels <- df |>
    dplyr::group_by(.data$gene.set) |>
    dplyr::summarise(max_log10 = max(.data$log10P_adj, na.rm = TRUE), .groups = "drop") |>
    dplyr::arrange(dplyr::desc(.data$max_log10)) |>
    dplyr::pull(.data$gene.set)
  
  df$gene.set    <- factor(df$gene.set, levels = y_levels)
  df$interaction <- factor(df$interaction)
  
  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(x = .data$interaction,
                 y = .data$gene.set,
                 size = .data$n_leading,
                 fill = .data$log10P_adj)
  ) +
    ggplot2::geom_point(shape = 21, stroke = 0.2) +
    ggplot2::scale_size_area(max_size = 9, name = "Leading Edge") +
    ggplot2::scale_fill_gradient(low = "blue", high = "red", name = "-log10(FDR)") +
    ggplot2::labs(x = "Interaction", y = "Pathway") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid.major.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      legend.position = "right"
    )
  return(p)
}

#' Average & detection by group for a gene set (NA-safe on meta_col)
#'
#' @param counts   Matrix or dgCMatrix: genes (rows) × cells (cols).
#' @param metadata Data frame with per-cell metadata. Must include `meta_col`
#'                 and either have barcodes as rownames or a column named in `barcode_col`.
#' @param features Character vector of gene symbols to include.
#' @param meta_col Name of the metadata column used to group cells (e.g., "celltype").
#' @param barcode_col Optional column name in `metadata` holding cell barcodes.
#'                    If NULL, rownames(metadata) are used as barcodes.
#'
#' @return data.frame with columns: {meta_col}, gene_symbol, pct.expr, avg.expr
#'
#' @examples
#' \dontrun{
#' # Example: average expression by SpaceMarkers interaction region
#'
#' # Suppose:
#' #   - spCounts is a genes x barcodes count matrix
#' #   - spHotspots_all is a metadata data.frame with a "barcode" column
#' #     and an interaction-region column "Epi_Plasma"
#' #   - df is a data.frame with a column "Gene" containing gene symbols of interest
#'
#' meta_col <- "Epi_Plasma"
#'
#' exp_df <- create_avg_exp_df(
#'   counts      = spCounts,
#'   metadata    = spHotspots_all,
#'   features    = df$Gene,
#'   meta_col    = meta_col,
#'   barcode_col = "barcode"
#' )
#'
#' head(exp_df)
#' }
create_avg_exp_df <- function(counts, metadata, features, meta_col, barcode_col = NULL) {
  # --- Validate inputs ---
  if (!is.matrix(counts) && !inherits(counts, "dgCMatrix")) {
    stop("`counts` must be a matrix or dgCMatrix with genes as rows and cells as columns.")
  }
  if (!is.data.frame(metadata)) stop("`metadata` must be a data.frame.")
  if (!meta_col %in% colnames(metadata)) {
    stop(sprintf("`meta_col` ('%s') not found in metadata.", meta_col))
  }
  if (is.null(features) || !length(features)) stop("`features` must be a non-empty character vector.")
  features <- unique(as.character(features))
  
  # --- Identify barcodes in metadata ---
  if (is.null(barcode_col)) {
    if (is.null(rownames(metadata))) {
      stop("`barcode_col` is NULL and metadata has no rownames; provide `barcode_col`.")
    }
    md_barcodes <- rownames(metadata)
  } else {
    if (!barcode_col %in% names(metadata)) {
      stop(sprintf("`barcode_col` ('%s') not found in metadata.", barcode_col))
    }
    md_barcodes <- as.character(metadata[[barcode_col]])
    rownames(metadata) <- md_barcodes
  }
  
  # --- Initial overlap & alignment by barcode ---
  cell_barcodes <- colnames(counts)
  if (is.null(cell_barcodes)) stop("`counts` must have column names (cell barcodes).")
  common_bc <- intersect(cell_barcodes, md_barcodes)
  if (!length(common_bc)) stop("No overlapping barcodes between counts and metadata.")
  metadata <- metadata[common_bc, , drop = FALSE]
  counts   <- counts[, common_bc, drop = FALSE]
  
  # --- NA handling ---
  # 1) Drop metadata rows ONLY when meta_col is NA
  keep_meta <- !is.na(metadata[[meta_col]])
  if (any(!keep_meta)) {
    metadata <- metadata[keep_meta, , drop = FALSE]
    counts   <- counts[, rownames(metadata), drop = FALSE]
  }
  # 2) Drop counts columns that contain any NA
  has_na_col <- colSums(is.na(counts)) > 0
  if (any(has_na_col)) {
    counts   <- counts[, !has_na_col, drop = FALSE]
    metadata <- metadata[colnames(counts), , drop = FALSE]
  }
  # Final sanity check
  if (ncol(counts) == 0L) stop("All cells removed after NA filtering.")
  if (!identical(colnames(counts), rownames(metadata))) {
    stop("Post-NA filtering alignment failed; barcodes are misaligned.")
  }
  
  # --- Subset genes (features) ---
  keep_features <- intersect(features, rownames(counts))
  if (!length(keep_features)) stop("None of the requested `features` are present in `counts` rownames.")
  if (length(keep_features) < length(features)) {
    missing <- setdiff(features, keep_features)
    warning(sprintf("Dropping %d missing feature(s): %s",
                    length(missing), paste(head(missing, 10), collapse = ", ")),
            call. = FALSE)
  }
  counts_sub <- counts[keep_features, , drop = FALSE]
  
  # --- Transpose to cells × genes and attach group column ---
  if (inherits(counts_sub, "dgCMatrix")) {
    if (requireNamespace("Matrix", quietly = TRUE)) {
      counts_t <- as.matrix(Matrix::t(counts_sub))
    } else {
      counts_t <- t(as.matrix(counts_sub))
    }
  } else {
    counts_t <- t(counts_sub)
  }
  df_expr <- data.frame(counts_t, check.names = FALSE)
  stopifnot(identical(rownames(df_expr), rownames(metadata)))
  df_expr[[meta_col]] <- as.character(metadata[[meta_col]])
  
  # --- Wide -> long (reshape2), then aggregate ---
  long_df <- reshape2::melt(
    df_expr,
    id.vars       = meta_col,
    variable.name = "gene_symbol",
    value.name    = "expr"
  )
  
  out <- long_df |>
    dplyr::group_by(.data[[meta_col]], .data$gene_symbol) |>
    dplyr::summarise(
      pct.expr = sum(.data$expr > 0, na.rm = TRUE) / dplyr::n() * 100,
      avg.expr = mean(.data$expr, na.rm = TRUE),
      .groups  = "drop"
    ) |>
    as.data.frame()
  
  names(out)[names(out) == meta_col] <- meta_col
  return(out)
} 

#' Dot plot of gene expression summary across interaction groups
#'
#' @description
#' Create a Seurat-style dot plot where:
#' - the x-axis represents interaction / group levels (e.g. `"Stom"`, `"Interacting"`, `"B"`),
#' - the y-axis represents genes (or any feature),
#' - dot size encodes percentage of cells expressing the gene,
#' - dot fill encodes average expression.
#'
#' The function can optionally cluster rows/columns using hierarchical clustering
#' on the average-expression matrix, or respect a user-specified ordering. If
#' an `"Interacting"` level is present on the x-axis, it is automatically
#' centered between the other two levels (if there are exactly two others).
#'
#' @param meta_df Data frame containing per-group gene summaries. Must include
#'   at least the columns:
#'   - `meta_col` (default `NULL`),
#'   - `y_var` (default `"gene_symbol"`),
#'   - `"pct.expr"` and `"avg.expr"`.
#' @param meta_col Character, metadata column representing the main grouping
#'   (usually the same as `x_var`). Default: `NULL`.
#' @param x_var Character, column to use on the x-axis (group labels).
#'   Default: `meta_col`.
#' @param y_var Character, column to use on the y-axis (gene labels).
#'   Default: `"gene_symbol"`.
#' @param title,xlab,ylab Plot title and axis labels.
#' @param x_text_angle,x_text_hjust,x_text_size Formatting of x-axis text.
#' @param y_text_size,title_size Font sizes for y-axis text and title.
#' @param size_name,fill_name Legend titles for size (percentage) and fill
#'   (average expression).
#' @param fill_colors Optional vector of colors for a continuous gradient.
#'   If `NULL`, a white→red gradient is generated with \pkg{circlize}.
#' @param size_var,fill_var Column names in `meta_df` for size and fill
#'   aesthetics. Defaults: `"pct.expr"` and `"avg.expr"`.
#' @param order Character; ordering strategy for rows/cols. One of
#'   `"both"`, `"row"`, `"column"`, or `NULL`. If `NULL`, user must supply
#'   both `user_order_rows` and `user_order_cols`.
#' @param user_order_rows,user_order_cols Optional explicit orderings for the
#'   y- and x-axis when `order = NULL` or partially overridden.
#' @param dot_shape,border_color,stroke Aesthetics for the points.
#' @param clust_method Clustering method passed to \code{hclust} when rows/cols
#'   are ordered by hierarchical clustering. Default: `"ward.D2"`.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' \dontrun{
#' # Suppose `avg_df` was created by create_avg_exp_df() and then subset:
#' # avg_df has columns: NULL, gene_symbol, pct.expr, avg.expr
#'
#' p <- plot_interaction_dotplot(
#'   meta_df = avg_df,
#'   meta_col = NULL,
#'   x_var    = NULL,
#'   y_var    = "gene_symbol",
#'   title    = "Senescence gene panel by subClone",
#'   order    = "both"
#' )
#' print(p)
#' }
#'
#' @importFrom ggplot2 ggplot aes_string geom_point scale_size_continuous
#'   scale_fill_gradientn scale_x_discrete scale_y_discrete theme_classic
#'   labs theme element_text guides guide_legend
#' @importFrom stats dist hclust as.formula
#' @importFrom reshape2 acast
#' @export
plot_interaction_dotplot <- function(
    meta_df,
    meta_col = NULL,
    x_var = meta_col,              # genes on y; groups on x
    y_var = "gene_symbol",
    title = "Dot Plot",
    xlab = meta_col,
    ylab = "Gene Symbol",
    x_text_angle = 90,
    x_text_hjust = 1,
    x_text_size = 12,
    y_text_size = 12,
    title_size = 14,
    size_name = "Percentage Expressed",
    fill_name = "Avg Expression",
    fill_colors = NULL,
    size_var = "pct.expr",
    fill_var = "avg.expr",
    order = "both",
    user_order_rows = NULL,
    user_order_cols = NULL,
    dot_shape = 21,
    border_color = "black",
    stroke = 0,
    clust_method = "ward.D2"
) {
  if (!requireNamespace("circlize", quietly = TRUE)) {
    stop("Package 'circlize' is required for the color palette. Please install it: install.packages('circlize')")
  }
  req_cols <- c(meta_col, x_var, "gene_symbol", "pct.expr", "avg.expr")
  miss <- setdiff(req_cols, names(meta_df))
  if (length(miss)) stop("Missing required columns in meta_df: ", paste(miss, collapse = ", "))
  
  form <- stats::as.formula(paste(y_var, "~", x_var))
  wide_df <- reshape2::acast(meta_df[, c(y_var, x_var, fill_var)], form, value.var = fill_var)
  if (anyNA(wide_df)) wide_df[is.na(wide_df)] <- 0
  result_matrix <- as.matrix(wide_df)
  
  if (!is.null(order)) {
    if (order == "both") {
      row_dist <- stats::dist(result_matrix)
      col_dist <- stats::dist(t(result_matrix))
      hc_row <- stats::hclust(row_dist, method = clust_method)
      hc_col <- stats::hclust(col_dist, method = clust_method)
      order_rows <- rownames(result_matrix)[hc_row$order]
      order_cols <- colnames(result_matrix)[hc_col$order]
    } else if (order == "row") {
      row_dist <- stats::dist(result_matrix)
      hc_row <- stats::hclust(row_dist, method = clust_method)
      order_rows <- rownames(result_matrix)[hc_row$order]
      order_cols <- if (!is.null(user_order_cols)) user_order_cols else sort(unique(meta_df[[x_var]]))
    } else if (order == "column") {
      col_dist <- stats::dist(t(result_matrix))
      hc_col <- stats::hclust(col_dist, method = clust_method)
      order_cols <- colnames(result_matrix)[hc_col$order]
      order_rows <- if (!is.null(user_order_rows)) user_order_rows else sort(unique(meta_df[[y_var]]))
    } else {
      stop("Invalid 'order'. Choose 'row', 'column', 'both', or NULL.")
    }
  } else {
    if (is.null(user_order_rows) || is.null(user_order_cols)) {
      stop("For order = NULL, provide both user_order_rows and user_order_cols.")
    }
    order_rows <- user_order_rows
    order_cols <- user_order_cols
  }
  
  # Enforce x-axis with 'Interacting' centered and the other two groups on either side
  ux <- unique(meta_df[[x_var]])
  others <- sort(setdiff(ux, "Interacting"))
  if (length(others) == 2 && "Interacting" %in% ux) {
    order_cols <- c(others[1], "Interacting", others[2])
  } else {
    order_cols <- ux
  }
  
  if (is.null(fill_colors)) {
    v <- meta_df[[fill_var]]
    vmin <- min(v, na.rm = TRUE); vmax <- max(v, na.rm = TRUE)
    if (vmin == vmax) vmax <- vmin + 1e-8
    pal_fun <- circlize::colorRamp2(c(vmin, vmax), c("#FFFFFF", "#B30000"))
    fill_colors <- pal_fun(seq(vmin, vmax, length.out = 100))
  }
  
  p <- ggplot2::ggplot(
    meta_df,
    ggplot2::aes_string(x = x_var, y = y_var, size = size_var, fill = fill_var)
  ) +
    ggplot2::geom_point(shape = dot_shape, stroke = stroke, colour = border_color) +
    ggplot2::scale_size_continuous(range = c(1, 10), name = size_name) +
    ggplot2::scale_fill_gradientn(colors = fill_colors, name = fill_name) +
    ggplot2::scale_x_discrete(limits = order_cols) +
    ggplot2::scale_y_discrete(limits = order_rows) +
    ggplot2::theme_classic() +
    ggplot2::labs(x = xlab, y = ylab, title = title) +
    ggplot2::theme(
      axis.text.x  = ggplot2::element_text(angle = x_text_angle, hjust = x_text_hjust, size = x_text_size),
      axis.text.y  = ggplot2::element_text(size = y_text_size),
      plot.title   = ggplot2::element_text(face = "bold", size = title_size),
      axis.title.x = ggplot2::element_text(size = x_text_size),
      axis.title.y = ggplot2::element_text(size = y_text_size),
      legend.title = ggplot2::element_text(size = x_text_size),
      legend.text  = ggplot2::element_text(size = y_text_size)
    ) +
    ggplot2::guides(
      size = ggplot2::guide_legend(
        override.aes = list(shape = dot_shape, stroke = 0.7, colour = border_color)
      )
    )
  return(p)
}
          
                                     

#' @title calculate_gene_set_score
#' @description Calculate the mean interaction score for a set of genes
#' @param IMscores A matrix of interaction scores
#' @param gene_sets A list of gene sets, where each set is a vector of gene names
#' @param weighted Logical; if TRUE, gene scores are weighted by their occurrence in multiple gene sets
#' @param method Character; specifies the aggregation method for gene set scores. Options are "geometric_mean" or "arithmetic_mean"
#' @details This function computes mean interaction scores for given gene sets across cell interactions.
#' It supports both geometric and arithmetic means, and can weight gene contributions based on their presence in multiple gene sets.
#' @return A matrix of mean interaction scores for genes in each gene set, with 
#' attributes for log p-value sums and number of genes for later fisher combination
#' @export
calculate_gene_set_score <- function(IMscores, gene_sets, weighted = TRUE, method = c("geometric_mean", "arithmetic_mean")) {
    method <- match.arg(method[1], choices = c("geometric_mean", "arithmetic_mean"))

    IMS <- split(IMscores, IMscores$cell_interaction)

    # Pre-calculate gene overlap counts
    # How many complexes does each gene appear in?
    all_genes <- unique(unlist(gene_sets))
    gene_complex_count <- sapply(all_genes, function(g) {
        sum(sapply(gene_sets, function(gs) g %in% gs))
    })

    scores <- matrix(NA, nrow=length(gene_sets), ncol=length(IMS))
    rownames(scores) <- names(gene_sets)
    colnames(scores) <- names(IMS)

    for (i in seq_along(IMS)) {
        temp <- lapply(gene_sets, function(geneSet) {
            if (length(geneSet) == 0) return(NA)
            ims <- IMS[[i]]
            geneSet <- intersect(geneSet, ims$gene)
            if (length(geneSet) == 0) return(NA)

            subset_data <- ims[ims$gene %in% geneSet,]
            
                        # Calculate uniqueness weight for each gene
            # Weight = 1 / (number of complexes containing this gene)
            if (!weighted) {
                gene_weights <- rep(1, nrow(subset_data))
            } else {
                gene_weights <- 1 / gene_complex_count[subset_data$gene]
                # Normalize weights to sum to number of genes
                # (preserves interpretation of the score)
                gene_weights <- gene_weights * length(gene_weights) / sum(gene_weights)
            }

            # Calculate weighted effect size
            weighted_effect <- switch(method,
                "geometric_mean" = exp(weighted.mean(log(subset_data$effect_size),gene_weights)),
                "arithmetic_mean" = weighted.mean(subset_data$effect_size, gene_weights),
                stop("Unknown method for gene set score calculation")
                )

            weighted_effect

        })
        scores[,i] <- unlist(temp)
    }
    
    scores <- .call(scores, 2, as.numeric)
    colnames(scores) <- names(IMS)
    rownames(scores) <- names(gene_sets)
        # Add gene overlap information as an attribute
    attr(scores, "gene_complex_count") <- gene_complex_count
    
    return(scores)
}


#' @title .pick_image
#' @description The function picks the appropriate histology image file from the
#' spatial directory based on the specified resolution.
#' @param sp_dir path to the spatial directory
#' @param res a character string specifying the resolution of the image
#' @return a character string of the image file name
.pick_image <- function(sp_dir, res) {
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
#' @importFrom ggplot2 scale_fill_gradientn ggplot annotation_raster 
#' geom_point aes coord_fixed labs theme_bw theme element_rect element_blank
#' @importFrom dplyr mutate rename

plot_spatial_data_over_image <- function(
    visiumDir, df, feature_col, barcode_col="barcode",
    resolution=c("lowres","hires","fullres"), version=NULL, colors=NULL,
    point_size=2.5, stroke=0.05, alpha=0.5, title="Spatial Heatmap",
    bg_color=NULL, crop=TRUE, text_size = 15) {
  resolution <- match.arg(resolution)
  if (is.null(barcode_col)) {
    if (!is.null(rownames(df))) df$barcode <- rownames(df) else
      stop("barcode_col is NULL and df has no rownames.")
  } else if (barcode_col != "barcode") 
    df <- dplyr::rename(
    df, barcode = !!rlang::sym(barcode_col))
  pos <- SpaceMarkers::load10XCoords(visiumDir, resolution, version)
  names(pos) <- c("barcode","y","x"); df <- merge(
    df[, c("barcode", feature_col)], pos, by="barcode")
  img <- readbitmap::read.bitmap(.pick_image(
    file.path(visiumDir, "spatial"), resolution))
  if (crop) {xr <- range(df$x, na.rm=TRUE); yr <- range(df$y, na.rm=TRUE)
  xmin <- max(1L, floor(xr[1])); xmax <- min(ncol(img), ceiling(xr[2]))
  ymin <- max(1L, floor(yr[1])); ymax <- min(nrow(img), ceiling(yr[2]))
  img <- img[ymin:ymax, xmin:xmax,, drop=FALSE]
  df <- dplyr::mutate(df, x_c=x-xmin+1L, y_c=y-ymin+1L)
  xl <- c(0, xmax-xmin+1L); yl <- c(0, ymax-ymin+1L)
  } else {
    df <- dplyr::mutate(df, x_c=x, y_c=y);xl<-c(0, ncol(img));yl<-c(0, nrow(img))}
    p <- ggplot2::ggplot() +
      ggplot2::annotation_raster(as.raster(img), 0, diff(xl), 0, diff(yl)) +
      ggplot2::geom_point(
        data = df,
        ggplot2::aes(x_c, yl[2] - y_c, fill = .data[[feature_col]]),
        shape = 21, color = "black", size = point_size,
        stroke = stroke, alpha = alpha
      ) +
      ggplot2::coord_fixed(xlim = xl, ylim = yl, expand = FALSE) +
      ggplot2::labs(fill = feature_col, x = NULL, y = NULL, title = title) +
      ggplot2::theme_bw(text_size) +
      ggplot2::theme(
        plot.background = if (!is.null(bg_color))
          ggplot2::element_rect(fill = bg_color, color = NA)
        else ggplot2::element_blank()
      )
  if (is.numeric(df[[feature_col]])) {
    p <- p + (if (!is.null(colors)) ggplot2::scale_fill_gradientn(colours=colors)
              else viridis::scale_fill_viridis())
  } else {
    vals <- unique(stats::na.omit(df[[feature_col]]))
    p <- p + if (!is.null(colors)) ggplot2::scale_fill_manual(values=colors) else
      if (length(vals) > 1) ggplot2::scale_fill_manual(values = stats::setNames(
        RColorBrewer::brewer.pal(max(3, min(length(vals), 9)), "Set1")
        [seq_along(vals)], vals)) else ggplot2::scale_fill_manual(values="red")
  }
  return(p)
}



#' @title calculate_lr_scores
#' @description Calculate L-R pair scores using Fisher's method
#' @param ligand_scores Output from getGeneSetScore for ligands
#' @param receptor_scores Output from getGeneSetScore for receptors
#' @param lr_pairs Data frame with columns 'ligand' and 'receptor'
#' @param ligand_test Character; specifies the type of test for ligand overexpression. Options are "greater" or "two.sided"
#' @param method Character; specifies the aggregation method for L-R scores. Options are "
#' geometric_mean" or "arithmetic_mean"
#' @param weighted Logical; if TRUE, L-R scores are weighted by their occurrence in
#' multiple L-R pairs
#' @details This function computes L-R pair scores by combining ligand and receptor
#' overexpression scores using either geometric or arithmetic mean. It can also
#' weight L-R pairs based on their presence in multiple pairs to reduce bias from
#' promiscuous ligands or receptors.
#' @return Data frame with L-R scores and p-values
#' @export
calculate_lr_scores <- function(ligand_scores, receptor_scores, lr_pairs,
                              ligand_test = c("greater", "two.sided"), method = c("geometric_mean", "arithmetic_mean"),
                              weighted = TRUE) {

    method <- match.arg(method[1], choices = c("geometric_mean", "arithmetic_mean"))
       # parameter checks
    ligand_test <- match.arg(ligand_test[1], choices = c("greater", "two.sided"))
    
    ## Scoring ligand overexpression near target cell type and receptor overexpression in target cell type
    if (!all(grepl("near",colnames(receptor_scores)))) {
        target_cells <- gsub("^.*_","",colnames(ligand_scores))
        if (!any(target_cells %in% colnames(receptor_scores))) {
            stop("Receptor scores do not have expected cell type names in column names.")
        }
        ligand_cols <- colnames(ligand_scores)[grepl(paste0("near_(", paste(target_cells, collapse="|"), ")"), colnames(ligand_scores))]
        names(ligand_cols) <- ligand_cols
        # Map ligand columns to receptor columns
        mapped_receptor_cols <- sapply(ligand_cols, function(lc) gsub("^.*_","",lc))
        names(mapped_receptor_cols) <- ligand_cols
        lr_scores <- matrix(0, nrow=nrow(lr_pairs), ncol=length(ligand_cols))
        if (weighted) {
            ligand_counts <- table(lr_pairs$ligand.symbol)
            receptor_counts <- table(lr_pairs$receptor.symbol)
            ligand_weights <- 1 / ligand_counts[lr_pairs$ligand.symbol]
            receptor_weights <- 1 / receptor_counts[lr_pairs$receptor.symbol]
            ligand_weights <- matrix(ligand_weights * length(ligand_weights) / sum(ligand_weights), nrow=nrow(lr_pairs), ncol=1)
            receptor_weights <- matrix(receptor_weights * length(receptor_weights) / sum(receptor_weights), nrow=nrow(lr_pairs), ncol=1)
            rownames(ligand_weights) <- rownames(receptor_weights) <- rownames(lr_pairs)
            weights <- cbind(ligand_weights, receptor_weights)
            colnames(weights) <- c("ligand_weights", "receptor_weights")
            rownames(weights) <- rownames(lr_pairs)
        } else {
            weights <- matrix(1, nrow=nrow(lr_pairs), ncol=2)
            colnames(weights) <- c("ligand_weights", "receptor_weights")
            rownames(weights) <- rownames(lr_pairs)
        }

        lr_scores <- matrix(NA, nrow=nrow(lr_pairs), ncol=length(ligand_cols))
        rownames(lr_scores) <- rownames(lr_pairs)
        mapped_lr_cols <- gsub("near_","to_",ligand_cols)
        names(mapped_lr_cols) <- ligand_cols
        colnames(lr_scores) <- mapped_lr_cols
            
        for (i in seq_along(ligand_cols)) {
            lc <- ligand_cols[i]
            rc <- mapped_receptor_cols[i]
            scores <- cbind(ligand_scores[,lc], receptor_scores[,rc])
            colnames(scores) <- c("ligand_score", "receptor_score")
            rownames(scores) <- rownames(lr_pairs)
            lr_scores[,i] <- switch(method,
                "geometric_mean" = exp(sapply(rownames(scores), function(r) {weighted.mean(log(scores[r, ]), w=weights[r,])})),
                "arithmetic_mean" = sapply(rownames(scores), function(r) {weighted.mean(scores[r, ], w=weights[r,])}),
                stop("Unknown method for L-R score calculation")
            )
        }
        colnames(lr_scores) <- mapped_lr_cols
        rownames(lr_scores) <- rownames(lr_pairs)
        return(lr_scores)
    } else {
        stop("Receptor scores have been calculated using calculate_gene_set_score (overexpression) of receptors near source cell; use calculate_gene_set_specificity for receptors instead.")
    }

}

#' Calculate Gene Set Specificity Scores
#' @title calculate_gene_set_specificity
#' @description
#' This function computes specificity scores for given gene sets across cell types or spatial patterns.
#' It uses fold-change scores and p-values to weight gene contributions, and supports both geometric and arithmetic means.
#'
#' @param data A numeric matrix or data frame of gene expression values (genes x samples).
#' @param spPatterns A data frame or matrix containing spatial pattern information, with columns for cell types and optionally "x", "y", "barcode".
#' @param gene_sets A named list of character vectors, where each vector contains gene names for a gene set.
#' @param weighted Logical; if TRUE, gene scores are weighted by their occurrence in multiple gene sets.
#' @param method Character; specifies the aggregation method for gene set scores. Options are "geometric_mean" or "arithmetic_mean".
#'
#' @return A numeric matrix of gene set specificity scores (gene sets x cell types).
#'
#' @details
#' - Genes not present in the data are excluded.
#' - Genes with all zero expression are removed.
#' - Fold-change scores and p-values are calculated using `.calculate_all_fc_scores`.
#' - Scores are normalized and weighted by p-value significance.
#' - For each gene set, scores are aggregated using the specified method and gene weights.
#'
#' @examples
#' # Example usage:
#' # gene_set_scores <- calculate_gene_set_specificity(expr_matrix, spPatterns, gene_sets)
#'
#' @export
calculate_gene_set_specificity <- function(data, spPatterns, gene_sets, weighted = TRUE, method = c("geometric_mean", "arithmetic_mean")) {
    method <- match.arg(method[1], choices = c("geometric_mean", "arithmetic_mean"))
    genes <- unique(unlist(gene_sets))
    genes <- intersect(genes, rownames(data))
    if (length(genes) == 0) {
        stop("No genes from gene sets found in data")
    }
    expr <- data[genes, , drop = FALSE]
    expr[is.na(expr)] <- 0

    expr <- expr[apply(expr, 1, function(x) any(x > 0)), , drop = FALSE]
    
    genes <- rownames(expr)

    gene_complex_count <- sapply(genes, function(g) {
        sum(sapply(gene_sets, function(gs) g %in% gs))
    })
    if (weighted) {
        gene_weights <- 1 / gene_complex_count[genes]
            # (preserves interpretation of the score)
        gene_weights <- gene_weights * length(gene_weights) / sum(gene_weights)
    } else {
        gene_weights <- rep(1, length(genes))
    }
    names(gene_weights) <- genes
            # (preserves interpretation of the score)

    if (nrow(expr) == 0) {
        return(matrix(0, nrow=length(gene_sets), ncol=ncol(spPatterns)))
    }

    # Group by cell type
    patnames <- setdiff(colnames(spPatterns),c("x","y","barcode"))
    cell_types <- unique(patnames)

    # Initialize result matrix
    gene_set_scores <- matrix(0, nrow = length(gene_sets), ncol = length(cell_types))
    rownames(gene_set_scores) <- names(gene_sets)
    colnames(gene_set_scores) <- cell_types

    gene_scores <- .calculate_all_fc_scores(expr, spPatterns, low_thr = 0.2, high_thr = 0.8)


    for (gene_set in names(gene_sets)) {
        genes <- gene_sets[[gene_set]]
        valid_genes <- intersect(genes, rownames(gene_scores))
        if (length(valid_genes) > 1) {
            
            # Calculate weighted scores
            gene_set_scores[gene_set, ] <- switch(method,
                "geometric_mean" = exp(apply(log(gene_scores[valid_genes, ]), 2, weighted.mean, w = gene_weights[valid_genes])),
                "arithmetic_mean" = apply(gene_scores[valid_genes, ], 2, weighted.mean, w = gene_weights[valid_genes]),
                stop("Unknown method for gene set score calculation")
            )
        } else if (length(valid_genes) == 1) {
        gene_set_scores[gene_set, ] <- gene_scores[valid_genes, ]
        } else {
        gene_set_scores[gene_set, ] <- 0
        }
    }

    return(gene_set_scores)
}

.calculate_fc_score <- function(expr, spPatterns, gene, ct,
                              low_thr = 0.2, high_thr = 0.8) {

  # Get thresholds
  low_thr <- quantile(spPatterns[, ct], probs = low_thr, na.rm = TRUE)
  high_thr <- quantile(spPatterns[, ct], probs = high_thr, na.rm = TRUE)

  # Get high and low bins
  high_bins <- which(spPatterns[, ct] > high_thr)
  low_bins <- which(spPatterns[, ct] < low_thr)
  
  # Calculate means
  mean_high <- mean(expr[gene, high_bins])
  mean_low <- mean(expr[gene, low_bins])
  
  # Log fold change
  lfc <- mean_high - mean_low
  
  # Get p-value
  if (all(expr[gene, c(high_bins, low_bins)] == 0)) {
    score <- NA
    attr(score, "p_value") <- 1
    return(score)
  }
  w_test <- wilcox.test(as.matrix(expr[gene, high_bins]), as.matrix(expr[gene, low_bins]))

  score <- lfc
  attr(score, "p_value") <- w_test$p.value
  return(score)
}

# Wrapper for all genes and cell types
# importFrom BiocParallel bplapply
.calculate_all_fc_scores <- function(expr, spPatterns, 
                                    low_thr = 0.2, high_thr = 0.9, ..., workers = NULL) {

  genes <- rownames(expr)
  cell_types <- setdiff(colnames(spPatterns), c("x", "y", "barcode"))
  
  # Initialize results
  score_matrix <- matrix(NA, nrow = length(genes), ncol = length(cell_types))
  rownames(score_matrix) <- genes
  colnames(score_matrix) <- cell_types
  
  gene_scores <- lfc_matrix <- p_values <- score_matrix
  use_biocparallel <- requireNamespace("BiocParallel", quietly = TRUE)
  # Calculate for each combination
  if (!use_biocparallel) {
    for (i in seq_along(genes)) {
        for (ct in cell_types) {
            result <- .calculate_fc_score(expr, spPatterns, 
                                    genes[i], ct,
                                    low_thr, high_thr)
            lfc_matrix[i, ct] <- result
            p_values[i, ct] <- attr(result, "p_value")
        }
    }
  } else {
    if (is.null(workers) || workers < 1) {
        workers <- BiocParallel::bpworkers(BiocParallel::bpparam())
    }
    BiocParallel::register(BiocParallel::MulticoreParam(workers))
    lfc <- BiocParallel::bplapply(cell_types, function(ct) {
        lfc_col <- p_val_col <- numeric(length(genes))
        names(lfc_col) <- names(p_val_col) <- genes
        for (i in seq_along(genes)) {
            result <- .calculate_fc_score(expr, spPatterns, 
                                genes[i], ct,
                                low_thr, high_thr)
            lfc_col[i] <- result
            p_val_col[i] <- attr(result, "p_value")
        }
        return(list(lfc = lfc_col, p_value = p_val_col))
    })
    for (i in seq_along(cell_types)) {
        ct <- cell_types[i]
        lfc_matrix[, ct] <- lfc[[i]]$lfc
        p_values[, ct] <- lfc[[i]]$p_value
    }
  }
    p_adj <- apply(p_values,2,p.adjust,method="BH")
    p_weights <- 1 / (1 + exp(2 * (log10(p_adj) - log10(0.05))))

    norm_scores <- 1/(1 + exp(-40 * (lfc_matrix-median(lfc_matrix[lfc_matrix>0], na.rm = TRUE))))
    gene_scores <- norm_scores * p_weights
    gene_scores[is.na(gene_scores)] <- 0
    attr(gene_scores,"p_values") <- p_values
    attr(gene_scores,"p_adj") <- p_adj
    attr(gene_scores,"p_weights") <- p_weights
    attr(gene_scores,"lfc") <- lfc_matrix
  return(gene_scores)
}

.call <- function(x, MARGIN = NULL, FUN, ...) {
    # Save attributes
    saved_attrs <- attributes(x)

    # Apply function
    if (is.null(MARGIN)) {
        result <- FUN(x, ...)
    } else {
        result <- apply(x, MARGIN, FUN, ...)
    }

    # Restore custom attributes
    for (attr_name in names(saved_attrs)) {
        attr(result, attr_name) <- saved_attrs[[attr_name]]
    }
    
    return(result)
}
