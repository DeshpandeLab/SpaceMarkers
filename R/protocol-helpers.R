# Protocol helper functions imported from the \`protocol\` branch to
# enable the SpaceMarkers_Protocol notebooks to run on top of the
# SpaceMarkersExperiment (SME) workflow. Function bodies are copied
# unchanged from the `protocol` branch to preserve behavior.

#--- from origin/protocol R/hdutils.R -----------------------------

#' @rdname classify_spots
#' @title Classify spots into interacting / non-interacting pattern regions
#' @param pat_hotspots Data frame of pattern hotspots.
#' @param influence_hotspots Data frame of influence hotspots. If NULL, symmetric classification.
#' @param patternpair Length-2 character vector of pattern names to classify.
#' @return A data frame of region labels (see details in the protocol vignette).
#' @export
classify_spots <- function(pat_hotspots, influence_hotspots = NULL, patternpair = NULL) {
    .classify_spots(pat_hotspots = pat_hotspots,
                    influence_hotspots = influence_hotspots,
                    patternpair = patternpair)
}

#' @title Classify all pattern pairs for each spot
#' @description
#' For all pairwise combinations of patterns in a hotspot matrix, classify
#' each spot as interacting or belonging to one pattern, optionally using
#' both pattern and influence hotspots.
#'
#' @param pat_hotspots Data frame with columns barcode, x, y and one or
#' more pattern hotspot columns.
#' @param influence_hotspots Optional data frame with the same structure
#' as pat_hotspots, containing influence hotspots. If NULL, only pattern
#' hotspots are used.
#'
#' @return A data.frame with one column per pattern pair and one row per
#' spot containing classification labels.
#' @export
classify_allspots <- function(pat_hotspots,
                              influence_hotspots = NULL) {
  key_cols <- c("barcode", "x", "y")
  patnames <- setdiff(colnames(pat_hotspots), key_cols)
  
  if (length(patnames) < 2L) {
    stop("Need at least 2 pattern columns in pat_hotspots.")
  }
  
  # All pairwise combinations of patterns
  patternpairs <- utils::combn(patnames, 2, simplify = FALSE)
  
  # Run classify_spots for every pair
  class_list <- lapply(patternpairs, function(pp) {
    classify_spots(
      pat_hotspots       = pat_hotspots,
      influence_hotspots = influence_hotspots,
      patternpair        = pp
    )
  })
  
  # cbind all interacting hotspot columns together
  class_df <- do.call(cbind, class_list)
  
  class_df
}

#' @title Create a long-format ligand–receptor score data frame
#'
#' @description Takes LR, ligand, and receptor score matrices and combines them into a single
#' long data frame with one row per interaction and source/target cell-type pair.
#' Ligand and receptor gene names are added from \code{lrpairs}. Any missing LR
#' scores are set to 0. Column names in \code{ligand_scores} that start with
#' \code{near_} are converted to \code{to_} to match LR columns.
#'
#' @param lrscores Matrix of LR scores (rows = interactions, columns = source-to-target pairs).
#' @param ligand_scores Matrix of ligand scores (same shape as \code{lrscores}).
#' @param receptor_scores Matrix of receptor scores (rows = interactions, columns = target cell types).
#' @param lrpairs Data frame with columns \code{ligand} and \code{receptor}; rownames are interaction IDs.
#' @param complex_sep Separator used for multi-gene complexes in \code{lrpairs}; replaced with \code{"_"}.
#'
#' @return A data frame with columns:
#' \code{source_cell_type}, \code{ligand}, \code{receptor}, \code{target_cell_type},
#' \code{ligand_score}, \code{receptor_score}, \code{score}, \code{interaction}.
#' The interaction and cell types are parsed from the combined \code{lr_cell} key.
#' @export
create_lr_dataframe <- function(lrscores, ligand_scores, receptor_scores,lrpairs,complex_sep = ", "){
  # Lr scores
  lrscores[is.na(lrscores)] <- 0
  lr_long <- reshape2::melt(as.matrix(lrscores),
                            varnames = c("interaction", "source_to_target"),
                            value.name = "score")
  lr_long$lr_cell <- paste0(
    lr_long$interaction,"_and_",lr_long$source_to_target) # merge on this
  # Ligand scores
  ligand_long <- reshape2::melt(as.matrix(ligand_scores),
                                varnames = c("interaction", "source_to_target"),
                                value.name = "ligand_score")
  ligand_long$source_to_target <- gsub(
    "near_", "to_", ligand_long$source_to_target)
  ligand_long$lr_cell <- paste0(
    ligand_long$interaction,"_and_",ligand_long$source_to_target) # merge on this
  # LR scores + ligand scores
  lr_ligand <- dplyr::inner_join(
    lr_long %>% dplyr::select(c(-source_to_target,-interaction)),
    ligand_long  %>% dplyr::select(c(-source_to_target, -interaction)),
    by = "lr_cell")
  ##  extract  relevant columns
  lr_ligand$interaction <- sub("_and_.*$", "", lr_ligand$lr_cell)
  lr_ligand$source_cell_type <- sub("^.*_and_(.*)_to_.*$", "\\1",
                                    lr_ligand$lr_cell)
  lr_ligand$target_cell_type <- sub("^.*_to_", "", lr_ligand$lr_cell)
  # Receptor scores (merge on targets))
  receptor_long <- reshape2::melt(
    as.matrix(receptor_scores),
    varnames = c("interaction", "target_cell_type"),
    value.name = "receptor_score")
  receptor_long$lr_target <- paste0(receptor_long$interaction,
                                    "_and_",receptor_long$target_cell_type)
  # LR scores + receptor scores
  lr_ligand$lr_target <- paste0(
    lr_ligand$interaction,"_and_",lr_ligand$target_cell_type)
  lr_receptor <- dplyr::inner_join(
    lr_ligand,receptor_long %>%
      dplyr::select(c(-target_cell_type,-interaction)) , by = "lr_target")
  # Get ligands and receptors
  lrpairs$interaction <- rownames(lrpairs)
  final_lrscores <- dplyr::inner_join(lr_receptor,lrpairs,by = "interaction")
  # format ligand and receptors
  final_lrscores$receptor <- gsub(pattern = complex_sep,
                                  replacement = "_",final_lrscores$receptor)
  final_lrscores$ligand <- gsub(pattern = complex_sep,
                                replacement = "_",final_lrscores$ligand)
  final_lrscores <- final_lrscores %>%
    dplyr::select(source_cell_type,ligand,receptor,target_cell_type,
                  ligand_score,receptor_score,score,interaction)
  final_lrscores$source_to_target <- paste0(final_lrscores$source_cell_type,
                                            "_to_",
                                            final_lrscores$target_cell_type)
  return(final_lrscores)
}



#--- from origin/protocol R/utils.R -------------------------------

#' @title calc_overlap_scores
#' @description Calculate the overlap scores between patterns in hotspots
#' @param hotspots A data frame with columns x, y, barcode and pattern names
#' @param patternList A character vector of pattern names to calculate overlap
#' scores for
#' @param method The method to calculate overlap scores. Options are
#' "Szymkiewicz-Simpson", "Jaccard", "Sorensen-Dice", "Ochiai", "absolute",
#' and "NPMI".
#' @details The function calculates the overlap scores between patterns hotspots
#' using the specified method. The default method is "Szymkiewicz-Simpson"
#' overlap coefficient.
#'
#' The "NPMI" method computes the Normalized Pointwise Mutual Information
#' (Bouma 2009), defined as:
#'   NPMI(x,y) = PMI(x,y) / -log(p(x,y))
#' where PMI(x,y) = log(p(x,y) / (p(x) * p(y))), and probabilities are
#' estimated from the fraction of spots in each pattern. NPMI ranges from
#' -1 (patterns never co-occur) through 0 (independent) to +1 (one pattern
#' is a perfect subset of the other). Pairs with zero intersection are
#' assigned -1 by convention.
#' @return A matrix of overlap scores (patterns x patterns)
#' @export
#' @examples
#' hotspots <- data.frame(x = c(1,2,3,4,5),
#'                         y = c(1,2,3,4,5),
#'                         barcode = c("A","B","C","D","E"),
#'                         pattern1 = c(1,0,1,0,1),
#'                         pattern2 = c(1,1,0,0,1))
#' calc_overlap_scores(hotspots = hotspots)
#' calc_overlap_scores(hotspots = hotspots, method = "NPMI")
#' @importFrom ggplot2 ggplot geom_tile geom_text theme_minimal
#' @importFrom reshape2 melt
#' @importFrom stats complete.cases
calc_overlap_scores <- function(hotspots, patternList = NULL,
                                method = c("Szymkiewicz-Simpson",
                                           "Jaccard", "Sorensen-Dice",
                                           "Ochiai", "absolute", "NPMI")) {
  # stop if more than one method is supplied, do not warn by default
  if (length(method) > 1) {
    method <- method[1]
    message("Only one method can be used at a time. Using ", method)
  }
  
  if (is.null(patternList)) {
    patternList <- setdiff(colnames(hotspots), c("x", "y", "barcode"))
  } else if (!all(patternList %in% colnames(hotspots))) {
    stop("Pattern names not found in hotspots")
  }
  
  binarized  <- (!is.na(hotspots[, patternList])) * 1
  intersects <- t(binarized) %*% binarized
  nHotspots  <- colSums(binarized)
  nHotsP1    <- t(t(nHotspots)) %*% array(1, length(patternList))
  nHotsP2    <- t(nHotsP1)
  
  overlapScore <- switch(
    method,
    "Szymkiewicz-Simpson" = intersects / pmin(nHotsP1, nHotsP2),
    "Jaccard"             = intersects / (nHotsP1 + nHotsP2 - intersects),
    "Sorensen-Dice"       = 2 * intersects / (nHotsP1 + nHotsP2),
    "Ochiai"              = intersects / sqrt(nHotsP1 * nHotsP2),
    "absolute"            = intersects,
    "NPMI"                = {
      N <- nrow(binarized)
      
      # Work in log space; guard against log(0) by using NA for zero intersects
      log_intersects        <- matrix(NA_real_, nrow = nrow(intersects),
                                      ncol = ncol(intersects),
                                      dimnames = dimnames(intersects))
      has_overlap           <- intersects > 0
      log_intersects[has_overlap] <- log(intersects[has_overlap])
      
      # PMI  = log(|A∩B|) + log(N) - log(|A|) - log(|B|)
      pmi_matrix <- log_intersects + log(N) - log(nHotsP1) - log(nHotsP2)
      
      # Normalisation denominator: -log( p(x,y) ) = log(N) - log(|A∩B|)
      norm_denom <- log(N) - log_intersects
      
      npmi_matrix <- pmi_matrix / norm_denom
      
      # Pairs with zero intersection → NPMI = -1 by convention
      npmi_matrix[!has_overlap] <- -1
      
      npmi_matrix
    },
    stop("Method not supported")
  )
  
  return(overlapScore)
}

#' @title Compute pairwise overlap scores between patterns
#' @description This function computes overlap scores between pattern hotspots,
#' either symmetrically within pat_hotspots or between pat_hotspots and in_hotspots,
#' using a specified overlap metric.
#' @param pat_hotspots A data frame with x, y, barcode columns and additional
#' columns for pattern-specific hotspot values. If in_hotspots is NULL, symmetric
#' overlaps are computed between these patterns.
#' @param in_hotspots An optional data frame with x, y, barcode columns and
#' pattern-specific hotspot values representing neighborhood or influence patterns.
#' Must share barcodes with pat_hotspots when provided.
#' @param patternList An optional character vector specifying a subset of patterns
#' to use when computing symmetric overlap scores.
#' @param method A character vector specifying the overlap metric to use.
#' Options include "Szymkiewicz-Simpson", "Jaccard", "Sorensen-Dice", "Ochiai",
#' "absolute", and "NPMI". Only the first value is used.
#' @return A data frame with columns pattern1, pattern2, and overlapScore,
#' containing pairwise overlap scores for the selected patterns.
#' @export
get_overlap_scores <- function(pat_hotspots = NULL, in_hotspots = NULL,
                               patternList = NULL,
                               method = c("Szymkiewicz-Simpson", "Jaccard",
                                          "Sorensen-Dice", "Ochiai",
                                          "absolute", "NPMI")) {
  if (!is.null(pat_hotspots) & is.null(in_hotspots)) {
    message("Assuming symmetric overlap scores. Setting upper triangle and diagonal to NA.")
    overlapScore <- calc_overlap_scores(hotspots = pat_hotspots,
                                        patternList = patternList, method = method)
    overlapScore[upper.tri(overlapScore, diag = TRUE)] <- NA
    dfOverlap <- reshape2::melt(overlapScore)
    dfOverlap <- dfOverlap[stats::complete.cases(dfOverlap), ]
    colnames(dfOverlap) <- c("pattern2", "pattern1", "overlapScore")
    dfOverlap <- dfOverlap[, c(2, 1, 3)]
  } else {
    patternList_in <- setdiff(colnames(in_hotspots), c("x", "y", "barcode"))
    rownames(in_hotspots) <- in_hotspots$barcode
    in_hotspots <- in_hotspots[, patternList_in, drop = FALSE]
    colnames(in_hotspots) <- paste0("near_", patternList_in)
    
    patternList_pat <- setdiff(colnames(pat_hotspots), c("x", "y", "barcode"))
    rownames(pat_hotspots) <- pat_hotspots$barcode
    hotspots <- cbind(pat_hotspots[, patternList_pat, drop = FALSE],
                      in_hotspots[rownames(pat_hotspots), , drop = FALSE])
    
    overlapScore <- calc_overlap_scores(hotspots = hotspots, method = method)
    dfOverlap <- reshape2::melt(overlapScore)
    dfOverlap <- dfOverlap[stats::complete.cases(dfOverlap), ]
    colnames(dfOverlap) <- c("pattern2", "pattern1", "overlapScore")
    dfOverlap <- dfOverlap[, c(2, 1, 3)]
    dfOverlap <- dfOverlap[!grepl("^near_", dfOverlap$pattern1), ]
    dfOverlap <- dfOverlap[grepl("^near_", dfOverlap$pattern2), ]
  }
  
  return(dfOverlap)
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
#' @importFrom fgsea fgsea
#' @importFrom dplyr rename arrange
#' 
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

#' Plot enriched pathways by interaction
#' @title plot_enriched_results
#' @description
#' Create a dot plot of pathways (gene sets) versus ligand–receptor
#' interactions from a list of enrichment result tables.
#' Dot size shows the size of the leading edge, and dot fill shows
#' minus log10 adjusted p-values.
#'
#' @param enr_list
#' A list of enrichment result data frames.
#' Each non-empty data frame must contain the columns:
#' gene.set, interaction, leadingEdge, and log10P_adj.
#' @param top_n
#' Integer. Number of top pathways to keep per interaction,
#' ranked by log10P_adj. Default is 5.
#'
#' @return
#' A ggplot object with pathways on the y-axis, interactions on
#' the x-axis, dot size for number of leading-edge genes, and
#' dot fill for minus log10(FDR).
#' @importFrom dplyr bind_rows group_by filter slice_max ungroup pull mutate
#'   summarise arrange
#' @importFrom rlang .data
#' @importFrom ggplot2 ggplot aes geom_point scale_size_area scale_fill_gradient
#'   labs theme_minimal element_blank element_text
#' @export
#' 
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

#' Average and detection by group for a gene set
#'
#' @description
#' Compute average expression and percent of cells expressing each gene in a
#' set, stratified by a metadata grouping column. Barcodes in metadata and
#' counts are aligned, rows with NA in the grouping column are dropped, and
#' columns with any NA values in the counts matrix are removed.
#'
#' @param counts Matrix or dgCMatrix with genes as rows and cells as columns.
#' @param metadata Data frame with per cell metadata. Must include meta_col and
#'   either have barcodes as row names or a column named in barcode_col.
#' @param features Character vector of gene symbols to include.
#' @param meta_col Character. Name of the metadata column used to group cells
#'   (for example "celltype").
#' @param barcode_col Optional character. Column name in metadata holding cell
#'   barcodes. If NULL, row names of metadata are used as barcodes.
#'
#' @return A data.frame with columns:
#'   meta_col, gene_symbol, pct.expr, avg.expr.
#'
#' @examples
#' \dontrun{
#' # Example: average expression by SpaceMarkers interaction region
#'
#' # Suppose:
#' #   - spCounts is a genes x barcodes count matrix
#' #   - spHotspots_all is a metadata data.frame with a "barcode" column
#' #     and an interaction region column "Epi_Plasma"
#' #   - df is a data.frame with a column "Gene" containing gene symbols
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
#' 
#' 
#' @importFrom reshape2 melt
#' @importFrom dplyr group_by summarise
#' @importFrom rlang .data
#' @export
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
  
  # --- Initial overlap and alignment by barcode ---
  cell_barcodes <- colnames(counts)
  if (is.null(cell_barcodes)) stop("`counts` must have column names (cell barcodes).")
  common_bc <- intersect(cell_barcodes, md_barcodes)
  if (!length(common_bc)) stop("No overlapping barcodes between counts and metadata.")
  metadata <- metadata[common_bc, , drop = FALSE]
  counts   <- counts[, common_bc, drop = FALSE]
  
  # --- NA handling ---
  # 1) Drop metadata rows where meta_col is NA
  keep_meta <- !is.na(metadata[[meta_col]])
  if (any(!keep_meta)) {
    metadata <- metadata[keep_meta, , drop = FALSE]
    counts   <- counts[, rownames(metadata), drop = FALSE]
  }
  # 2) Drop counts columns that contain any NA
  if (anyNA(counts)) {
    if (inherits(counts, "dgCMatrix")) {
      na_mat <- is.na(as.matrix(counts))
    } else {
      na_mat <- is.na(counts)
    }
    has_na_col <- colSums(na_mat) > 0
    if (any(has_na_col)) {
      counts   <- counts[, !has_na_col, drop = FALSE]
      metadata <- metadata[colnames(counts), , drop = FALSE]
    }
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
  
  # --- Transpose to cells x genes and attach group column ---
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
  
  # --- Wide to long and aggregate ---
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
  out
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
    ggplot2::aes(
      x     = !!rlang::sym(x_var),
      y     = !!rlang::sym(y_var),
      size  = !!rlang::sym(size_var),
      fill  = !!rlang::sym(fill_var)
    )
  ) +
    ggplot2::geom_point(shape = dot_shape,
                        stroke = stroke,
                        colour = border_color) +
    ggplot2::scale_size_continuous(range = c(1, 10), name = size_name) +
    ggplot2::scale_fill_gradientn(colors = fill_colors, name = fill_name) +
    ggplot2::scale_x_discrete(limits = order_cols) +
    ggplot2::scale_y_discrete(limits = order_rows) +
    ggplot2::theme_classic() +
    ggplot2::labs(x = xlab, y = ylab, title = title) +
    ggplot2::theme(
      axis.text.x  = ggplot2::element_text(angle = x_text_angle,
                                           hjust = x_text_hjust,
                                           size  = x_text_size),
      axis.text.y  = ggplot2::element_text(size = y_text_size),
      plot.title   = ggplot2::element_text(face = "bold", size = title_size),
      axis.title.x = ggplot2::element_text(size = x_text_size),
      axis.title.y = ggplot2::element_text(size = y_text_size),
      legend.title = ggplot2::element_text(size = x_text_size),
      legend.text  = ggplot2::element_text(size = y_text_size)
    ) +
    ggplot2::guides(
      size = ggplot2::guide_legend(
        override.aes = list(shape = dot_shape,
                            stroke = 0.7,
                            colour = border_color)
      )
    )
  
  return(p)
}

