#' @title Align counts and coordinates
#' @description
#' Align a sparse counts matrix and a coordinates data frame so that they use
#' the same set of barcodes in the same order. The coordinates are reordered
#' to match the column order of the counts matrix.
#'
#' @param counts A sparse genes x spots matrix with column names as barcodes.
#' @param coords_df A data.frame with a column called barcode.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{counts}{Filtered counts matrix with columns restricted to common barcodes.}
#'   \item{coords}{Filtered coordinate data.frame with rows in the same order as the counts columns.}
#' }
#'
#' @examples
#' aligned <- align_counts_coords(counts, coords)
#' aligned$counts[1:5, 1:5]
#' aligned$coords[1:5, ]
#'
#' @keywords internal
align_counts_coords <- function(counts, coords_df) {
  # --- sanity checks ---
  if (is.null(colnames(counts))) {
    stop("`counts` must have column names corresponding to barcodes.")
  }
  if (!("barcode" %in% colnames(coords_df))) {
    stop("`coords_df` must contain a 'barcode' column.")
  }
  
  # --- find common barcodes ---
  common <- intersect(colnames(counts), coords_df$barcode)
  if (length(common) == 0L) {
    stop("No overlapping barcodes between counts and coords.")
  }
  
  # --- subset and reorder ---
  counts    <- counts[, common, drop = FALSE]
  coords_df <- coords_df[match(common, coords_df$barcode), , drop = FALSE]
  rownames(coords_df) <- coords_df$barcode
  
  # --- return aligned objects ---
  list(counts = counts, coords = coords_df)
}



#' @title Process multiple Visium or VisiumHD samples
#' @description
#' Load counts and spatial coordinates for multiple samples, filter genes by
#' a per-sample nonzero threshold, align barcodes between counts and coordinates,
#' and save per-sample and combined RDS files.
#'
#' @param data_dir Character. Root directory containing sample folders.
#' @param samples Character vector of sample folder names under data_dir.
#' @param counts_file Character. H5 file name in each sample bin directory.
#' @param good_gene_threshold Integer >= 1. Keep genes detected in at least
#' this many spots per sample.
#' @param bin_subdir Character. Relative path under each sample directory
#' where the 10X files live.
#' @param out_dir Character or NULL. Output directory for RDS files. If NULL,
#' file.path(data_dir, "outputs") is used.
#' @param coords_resolution Character. Resolution argument passed to load10XCoords.
#' @param verbose Logical. If TRUE, print progress messages.
#'
#' @details
#' This function assumes load10XExpr returns a dgCMatrix with log-transformed
#' counts, and load10XCoords returns a data.frame with a barcode column.
#'
#' @return A list with three elements:
#' \describe{
#'   \item{counts}{Named list of filtered dgCMatrix objects per sample.}
#'   \item{coords}{Named list of aligned coordinate data.frames per sample.}
#'   \item{paths}{List of per-sample RDS paths and a combined RDS path.}
#' }
#'
#' @export
process_visium_samples <- function(
    data_dir,
    samples,
    counts_file         = "filtered_feature_bc_matrix.h5",
    good_gene_threshold = 10L,
    bin_subdir          = "binned_outputs/square_016um",
    out_dir             = NULL,
    coords_resolution   = "fullres",
    verbose             = TRUE
) {
  # ---- validate inputs ----
  if (!is.character(data_dir) || length(data_dir) != 1L)
    stop("`data_dir` must be a single character path.")
  if (!dir.exists(data_dir))
    stop("`data_dir` does not exist: ", data_dir)
  
  if (!is.character(samples) || length(samples) == 0L)
    stop("`samples` must be a non-empty character vector of folder names.")
  
  if (!is.character(counts_file) || length(counts_file) != 1L)
    stop("`counts_file` must be a single character filename.")
  
  if (!is.numeric(good_gene_threshold) || length(good_gene_threshold) != 1L || good_gene_threshold < 1)
    stop("`good_gene_threshold` must be a single integer >= 1.")
  
  if (is.null(out_dir)) {
    out_dir <- file.path(data_dir, "outputs")
  }
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  # containers
  spCounts_list <- vector("list", length(samples)); names(spCounts_list) <- samples
  spCoords_list <- vector("list", length(samples)); names(spCoords_list) <- samples
  path_counts   <- setNames(vector("character", length(samples)), samples)
  path_coords   <- setNames(vector("character", length(samples)), samples)
  
  # ---- process each sample ----
  for (s in samples) {
    if (verbose) cat("\n=== Processing:", s, "===\n")
    visium_dir <- file.path(data_dir, s, bin_subdir)
    
    if (!dir.exists(visium_dir)) {
      stop("Expected directory not found for sample '", s, "': ", visium_dir)
    }
    
    # 1) Load counts & coords
    counts <- load10XExpr(visiumDir = visium_dir, h5filename = counts_file)
    coords <- load10XCoords(visiumDir = visium_dir, resolution = coords_resolution)
    
    if (!("barcode" %in% colnames(coords))) {
      stop("`coords` must contain a 'barcode' column for sample '", s, "'.")
    }
    rownames(coords) <- coords$barcode
    
    # Ensure sparse type
    if (!inherits(counts, "dgCMatrix")) {
      if (is.matrix(counts)) {
        counts <- Matrix::Matrix(counts, sparse = TRUE)
        counts <- methods::as(counts, "dgCMatrix")
      } else {
        stop("`counts` is not a dgCMatrix and could not be coerced (sample '", s, "').")
      }
    }
    
    # 2) Per-sample gene filter: nnz per gene (row) from CSC @i
    nnz_by_row <- tabulate(counts@i + 1L, nbins = nrow(counts))
    keep <- nnz_by_row >= as.integer(good_gene_threshold)
    if (!any(keep)) {
      stop("After filtering, zero genes remain for sample '", s, "'. ",
           "Consider lowering `good_gene_threshold`.")
    }
    counts <- counts[keep, , drop = FALSE]
    
    # 3) Align barcodes between counts/coords
    aligned <- align_counts_coords(counts, coords)
    counts  <- aligned$counts
    coords  <- aligned$coords
    
    # 4) Save per-sample .rds
    counts_path <- file.path(out_dir, paste0(s, "_counts.rds"))
    coords_path <- file.path(out_dir, paste0(s, "_coords.rds"))
    saveRDS(counts, counts_path)
    saveRDS(coords,  coords_path)
    
    # 5) Keep in memory
    spCounts_list[[s]] <- counts
    spCoords_list[[s]] <- coords
    path_counts[[s]]   <- counts_path
    path_coords[[s]]   <- coords_path
    
    if (verbose) {
      cat(sprintf("Saved: %s [%d genes x %d spots]\n",
                  basename(counts_path), nrow(counts), ncol(counts)))
      cat(sprintf("Saved: %s [%d spots]\n",
                  basename(coords_path), nrow(coords)))
    }
  }
  
  # 6) Save combined object
  all_path <- file.path(out_dir, "all_samples_counts_coords.rds")
  saveRDS(list(counts = spCounts_list, coords = spCoords_list), all_path)
  if (verbose) cat("\nAll samples saved to:", all_path, "\n")
  
  # return structured result
  list(
    counts = spCounts_list,
    coords = spCoords_list,
    paths  = list(counts = path_counts, coords = path_coords, all = all_path)
  )
}


#' @title Build a per-spot spatial feature table
#' @description
#' Combine spatial coordinates and spatial features into a single data.frame
#' for one sample. Barcodes are intersected across counts, coordinates, and
#' features and aligned without merges.
#'
#' @param spCounts A genes x barcodes matrix (dense or sparse). Column names
#' must be barcodes.
#' @param spCoords A coordinates table (data.frame or matrix) with numeric
#' columns x and y. If barcode is missing, row names are used.
#' @param spFeatures A per-barcode feature table. Can be a matrix or sparse
#' Matrix with row names as barcodes, or a data.frame with or without a
#' barcode column.
#'
#' @return A data.frame with barcodes as row names and columns:
#' barcode, x, y, and feature columns.
#'
#' @export
make_sp_patterns <- function(spCounts, spCoords, spFeatures) {
  # --- coords: ensure barcode column & de-dupe; keep as data.frame ---
  if (!"barcode" %in% names(spCoords)) {
    spCoords <- data.frame(barcode = rownames(spCoords), spCoords, check.names = FALSE)
  } else {
    spCoords <- as.data.frame(spCoords, stringsAsFactors = FALSE)
  }
  spCoords <- spCoords[!duplicated(spCoords$barcode), , drop = FALSE]
  
  # --- features: get barcode index and a matrix view without forcing a data.frame yet ---
  if (is.matrix(spFeatures) || inherits(spFeatures, "Matrix")) {
    feat_barcodes <- rownames(spFeatures)
    feat_mat <- spFeatures
  } else if (is.data.frame(spFeatures)) {
    if ("barcode" %in% names(spFeatures)) {
      feat_barcodes <- spFeatures$barcode
      rownames(spFeatures) <- feat_barcodes
      feat_mat <- as.matrix(spFeatures[, setdiff(names(spFeatures), "barcode"), drop = FALSE])
    } else {
      feat_barcodes <- rownames(spFeatures)
      feat_mat <- as.matrix(spFeatures)
    }
  } else {
    stop("Unsupported spFeatures type: ", paste(class(spFeatures), collapse = ", "))
  }
  
  # --- intersect barcodes across counts/coords/features ---
  if (is.null(colnames(spCounts))) stop("spCounts must have colnames (barcodes).")
  barcodes <- Reduce(intersect, list(colnames(spCounts), spCoords$barcode, feat_barcodes))
  if (length(barcodes) == 0L) stop("No overlapping barcodes among counts/coords/features.")
  
  # --- align (NO merge): keep order stable with match() ---
  coords_ord <- spCoords[match(barcodes, spCoords$barcode), , drop = FALSE]
  if (!all(c("x","y") %in% names(coords_ord)))
    stop("spCoords must contain numeric columns 'x' and 'y'.")
  
  feat_ord <- feat_mat[barcodes, , drop = FALSE]
  
  out <- cbind(coords_ord, as.data.frame(as.matrix(feat_ord), check.names = FALSE))
  rownames(out) <- out$barcode
  out
}

#' @title Compute overlap scores across multiple samples
#' @description This function applies get_overlap_scores to multiple samples and
#' reshapes the results into a wide matrix with interactions as rows and datasets
#' as columns.
#' @param pattern_hotspots A named list of per-sample hotspot data frames to be
#' passed as pat_hotspots to get_overlap_scores.
#' @param influence_hotspots An optional named list of per-sample hotspot data
#' frames to be passed as in_hotspots to get_overlap_scores. If NULL, symmetric
#' overlaps are computed within each element of pattern_hotspots.
#' @param datasets A character vector of dataset names corresponding to the
#' elements of pattern_hotspots. These become the column names in the wide output.
#' @param ... Additional arguments passed to get_overlap_scores, such as
#' patternList or method.
#' @return A data frame in wide format with one row per interaction and one
#' column per dataset containing overlap scores.
#' @export
get_multi_sample_overlaps <- function(pattern_hotspots,
                                      influence_hotspots = NULL,
                                      datasets,
                                      ...) {
  # 1) If influence_hotspots is NULL, make a list of NULLs of the same length
  if (is.null(influence_hotspots)) {
    influence_hotspots <- replicate(length(pattern_hotspots), NULL, simplify = FALSE)
  }
  
  combined_df <- do.call(
    "rbind",
    Map(
      function(pat_hotspots, inf_hotspots, dataset, ...) {
        # 2) Fix typo in argument names passed to get_overlap_scores
        df <- get_overlap_scores(
          pat_hotspots = pat_hotspots,
          in_hotspots  = inf_hotspots,
          ...
        )
        
        df <- df[stats::complete.cases(df), , drop = FALSE]
        
        df <- transform(
          df,
          interaction = paste0(pattern1, "_", pattern2),
          dataset     = dataset
        )[, c("interaction", "dataset", "overlapScore"), drop = FALSE]
        
        df
      },
      pattern_hotspots,
      influence_hotspots,
      datasets,
      ...
    )
  )
  
  wide_df <- reshape2::dcast(combined_df, interaction ~ dataset, value.var = "overlapScore")
  return(wide_df)
}


#' @title Plot overlap scores across samples as a heatmap
#' @description
#' Visualize sample-level overlap scores in a heatmap using a wide overlap
#' data frame, with optional scaling, clustering, column ordering, and
#' numeric labels.
#'
#' @param wide_df Data frame with column interaction and one or more
#' additional columns for datasets containing overlap scores.
#' @param scale Character indicating if values should be scaled by row,
#' column, or not at all. One of "none", "row", or "column".
#' @param cluster_rows Logical indicating whether to cluster rows.
#' @param cluster_cols Logical indicating whether to cluster columns.
#' @param color Optional vector of colors to use for the heatmap.
#' @param na_color Color for missing values.
#' @param display_numbers Logical indicating whether to overlay numeric
#' values on the heatmap.
#' @param number_format Format string for numeric labels when
#' display_numbers is TRUE.
#' @param fontsize Numeric base font size for the heatmap.
#' @param fontsize_row Numeric font size for row labels.
#' @param fontsize_col Numeric font size for column labels.
#' @param main Character title for the heatmap.
#' @param filename Character file name for the output image.
#' @param col_order Optional character vector specifying the order of
#' dataset columns to plot. Columns not listed are omitted.
#' @param ... Additional arguments passed to pheatmap.
#'
#' @return Invisibly returns the input wide_df.
#' @export
plot_sample_overlaps <- function(wide_df,
                                 scale = c("none", "row", "column"),
                                 cluster_rows = TRUE,
                                 cluster_cols = TRUE,
                                 color = NULL,
                                 na_color = "grey90",
                                 display_numbers = FALSE,
                                 number_format = "%.2f",
                                 fontsize = 9,
                                 fontsize_row = 7,
                                 fontsize_col = 9,
                                 main = "Overlap scores across samples",
                                 filename = "SampleOverlapScores.png",
                                 col_order = NULL,
                                 ...) {
  scale <- match.arg(scale)
  
  if (ncol(wide_df) < 2L)
    stop("wide_df must have ≥2 columns: interaction + ≥1 dataset.")
  if (!"interaction" %in% names(wide_df))
    stop("wide_df must contain column 'interaction'.")
  if (nrow(wide_df) == 0L) {
    warning("No interactions left after filtering.")
    return(wide_df)
  }
  if (anyDuplicated(wide_df$interaction)) {
    stop("Duplicate 'interaction' rows remain after filtering; make them unique before plotting.")
  }
  
  plot_df <- wide_df
  rownames(plot_df) <- plot_df$interaction
  plot_df$interaction <- NULL
  
  # Coerce to numeric
  for (j in seq_len(ncol(plot_df))) {
    plot_df[[j]] <- as.numeric(plot_df[[j]])
  }
  
  # ---- Optional column order override ----
  if (!is.null(col_order)) {
    # Check that requested columns exist
    missing_cols <- setdiff(col_order, colnames(plot_df))
    if (length(missing_cols) > 0L) {
      stop("These col_order columns are not in wide_df: ",
           paste(missing_cols, collapse = ", "))
    }
    
    # Warn if some columns will be dropped (not included in col_order)
    extra_cols <- setdiff(colnames(plot_df), col_order)
    if (length(extra_cols) > 0L) {
      warning("These columns in wide_df are not in col_order and will be omitted: ",
              paste(extra_cols, collapse = ", "))
    }
    
    # Reorder / subset columns
    plot_df <- plot_df[, col_order, drop = FALSE]
    
    # If user is explicitly specifying order, disable column clustering
    cluster_cols <- FALSE
  }
  
  # ---- Colors ----
  if (is.null(color)) {
    if (scale == "none") {
      color <- c("#FFF7EC", "#FDBB84", "#D7301F")
    } else {
      color <- c("#2166AC", "#F7F7F7", "#B2182B")
    }
  }
  
  # ---- Optional numbers overlay ----
  num_mat <- NULL
  if (display_numbers) {
    m <- as.matrix(plot_df)
    num_mat <- matrix(
      sprintf(number_format, m),
      nrow = nrow(m),
      ncol = ncol(m),
      dimnames = dimnames(m)
    )
  }
  
  # ---- Heatmap ----
  pheatmap::pheatmap(
    plot_df,
    scale            = scale,
    cluster_rows     = cluster_rows,
    cluster_cols     = cluster_cols,
    color            = color,
    na_col           = na_color,
    display_numbers  = num_mat,
    fontsize         = fontsize,
    fontsize_row     = fontsize_row,
    fontsize_col     = fontsize_col,
    main             = main,
    filename         = filename,
    ...
  )
  
  invisible(wide_df)
}



#' @title Compute pairwise overlap differences
#' @description
#' Given a wide overlap matrix with one column for \code{interaction} and the
#' remaining columns as datasets, compute all pairwise differences between
#' datasets in long format. Optionally, map datasets to higher-level groups
#' and compute group-level summaries.
#'
#' @param mtx A data frame where the column given by \code{feature_col}
#'   contains interaction names and the remaining columns are datasets with
#'   numeric overlap scores.
#' @param sample_groups Optional named character vector mapping dataset
#'   names (column names of \code{mtx} excluding \code{feature_col}) to
#'   user-defined groups or conditions. If \code{NULL}, only dataset-level
#'   comparisons are returned.
#' @param feature_col Character string giving the name of the interaction
#'   column in \code{mtx}. Default is \code{"interaction"}.
#'
#' @return A data frame with one row per dataset pair and interaction,
#'   including columns for dataset-level differences. If \code{sample_groups}
#'   is not \code{NULL}, additional columns \code{comparison_group},
#'   \code{group1}, \code{group2}, \code{median_diff}, and
#'   \code{median_overlapScore} are added.
#' @export
compare_scores <- function(mtx,sample_groups = NULL,feature_col = "interaction") {
  # wide matrix from df
  rows     <- mtx[[feature_col ]]
  wide_mat <- as.matrix(mtx[, -1, drop = FALSE])
  rownames(wide_mat) <- rows
  
  # all pairwise dataset combinations
  datasets <- colnames(wide_mat)
  combs    <- t(combn(datasets, 2))
  
  # build two wide matrices for each side of the pair
  wide_mat1 <- t(wide_mat[, combs[, 1], drop = FALSE])
  wide_mat2 <- t(wide_mat[, combs[, 2], drop = FALSE])
  # long format for dataset-level differences
  overlap1 <- reshape2::melt(
    wide_mat1,
    varnames  = c("dataset1", "interaction"),
    value.name = "dataset1_overlapScore"
  ) %>% dplyr::pull(dataset1_overlapScore)
  overlap2 <- reshape2::melt(
    wide_mat2,
    varnames  = c("dataset2", "interaction"),
    value.name = "dataset2_overlapScore"
  ) %>% dplyr::pull(dataset2_overlapScore)
  
  # raw differences
  diff_mat <- wide_mat1 - wide_mat2
  
  # rownames encode dataset1_minus_dataset2
  diff_rows <- paste0(rownames(wide_mat1), "_minus_", rownames(wide_mat2))
  rownames(diff_mat) <- diff_rows
  
  # long format for dataset-level differences
  diff_df <- reshape2::melt(
    diff_mat,
    varnames  = c("comparison", "interaction"),
    value.name = "diff_overlapScore"
  )
  diff_df$comparison  <- as.character(diff_df$comparison)
  diff_df$interaction <- as.character(diff_df$interaction)
  
  # split comparison into dataset1 / dataset2
  split_data <- do.call("rbind", strsplit(diff_df$comparison, "_minus_"))
  colnames(split_data) <- c("dataset1", "dataset2")
  
  final_df <- cbind(cbind(diff_df, split_data, stringsAsFactors = FALSE),
                    dataset1_overlapScore = overlap1,
                    dataset2_overlapScore = overlap2)
  
  # if sample_groups is provided, compute group-level comparison labels
  if (!is.null(sample_groups)) {
    # align sample_groups to datasets used in wide_mat1 / wide_mat2
    sample_groups <- sample_groups[rownames(wide_mat1)]
    
    diff_rows_group <- paste0(
      sample_groups[rownames(wide_mat1)],
      "_minus_",
      sample_groups[rownames(wide_mat2)]
    )
    rownames(diff_mat) <- diff_rows_group
    
    # long format for group-level comparisons
    diff_df2 <- reshape2::melt(
      diff_mat,
      varnames  = c("comparison_group", "interaction"),
      value.name = "diff_overlapScore"
    )
    diff_df2$comparison_group <- as.character(diff_df2$comparison_group)
    diff_df2$interaction      <- as.character(diff_df2$interaction)
    
    split_data2 <- do.call("rbind", strsplit(diff_df2$comparison_group, "_minus_"))
    colnames(split_data2) <- c("group1", "group2")
    
    diff_df2 <- cbind(
      diff_df2[, c("comparison_group", "interaction")],
      split_data2,
      stringsAsFactors = FALSE
    )
    
    # combine dataset-level and group-level info (same ordering from melt)
    final_df <- cbind(final_df, diff_df2[, c("comparison_group", "group1",
                                             "group2")])
    
    # add median_diff per (group1, interaction)
    final_df <- final_df %>%
      dplyr::group_by(group1, interaction) %>%
      dplyr::mutate(median_diff = median(diff_overlapScore, na.rm = TRUE)) %>%
      dplyr::ungroup()
    final_df <- final_df %>%
      dplyr::group_by(group1, interaction) %>%
      dplyr::mutate(median_overlapScore = median(dataset1_overlapScore, na.rm = TRUE)) %>%
      dplyr::ungroup()
  }
  
  final_df
}


#' @title Plot group-level median overlap scores
#' @description
#' Create a horizontal bar plot of median overlap scores for the top
#' interactions, typically using the output of compare_scores. Interactions
#' are ranked by the absolute median difference between groups.
#'
#' @param df Data frame containing at least the columns interaction, group1,
#' median_diff, and median_overlapScore.
#' @param top Integer specifying the number of top interactions to display.
#' @param out Character file name for the saved TIFF figure.
#'
#' @return An invisible list with components data (the plotted subset) and
#' plot (the ggplot object).
#' @export
plot_overlap_scores_bar <- function(df,
                                    top = 10,
                                    out = "median_overlap_scores.tiff") {
  
  df <- df %>%
    dplyr::select(interaction, group1, median_diff, median_overlapScore) %>%
    dplyr::distinct()
  
  # helper: remove self interactions like B_near_B, Strom_to_Strom, B_B, Strom_Strom
  drop_self <- function(x) {
    x[
      {
        left1  <- sub("^(.*)_(near|to)_.*$", "\\1", x)
        right1 <- sub("^.*_(near|to)_", "", x)
        is_near_to <- grepl("_(near|to)_", x)
        
        left2  <- sub("_.*$", "", x)
        right2 <- sub("^.*_", "", x)
        is_uscore <- grepl("_", x) & !is_near_to
        
        keep1 <- !is_near_to | (left1 != right1)
        keep2 <- !is_uscore  | (left2 != right2)
        
        keep1 & keep2
      }
    ]
  }
  
  # rank interactions per group
  ranked_by_group <- df %>%
    dplyr::arrange(group1, dplyr::desc(abs(median_diff))) %>%
    dplyr::group_by(group1) %>%
    dplyr::summarise(ranked = list(drop_self(unique(interaction))), .groups = "drop")
  
  # build ordered_levels so each group contributes up to `top` *new* interactions
  used <- character(0)
  ordered_levels <- character(0)
  top_table_list <- vector("list", nrow(ranked_by_group))
  
  for (i in seq_len(nrow(ranked_by_group))) {
    g <- ranked_by_group$group1[i]
    cand <- ranked_by_group$ranked[[i]]
    new_for_group <- cand[!cand %in% used]
    
    # take top new ones; if not enough exist, you'll get fewer than `top`
    pick <- head(new_for_group, top)
    
    ordered_levels <- c(ordered_levels, pick)
    used <- c(used, pick)
    
    top_table_list[[i]] <- data.frame(
      group1 = g,
      top_interactions = I(list(pick)),
      group_order = i,
      stringsAsFactors = FALSE
    )
  }
  
  top_tbl <- dplyr::bind_rows(top_table_list)
  
  # keep only those interactions (union across groups, but "filled" to preserve top-per-group intent)
  median_df <- df %>%
    dplyr::mutate(group = group1) %>%
    dplyr::select(group, interaction, median_overlapScore) %>%
    dplyr::distinct() %>%
    dplyr::filter(interaction %in% unique(ordered_levels))
  
  # y-axis order: group1 block, then group2 block, etc.
  median_df$interaction <- factor(
    median_df$interaction,
    levels = rev(unique(ordered_levels))
  )
  
  p <- ggplot2::ggplot(
    median_df,
    ggplot2::aes(x = median_overlapScore, y = interaction, fill = group)
  ) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.8)) +
    ggplot2::labs(x = "Median overlap score", y = "Interaction", fill = "Group") +
    ggplot2::theme_classic()
  
  ggplot2::ggsave(
    filename    = out,
    plot        = p,
    device      = "tiff",
    width       = 8,
    height      = 6,
    units       = "in",
    dpi         = 300,
    compression = "lzw"
  )
  
  invisible(list(data = median_df, plot = p, top_table = top_tbl, ordered_levels = ordered_levels))
}



#' @title Run directed IMscores for multiple samples
#' @description
#' Run SpaceMarkers::calculate_gene_scores_directed for each sample after
#' aligning sample names and enforcing a global feature intersection across
#' spPatternsList. Optionally save IMscores and timing information.
#'
#' @param spCountsList Named list of count matrices (genes x barcodes).
#' @param spPatternsList Named list of pattern tables with barcode, x, y,
#' and feature columns.
#' @param spHotspotsList Named list of pattern hotspot objects.
#' @param spHotspotsInfluenceList Named list of influence hotspot objects.
#' @param sampleNames Optional character vector of sample names to run. If NULL,
#' use the intersection of names across all lists.
#' @param outDir Optional output directory for RDS files. If NULL, nothing
#' is written to disk.
#' @param imscoresSuffix Character. Suffix for IMscores RDS file names.
#' @param elapsedFilename Character. File name for the elapsed seconds RDS.
#' @param verbose Logical. If TRUE, print progress and timing information.
#'
#' @return A list with:
#' \describe{
#'   \item{imscores}{Named list of IMscores objects per sample.}
#'   \item{elapsed_sec}{Named numeric vector of elapsed seconds per sample.}
#'   \item{paths}{Optional list of file paths if outDir is provided.}
#' }
#'
#' @importFrom utils combn
#' @export
run_imscores_one <- function(
    spCountsList,
    spPatternsList,
    spHotspotsList,
    spHotspotsInfluenceList,
    sampleNames      = NULL,
    outDir           = NULL,
    imscoresSuffix   = "_IMscores_directed.rds",
    elapsedFilename  = "IMscoresHD_elapsed_seconds.rds",
    verbose          = TRUE
) {
  common_all <- Reduce(intersect, list(
    names(spCountsList),
    names(spPatternsList),
    names(spHotspotsList),
    names(spHotspotsInfluenceList)
  ))
  if (!length(common_all)) stop("No common sample names across input lists.")
  
  if (is.null(sampleNames)) {
    common <- common_all
  } else {
    common <- intersect(sampleNames, common_all)
    if (!length(common))
      stop("None of `sampleNames` found across all input lists.")
  }
  
  spCountsList            <- spCountsList[common]
  spPatternsList          <- spPatternsList[common]
  spHotspotsList          <- spHotspotsList[common]
  spHotspotsInfluenceList <- spHotspotsInfluenceList[common]
  
  global_features <- Reduce(
    intersect,
    lapply(spPatternsList, function(df) setdiff(colnames(df), c("barcode","x","y")))
  )
  if (length(global_features) < 2) {
    stop("Global feature intersection resulted in <2 features; cannot compute IMscores.")
  }
  if (verbose) message("Using global feature set (n = ", length(global_features), ").")
  
  if (!is.null(outDir) && !dir.exists(outDir)) {
    dir.create(outDir, recursive = TRUE, showWarnings = FALSE)
  }
  
  imscores_list <- setNames(vector("list", length(common)), common)
  imscores_paths <- if (!is.null(outDir)) setNames(character(length(common)), common) else NULL
  elapsed_sec <- setNames(numeric(length(common)), common)
  
  for (s in common) {
    if (verbose) message("\n=== IMscoresHD: ", s, " ===")
    
    feat_cols <- intersect(
      setdiff(colnames(spPatternsList[[s]]), c("barcode","x","y")),
      global_features
    )
    if (length(feat_cols) < 2) {
      stop(sprintf("[%s] <2 global features available; cannot form pairs.", s))
    }
    
    spPat <- spPatternsList[[s]]
    if (!"barcode" %in% names(spPat)) stop(sprintf("[%s] spPatterns must have a 'barcode' column.", s))
    rownames(spPat) <- spPat$barcode
    barcodes <- intersect(colnames(spCountsList[[s]]), rownames(spPat))
    if (!length(barcodes)) stop(sprintf("[%s] No overlapping barcodes between counts and patterns.", s))
    data_counts <- spCountsList[[s]][, barcodes, drop = FALSE]
    
    patternPairs <- t(utils::combn(feat_cols, 2))
    
    t0 <- proc.time()
    ims <- SpaceMarkers::calculate_gene_scores_directed(
      data         = data_counts,
      pat_hotspots = spHotspotsList[[s]],
      influence_hotspots = spHotspotsInfluenceList[[s]],
      pattern_pairs = patternPairs
    )
    t1 <- proc.time()
    elapsed_sec[[s]] <- as.numeric((t1 - t0)["elapsed"])
    
    imscores_list[[s]] <- ims
    
    if (!is.null(outDir)) {
      ims_path <- file.path(outDir, paste0(s, imscoresSuffix))
      saveRDS(ims, ims_path)
      imscores_paths[[s]] <- ims_path
      if (verbose) message("Saved: ", basename(ims_path))
    }
  }
  
  elapsed_path <- NULL
  if (!is.null(outDir)) {
    elapsed_path <- file.path(outDir, elapsedFilename)
    saveRDS(elapsed_sec, elapsed_path)
    if (verbose) message("\nSaved elapsed seconds vector: ", basename(elapsed_path))
  }
  
  out <- list(
    imscores    = imscores_list,
    elapsed_sec = elapsed_sec
  )
  if (!is.null(outDir)) {
    out$paths <- list(
      imscores     = imscores_paths,
      elapsed_path = elapsed_path
    )
  }
  out
}


# ---- top-level worker (ONE sample) ----
.lr_worker <- function(test, cnt, pat,
                       ligands_list, receptors_list, lrpairs,
                       weighted, ligand_method, receptor_method, lr_method) {
  ligand_scores <- calculate_gene_set_score(
    test, gene_sets = ligands_list, weighted = weighted, method = ligand_method
  )
  receptor_scores <- calculate_gene_set_specificity(
    cnt, pat, gene_sets = receptors_list, weighted = weighted, method = receptor_method
  )
  lr_scores <- calculate_lr_scores(
    ligand_scores, receptor_scores, lr_pairs = lrpairs,
    method = lr_method, weighted = weighted
  )
  # return all three
  list(
    lr_scores       = lr_scores,
    ligand_scores   = ligand_scores,
    receptor_scores = receptor_scores
  )
}

#' @title Parallel LR scoring over samples
#' @description
#' Compute ligand, receptor, and ligand-receptor scores for each sample in
#' parallel, given IMscores, counts, patterns, and ligand/receptor gene sets.
#'
#' @param IMscores_directed_list Named list of matrices used as input for
#' ligand scoring.
#' @param spCounts_list Named list of per-sample count matrices.
#' @param spPatterns_list Named list of per-sample pattern tables.
#' @param ligands_list Named list of ligand gene sets.
#' @param receptors_list Named list of receptor gene sets.
#' @param lrpairs Data.frame of ligand-receptor pairs.
#' @param workers Integer. Number of parallel workers.
#' @param weighted Logical. If TRUE, use weighted scoring.
#' @param ligand_method Character. Method used to aggregate ligand scores.
#' @param receptor_method Character. Method used to aggregate receptor scores.
#' @param lr_method Character. Method used to combine ligand and receptor
#' scores.
#'
#' @return A list with three named lists:
#' res_list (LR score matrices), ligand_scores_list, and receptor_scores_list.
#'
#' @export
generate_sample_lr_scores <- function(
    IMscores_directed_list,
    spCounts_list,
    spPatterns_list,
    ligands_list,
    receptors_list,
    lrpairs,
    workers        = 5,
    weighted       = TRUE,
    ligand_method  = "arithmetic_mean",
    receptor_method= "arithmetic_mean",
    lr_method      = "geometric_mean"
) {
  if (!requireNamespace("future.apply", quietly = TRUE)) {
    stop("Please install 'future.apply' (install.packages('future.apply')).")
  }
  n <- length(IMscores_directed_list)
  if (length(spCounts_list) != n || length(spPatterns_list) != n) {
    stop("spCounts_list, spPatterns_list, and IMscores_directed_list must have the same length.")
  }
  samp_names <- names(IMscores_directed_list)
  if (is.null(samp_names) || any(samp_names == "")) {
    samp_names <- paste0("sample", seq_len(n))
  }
  
  old_plan <- future::plan(); on.exit(future::plan(old_plan), add = TRUE)
  future::plan(future::multisession, workers = workers)
  
  res_per_sample <- future.apply::future_mapply(
    FUN  = .lr_worker,
    test = IMscores_directed_list,
    cnt  = spCounts_list,
    pat  = spPatterns_list,
    MoreArgs = list(
      ligands_list    = ligands_list,
      receptors_list  = receptors_list,
      lrpairs         = lrpairs,
      weighted        = weighted,
      ligand_method   = ligand_method,
      receptor_method = receptor_method,
      lr_method       = lr_method
    ),
    SIMPLIFY   = FALSE,
    future.seed = TRUE
  )
  
  res_list <- lapply(res_per_sample, `[[`, "lr_scores")
  ligand_scores_list <- lapply(res_per_sample, `[[`, "ligand_scores")
  receptor_scores_list <- lapply(res_per_sample, `[[`, "receptor_scores")
  
  names(res_list) <- samp_names
  names(ligand_scores_list) <- samp_names
  names(receptor_scores_list) <- samp_names
  
  list(
    res_list             = res_list,
    ligand_scores_list   = ligand_scores_list,
    receptor_scores_list = receptor_scores_list
  )
}


#' @title Create ligand–receptor score table across multiple datasets
#'
#' @description Applies \code{create_lr_dataframe} to each dataset in a score-list object, then
#' combines results into one data frame and adds dataset-level metadata.
#'
#' @param lr_scores_list A list containing \code{res_list} (LR score matrices),
#' \code{ligand_scores_list} (ligand score matrices), and \code{receptor_scores_list}
#' (receptor score matrices). Each must be a named list with matching names.
#' @param lrpairs A data frame with columns \code{ligand} and \code{receptor} where
#' row names are the ligand–receptor interaction IDs used to join to scores.
#' @param sample_groups Optional named character vector mapping dataset names
#' (names of \code{lr_scores_list$res_list}) to a condition label. If \code{NULL},
#' \code{condition} is set to \code{NA}.
#' @param complex_sep Character separator used for ligand/receptor complexes, passed
#' to \code{create_lr_dataframe}.
#'
#' @return A data frame combining all per-dataset outputs of \code{create_lr_dataframe},
#' with two additional columns: \code{dataset} (dataset name) and \code{condition}
#' (mapped from \code{sample_groups}, or \code{NA}).
#' 
#' @export
create_lr_dataframe_multi <- function(lr_scores_list, lrpairs, sample_groups = NULL, complex_sep = ", ") {
  lrscores_list        <- lr_scores_list$res_list
  ligand_scores_list   <- lr_scores_list$ligand_scores_list
  receptor_scores_list <- lr_scores_list$receptor_scores_list
  
  # basic checks
  stopifnot(is.list(lrscores_list), is.list(ligand_scores_list), is.list(receptor_scores_list))
  stopifnot(length(lrscores_list) == length(ligand_scores_list),
            length(lrscores_list) == length(receptor_scores_list))
  
  datasets <- names(lrscores_list)
  if (is.null(datasets) || any(datasets == "")) {
    stop("lr_scores_list$res_list must be a *named* list so we can populate the 'dataset' column.")
  }
  
  out_list <- Map(function(lrs, lig, rec, ds) {
    df <- create_lr_dataframe(
      lrscores        = lrs,
      ligand_scores   = lig,
      receptor_scores = rec,
      lrpairs         = lrpairs,
      complex_sep     = complex_sep
    )
    
    df$dataset <- ds
    df$condition <- if (!is.null(sample_groups)) unname(sample_groups[ds]) else NA_character_
    df
  }, lrscores_list, ligand_scores_list, receptor_scores_list, datasets)
  
  final_lrscores<- dplyr::bind_rows(out_list)
  final_lrscores$source_to_target <- paste0(final_lrscores$source_cell_type,
                                    "_to_",final_lrscores$target_cell_type)
  return(final_lrscores)
}


#' @title Get top ligand–receptor interactions per condition
#'
#' @description Identifies the highest-scoring ligand–receptor interactions within each condition
#' by ranking interactions based on the median LR score across datasets. The function first removes
#' self-interactions where ligand and receptor are identical, then computes a per-condition median
#' score for each (source–target, interaction) pair. The top interactions per condition are selected
#' and returned with all original per-dataset rows retained, along with the computed median score.
#'
#' @param df A data frame of ligand–receptor scores containing at least the columns
#' \code{ligand}, \code{receptor}, \code{score}, \code{dataset}, \code{condition},
#' \code{interaction}, and \code{source_to_target}.
#' @param top Integer specifying the number of top interactions to return per condition
#' based on median score.
#'
#' @return A data frame containing all rows from the original input corresponding to the
#' top-ranked interactions per condition, with two additional columns:
#' \code{median_score} (median LR score across datasets) and \code{n_datasets}
#' (number of datasets contributing to the interaction).
#'
#' @export
get_top_lr_interactions <- function(df, top = 10) {
  # 1) drop ligand == receptor
  df2 <- df %>%
    dplyr::filter(is.na(ligand) | is.na(receptor) | ligand != receptor)
  
  # 2) median score per (condition, source_to_target, interaction)
  med_tbl <- df2 %>%
    dplyr::group_by(condition, source_to_target, interaction) %>%
    dplyr::summarise(
      median_score = median(score, na.rm = TRUE),
      n_datasets   = dplyr::n_distinct(dataset),
      .groups = "drop"
    )
  
  # 3) pick top N per condition (tie-broken by n_datasets)
  top_keys <- med_tbl %>%
    dplyr::arrange(condition, dplyr::desc(median_score), dplyr::desc(n_datasets)) %>%
    dplyr::group_by(condition) %>%
    dplyr::slice_head(n = top) %>%
    dplyr::ungroup()
  
  # 4) return full rows from original df + median_score (+ n_datasets)
  out <- df2 %>%
    dplyr::inner_join(
      top_keys,
      by = c("condition", "source_to_target", "interaction")
    ) %>%
    dplyr::arrange(condition, dplyr::desc(median_score), dplyr::desc(score))
  
  out
}

#' @title Plot ligand receptor alluvial diagram
#'
#' @description
#' Creates an alluvial (Sankey-style) plot showing flows from source cell type
#' to ligand, receptor, and target cell type. Flow width is proportional to an
#' interaction score, and flows are colored by the source cell type.
#'
#' @param df A data frame containing ligand–receptor scores and the required
#' columns for source, ligand, receptor, target, and score.
#' @param score_col Character. Name of the score column. Default is \code{"score"}.
#' @param source_col Character. Name of the source cell type column. Default is
#' \code{"source_cell_type"}.
#' @param target_col Character. Name of the target cell type column. Default is
#' \code{"target_cell_type"}.
#' @param ligand_col Character. Name of the ligand column. Default is \code{"ligand"}.
#' @param receptor_col Character. Name of the receptor column. Default is \code{"receptor"}.
#' @param title Optional plot title. If \code{NULL}, a default title is used.
#' @param min_score Numeric. Minimum score to keep before plotting.
#' @param alpha Numeric between 0 and 1. Alpha transparency for flows.
#' @param knot_pos Numeric between 0 and 1. Knot position along the flow.
#' @param label_size Numeric. Text size for node labels.
#' @param wrap_width Integer or \code{NULL}. If not \code{NULL}, wraps labels to
#' this width.
#'
#' @return A ggplot object.
#'
#' @export
plot_lr_alluvial <- function(
    df,
    score_col    = "score",
    source_col   = "source_cell_type",
    target_col   = "target_cell_type",
    ligand_col   = "ligand",
    receptor_col = "receptor",
    title        = NULL,
    min_score    = 0,
    alpha        = 0.8,
    knot_pos     = 0.4,
    label_size   = 2.7,
    wrap_width   = 12
) {
  for (pkg in c("ggplot2","ggalluvial","RColorBrewer","circlize")) {
    if (!requireNamespace(pkg, quietly = TRUE)) stop(sprintf("Please install '%s'.", pkg))
  }
  
  # resolve columns (case tolerant)
  .res <- function(nm) {
    hit <- which(tolower(names(df)) == tolower(nm))
    if (length(hit)) names(df)[hit[1]] else nm
  }
  score_col    <- .res(score_col)
  source_col   <- .res(source_col)
  target_col   <- .res(target_col)
  ligand_col   <- .res(ligand_col)
  receptor_col <- .res(receptor_col)
  
  need <- c(score_col, source_col, target_col, ligand_col, receptor_col)
  miss <- setdiff(need, names(df))
  if (length(miss)) stop("Missing columns: ", paste(miss, collapse = ", "))
  
  d <- df
  d[[score_col]] <- as.numeric(d[[score_col]])
  d <- d[is.finite(d[[score_col]]) & d[[score_col]] >= min_score, , drop = FALSE]
  if (!nrow(d)) stop("No rows to plot after filtering.")
  
  # wrap labels only (no abbreviation)
  wrap_names <- function(v) {
    if (is.null(wrap_width)) return(v)
    vapply(v, function(s) paste(strwrap(s, width = wrap_width), collapse = "\n"), "")
  }
  d[[source_col]]   <- wrap_names(d[[source_col]])
  d[[ligand_col]]   <- wrap_names(d[[ligand_col]])
  d[[receptor_col]] <- wrap_names(d[[receptor_col]])
  d[[target_col]]   <- wrap_names(d[[target_col]])
  
  # color by source cell type
  d$fill_key <- factor(d[[source_col]])
  
  .discrete_palette <- function(n) {
    if (n <= 8) {
      RColorBrewer::brewer.pal(max(3, n), "Set2")[1:n]
    } else if (n <= 12) {
      RColorBrewer::brewer.pal(max(3, n), "Set3")[1:n]
    } else {
      base <- RColorBrewer::brewer.pal(12, "Set3")
      colfun <- circlize::colorRamp2(seq(0, 1, length.out = length(base)), base)
      colfun(seq(0, 1, length.out = n))
    }
  }
  pal_vals <- .discrete_palette(nlevels(d$fill_key))
  names(pal_vals) <- levels(d$fill_key)
  
  dp <- data.frame(
    source   = d[[source_col]],
    ligand   = d[[ligand_col]],
    receptor = d[[receptor_col]],
    target   = d[[target_col]],
    score    = d[[score_col]],
    fill_key = d$fill_key,
    stringsAsFactors = FALSE
  )
  
  default_title <- "LR Alluvial: Source -> Ligand -> Receptor -> Target"
  
  ggplot2::ggplot(
    dp,
    ggplot2::aes(axis1 = source, axis2 = ligand, axis3 = receptor, axis4 = target, y = score)
  ) +
    ggplot2::scale_x_discrete(
      limits = c("Source", "Ligand", "Receptor", "Target"),
      expand = c(.08, .05)
    ) +
    ggalluvial::geom_alluvium(
      ggplot2::aes(fill = fill_key),
      alpha = alpha,
      knot.pos = knot_pos,
      show.legend = TRUE
    ) +
    ggalluvial::geom_stratum(width = 0.25, color = "grey30", fill = "grey85") +
    ggalluvial::stat_stratum(
      geom = "text",
      ggplot2::aes(label = after_stat(stratum)),
      size = label_size,
      check_overlap = TRUE
    ) +
    ggplot2::scale_fill_manual(values = pal_vals, name = "Source cell type", drop = FALSE) +
    ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(alpha = 1))) +
    ggplot2::labs(
      title    = if (is.null(title)) default_title else title,
      subtitle = "Flow width ~ score; color ~ source_cell_type",
      y = "Score",
      x = NULL
    ) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      legend.position = "right",
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(size = 11),
      plot.margin = ggplot2::margin(10, 30, 10, 30)
    )
}



#' @title Plot ligand vs receptor scores scatter
#'
#' @description Creates a scatter plot of ligand scores versus receptor scores and
#' saves it to a TIFF file. Points are colored by a user-specified column and shaped
#' by a user-specified column, using ggplot2 default discrete color scales.
#'
#' @param df A data frame containing at least \code{ligand_score} and
#' \code{receptor_score}, plus the columns specified by \code{color_by} and
#' \code{shape_by}.
#' @param out Output file path for the saved plot (should end in \code{.tiff}).
#' @param color_by Character. Column name used for point color. Default is
#' \code{"interaction"}.
#' @param shape_by Character. Column name used for point shape. Default is
#' \code{"source_to_target"}.
#' @param point_size Numeric. Point size.
#' @param point_alpha Numeric between 0 and 1. Point transparency.
#' @param width Numeric. Plot width in inches.
#' @param height Numeric. Plot height in inches.
#' @param dpi Integer. Plot resolution.
#'
#' @return Invisibly returns the output file path \code{out}.
#'
#' @export
plot_lr_scores_scatter <- function(df,
                                   out = "lr_scores_scatter.tiff",
                                   color_by = "interaction",
                                   shape_by = "source_to_target",
                                   point_size = 2,
                                   point_alpha = 0.8,
                                   width = 8,
                                   height = 6,
                                   dpi = 300) {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Please install 'ggplot2'.")
  
  needed_cols <- c("ligand_score", "receptor_score", color_by, shape_by)
  missing_cols <- setdiff(needed_cols, colnames(df))
  if (length(missing_cols)) {
    stop("Missing required columns in 'df': ",
         paste(missing_cols, collapse = ", "))
  }
  
  out_dir <- dirname(out)
  if (!dir.exists(out_dir))
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  p <- ggplot2::ggplot(
    df,
    ggplot2::aes_string(
      x = "ligand_score",
      y = "receptor_score",
      color = color_by,
      shape = shape_by
    )
  ) +
    ggplot2::geom_point(size = point_size, alpha = point_alpha) +
    ggplot2::labs(
      title = "Ligand vs Receptor scores",
      x = "Ligand score",
      y = "Receptor score",
      color = color_by,
      shape = shape_by
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      legend.position = "right",
      panel.grid.minor = ggplot2::element_blank()
    )
  
  ggplot2::ggsave(
    filename = out,
    plot = p,
    device = "tiff",
    width = width,
    height = height,
    dpi = dpi,
    units = "in",
    compress = "lzw"
  )
  
  invisible(out)
}


