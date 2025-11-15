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
#' @keywords internal
align_counts_coords <- function(counts, coords_df) {
  common <- intersect(colnames(counts), coords_df$barcode)
  if (length(common) == 0L) stop("No overlapping barcodes between counts and coords.")
  counts    <- counts[, common, drop = FALSE]
  coords_df <- coords_df[match(common, coords_df$barcode), , drop = FALSE]
  rownames(coords_df) <- coords_df$barcode
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

#' @title Compute multi-sample overlap scores
#' @description
#' Compute pairwise overlap scores per sample and assemble long and wide
#' summaries across samples.
#'
#' @param pat_hotspots_list Named list of hotspot tables for each sample.
#' @param in_hotspots_list Optional named list of influence hotspot tables
#' for each sample. If NULL, only pat_hotspots_list is used.
#' @param patternList Optional character vector of pattern names to include.
#' @param method Character. Overlap metric to use. One of
#' "Szymkiewicz-Simpson", "Jaccard", "Sorensen-Dice", "Ochiai", or "absolute".
#'
#' @return A list with two data.frames:
#' \describe{
#'   \item{long}{Long-format table with one row per dataset and interaction.}
#'   \item{wide}{Wide-format table with interactions as rows and datasets as columns.}
#' }
#'
#' @export
get_multi_sample_overlaps <- function(pat_hotspots_list,
                                     in_hotspots_list = NULL,
                                     patternList = NULL,
                                     method = c("Szymkiewicz-Simpson",
                                                "Jaccard",
                                                "Sorensen-Dice",
                                                "Ochiai",
                                                "absolute")) {
  # ---- method guard ----
  if (length(method) > 1) {
    method <- method[1]
    message("Only one method can be used at a time. Using ", method)
  }
  
  if (is.null(pat_hotspots_list) || !length(pat_hotspots_list)) {
    stop("pat_hotspots_list must be a non-empty named list.")
  }
  if (is.null(names(pat_hotspots_list)) || any(!nzchar(names(pat_hotspots_list)))) {
    stop("pat_hotspots_list must be a *named* list (sample names as names()).")
  }
  
  if (!is.null(in_hotspots_list)) {
    common_names <- intersect(names(pat_hotspots_list), names(in_hotspots_list))
    if (!length(common_names)) {
      stop("No overlapping sample names between pat_hotspots_list and in_hotspots_list.")
    }
    pat_hotspots_list <- pat_hotspots_list[common_names]
    in_hotspots_list  <- in_hotspots_list [common_names]
  } else {
    common_names <- names(pat_hotspots_list)
  }
  
  per_sample <- lapply(common_names, function(s) {
    if (!is.null(in_hotspots_list)) {
      df <- get_overlap_scores(
        hotspots     = NULL,
        in_hotspots  = in_hotspots_list[[s]],
        pat_hotspots = pat_hotspots_list[[s]],
        patternList  = patternList,
        method       = method
      )
    } else {
      df <- get_overlap_scores(
        hotspots    = pat_hotspots_list[[s]],
        patternList = patternList,
        method      = method
      )
    }
    
    if (!is.null(df) && nrow(df)) {
      df <- df[stats::complete.cases(df), , drop = FALSE]
      df$interaction <- paste0(df$pattern1, "_", df$pattern2)
      df$dataset     <- s
      df <- df[, c("dataset", "pattern1", "pattern2", "interaction", "overlapScore")]
    } else {
      df <- data.frame(
        dataset      = character(0),
        pattern1     = character(0),
        pattern2     = character(0),
        interaction  = character(0),
        overlapScore = numeric(0),
        stringsAsFactors = FALSE
      )
    }
    df
  })
  names(per_sample) <- common_names
  
  long_df <- do.call(rbind, per_sample)
  if (nrow(long_df)) {
    p1 <- as.character(long_df$pattern1)
    p2 <- as.character(long_df$pattern2)
    long_df$._rm_self <- mapply(function(a, b) grepl(paste0("^near_", a, "$"), b),
                                a = p1, b = p2)
    long_df <- long_df[!long_df$._rm_self, , drop = FALSE]
    long_df$._rm_self <- NULL
  }
  
  if (nrow(long_df)) {
    wide_df <- reshape2::dcast(long_df, interaction ~ dataset, value.var = "overlapScore")
  } else {
    wide_df <- data.frame(interaction = character(0), stringsAsFactors = FALSE)
  }
  
  list(long = long_df, wide = wide_df)
}


#' @title Plot overlap scores across samples
#' @description
#' Plot a heatmap of overlap scores where rows are interactions and columns
#' are samples, using pheatmap.
#'
#' @param wide_df A data.frame with column interaction and additional columns
#' for each dataset containing numeric overlap scores.
#' @param scale Character. One of "none", "row", or "column". Passed to pheatmap.
#' @param cluster_rows Logical. If TRUE, cluster rows.
#' @param cluster_cols Logical. If TRUE, cluster columns.
#' @param color Optional color palette. If NULL, a default palette is used.
#' @param na_color Color for NA tiles.
#' @param display_numbers Logical. If TRUE, print values inside tiles.
#' @param number_format Character. sprintf format for numbers, for example "%.2f".
#' @param fontsize Numeric. Base font size.
#' @param fontsize_row Numeric. Font size for row labels.
#' @param fontsize_col Numeric. Font size for column labels.
#' @param main Character. Plot title.
#' @param filename Character. Output file name for pheatmap. If NULL, the plot
#' is returned but not written.
#' @param ... Additional arguments passed to pheatmap::pheatmap.
#'
#' @return The filtered wide_df that was used for plotting.
#'
#' @importFrom pheatmap pheatmap
#' @export
plot_sample_overlaps <- function(wide_df,
                                 scale = c("none","row","column"),
                                 cluster_rows = TRUE,
                                 cluster_cols = TRUE,
                                 color = NULL,
                                 na_color = "grey90",
                                 display_numbers = FALSE,
                                 number_format = "%.2f",
                                 fontsize = 9,
                                 fontsize_row = 7 ,
                                 fontsize_col = 9,
                                 main = "Overlap scores across samples",
                                 filename = "SampleOverlapScores.png",
                                 ...) {
  scale <- match.arg(scale)
  if (ncol(wide_df) < 2L) stop("wide_df must have \u22652 columns: interaction + \u22651 dataset.")
  if (!"interaction" %in% names(wide_df)) stop("wide_df must contain column 'interaction'.")
  
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
  for (j in seq_len(ncol(plot_df))) plot_df[[j]] <- as.numeric(plot_df[[j]])
  
  if (is.null(color)) {
    if (scale == "none") {
      color <- c("#FFF7EC", "#FDBB84", "#D7301F")
    } else {
      color <- c("#2166AC", "#F7F7F7", "#B2182B")
    }
  }
  
  num_mat <- NULL
  if (display_numbers) {
    m <- as.matrix(plot_df)
    num_mat <- matrix(
      sprintf(number_format, m),
      nrow = nrow(m), ncol = ncol(m),
      dimnames = dimnames(m)
    )
  }
  
  pheatmap::pheatmap(
    plot_df,
    scale           = scale,
    cluster_rows    = cluster_rows,
    cluster_cols    = cluster_cols,
    color           = color,
    na_col          = na_color,
    display_numbers = num_mat,
    fontsize        = fontsize,
    fontsize_row    = fontsize_row,
    fontsize_col    = fontsize_col,
    main            = main,
    filename        = filename,
    ...
  )
  
  wide_df
}

#' @title Compare scores between a reference and a condition
#' @description
#' For each feature, compare values between a reference condition and a
#' second condition using a Wilcoxon rank-sum test. Add p values and
#' adjusted p values to the input data.
#'
#' @param df Data.frame containing at least the group, feature, and value columns.
#' @param ref_condition Character. Name of the reference condition.
#' @param group_col Character. Column name for the group variable.
#' @param feature_col Character. Column name for the feature identifier.
#' @param value_col Character. Column name for the numeric values to compare.
#'
#' @return A data.frame with the same rows as df, grouped by feature, and
#' additional columns p_value, fdr, ref_median, and condition_median.
#'
#' @importFrom dplyr filter pull
#' @importFrom rlang sym
#' @importFrom stats wilcox.test p.adjust
#' @export
compare_scores <- function(
    df,
    ref_condition   = "Normal",
    group_col       = "group",
    feature_col     = "interaction",
    value_col       = "overlapScore"
) {
  overlaps_list <- split(x = df, f = df[[feature_col]])
  both_conditions <- unique(df[[group_col]])
  conditions <- factor(
    df[[group_col]],
    levels = c(ref_condition, setdiff(both_conditions, c(ref_condition)))
  )
  other <- levels(conditions)[2]
  df[[group_col]] <- conditions
  
  test_wrapper <- function(odf){
    cd1_df <- odf %>% dplyr::filter(!!rlang::sym(group_col) == ref_condition)
    cd2_df <- odf %>% dplyr::filter(!!rlang::sym(group_col) == other)
    
    tt <- stats::wilcox.test(
      x = cd1_df %>% dplyr::pull(!!rlang::sym(value_col)),
      y = cd2_df %>% dplyr::pull(!!rlang::sym(value_col)),
      alternative = "two.sided"
    )
    odf$p_value <- tt$p.value
    odf$fdr <- stats::p.adjust(odf$p_value, method = "BH")
    cd1_median <- median(cd1_df[[value_col]])
    cd2_median <- median(cd2_df[[value_col]])
    
    odf$ref_median <- cd1_median
    odf$condition_median <- cd2_median
    odf
  }
  
  out <- do.call('rbind', lapply(overlaps_list, test_wrapper))
  rownames(out) <- NULL
  out
}


#' @title Bar plot of median overlap scores
#' @description
#' Plot horizontal bars for median overlap scores in a reference group
#' and a condition group for each feature. Features are ordered by the
#' difference between condition and reference medians.
#'
#' @param df Data.frame with feature, overlap, and group columns.
#' @param condition_group Character. Name of the condition group.
#' @param feature_col Character. Column name for the feature identifier.
#' @param overlap_col Character. Column name for the overlap score.
#' @param group_col Character. Column name for the group variable.
#' @param n_table Integer or NULL. If not NULL, keep only the bottom and top
#' n_table features by difference in medians.
#' @param reference_label Character. Legend label for the reference bar.
#' @param title Character or NULL. Plot title.
#' @param width Numeric. Width in inches for saved output.
#' @param height Numeric. Height in inches for saved output.
#' @param save_path Character or NULL. If not NULL, path to save a TIFF file.
#' @param dpi Integer. Resolution in dots per inch for saved output.
#' @param transparent Logical. If TRUE, use a transparent background.
#' @param ref_color Character. Color for the reference bars.
#' @param cond_color Character. Color for the condition bars.
#'
#' @return A ggplot object.
#'
#' @export
plot_overlap_scores_bar <- function(
    df,
    condition_group,
    feature_col = "interaction",
    overlap_col = "overlapScore",
    group_col   = "group",
    n_table = NULL,
    reference_label = "reference",
    title = NULL,
    width = 7, height = 10,
    save_path = NULL,
    dpi = 300,
    transparent = FALSE,
    ref_color  = "#1F77B4",
    cond_color = "#D62728"
) {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Please install ggplot2.")
  if (missing(condition_group) || is.null(condition_group) || !nzchar(condition_group))
    stop("`condition_group` must be provided to define the condition.")
  
  ref_levels <- setdiff(unique(df[[group_col]]), condition_group)
  if (length(ref_levels) < 1)
    stop("Could not infer a reference group distinct from `condition_group`.")
  
  by_feat <- split(df, df[[feature_col]])
  agg <- lapply(by_feat, function(x) {
    cond_vals <- x[x[[group_col]] == condition_group, overlap_col, drop = TRUE]
    ref_vals  <- x[x[[group_col]] %in% ref_levels, overlap_col, drop = TRUE]
    med_cond <- if (length(cond_vals)) median(cond_vals, na.rm = TRUE) else NA_real_
    med_ref  <- if (length(ref_vals))  median(ref_vals,  na.rm = TRUE) else NA_real_
    diff_val <- med_cond - med_ref
    data.frame(
      feature     = unique(x[[feature_col]])[1],
      median_ref  = med_ref,
      median_cond = med_cond,
      diff        = diff_val,
      stringsAsFactors = FALSE
    )
  })
  agg <- do.call(rbind, agg)
  names(agg)[names(agg) == "feature"] <- feature_col
  
  agg <- agg[is.finite(agg$median_ref) & is.finite(agg$median_cond), , drop = FALSE]
  if (nrow(agg) == 0) stop("No features with finite medians to plot.")
  
  agg <- agg[order(agg$diff), , drop = FALSE]
  
  if (!is.null(n_table) && is.numeric(n_table) && n_table > 0) {
    n_table <- as.integer(min(n_table, floor(nrow(agg) / 2)))
    if (n_table > 0) {
      keep_idx <- unique(c(seq_len(n_table), seq.int(from = nrow(agg) - n_table + 1, to = nrow(agg))))
      agg <- agg[keep_idx, , drop = FALSE]
    }
  }
  
  agg[[feature_col]] <- factor(agg[[feature_col]], levels = agg[[feature_col]])
  
  plot_df <- rbind(
    data.frame(feature = agg[[feature_col]], group = "reference", value = agg$median_ref, stringsAsFactors = FALSE),
    data.frame(feature = agg[[feature_col]], group = condition_group, value = agg$median_cond, stringsAsFactors = FALSE)
  )
  names(plot_df)[names(plot_df) == "feature"] <- feature_col
  plot_df$group <- factor(plot_df$group, levels = c("reference", condition_group))
  
  cols <- c("reference" = ref_color, setNames(cond_color, condition_group))
  
  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = .data[["value"]], y = .data[[feature_col]], fill = .data[["group"]])
  ) +
    ggplot2::geom_col(position = "dodge", width = 0.75) +
    ggplot2::scale_fill_manual(
      values = cols,
      breaks = c("reference", condition_group),
      labels = c(reference_label, condition_group),
      name = NULL
    ) +
    ggplot2::labs(
      x = "Median Overlap Score",
      y = "Feature",
      title = if (is.null(title)) "" else title
    ) +
    ggplot2::theme_bw(base_size = 14) +
    ggplot2::theme(
      legend.position = "bottom",
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()
    )
  
  if (!is.null(save_path)) {
    dir.create(dirname(save_path), recursive = TRUE, showWarnings = FALSE)
    ggplot2::ggsave(
      filename = save_path, plot = p, width = width, height = height,
      dpi = dpi, device = "tiff",
      bg = ifelse(transparent, "transparent", "white")
    )
  }
  
  return(p)
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
      patHotspots  = spHotspotsList[[s]],
      infHotspots  = spHotspotsInfluenceList[[s]],
      patternPairs = patternPairs
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


#' @title Assemble ligand receptor interaction table
#' @description
#' Combine LR scores, ligand-only scores, and receptor-only scores into a single
#' long-format data.frame for downstream plotting. Works with either single
#' matrices or named lists of matrices (one per sample).
#'
#' @param lr_scores Matrix or named list of matrices with LR scores. Rows are
#' interaction IDs and columns are Source_to_Target combinations.
#' @param ligand_scores Matrix or named list of matrices with ligand-only
#' scores. Rows are interaction IDs and columns are Source_near_Target labels.
#' @param receptor_scores Matrix or named list of matrices with receptor-only
#' scores. Rows are interaction IDs and columns are target cell types.
#' @param lrpairs Data.frame with row names equal to interaction IDs and
#' annotation columns for ligand and receptor genes.
#' @param ligand_col Character. Column name in lrpairs containing ligand genes.
#' @param receptor_col Character. Column name in lrpairs containing receptor genes.
#' @param name_split_token Character. Token used to split Source_to_Target
#' column names in lr_scores.
#' @param ligand_near_token Character. Token used to construct Source_near_Target
#' column names in ligand_scores.
#' @param na_replace Numeric. Value used to replace NA scores.
#'
#' @return A data.frame where each row is one interaction for one
#' Source_to_Target column and one sample. Columns include ligand, ligand_score,
#' receptor, receptor_score, interaction, score, source_cell_type,
#' target_cell_type, source_to_target, and sample.
#'
#' @examples
#' \dontrun{
#' # Toy example with a single "sample"
#' set.seed(1)
#'
#' # interaction IDs
#' ints <- c("INT1","INT2")
#'
#' # LR scores: rows = interactions, cols = "Source_to_Target"
#' lr_mat <- matrix(runif(4), nrow = 2,
#'                  dimnames = list(ints, c("Epi_to_Tcell","Fibro_to_Tcell")))
#'
#' # ligand scores: rows = interactions, cols = "Source_near_Target"
#' lig_mat <- matrix(runif(4), nrow = 2,
#'                   dimnames = list(ints, c("Epi_near_Tcell","Fibro_near_Tcell")))
#'
#' # receptor scores: rows = interactions, cols = target cell types
#' rec_mat <- matrix(runif(2), nrow = 2,
#'                   dimnames = list(ints, "Tcell"))
#'
#' # ligand/receptor annotation
#' lrpairs <- data.frame(
#'   ligand   = c("LIG1","LIG2"),
#'   receptor = c("REC1","REC2"),
#'   row.names = ints
#' )
#'
#' df_long <- get_all_lr_interactions(
#'   lr_scores       = lr_mat,
#'   ligand_scores   = lig_mat,
#'   receptor_scores = rec_mat,
#'   lrpairs         = lrpairs
#' )
#' head(df_long)
#' }
#'
#' @importFrom dplyr bind_rows arrange
#' @export
get_all_lr_interactions <- function(
    lr_scores,                 # matrix OR *named list* of matrices (rows=interaction IDs, cols="Source_to_Target")
    ligand_scores,             # matrix OR *named list* (rows=interaction IDs, cols like "Source_near_Target")
    receptor_scores,           # matrix OR *named list* (rows=interaction IDs, cols are target cell types)
    lrpairs,                   # data.frame with rownames = interaction IDs; has ligand/receptor cols
    ligand_col       = "ligand",
    receptor_col     = "receptor",
    name_split_token = "_to_",
    ligand_near_token= "_near_",
    na_replace       = 0
){
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Please install 'dplyr'.")
  
  # --- basic checks ---
  if (!is.data.frame(lrpairs) || is.null(rownames(lrpairs)))
    stop("`lrpairs` must be a data.frame with rownames = interaction IDs.")
  if (!all(c(ligand_col, receptor_col) %in% colnames(lrpairs)))
    stop(sprintf("`lrpairs` must contain columns '%s' and '%s'.", ligand_col, receptor_col))
  
  .chk_same_type <- function(a, b, nmA, nmB){
    if (is.list(a) != is.list(b)) {
      stop(sprintf("`%s` and `%s` must both be matrices or both be *named* lists.", nmA, nmB))
    }
  }
  .chk_same_type(lr_scores, ligand_scores, "lr_scores", "ligand_scores")
  .chk_same_type(lr_scores, receptor_scores, "lr_scores", "receptor_scores")
  
  # --- coerce to unified iterator of (sample_name, lr_mat, lig_mat, rec_mat) ---
  .as_iter <- function(lr, lig, rec){
    if (is.list(lr)) {
      if (is.null(names(lr)) || any(names(lr) == "")) {
        stop("`lr_scores` list must be *named* (sample IDs).")
      }
      if (!all(names(lr) %in% names(lig))) {
        stop("`ligand_scores` missing samples found in `lr_scores`.")
      }
      if (!all(names(lr) %in% names(rec))) {
        stop("`receptor_scores` missing samples found in `lr_scores`.")
      }
      lapply(
        names(lr),
        function(s) list(sample = s, lr = lr[[s]], lig = lig[[s]], rec = rec[[s]])
      )
    } else {
      list(list(sample = "sample", lr = lr, lig = lig, rec = rec))
    }
  }
  iter <- .as_iter(lr_scores, ligand_scores, receptor_scores)
  
  # --- helpers ---
  .check_mat <- function(m, nm){
    if (!is.matrix(m) || is.null(rownames(m))) {
      stop(sprintf("`%s` must be a matrix with rownames=interaction IDs.", nm))
    }
  }
  
  # Build long DF
  out_list <- vector("list", length(iter))
  for (i in seq_along(iter)) {
    samp <- iter[[i]]$sample
    mLR  <- iter[[i]]$lr
    mL   <- iter[[i]]$lig
    mR   <- iter[[i]]$rec
    
    .check_mat(mLR, paste0("lr_scores[['", samp, "']]"))
    .check_mat(mL,  paste0("ligand_scores[['", samp, "']]"))
    .check_mat(mR,  paste0("receptor_scores[['", samp, "']]"))
    
    # keep only interactions present in lrpairs
    keep_rown <- intersect(rownames(mLR), rownames(lrpairs))
    if (!length(keep_rown)) next
    mLR <- mLR[keep_rown, , drop = FALSE]
    mL  <- mL [keep_rown, , drop = FALSE]
    mR  <- mR [keep_rown, , drop = FALSE]
    
    # long form for LR
    dfLR <- as.data.frame(as.table(mLR), stringsAsFactors = FALSE)
    colnames(dfLR) <- c("interaction","colname","score")
    dfLR$score <- suppressWarnings(as.numeric(dfLR$score))
    dfLR$score[is.na(dfLR$score)] <- na_replace
    
    # parse source/target from "Source_to_Target"
    parts <- strsplit(dfLR$colname, split = name_split_token, fixed = TRUE)
    dfLR$source_cell_type <- vapply(
      parts, function(x) if (length(x)) x[[1]] else NA_character_, character(1)
    )
    dfLR$target_cell_type <- vapply(
      parts, function(x) if (length(x) >= 2) x[[2]] else NA_character_, character(1)
    )
    dfLR$source_to_target <- paste0(dfLR$source_cell_type, name_split_token, dfLR$target_cell_type)
    dfLR$sample <- samp
    
    # locate corresponding ligand & receptor columns
    lig_cols <- paste0(dfLR$source_cell_type, ligand_near_token, dfLR$target_cell_type)
    rec_cols <- dfLR$target_cell_type
    
    # pull ligand/receptor scores by matching columns; use NA if col missing
    dfLR$ligand_score <- mapply(function(inter, lc){
      if (!is.null(lc) && lc %in% colnames(mL)) {
        as.numeric(mL[inter, lc])
      } else NA_real_
    }, dfLR$interaction, lig_cols)
    
    dfLR$receptor_score <- mapply(function(inter, rc){
      if (!is.null(rc) && rc %in% colnames(mR)) {
        as.numeric(mR[inter, rc])
      } else NA_real_
    }, dfLR$interaction, rec_cols)
    
    # attach ligand/receptor names from lrpairs
    lrsub <- lrpairs[dfLR$interaction, , drop = FALSE]
    dfLR$ligand   <- gsub(", ", "|", lrsub[[ligand_col]],   fixed = TRUE)
    dfLR$receptor <- gsub(", ", "|", lrsub[[receptor_col]], fixed = TRUE)
    
    out_list[[i]] <- dfLR[, c(
      "ligand","ligand_score",
      "receptor","receptor_score",
      "interaction","score",
      "source_cell_type","target_cell_type","source_to_target",
      "sample"
    ), drop = FALSE]
  }
  
  out <- dplyr::bind_rows(out_list)
  out$ligand_score[is.na(out$ligand_score)]     <- na_replace
  out$receptor_score[is.na(out$receptor_score)] <- na_replace
  
  out <- out |>
    dplyr::arrange(.data$source_to_target, .data$sample, dplyr::desc(.data$score), .data$interaction)
  
  as.data.frame(out)
}


#' @title Plot ligand receptor alluvial diagram
#' @description
#' Make an alluvial (Sankey style) plot showing flows from source cell type
#' to ligand, receptor, and target cell type. Flow width is proportional to an
#' interaction score. Flows can be colored by sample group or by the
#' source to target pair.
#'
#' @param df Data.frame with score, source, target, ligand, receptor, and
#' optionally sample columns.
#' @param sample_groups Optional named vector that maps sample IDs to group
#' labels. Names must match the sample column values.
#' @param col_score Character. Column name for the interaction score.
#' @param col_source Character. Column name for the source cell type.
#' @param col_target Character. Column name for the target cell type.
#' @param col_ligand Character. Column name for the ligand.
#' @param col_receptor Character. Column name for the receptor.
#' @param col_sample Character. Column name for the sample ID.
#' @param title Optional plot title. If NULL, a default title is used.
#' @param min_score Numeric. Minimum score to keep.
#' @param top_k_per_pair Optional integer. If not NULL, keep at most this many
#' top rows per source and target pair.
#' @param alpha Numeric between 0 and 1. Alpha transparency for flows.
#' @param knot_pos Numeric between 0 and 1. Knot position along the flow.
#' @param label_size Numeric. Text size for node labels.
#' @param wrap_width Integer or NULL. If not NULL, wrap long labels to this width.
#'
#' @return A ggplot object.
#'
#' @examples
#' \dontrun{
#' if (requireNamespace("ggplot2", quietly = TRUE) &&
#'     requireNamespace("ggalluvial", quietly = TRUE) &&
#'     requireNamespace("RColorBrewer", quietly = TRUE) &&
#'     requireNamespace("circlize", quietly = TRUE)) {
#'
#'   set.seed(123)
#'   top_normal <- data.frame(
#'     score            = runif(20, 0.1, 1),
#'     source_cell_type = sample(c("Keratinocyte", "Fibroblast"), 20, TRUE),
#'     target_cell_type = sample(c("T cell", "Macrophage"), 20, TRUE),
#'     ligand           = sample(c("CXCL12", "CCL5", "TNF"), 20, TRUE),
#'     receptor         = sample(c("CXCR4", "CCR5", "TNFRSF1A"), 20, TRUE),
#'     sample           = sample(c("NS_1", "NS_2", "VU_1", "VU_2"), 20, TRUE),
#'     stringsAsFactors = FALSE
#'   )
#'
#'   # Sample to group mapping
#'   sample_groups <- c(
#'     NS_1 = "Normal skin",
#'     NS_2 = "Normal skin",
#'     VU_1 = "Venous ulcer",
#'     VU_2 = "Venous ulcer"
#'   )
#'
#'   # Basic usage with grouping and custom title
#'   p <- plot_lr_alluvial(
#'     df             = top_normal,
#'     sample_groups  = sample_groups,
#'     top_k_per_pair = NULL,
#'     title = "Top Normal LR Interactions"
#'   )
#'   print(p)
#' }
#' }
#'
#' @import ggplot2
#' @import ggalluvial
#' @import RColorBrewer
#' @import circlize
#' @export
plot_lr_alluvial <- function(
    df,
    sample_groups = NULL,      # named vector: names = sample IDs, values = group labels (optional)
    col_score       = "score",
    col_source      = "source_cell_type",
    col_target      = "target_cell_type",
    col_ligand      = "ligand",
    col_receptor    = "receptor",
    col_sample      = "sample",
    title = NULL,
    min_score       = 0,       # drop tiny flows
    top_k_per_pair  = NULL,    # optional: keep top-K per (source,target) by score
    alpha           = 0.8,
    knot_pos        = 0.4,
    label_size      = 2.7,     # slightly smaller for dense plots
    wrap_width      = 12       # set NULL to disable wrapping
) {
  for (pkg in c("ggplot2","ggalluvial","RColorBrewer","circlize")) {
    if (!requireNamespace(pkg, quietly = TRUE)) stop(sprintf("Please install '%s'.", pkg))
  }
  
  # resolve columns (case tolerant)
  .res <- function(nm) { hit <- which(tolower(names(df)) == tolower(nm)); if (length(hit)) names(df)[hit[1]] else nm }
  col_score    <- .res(col_score);  col_source <- .res(col_source); col_target <- .res(col_target)
  col_ligand   <- .res(col_ligand); col_receptor <- .res(col_receptor); col_sample <- .res(col_sample)
  
  need <- c(col_score, col_source, col_target, col_ligand, col_receptor)
  miss <- setdiff(need, names(df)); if (length(miss)) stop("Missing columns: ", paste(miss, collapse = ", "))
  
  d <- df
  d[[col_score]] <- as.numeric(d[[col_score]])
  d <- d[is.finite(d[[col_score]]) & d[[col_score]] >= min_score, , drop = FALSE]
  if (!nrow(d)) stop("No rows to plot after filtering.")
  
  if (!is.null(top_k_per_pair) && top_k_per_pair > 0) {
    ord <- order(d[[col_source]], d[[col_target]], -d[[col_score]]); d <- d[ord, ]
    idx <- ave(d[[col_score]], paste(d[[col_source]], d[[col_target]], sep = "->"),
               FUN = function(x) seq_along(x) <= top_k_per_pair)
    d <- d[idx, , drop = FALSE]
  }
  
  # map fill: groups if provided, otherwise Source->Target interaction
  if (!is.null(sample_groups)) {
    if (!(col_sample %in% names(d))) stop("`sample_groups` provided but no 'sample' column in data.")
    d$.__group__ <- unname(sample_groups[ as.character(d[[col_sample]]) ])
    d <- d[!is.na(d$.__group__), , drop = FALSE]
    if (!nrow(d)) stop("No rows left after mapping samples to `sample_groups`.")
    d$fill_key <- factor(d$.__group__)
    fill_name  <- "Group"
  } else {
    d$fill_key <- factor(paste(d[[col_source]], d[[col_target]], sep = "->"))
    fill_name  <- "Cell interaction"
  }
  
  # wrapping only (no abbreviation)
  wrap_names <- function(v) if (is.null(wrap_width)) v else vapply(v, function(s) paste(strwrap(s, width = wrap_width), collapse = "\n"), "")
  d[[col_source]]   <- wrap_names(d[[col_source]])
  d[[col_ligand]]   <- wrap_names(d[[col_ligand]])
  d[[col_receptor]] <- wrap_names(d[[col_receptor]])
  d[[col_target]]   <- wrap_names(d[[col_target]])
  
  # palette: RColorBrewer first, extend with circlize for many levels
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
    source   = d[[col_source]],
    ligand   = d[[col_ligand]],
    receptor = d[[col_receptor]],
    target   = d[[col_target]],
    score    = d[[col_score]],
    fill_key = d$fill_key,
    stringsAsFactors = FALSE
  )
  
  p <- ggplot2::ggplot(
    dp,
    ggplot2::aes(axis1 = source, axis2 = ligand, axis3 = receptor, axis4 = target, y = score)
  ) +
    ggplot2::scale_x_discrete(limits = c("Source", "Ligand", "Receptor", "Target"),
                              expand = c(.08, .05)) +
    ggalluvial::geom_alluvium(
      ggplot2::aes(fill = fill_key),
      alpha = alpha, knot.pos = knot_pos, show.legend = TRUE
    ) +
    ggalluvial::geom_stratum(width = 0.25, color = "grey30", fill = "grey85") +
    ggalluvial::stat_stratum(
      geom = "text",
      ggplot2::aes(label = after_stat(stratum)),
      size = label_size,
      check_overlap = TRUE
    ) +
    ggplot2::scale_fill_manual(values = pal_vals, name = fill_name, drop = FALSE) +
    ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(alpha = 1))) +
    ggplot2::labs(
      title = "LR Alluvial: Source -> Ligand -> Receptor -> Target",
      subtitle = if (is.null(sample_groups))
        "Flow width ~ score; color ~ Cell interaction (Source->Target)"
      else
        "Flow width ~ score; color ~ Group",
      y = "Score", x = NULL
    ) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      legend.position = "right",
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(size = 11),
      plot.margin = ggplot2::margin(10, 30, 10, 30)
    )
  
  default_title <- "LR Alluvial: Source -> Ligand -> Receptor -> Target"
  p <- p +
    ggplot2::labs(
      title    = if (is.null(title)) default_title else title,
      y = "Score",
      x = NULL
    )
  return(p)
}

#' @title Plot ligand receptor scores by interaction
#' @description
#' Make scatter plots of ligand score versus receptor score, colored by
#' interaction and shaped by sample, and save one TIFF file per
#' source_to_target combination.
#'
#' @param df_long Long-format data.frame, usually from get_all_lr_interactions,
#' containing ligand_score, receptor_score, interaction, sample, and
#' source_to_target columns.
#' @param out_dir Character. Output directory for TIFF files. Created if needed.
#' @param filename_prefix Character. Prefix for output file names.
#' @param source_to_target Character vector of source_to_target values to plot.
#' If NULL, all unique values are used.
#' @param point_size Numeric. Point size for geom_point.
#' @param point_alpha Numeric between 0 and 1. Alpha transparency for points.
#' @param width Numeric. Plot width in inches.
#' @param height Numeric. Plot height in inches.
#' @param dpi Numeric. Resolution in dots per inch for the saved TIFF files.
#'
#' @return Invisibly, a character vector of paths to the saved TIFF files.
#'
#' @examples
#' \dontrun{
#' if (requireNamespace("ggplot2", quietly = TRUE) &&
#'     requireNamespace("RColorBrewer", quietly = TRUE) &&
#'     requireNamespace("circlize", quietly = TRUE)) {
#'
#'   set.seed(123)
#'   df_long <- data.frame(
#'     ligand_score    = runif(50, 0, 1),
#'     receptor_score  = runif(50, 0, 1),
#'     interaction     = sample(c("LGF1_REC1", "LGF2_REC2", "LGF3_REC3"), 50, TRUE),
#'     sample          = sample(c("NS_1", "NS_2", "VU_1"), 50, TRUE),
#'     source_to_target = sample(c("Epi_to_Tcell", "Fibro_to_Myeloid"), 50, TRUE),
#'     stringsAsFactors = FALSE
#'   )
#'
#'   out_dir <- tempdir()
#'
#'   # Plot all source_to_target combinations
#'   files_all <- plot_lr_scores_scatter(
#'     df_long        = df_long,
#'     out_dir        = out_dir,
#'     filename_prefix = "lr_scatter",
#'     source_to_target = NULL,
#'     point_size      = 2,
#'     point_alpha     = 0.8
#'   )
#'
#'   # Plot only one specific source_to_target
#'   files_one <- plot_lr_scores_scatter(
#'     df_long        = df_long,
#'     out_dir        = out_dir,
#'     filename_prefix = "Epi_Tcell_lr_scatter",
#'     source_to_target = "Epi_to_Tcell"
#'   )
#' }
#' }
#'
#' @import ggplot2
#' @import RColorBrewer
#' @import circlize
#' @export
plot_lr_scores_scatter <- function(
    df_long,                       # from get_all_lr_interactions()
    out_dir,
    filename_prefix = "lr_scatter",
    source_to_target = NULL,       # if set, filter to this "Source_to_Target" and plot just one
    point_size      = 2,
    point_alpha     = 0.8,
    width           = 8,           # inches
    height          = 6,           # inches
    dpi             = 300
){
  if (!requireNamespace("ggplot2",      quietly = TRUE)) stop("Please install 'ggplot2'.")
  if (!requireNamespace("RColorBrewer", quietly = TRUE)) stop("Please install 'RColorBrewer'.")
  if (!requireNamespace("circlize",     quietly = TRUE)) stop("Please install 'circlize'.")
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # choose which source_to_target values to plot
  st_vals <- if (is.null(source_to_target)) unique(df_long$source_to_target) else source_to_target
  st_vals <- st_vals[st_vals %in% unique(df_long$source_to_target)]
  if (!length(st_vals)) stop("No matching 'source_to_target' in data.")
  
  # discrete palette using RColorBrewer + circlize (no 'scales' dependency)
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
  
  saved <- character(0)
  for (st in st_vals) {
    dd <- df_long[df_long$source_to_target == st, , drop = FALSE]
    if (!nrow(dd)) next
    
    # color by interaction, shape by sample
    n_inter <- length(unique(dd$interaction))
    cols <- .discrete_palette(n_inter)
    names(cols) <- unique(dd$interaction)
    
    p <- ggplot2::ggplot(
      dd,
      ggplot2::aes(x = ligand_score, y = receptor_score,
                   color = interaction, shape = sample)
    ) +
      ggplot2::geom_point(size = point_size, alpha = point_alpha) +
      ggplot2::scale_color_manual(values = cols, name = "Interaction") +
      ggplot2::labs(
        title = st,
        x = "Ligand score",
        y = "Receptor score"
      ) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(
        legend.position = "right",
        panel.grid.minor = ggplot2::element_blank()
      )
    
    out_file <- file.path(
      out_dir,
      paste0(
        filename_prefix, "_",
        gsub("[^A-Za-z0-9._-]+","_", st),
        ".tiff"
      )
    )
    ggplot2::ggsave(
      filename = out_file, plot = p, device = "tiff",
      width = width, height = height, dpi = dpi,
      units = "in", compress = "lzw"
    )
    saved <- c(saved, out_file)
  }
  
  invisible(saved)
}



