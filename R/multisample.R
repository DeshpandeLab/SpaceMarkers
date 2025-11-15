#' Align counts and coordinates by shared barcodes (same set, same order)
#'
#' @description
#' Ensures that a sparse counts matrix and a coordinates data.frame
#' contain the **same barcodes in the same order**. Rows of `coords_df`
#' are re-ordered to match the column order of `counts`. Returns both
#' aligned objects.
#'
#' @param counts A sparse gene-by-spot matrix (typically a \code{dgCMatrix})
#'   whose \code{colnames} are spot barcodes.
#' @param coords_df A data.frame of spot coordinates with a \code{barcode}
#'   column containing the spot barcodes.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{counts}: the input \code{counts} restricted to shared barcodes and
#'   with columns ordered to match \code{coords_df} after alignment
#'   \item \code{coords}: the input \code{coords_df} restricted and re-ordered so
#'   that \code{rownames(coords)} equals \code{colnames(counts)}
#' }
#'
#' @details
#' If there are no overlapping barcodes between the two inputs, the function
#' throws an error.
#'
#' @examples
#' \dontrun{
#' aligned <- align_counts_coords(counts, coords_df)
#' counts_al <- aligned$counts
#' coords_al <- aligned$coords
#' stopifnot(identical(colnames(counts_al), rownames(coords_al)))
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


#' Process multiple Visium/VisiumHD samples into filtered counts/coords RDS files
#'
#' @description
#' Loads counts and coordinates for each sample, applies a per-sample sparse
#' gene filter (keep genes detected in at least `good_gene_threshold` spots),
#' aligns barcodes between counts and coords via \code{align_counts_coords()},
#' saves per-sample `.rds` files, and writes a combined
#' `all_samples_counts_coords.rds`.
#'
#' @param data_dir Character. Root directory containing the sample folders.
#' @param samples Character vector. Sample folder names under `data_dir`.
#' @param counts_file Character. H5 filename in each sample's bin directory.
#'   Default: \code{"filtered_feature_bc_matrix.h5"}.
#' @param good_gene_threshold Integer (>=1). Keep genes detected (>0) in at least
#'   this many spots (per sample). Default: \code{10}.
#' @param bin_subdir Character. Relative path under each sample directory where
#'   the 10X files live (e.g., VisiumHD bins). Default:
#'   \code{"binned_outputs/square_016um"}.
#' @param out_dir Character or \code{NULL}. Directory to save outputs. If \code{NULL},
#'   uses \code{file.path(data_dir, "outputs")}. Created if it doesn't exist.
#' @param coords_resolution Character. Resolution hint passed to
#'   \code{load10XCoords()} (e.g., \code{"fullres"}). Default: \code{"fullres"}.
#' @param verbose Logical. If \code{TRUE}, prints progress. Default: \code{TRUE}.
#'
#' @details
#' Assumes \code{load10XExpr()} returns a \code{dgCMatrix} (log1p-transformed per your
#' pipeline) and \code{load10XCoords()} returns a data.frame with a \code{barcode}
#' column. Nonzero-per-gene is computed from the \code{@i} slot of the CSC matrix.
#'
#' @return A named list with:
#' \itemize{
#'   \item \code{counts}: list of filtered \code{dgCMatrix} objects per sample
#'   \item \code{coords}: list of aligned coordinate data.frames per sample
#'   \item \code{paths}: list with per-sample \code{counts} and \code{coords} RDS paths and \code{all} combined path
#' }
#'
#' @examples
#' \dontrun{
#' samples <- c(
#'   "visiumHD_crc_p1",
#'   "visiumHD_crc_p2",
#'   "visiumHD_crc_p5",
#'   "visiumHD_normal_colon_p3",
#'   "visiumHD_normal_colon_p5"
#' )
#' res <- process_visium_samples(
#'   data_dir            = "spacemarkers/data",
#'   samples             = samples,
#'   counts_file         = "filtered_feature_bc_matrix.h5",
#'   good_gene_threshold = 10,
#'   bin_subdir          = "binned_outputs/square_016um",
#'   out_dir             = NULL,
#'   coords_resolution   = "fullres",
#'   verbose             = TRUE
#' )
#' res$paths$all
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


#' Build a per-spot spatial features table (coords + features) for one sample
#'
#' @description
#' Robust helper that accepts coordinates as a data.frame (or matrix) and features
#' as a data.frame, base matrix, or \code{Matrix} sparse matrix. It:
#' \itemize{
#'   \item ensures a \code{barcode} column exists in \code{spCoords} (deduplicated),
#'   \item extracts/standardizes the feature matrix from \code{spFeatures},
#'   \item intersects barcodes across counts/coords/features,
#'   \item aligns by barcode (no merges) and returns a data.frame:
#'         \code{[barcode, ...coord columns..., <feature columns>]}.
#' }
#'
#' @param spCounts A genes x barcodes matrix (dense or \code{Matrix}) used only
#'   to define the barcode set via \code{colnames(spCounts)}.
#' @param spCoords A coordinates table (data.frame or matrix). Must contain
#'   numeric columns \code{x} and \code{y}. If it lacks a \code{barcode} column,
#'   rownames are used to construct one.
#' @param spFeatures A table of per-barcode features. Supported inputs:
#'   \itemize{
#'     \item matrix/\code{Matrix}: rownames are barcodes; columns are features
#'     \item data.frame with \code{barcode}: that column is used as row ID and
#'           removed from the feature matrix
#'     \item data.frame without \code{barcode}: rownames are treated as barcodes
#'   }
#'
#' @return A data.frame with barcodes as rownames containing:
#'   \code{barcode}, \code{x}, \code{y}, and all feature columns.
#'
#' @examples
#' \dontrun{
#' # Single sample
#' pat <- make_sp_patterns(spCounts, spCoords, spFeatures)
#' head(pat)
#'
#' # Multiple samples with parallel lists
#' common <- Reduce(intersect, list(names(spCountsList), names(spCoordsList), names(spFeaturesList)))
#' spPatternsList <- Map(
#'   make_sp_patterns,
#'   spCountsList [common],
#'   spCoordsList [common],
#'   spFeaturesList[common]
#' )
#' }
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

#' @title get_multi_sample_overlaps
#' @description Compute overlap scores per sample and assemble long + wide summaries
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


#' @title plot_sample_overlaps
#' @description pheatmap of multi-sample overlap scores (rows = interactions, cols = samples)
#'              after removing "similar" interactions by default (e.g., keep Epi_Epi over Epi_Epi_in).
#' @param wide_df data.frame with first column 'interaction' and remaining columns = datasets
#' @param scale one of "none","row","column" (passed to pheatmap)
#' @param cluster_rows,cluster_cols logical; cluster rows/cols
#' @param color optional color palette; if NULL, a sensible palette is chosen
#' @param na_color tile color for NAs
#' @param display_numbers logical; print values inside tiles
#' @param number_format sprintf format for numbers (e.g. "%.2f")
#' @param fontsize base font size; fontsize_row, fontsize_col for axes
#' @param main title
#' @param filename optional output path (pheatmap saves if provided)
#' @param ... passed through to pheatmap::pheatmap
#' @return Filtered wide data.frame (what was plotted)
#' @importFrom pheatmap pheatmap
#' @export
#' @title plotSampleOverlaps
#' @description pheatmap of multi-sample overlap scores (rows = interactions, cols = samples)
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


#' Compare scores between a reference condition and a comparator across features
#'
#' @description
#' For each \code{feature_col} (e.g., ligand–receptor interaction), split the
#' data and perform an unpaired two-sided Welch t-test on \code{value_col}
#' between the reference condition (\code{ref_condition}) and the comparator
#' condition (taken as the first non-reference level observed in \code{group_col}).
#' Returns per-row annotations including raw \code{p_value}, BH-adjusted \code{fdr},
#' median-based \code{log2fc}, and the direction (\code{fc_direction}).
#'
#' @param df A data.frame containing at least \code{group_col}, \code{feature_col},
#'   and \code{value_col}.
#' @param ref_condition Character scalar. Name of the reference condition
#'   (baseline) in \code{group_col}. Default \code{"Normal"}.
#' @param group_col Character scalar. Column with two conditions to compare.
#'   Default \code{"group"}.
#' @param feature_col Character scalar. Column that defines the feature split
#'   (e.g., interaction ID). Default \code{"interaction"}.
#' @param value_col Character scalar. Column with numeric scores to compare
#'   (e.g., overlap score). Default \code{"overlapScore"}.
#'
#' @details
#' The comparator condition is chosen as the first non-reference level present
#' in \code{df[[group_col]]}. If more than two levels are present, only the
#' reference and the first other level are used. Fold-change is computed from
#' medians as \code{log2(median(comparator) / (median(reference) + 1e-5))}.
#'
#' @return A data.frame matching the input rows (stacked across features) with
#'   added columns: \code{p_value}, \code{fdr}, \code{log2fc}, \code{fc_direction}.
#'
#' @examples
#' \dontrun{
#' res <- compare_scores(
#'   df,
#'   ref_condition = "Normal",
#'   group_col     = "group",
#'   feature_col   = "interaction",
#'   value_col     = "overlapScore"
#' )
#' }
#'
#' @importFrom dplyr filter pull
#' @importFrom rlang sym
#' @importFrom magrittr %>%
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


#' Bar plot of median overlap scores for condition and reference per feature
#'
#' @description
#' For each feature, computes the median \code{overlap_col} for the
#' \emph{condition} (\code{condition_group}) and for the \emph{reference}
#' (all other levels of \code{group_col}), then plots both medians as
#' side-by-side horizontal bars. Features are ordered by the signed
#' difference \code{median(condition) - median(reference)} (ascending).
#'
#' @param df Data.frame with columns \code{feature_col}, \code{overlap_col}, and \code{group_col}.
#' @param condition_group Character. Name of the condition group (the "condition" bar).
#' @param feature_col Character. Feature column (y-axis). Default \code{"interaction"}.
#' @param overlap_col Character. Column with overlap scores. Default \code{"overlapScore"}.
#' @param group_col Character. Group column. Default \code{"group"}.
#' @param n_table Integer or \code{NULL}. If provided, keeps only the bottom
#'   \code{n_table} and top \code{n_table} features by signed difference
#'   (condition − reference). Default \code{NULL}.
#' @param reference_label Character. How to display the reference group in the legend.
#'   Default \code{"reference"}.
#' @param title Character or \code{NULL}. Plot title; default \code{NULL} (empty).
#' @param width,height Numeric. Size for saved output in inches. Defaults \code{7}, \code{10}.
#' @param save_path Character or \code{NULL}. If provided, saves a TIFF to this path.
#' @param dpi Integer. Resolution for saved figure. Default \code{300}.
#' @param transparent Logical. Saved background transparency. Default \code{FALSE}.
#' @param ref_color,cond_color Character. Hex colors for reference and condition bars.
#'
#' @return A \code{ggplot} object.
#'
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



#' Run SpaceMarkers IMscores (HD) across samples with global feature alignment
#'
#' @description
#' Computes directed IMscores for each sample using
#' \code{SpaceMarkers::calculate_gene_scores_directed()} after:
#' \enumerate{
#'   \item aligning sample names across the four input lists, and
#'   \item enforcing a \strong{global feature intersection} across all samples
#'         (from \code{spPatternsList} columns excluding \code{barcode}, \code{x}, \code{y}).
#' }
#' Runtime is measured per sample (via \code{base::system.time}) and returned as
#' \code{elapsed_sec}. Optionally writes per-sample IMscores and the elapsed
#' vector to disk.
#'
#' @param spCountsList Named list of genes x barcodes matrices (typically \code{dgCMatrix}).
#' @param spPatternsList Named list of data.frames with columns \code{barcode}, \code{x}, \code{y},
#'   and \emph{feature columns} (must yield \eqn{\ge} 2 shared features globally).
#' @param spHotspotsList Named list of pattern hotspot objects (from \code{findAllHotspots.value()}).
#' @param spHotspotsInfluenceList Named list of influence hotspot objects.
#' @param sampleNames \code{NULL} or character vector of sample names to run. If \code{NULL},
#'   uses the intersection of names across all four lists; otherwise intersects with those names.
#' @param outDir \code{NULL} or character path. If provided, saves per-sample IMscores as RDS files
#'   and also saves an RDS with the named numeric vector \code{elapsed_sec}. The directory is
#'   created if it does not exist. Default \code{NULL}.
#' @param imscoresSuffix Filename suffix for IMscores RDS when saving. Default \code{"_IMscores_directed.rds"}.
#' @param elapsedFilename Filename (no sample prefix) for the elapsed vector RDS. Default
#'   \code{"IMscoresHD_elapsed_seconds.rds"}.
#' @param verbose Logical; print progress. Default \code{TRUE}.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{imscores}: named list of per-sample IMscores tables
#'   \item \code{elapsed_sec}: named numeric vector of elapsed seconds per sample
#'   \item \code{paths}: (only when \code{outDir} provided) a list with
#'         \code{imscores} (per-sample paths) and \code{elapsed_path}
#' }
#'
#' @examples
#' \dontrun{
#' res <- run_imscores_one(
#'   spCountsList            = spCounts_list,
#'   spPatternsList          = spPatterns_list,
#'   spHotspotsList          = spHotspots_list,
#'   spHotspotsInfluenceList = spHotspots_influence_list,
#'   sampleNames             = c("visiumHD_crc_p1","visiumHD_crc_p2"),
#'   outDir                  = file.path(data_dir, "outputs_spacemarkers")
#' )
#' res$elapsed_sec
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

#' Parallel LR-scoring over samples
#'
#' @param IMscores_directed_list list of per-sample matrices for ligand scoring
#' @param spCounts_list          list of per-sample count matrices
#' @param spPatterns_list        list of per-sample pattern matrices
#' @param ligands_list           named list of ligand gene sets
#' @param receptors_list         named list of receptor gene sets
#' @param lrpairs                data.frame of ligand-receptor pairs
#' @param workers                number of parallel workers (default 5)
#' @param weighted               logical, use weighted scoring (default TRUE)
#' @param ligand_method          aggregation method for ligands (default "arithmetic_mean")
#' @param receptor_method        aggregation method for receptors (default "arithmetic_mean")
#' @param lr_method              method to combine L/R (default "geometric_mean")
#' @return A list with three named lists (one element per sample):
#' \itemize{
#'   \item \code{$res_list}             — LR score matrices
#'   \item \code{$ligand_scores_list}   — ligand score matrices
#'   \item \code{$receptor_scores_list} — receptor score matrices
#' }
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


#' Assemble ligand–receptor interaction table across samples
#'
#' @description
#' Combines ligand–receptor (LR) scores, ligand-only scores, and receptor-only
#' scores into a single long-format data.frame suitable for downstream plotting
#' (e.g., scatter plots, alluvial diagrams). Supports both:
#' \itemize{
#'   \item a single matrix per object (\code{lr_scores}, \code{ligand_scores},
#'         \code{receptor_scores}), or
#'   \item parallel *named lists* of matrices (sample-wise), where each list
#'         element name is treated as a sample ID.
#' }
#'
#' Rows of the returned table correspond to individual ligand–receptor
#' \emph{interactions} for specific Source→Target pairs and samples.
#'
#' @param lr_scores A numeric matrix or a *named list* of matrices with
#'   LR scores. Each matrix must have:
#'   \itemize{
#'     \item \strong{rows} = interaction IDs (matching \code{rownames(lrpairs)})
#'     \item \strong{columns} = Source→Target combinations, e.g.
#'           \code{"Epi_to_Tcell"} (see \code{name_split_token}).
#'   }
#' @param ligand_scores A numeric matrix or a *named list* of matrices
#'   (same structure as \code{lr_scores}) containing ligand-only scores.
#'   Each matrix must have:
#'   \itemize{
#'     \item \strong{rows} = interaction IDs (matching \code{lr_scores} rows)
#'     \item \strong{columns} = Source\code{_near_}Target-style labels, e.g.
#'           \code{"Epi_near_Tcell"} by default (see \code{ligand_near_token}).
#'   }
#' @param receptor_scores A numeric matrix or a *named list* of matrices
#'   containing receptor-only scores. Each matrix must have:
#'   \itemize{
#'     \item \strong{rows} = interaction IDs (matching \code{lr_scores} rows)
#'     \item \strong{columns} = target cell types, e.g. \code{"Tcell"}.
#'   }
#' @param lrpairs A data.frame with rownames equal to interaction IDs
#'   (matching rownames of the matrices) and at least the ligand and receptor
#'   annotation columns specified by \code{ligand_col} and \code{receptor_col},
#'   e.g. \code{"ligand"} and \code{"receptor"}. These are used to attach
#'   gene names to each interaction.
#' @param ligand_col Character scalar. Name of the ligand annotation column
#'   in \code{lrpairs}. Default \code{"ligand"}.
#' @param receptor_col Character scalar. Name of the receptor annotation
#'   column in \code{lrpairs}. Default \code{"receptor"}.
#' @param name_split_token Character scalar used to split Source→Target
#'   column names in \code{lr_scores}. For example, if column names look like
#'   \code{"Epi_to_Tcell"}, then \code{name_split_token = "_to_"} (default).
#' @param ligand_near_token Character scalar used to construct the expected
#'   column names in \code{ligand_scores}. For a given source and target,
#'   the corresponding ligand column is assumed to be:
#'   \code{paste0(source, ligand_near_token, target)}, e.g.
#'   \code{"Epi_near_Tcell"} by default.
#' @param na_replace Numeric scalar. Value used to replace any \code{NA}
#'   scores in the assembled table (both LR and ligand/receptor scores).
#'   Default \code{0}.
#'
#' @details
#' This function:
#' \enumerate{
#'   \item Validates that \code{lr_scores}, \code{ligand_scores},
#'         and \code{receptor_scores} are either all matrices or all *named*
#'         lists of matrices.
#'   \item Aligns rows (interaction IDs) to those present in \code{lrpairs}.
#'   \item Converts \code{lr_scores} to long format; parses \code{source} and
#'         \code{target} cell types per column via \code{name_split_token}.
#'   \item Looks up the matching ligand and receptor columns in
#'         \code{ligand_scores} and \code{receptor_scores} using:
#'         \itemize{
#'           \item ligand column: \code{paste0(source, ligand_near_token, target)}
#'           \item receptor column: \code{target}
#'         }
#'   \item Attaches ligand and receptor gene names from \code{lrpairs}
#'         (columns \code{ligand_col}, \code{receptor_col}), with any
#'         comma+space sequences replaced by \code{"|"} for downstream parsing.
#' }
#'
#' When lists are provided, each list element is treated as a separate
#' sample, and a \code{sample} column is added to the output. When a single
#' matrix is provided, \code{sample} is set to \code{"sample"}.
#'
#' @return
#' A data.frame with one row per (sample, interaction ID, Source→Target column)
#' containing the following columns:
#' \itemize{
#'   \item \code{ligand}         — ligand gene(s) (from \code{lrpairs})
#'   \item \code{ligand_score}   — ligand-only score
#'   \item \code{receptor}       — receptor gene(s) (from \code{lrpairs})
#'   \item \code{receptor_score} — receptor-only score
#'   \item \code{interaction}    — interaction ID (rownames of \code{lr_scores})
#'   \item \code{score}          — LR score (from \code{lr_scores})
#'   \item \code{source_cell_type} — parsed source cell type
#'   \item \code{target_cell_type} — parsed target cell type
#'   \item \code{source_to_target} — Source→Target label, e.g. \code{"Epi_to_Tcell"}
#'   \item \code{sample}         — sample ID (list element name) or
#'                                  \code{"sample"} for single matrices
#' }
#'
#' Scores with missing values are replaced by \code{na_replace} before
#' returning. The rows are arranged by \code{source_to_target}, then
#' \code{sample}, then descending \code{score}, and finally by
#' \code{interaction}.
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


#' Plot Ligand–Receptor Alluvial Diagram
#'
#' Creates an alluvial (Sankey-style) plot showing flows from source cell type
#' to ligand, receptor, and target cell type, with flow width proportional to
#' an interaction score. Optionally, flows can be colored by sample groups
#' (e.g., conditions) or by the source–target interaction.
#'
#' @param df A data.frame containing at least the columns specified by
#'   \code{col_score}, \code{col_source}, \code{col_target},
#'   \code{col_ligand}, and \code{col_receptor}. Optionally a sample column
#'   (\code{col_sample}) when \code{sample_groups} is used.
#' @param sample_groups Optional named vector mapping sample IDs to group labels.
#'   Names must correspond to values in \code{df[[col_sample]]}. When provided,
#'   flows are colored by group; otherwise, flows are colored by
#'   Source→Target interaction.
#' @param col_score Character. Column name in \code{df} containing the numeric
#'   interaction scores used as flow widths. Default is \code{"score"}.
#' @param col_source Character. Column name for source cell type. Default
#'   \code{"source_cell_type"}.
#' @param col_target Character. Column name for target cell type. Default
#'   \code{"target_cell_type"}.
#' @param col_ligand Character. Column name for ligand. Default \code{"ligand"}.
#' @param col_receptor Character. Column name for receptor.
#'   Default \code{"receptor"}.
#' @param col_sample Character. Column name for sample ID (used only if
#'   \code{sample_groups} is not \code{NULL}). Default \code{"sample"}.
#' @param title Optional character string for the plot title. If \code{NULL},
#'   a default title ("LR Alluvial: Source → Ligand → Receptor → Target")
#'   is used.
#' @param min_score Numeric. Minimum score threshold; rows with scores below
#'   this value are dropped prior to plotting. Default is \code{0}.
#' @param top_k_per_pair Optional integer. If not \code{NULL} and > 0, keeps at
#'   most this many top-scoring rows per Source→Target pair (based on
#'   \code{col_score}). Default \code{NULL} (no additional filtering).
#' @param alpha Numeric (0–1). Alpha transparency for alluvial flows.
#'   Default is \code{0.8}.
#' @param knot_pos Numeric (0–1). Position of the knot (curve) along the flow
#'   in \code{geom_alluvium}. Default is \code{0.4}.
#' @param label_size Numeric. Text size for stratum (node) labels.
#'   Default is \code{2.7}.
#' @param wrap_width Integer or \code{NULL}. Width (in characters) used to wrap
#'   long labels (source, ligand, receptor, target) onto multiple lines.
#'   Set to \code{NULL} to disable wrapping. Default is \code{12}.
#'
#' @return A \code{ggplot} object representing the alluvial plot.
#'
#' @details
#' Column matching is case-tolerant: if, for example, \code{col_score = "score"},
#' the function will match any column whose name equals "score" ignoring case.
#' If \code{sample_groups} is provided, only rows whose sample IDs map to a
#' non-\code{NA} group are retained.
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
#'   # Sample → group mapping
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
#'     title          = "Top Normal Ligand–Receptor Interactions"
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
    idx <- ave(d[[col_score]], paste(d[[col_source]], d[[col_target]], sep = "→"),
               FUN = function(x) seq_along(x) <= top_k_per_pair)
    d <- d[idx, , drop = FALSE]
  }
  
  # map fill: groups if provided, otherwise Source→Target interaction
  if (!is.null(sample_groups)) {
    if (!(col_sample %in% names(d))) stop("`sample_groups` provided but no 'sample' column in data.")
    d$.__group__ <- unname(sample_groups[ as.character(d[[col_sample]]) ])
    d <- d[!is.na(d$.__group__), , drop = FALSE]
    if (!nrow(d)) stop("No rows left after mapping samples to `sample_groups`.")
    d$fill_key <- factor(d$.__group__)
    fill_name  <- "Group"
  } else {
    d$fill_key <- factor(paste(d[[col_source]], d[[col_target]], sep = "→"))
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
      title = "LR Alluvial: Source → Ligand → Receptor → Target",
      subtitle = if (is.null(sample_groups))
        "Flow width ~ score; color ~ Cell interaction (Source→Target)"
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
  
  default_title <- "LR Alluvial: Source → Ligand → Receptor → Target"
  p <- p +
    ggplot2::labs(
      title    = if (is.null(title)) default_title else title,
      y = "Score",
      x = NULL
    )
  return(p)
}

#' Plot Ligand–Receptor Scores Scatter Plots by Interaction
#'
#' Creates one or more scatter plots of ligand vs receptor scores, colored by
#' ligand–receptor interaction and shaped by sample. Each plot corresponds to a
#' single \code{source_to_target} pair and is saved as a TIFF file.
#'
#' @param df_long A data.frame in long format, typically produced by
#'   \code{get_all_lr_interactions()}, containing at least the columns:
#'   \itemize{
#'     \item \code{ligand_score} (numeric)
#'     \item \code{receptor_score} (numeric)
#'     \item \code{interaction} (factor or character; ligand–receptor ID)
#'     \item \code{sample} (factor or character; sample ID)
#'     \item \code{source_to_target} (factor or character; e.g. \code{"Epi_to_Tcell"})
#'   }
#' @param out_dir Character. Output directory where TIFF files will be saved.
#'   Created if it does not exist.
#' @param filename_prefix Character. Prefix for output filenames. Each file will
#'   be named as \code{<filename_prefix>_<source_to_target_lr_scatter>.tiff}
#'   after sanitizing \code{source_to_target}. Default is \code{"lr_scatter"}.
#' @param source_to_target Optional character vector of one or more
#'   \code{source_to_target} values to plot. If \code{NULL} (default), all
#'   unique \code{source_to_target} values in \code{df_long} are used.
#' @param point_size Numeric. Point size for \code{geom_point}. Default \code{2}.
#' @param point_alpha Numeric (0–1). Alpha transparency for points.
#'   Default \code{0.8}.
#' @param width Numeric. Width of the saved TIFF(s) in inches. Default \code{8}.
#' @param height Numeric. Height of the saved TIFF(s) in inches. Default \code{6}.
#' @param dpi Numeric. Resolution (dots per inch) for the saved TIFF(s).
#'   Default \code{300}.
#'
#' @return
#' Invisibly returns a character vector of file paths to the saved TIFFs.
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



