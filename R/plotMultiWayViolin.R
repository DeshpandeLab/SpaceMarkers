# Ensure required packages are installed (if not already)
# --- Helper: addInteractCol ---------------------------------------------------
#' Add an Interaction Label Column for a Pair of Hotspot Features
#'
#' Creates a factor label per row indicating whether neither, one, or both
#' hotspot features are present (non-`NA`). The new factor has levels in the
#' order: `<ref>`, `<interact_name>`, `<other>`.
#'
#' @param hotspots A `data.frame` with at least two columns that
#'   encode hotspot membership as `NA`/non-`NA`.
#' @param interaction_cols Character vector of length 2 giving the column names
#'  to compare (e.g., `c("refPattern", "OtherPattern")`). The first name is
#'  treated as the reference (`ref`), the second as the `other`.
#' @param interact_name Label to assign when both features are present
#'  in the same row. Default `"Interacting"`.
#' @param interaction_col_name Name of the new column to create.
#'
#' @return `hotspots` with one additional factor column `interaction_col_name`
#'   whose levels are `c(<ref>, <interact_name>, <other>)`.
#'
#' @details
#' A row is labeled:
#' - `<interact_name>` if both columns are non-`NA`
#' - `<ref>` if only the first column is non-`NA`
#' - `<other>` if only the second column is non-`NA`
#' - `NA` if both are `NA`
#'
#' A warning is emitted if all labels are `NA`.
#'
#' @examples
#' tmp <- data.frame(refPattern = c("a", NA, NA),
#'                   OtherPattern = c("b", "b", NA))
#' tmp2 <- addInteractCol(tmp,
#'   interaction_cols = c("refPattern", "OtherPattern"),
#'   interaction_col_name = "ref_vs_other"
#' )
#' print(tmp2)
#' @export
#' @seealso [addPairwiseInteractCol()], [plotMultiWayViolin()]
#' 
addInteractCol <- function(hotspots,
                           interaction_cols = c("refPattern", "OtherPattern"),
                           interact_name = "Interacting",
                           interaction_col_name = "interaction") {
  # Validate inputs
  miss <- setdiff(interaction_cols, colnames(hotspots))
  if (length(miss)) {
    stop(sprintf("Missing interaction columns in 'hotspots': %s",
                 paste(miss, collapse = ", ")))
  }
  
  # Compute region (formerly createInteractCol)
  refPattern <- hotspots[[interaction_cols[1]]]
  pattern    <- hotspots[[interaction_cols[2]]]
  
  ref   <- interaction_cols[1]
  other <- interaction_cols[2]
  
  region <- ifelse(!is.na(refPattern) & !is.na(pattern),
                   interact_name,
                   ifelse(!is.na(refPattern), ref,
                          ifelse(!is.na(pattern), other, NA)))
  names(region) <- rownames(hotspots)
  
  if (all(is.na(region))) {
    warning("No interaction found")
  }
  
  region <- factor(region, levels = c(ref, interact_name, other))
  
  
  hotspots[[interaction_col_name]] <- region
  return(hotspots)
}
# --- Helper: addPairwiseInteractCol ------------------------------------------
#' Add All Pairwise Interaction Columns Across Features
#'
#' For a set of hotspot feature columns, compute an interaction label
#' (see [addInteractCol()]) for every unique ordered pair `A_B`. Each resulting
#' column is a factor with levels `c("A", "Interacting", "B")`.
#'
#' @param hotspots A `data.frame` with hotspot feature columns and (optionally)
#'   an identifier column named `"barcode"`.
#' @param features Character vector of feature column names to pair. If `NULL`,
#'   all columns except `"barcode"` are used.
#' @param interact_name Label to assign when both features are present in the
#'   same row Default `"Interacting"`
#'
#' @return The input `hotspots` with one extra factor column per pair named
#'   `"A_B"`. Duplicate columns (by name) are deduplicated.
#'
#' @examples
#' df <- data.frame(barcode = paste0("s", 1:3),
#'                 A = c("a", NA, "a"),
#'                 B = c(NA, "B", "B"),
#'                 C = c(NA, NA,"C"))
#' out <- addPairwiseInteractCol(df, features = c("A", "B", "C"))
#' print(out)
#' @export
#' @seealso [addInteractCol()], [plotMultiWayViolin()]
#' @importFrom utils combn
#' 
addPairwiseInteractCol <- function(hotspots,
                                   features = NULL,
                                   interact_name = "Interacting") {
  # Choose feature columns
  if (is.null(features)) {
    features <- setdiff(colnames(hotspots), "barcode")
  }
  features <- unique(features)
  if (length(features) < 2L) {
    stop("Need at least two feature columns to form pairs.")
  }
  
  # Validate existence
  missing <- setdiff(features, colnames(hotspots))
  if (length(missing)) {
    stop(sprintf("These features are not in 'hotspots': %s",
                 paste(missing, collapse = ", ")))
  }
  
  # All unique forward pairs as rows: [cell1, cell2]
  patternPairs <- t(utils::combn(features, 2, simplify = TRUE))
  
  # If reverse=TRUE, append reversed pairs: [cell2, cell1]
  # if (isTRUE(reverse)) {
  #   patternPairs <- rbind(patternPairs, patternPairs[, c(2, 1), drop = FALSE])
  # }
  
  # Add one interaction column per (possibly reversed) pair
  for (i in seq_len(nrow(patternPairs))) {
    pair <- patternPairs[i, ]
    colname <- paste(pair[1], pair[2], sep = "_")
    hotspots <- addInteractCol(
      hotspots             = hotspots,
      interaction_cols     = pair,
      interact_name        = interact_name,
      interaction_col_name = colname
    )
  }
  hotspots <- hotspots[,!duplicated(colnames(hotspots))]
  return(hotspots)
}
#'
# --- Helper: mergeHotspotDfs ----------------------------------------------
#' Merge Pattern/Influence Tables Across Many Features
#'
#' Iterates over `features`, calling [mergeHotspotDfs()] for each, and returns a
#' single wide table containing pattern and influence columns per feature,
#' aligned by `cell_id`.
#' @import dplyr
#' @param pattern_hotspots `data.frame` with `cell_id` and feature columns.
#' @param influence_hotspots `data.frame` with `cell_id` and matching features.
#' @param pattern_cell Name of the feature column in `pattern_hotspots` to use.
#' @param influence_cell Name of the feature column in `influence_hotspots` 
#' @param cell_id Common identifier column (default `"barcode"`).
#' @param influence_suffix Suffix for influence columns (default `"-influence"`).
#' @return A wide `data.frame` with one row per `cell_id` and, for each feature
#'   `f`, two columns: `f` (pattern) and `paste0(f, influence_suffix)`.
#'
#' @examples
#' p <- data.frame(barcode = c("a","b"),
#'                 A = c(1, NA),
#'                 B = c(NA, 2))
#' i <- data.frame(barcode = c("a","b"),
#'                A = c(NA, "a"),
#'                 B = c("b", NA))
#' mergeHotspotDfs( pattern_hotspot = p, influence_hotspots = i,
#'                  pattern_cell ="A", influence_cell = "B")
#' @export
#' @seealso [mergeAllHotspotDfs()]
#' @importFrom dplyr select inner_join mutate
#' 
mergeHotspotDfs <- function(pattern_hotspots,influence_hotspots,
                            pattern_cell,influence_cell,cell_id = 'barcode',
                            influence_suffix = '-influence') {
  df1 <- pattern_hotspots %>% dplyr::select(.data[[cell_id]],
                                            .data[[pattern_cell]])
  df2 <- influence_hotspots %>% dplyr::select(.data[[cell_id]],
                                              .data[[influence_cell]] )
  influence_name <- paste0(influence_cell,influence_suffix)
  colnames(df2)[colnames(df2) == influence_cell] <- influence_name
  df <- dplyr::inner_join(df1,df2,by = cell_id )
  df <- df %>% dplyr::mutate(!!influence_name := ifelse(!is.na(.data[[influence_name]]),
                                                        influence_name,.data[[influence_name]]))
  
  return(df)
}
#' 
# --- Helper: mergeAllHotspotDfs ----------------------------------------------
#' Merge Pattern/Influence Tables Across Many Features
#'
#' Iterates over `features`, calling [mergeHotspotDfs()] for each, and returns a
#' single wide table containing pattern and influence columns per feature,
#' aligned by `cell_id`.
#'
#' @param pattern_hotspots `data.frame` with `cell_id` and feature columns.
#' @param influence_hotspots `data.frame` with `cell_id` and matching features.
#' @param features Character vector of features to merge. If `NULL`, the
#'   intersection of non-`cell_id` column names is used.
#' @param cell_id Common identifier column (default `"barcode"`).
#' @param influence_suffix Suffix for influence columns (default `"-influence"`).
#'
#' @return A wide `data.frame` with one row per `cell_id` and, for each feature
#'   `f`, two columns: `f` (pattern) and `paste0(f, influence_suffix)`.
#'
#' @examples
#' p <- data.frame(barcode = c("a","b"),
#'                 A = c(1, NA),
#'                 B = c(NA, 2))
#'i <- data.frame(barcode = c("a","b"),
#'                 A = c(NA, "a"),
#'                 B = c("b", NA))
#' mergeAllHotspotDfs(p, i, features = c("A","B"))
#' @export
#' @seealso [mergeHotspotDfs()]
#' @importFrom dplyr inner_join select distinct left_join
#' 
mergeAllHotspotDfs <- function(pattern_hotspots, influence_hotspots,
                               features = NULL,
                               cell_id = "barcode",
                               influence_suffix = "-influence") {
  # sanity
  if (!cell_id %in% names(pattern_hotspots) || !cell_id %in% names(influence_hotspots)) {
    stop(sprintf("Both data frames must contain '%s' column.", cell_id))
  }
  
  # default: use intersection of columns (minus id)
  if (is.null(features)) {
    features <- intersect(setdiff(names(pattern_hotspots), cell_id),
                          setdiff(names(influence_hotspots), cell_id))
  }
  missing <- setdiff(features, intersect(names(pattern_hotspots), names(influence_hotspots)))
  if (length(missing)) {
    stop(sprintf("Missing in inputs: %s", paste(missing, collapse = ", ")))
  }
  
  # verify identical ids
  base_ids <- dplyr::inner_join(
    pattern_hotspots %>% dplyr::select(.data[[cell_id]]),
    influence_hotspots %>% dplyr::select(.data[[cell_id]]),
    by = cell_id
  ) %>% dplyr::distinct()
  
  out <- base_ids
  
  for (f in features) {
    tmp <- mergeHotspotDfs(
      pattern_hotspots   = pattern_hotspots,
      influence_hotspots = influence_hotspots,
      pattern_cell       = f,
      influence_cell     = f,
      cell_id            = cell_id,
      influence_suffix   = influence_suffix
    )
    
    
    out <- dplyr::left_join(out, tmp, by = cell_id)
  }
  
  return(out)
}
#' 
# --- Helper: getHotSpotsGeneExpr ---------------------------------------------
#' Join Spot-Level Gene Expression to Hotspot Labels
#'
#' Builds a long table of `(gene_name, barcode, gene_exp)` from `spCounts` and
#' joins hotspot columns from `hotspots` by `barcode`. Intended as input for
#' plotting functions (e.g., [plotMultiWayViolin()]).
#'
#' @param hotspots `data.frame` with a `"barcode"` column (if absent, row names
#'   are used to construct one) and any number of hotspot label columns.
#' @param spCounts A numeric matrix or a sparse `Matrix` (genes x barcodes)
#'   containing expression values.
#' @param gene Character vector of gene names (row names of `spCounts`) to
#'   extract. If `NULL`, all genes are used.
#' @param ... Ignored; reserved for compatibility.
#'
#' @return A `data.frame` with columns `gene_name`, `barcode`, `gene_exp` and
#'   all columns from `hotspots`.
#'
#' @details
#' If `spCounts` is a sparse `Matrix`, only **non-zero** entries are returned
#' (via `Matrix::summary`). If you require explicit zeros, convert the matrix to
#' dense first (e.g., `as.matrix(spCounts)`), understanding this can be memory
#' intensive.
#'
#' @examples
#' set.seed(1)
#' m <- matrix(rpois(12, 1), nrow = 3,
#'             dimnames = list(c("GeneA","GeneB","GeneC"),
#'                             paste0("bc",1:4)))
#' hs <- data.frame(barcode = paste0("bc",1:4),
#'                  A = c("A", NA, "A", NA),
#'                  B = c(NA, "B", NA, "B"))
#' ge <- getHotSpotsGeneExpr(hotspots = hs, m, gene = "GeneA")
#' print(ge)
#' 
#' @export
#' @seealso [plotMultiWayViolin()]
#' @importFrom Matrix summary
#' @importFrom dplyr inner_join
#' 
getHotSpotsGeneExpr <- function(hotspots, spCounts, gene = NULL, ...) {
  
  # ensure a 'barcode' column exists for the join (fallback to rownames if needed)
  if (!"barcode" %in% names(hotspots)) {
    hotspots <- cbind(barcode = rownames(hotspots), hotspots)
    rownames(hotspots) <- NULL
  }
  
  # allow all genes by default (or a vector)
  if (is.null(gene)) gene <- rownames(spCounts)
  
  # build one long table for all selected genes
  mat_sub <- spCounts[gene, , drop = FALSE]
  
  if (inherits(mat_sub, "Matrix")) {
    sm <- Matrix::summary(mat_sub)  # non-zeros only (sparse fast path)
    spCounts_df <- data.frame(
      gene_name = rownames(mat_sub)[sm$i],
      barcode   = colnames(mat_sub)[sm$j],
      gene_exp  = sm$x,
      stringsAsFactors = FALSE
    )
  } else {
    tmp <- as.data.frame(as.table(as.matrix(mat_sub)), stringsAsFactors = FALSE)
    spCounts_df <- tmp[, c(3, 2, 1)]
    colnames(spCounts_df) <- c("gene_exp", "barcode", "gene_name")
  }
  
  # join only everything
  gene_df <- dplyr::inner_join(spCounts_df, hotspots, by = "barcode")
  return(gene_df)
}
# --- Main: plotMultiWayViolin -------------------------------------------------
#' Plot Multi-Way Violin for a Gene Across Pattern/Influence Hotspots
#'
#' Produces a violin + jitter plot of a single gene's expression across a set
#' of hotspot categories derived from `cell1`, `cell2`, and their *influence*
#' interactions. Intended to visualize how expression differs in the pattern
#' cells, the influenced counterparts, and (optionally) their legacy pair.
#'
#' @param hotspots A `data.frame` that already contains:
#'   - an identifier column given by `cell_id` (default `"barcode"`);
#'   - a gene identifier column `gene_id_col` (default `"gene_name"`);
#'   - a numeric expression column `counts_col` (default `"gene_exp"`);
#'   - hotspot label columns for `cell1`, `cell2`, and their corresponding
#'     influence columns (e.g., `"cell1-influence"`, `"cell2-influence"`), and
#'     optionally a legacy pair column (`"cell1_cell2"` or `"cell2_cell1"`).
#'   Typically this is the output of [getHotSpotsGeneExpr()].
#' @param cell1 Character scalar: the name of the first cell pattern column to
#'   include on the x-axis.
#' @param cell2 Character scalar: the name of the second cell pattern column to
#'   include on the x-axis.
#' @param gene Character scalar: the gene to plot (must match values in
#'   `gene_id_col`).
#' @param gene_id_col Column name holding gene identifiers (default `"gene_name"`).
#' @param counts_col Column name holding expression values (default `"gene_exp"`).
#' @param influence_id Substring used to identify influence columns. If `NULL`,
#'   it is auto-detected as `"influence"` (fallback `"near"`). This string is
#'   matched in column names to select e.g. `"cell1-influence"`.
#' @param cell_id Identifier column name (default `"barcode"`). If `NULL`,
#'   row names are treated as barcodes.
#' @param violin_order Optional character vector specifying the order of
#'   hotspot groups on the x-axis. If `NULL`, the order is inferred as
#'   `c(cell1, <influences>, cell2, [legacy_if_present])`.
#' @param colors Fill colors passed to `scale_fill_manual()`. The length should
#'   match the number of hotspot groups plotted. If a legacy `cell1_cell2` (or
#'   `cell2_cell1`) column is detected, `"darkgreen"` is appended to `colors`.
#' @param out File path to save the plot (e.g., `"SpaceMarkerViolin.png"`).
#'   Use `NULL` to skip saving and only return the `ggplot` object.
#' @param width,height Dimensions (inches) for the saved figure.
#'
#' @return A `ggplot2` object. If `out` is non-`NULL`, the plot is also written
#'   to disk via `ggsave()`.
#'
#' @details
#' The function:
#' 1. Detects and selects relevant hotspot columns containing `cell1`, `cell2`,
#'    and the `influence_id` substring, while avoiding "double-influence"
#'    columns that contain the substring twice.
#' 2. Optionally includes a legacy `cell1_cell2` (or `cell2_cell1`) column if
#'    present.
#' 3. Filters rows to the requested `gene`, pivots hotspot labels to long format
#'    (keeping only non-`NA` labels), and draws a violin + jitter plot with
#'    manual fills and angled x-axis labels.
#'
#' @examples
#' # Minimal reproducible sketch:
#' set.seed(42)
#' counts <- matrix(rpois(12, 2), nrow = 3,
#'                 dimnames = list(c("GeneA","GeneB","GeneC"), paste0("bc",1:4)))
#' #Hotspot labels (NA = absent, non-NA = present)
#' ## Pattern Hotspots
#' pat_hs <- data.frame(
#'  barcode = paste0("bc",1:4),
#'   Acell   = c("Acell", NA, "Acell", NA),
#'    Bcell   = c(NA, "Bcell", NA, "Bcell"))
#'    ## Influence hot
#' inf_hs <- data.frame(
#' barcode = paste0("bc",1:4),
#' Acell = c(NA, "Acell", "Acell", NA),
#' Bcell = c(NA, NA, "Bcell", NA),
#' check.names = FALSE)
#' # Combine pattern and influence hotspots
#' features <- setdiff(colnames(pat_hs),c("barcode","x","y"))
#' all_hs <- mergeAllHotspotDfs(pattern_hotspots = pat_hs,
#' influence_hotspots = inf_hs,features = features)
#' # Identify pairwise interactions between pattern hotspots and influence hotspots
#' pairwise_hs <- addPairwiseInteractCol(all_hs)
#' gene_df <- getHotSpotsGeneExpr(pairwise_hs, counts, gene = "GeneA")
#' # Plot 4 way (optional 5) of pattern-influence
#' p <- plotMultiWayViolin(
#' hotspots = gene_df,
#' cell1 = "Acell", cell2 = "Bcell",
#' gene = "GeneA",
#' influence_id = "influence",
#' out = NULL)
#' print(p)
#' @export
#' @seealso [getHotSpotsGeneExpr()], [addPairwiseInteractCol()], [addInteractCol()]
#' @importFrom dplyr filter
#' @importFrom reshape2 melt
plotMultiWayViolin <- function(hotspots, cell1, cell2,gene, gene_id_col = "gene_name", counts_col = "gene_exp",
                               influence_id = NULL,cell_id ="barcode",violin_order = NULL,
                               colors = c("blue",'yellow',"red","orange"),
                               out = "SpaceMarkerViolin.png",
                               width = 10, height = 10) {
  # check for influence columns
  if(is.null(influence_id)){
    if(any(grepl("influence",colnames(hotspots)))){
      influence_id <- "influence"
    } else if (any(grepl("near",colnames(hotspots)))) {
      influence_id <- "near"
    } else {
      stop("Please specify an 'influence_id' to distinguish between pattern
      hotspots and influence hotspots")
    }
  }
  # check for barcode in columns or rownames
  if (is.null(cell_id)){
    message("Assuming the cell ids are in the rows")
    hotspots <- data.frame("barcode" =rownames(hotspots),hotspots)
    cell_id <- "barcode"
  } 
  
  
  # remove dual influence eg cell1-influence_cell2-influence
  x <- colnames(hotspots)
  k <- lengths(regmatches(
    x, gregexpr(influence_id, x, fixed = TRUE)))
  x <- x[k < 2]
  # Get pattern and pattern influence interactions
  sel <- grepl(cell1,x, ignore.case = FALSE) &
    grepl(influence_id, x, ignore.case = TRUE) &
    grepl(cell2,x, ignore.case = FALSE)
  
  non_hotspot_columns <- c(cell_id,gene_id_col,counts_col)
  hotspot_columns <- c(cell1,x[sel],cell2)
  
  # Finally check if cell1_cell2 is in the column; if not continue
  # Construct both possible legacy names
  legacy1 <- paste0(cell1, "_", cell2)
  legacy2 <- paste0(cell2, "_", cell1)
  
  # Check which exists
  if (legacy1 %in% x) {
    legacy <- legacy1
    idx <- match(legacy1, x)   # index of the column
  } else if (legacy2 %in% x) {
    legacy <- legacy2
    idx <- match(legacy2, x)
  } else {
    legacy <- NA
    idx <- NA
  }
  
  # Update hotspot_columns and colors if a match was found
  if (!is.na(legacy)) {
    hotspot_columns <- c(hotspot_columns, legacy)
    colors <- c(colors, "darkgreen")
  }
  
  selected_columns <- c(non_hotspot_columns,hotspot_columns)
  hotspots <- hotspots[,selected_columns]
  if (is.null(violin_order)){
    violin_order <- hotspot_columns
  }
  
  gene_df <- hotspots %>% dplyr::filter(.data[[gene_id_col]] == gene)
  
  # keep only the columns you need first
  tmp <- gene_df[, selected_columns]
  
  plot_df_long <- reshape2::melt(
    tmp,
    id.vars      = setdiff(selected_columns, hotspot_columns),  
    measure.vars = hotspot_columns,
    variable.name = "hotspots",
    value.name    = "label",
    na.rm         = TRUE  
  )
  
  # match your final select() (drop the 'label' column)
  plot_df_long <- plot_df_long[, c("barcode", "gene_name", "gene_exp", "hotspots")]
  
  

  plot_df_long$hotspots <- as.character(plot_df_long$hotspots)
  
  plot_df_long$hotspots <- factor(plot_df_long$hotspots, levels = violin_order )
  #return(gene_df)
  # set plot title from gene_name
  plot_title <- plot_df_long$gene_name[1]
  
  
  # build violin + jitter plot
  p <- ggplot(plot_df_long, aes(x = hotspots, y = gene_exp, fill = hotspots)) +
    geom_violin(trim = TRUE) +
    geom_jitter(
      width = 0.2,
      size  = 0.5,
      shape = 21,
      color = "black"
    ) +
    scale_fill_manual(values = colors) +
    labs(
      title = plot_title,
      x     = "Hotspots",
      y     = "Gene Expression"
    ) + theme_minimal(base_size = 20) + theme(
      axis.text.x = element_text(angle = 45, hjust = 0.8))
  
  if (!is.null(out)) ggsave(out, plot   = p,width  = width,height = height,
                            units  = "in")
  
  return(p)
}