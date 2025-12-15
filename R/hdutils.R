#' @title Compute the spatial influence of a spatial feature
#' @description This function computes the spatial influence of a specified pattern
#' @param spPatterns A data frame containing x, y coordinates and pattern name
#' @param optParams A data frame with optimal parameters for the pattern
#' @param ... Additional parameters for the Smooth function
#' @return A data frame with the spatial influence of the specified pattern
#' @export
calculate_influence <- function(spPatterns, optParams,...) {
    patnames <- setdiff(colnames(spPatterns),
                       c("x", "y", "barcode"))

    allwin <- spatstat.geom::owin(
    range(spPatterns$x),
    range(spPatterns$y))
    X <- spatstat.geom::ppp(x = spPatterns$x, y = spPatterns$y,
                            window = allwin,
                            marks = spPatterns[,patnames[1]])

    spInfluence <- sapply(patnames, function(pat) {
    # Create a point pattern object for each pattern
    X <- spatstat.geom::ppp(x = spPatterns$x, y = spPatterns$y,
                            window = allwin,
                            marks = spPatterns[,pat])
    
    # Calculate the kernel for the specified pattern
    Kact1 <- spatstat.explore::Smooth(
      X, at = "points", sigma = optParams[1,pat],...)
    
    # Plot the K-function
    return(Kact1)
    })
    spInfluence <- as.data.frame(spInfluence)
    colnames(spInfluence) <- patnames
    spInfluence <- cbind(spPatterns[,c("barcode","x", "y")], spInfluence)

    return(spInfluence)
}

#' @title Compute the threshold for identifying outlier values or hotspots
#' @description This function computes the threshold for identifying outlier 
#' values or hotspots by fitting a normal mixture model to the data.
#' @param df A vector containing pattern values
#' @param minval Minimum value for quantile threshold
#' @param maxval Maximum value for quantile threshold
#' @param method Method to use for threshold calculation. Options are "abs" for absolute (default) and "pct" for percentile.
#' @return A list containing the computed thresholds
#' @importFrom mixtools normalmixEM
.calc_threshold <- function(df, minval = 0.01, maxval = 0.99, method=c("abs","pct")) {
if (method[1]=="pct"){
    minthresh <- quantile(df,minval)
    maxthresh <- quantile(df,maxval)
} else {
   minthresh <- minval
   maxthresh <- maxval
}
#Calculate the two compmonents of the normal mixture model
res <- try(mixtools::normalmixEM(df, k = 2, maxit = 1000, epsilon = 1e-8), silent = TRUE)
if (inherits(res, "try-error") || is.null(res) || !all(c("mu", "sigma") %in% names(res)) || res$ft == 1000) {
    warning("mixtools::normalmixEM failed or reached maxit; using minval as threshold")
    return(minthresh)
}
comps <- res

# Identify smaller component
small <- which.min(comps$mu)

thresh <- min(max(comps$mu[small] + (comps$sigma[small] * 4), minthresh),maxthresh)

return(thresh)
}
#' @title Compute the thresholds for all columns in a data frame
#' @description This function computes the thresholds for all columns in a 
#' data frame. The data frame could be an spPatterns object or an spInfluence
#' object.
#' @param df A data frame with pattern values (optionally 
#' with x, y, barcode columns)
#' @param minvals Minimum value for quantile threshold
#' @param maxvals Maximum value for quantile threshold
#' @param ... Additional parameters to pass to lower level functions
#' @return A list containing the computed thresholds for each pattern
#' @export 
calculate_thresholds <- function(df, minvals = 0.01, maxvals = 0.99,...) {
  
    patnames <- setdiff(colnames(df), c("x", "y", "barcode"))
    #Check if minvals and maxvals are vectors
    if (length(minvals) == 1) minvals <- rep(minvals, length(patnames))
    if (length(maxvals) == 1) maxvals <- rep(maxvals, length(patnames))
    # Check if minvals and maxvals are of the same length and match 
    # the number of columns in df
    if ((length(minvals) != length(maxvals)) || length(patnames) != length(minvals)) {
        stop("minvals and maxvals must be scalar or vectors of 
              the length ncol(df)")
    }
    names(minvals) <- patnames
    names(maxvals) <- patnames

    # Calculate thresholds for all patterns
    thresholds <- sapply(patnames, function(pat) {
        .calc_threshold(df[,pat], minval=minvals[pat], maxval=maxvals[pat],...)
    })
    return(thresholds)
}


#' @title Find hotspots for all patterns or influences based on values
#' @description Convenience function to find hotspots for all spatial patterns 
#' or influence dataframes based on provided thresholds
#' @inheritParams calculate_thresholds
#' @param threshold a scalar or vector of thresholds for each column in the data frame.
#'  Either user provided or the output of @calculate_thresholds
#' @return a data frame with the same dimensions as the input data frame.
#' @export 
find_hotspots_gmm <- function(df, threshold = 0.1,...){
    patnames <- setdiff(colnames(df),c("x","y","barcode"))
    if (length(threshold)==1){
        threshold <- rep(threshold,length(patnames))
    }
    if (length(threshold)!=length(patnames)){
        stop("Length of threshold must be 1 or equal to number of patterns.")
    }
    names(threshold) <- patnames
    
    hotspots <- matrix(NA, nrow=nrow(df), ncol=length(patnames))
    colnames(hotspots) <- patnames
    for (pat in patnames){
        hotspots[,pat] <- ifelse(
            df[,pat]>threshold[pat],pat,NA)
    }
    hotspots <- cbind(df[c("barcode","y","x")],hotspots)
    row.names(hotspots) <- hotspots$barcode
    hotspots <- as.data.frame(hotspots)
    return(hotspots)
}

#' Classify spots into interacting / non-interacting pattern regions
#'
#' @description
#' Given pattern hotspots and their corresponding influence hotspots,
#' classify each spot into:
#' - undirected interacting region (when \code{influence_hotspots = NULL}), or
#' - directed regions for each pattern pair (e.g. \code{"Epi_near_Plasma"},
#'   \code{"Plasma_near_Epi"}).
#'
#' The function expects identical dimensions and pattern columns for
#' \code{pat_hotspots} and \code{influence_hotspots}. When only
#' \code{pat_hotspots} is provided, a single undirected column is created.
#'
#' @param pat_hotspots Data frame of pattern hotspots with columns:
#'   \code{x}, \code{y}, \code{barcode}, and one or more pattern columns.
#' @param influence_hotspots Data frame of influence hotspots with the same
#'   dimensions and pattern columns as \code{pat_hotspots}. If \code{NULL},
#'   a symmetric (undirected) classification is returned. Default: \code{NULL}.
#' @param patternpair Character vector of length 2 giving the two pattern
#'   names to classify, e.g. \code{c("Epi", "Plasma")}. If \code{NULL} and
#'   more than 2 patterns are present, an error is thrown.
#'
#' @return
#' A data frame of region labels:
#' \itemize{
#'   \item If \code{influence_hotspots} is \code{NULL}: one column named
#'         \code{"<pattern1>_<pattern2>"} with values
#'         \code{"Interacting"}, \code{pattern1}, \code{pattern2}, or \code{NA}.
#'   \item Otherwise: two columns
#'         \code{"<pattern1>_near_<pattern2>"} and
#'         \code{"<pattern2>_near_<pattern1>"} with values
#'         \code{"Interacting"}, \code{pattern1} / \code{pattern2}, or \code{NA}.
#' }
#'
#' @examples
#' \dontrun{
#' pat_hotspots <- data.frame(
#'   x = 1:5,
#'   y = 1:5,
#'   barcode = paste0("Spot", 1:5),
#'   Epi    = c(1, NA, 1, NA, NA),
#'   Plasma = c(NA, 1, NA, 1, NA)
#' )
#'
#' influence_hotspots <- data.frame(
#'   x = 1:5,
#'   y = 1:5,
#'   barcode = paste0("Spot", 1:5),
#'   Epi    = c(NA, 1, 1, NA, NA),
#'   Plasma = c(1, NA, NA, 1, 1)
#' )
#'
#' # Directed classification
#' cls <- classify_spots(
#'   pat_hotspots       = pat_hotspots,
#'   influence_hotspots = influence_hotspots,
#'   patternpair        = c("Epi","Plasma")
#' )
#'
#' # Undirected (influence_hotspots = NULL)
#' und <- classify_spots(
#'   pat_hotspots       = pat_hotspots,
#'   influence_hotspots = NULL,
#'   patternpair        = c("Epi","Plasma")
#' )
#' }
#' @export 
classify_spots <- function(pat_hotspots,
                           influence_hotspots,
                           patternpair = NULL) {
  patnames <- setdiff(colnames(pat_hotspots), c("x", "y", "barcode"))
  infnames <- setdiff(colnames(influence_hotspots), c("x", "y", "barcode"))
  
  if (!all(dim(pat_hotspots) == dim(influence_hotspots))) {
    stop("pat_hotspots and influence_hotspots must have the same dimensions.")
  } else if (!all(patnames == infnames)) {
    stop("pat_hotspots and influence_hotspots must have the same column names.")
  }
  
  if (!is.null(patternpair) && !all(patternpair %in% patnames)) {
    stop("patternpair must be NULL or contained in patnames.")
  } else if (is.null(patternpair) && length(patnames) > 2) {
    stop("More than 2 patterns found. Please provide patternpair.")
  }
  
  ## Case 1: influence_hotspots is NULL  ----
  if (is.null(influence_hotspots)) {
    message("You only provided the pattern hotspots. Specify influence hotspots if present")
    
    region <- ifelse(
      !is.na(pat_hotspots[, patternpair[1]]) &
        !is.na(pat_hotspots[, patternpair[2]]),
      "Interacting",
      ifelse(
        !is.na(pat_hotspots[, patternpair[1]]),
        pat_hotspots[, patternpair[1]],
        ifelse(
          !is.na(pat_hotspots[, patternpair[2]]),
          pat_hotspots[, patternpair[2]],
          NA
        )
      )
    )
    
    df <- data.frame(x = region)
    colnames(df)[1] <- paste0(patternpair[1], "_", patternpair[2])
    return(df)
  }
  
  ## Case 2: pattern + influence hotspots ----
  region <- pat_hotspots[, patternpair[1]]
  pat1.inf2 <- ifelse(
    !is.na(region) & !is.na(influence_hotspots[, patternpair[2]]),
    "Interacting",
    ifelse(!is.na(region), region, NA)
  )
  
  region <- pat_hotspots[, patternpair[2]]
  pat2.inf1 <- ifelse(
    !is.na(region) & !is.na(influence_hotspots[, patternpair[1]]),
    "Interacting",
    ifelse(!is.na(region), region, NA)
  )
  
  df <- data.frame(
    pat1.inf2 = pat1.inf2,
    pat2.inf1 = pat2.inf1
  )
  colnames(df) <- c(
    paste0(patternpair[1], "_near_", patternpair[2]),
    paste0(patternpair[2], "_near_", patternpair[1])
  )
  
  return(df)
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

#' @title Perform row-wise t-tests from scratch
#' @description This function iterates over the rows of a matrix and performs a
#' t-test comparing two groups of columns. It calculates the t-statistic, p-value,
#' and sample sizes without relying on `stats::t.test()` for the core logic.
#' @param in.data A numeric matrix. Rows represent features, columns represent samples.
#' @param region A factor or vector indicating the group membership for each column of `in.data`.
#'                Must have exactly two levels/unique values. Its length must equal `ncol(in.data)`.
#' @param min_bins Minimum number of non-missing observations required in each group to perform the t-test.
#' @param ... Additional parameters to pass to the t-test function.
#' @return A matrix with rows corresponding to the features and columns:
#'          - `statistic`: The calculated t-statistic.
#'          - `p.value`: The calculated two-sided p-value.
#'         - `n1`: Number of non-missing observations in group 1 for that row.
#'         - `n2`: Number of non-missing observations in group 2 for that row.
#' @importFrom stats t.test
#' @importFrom effsize cohen.d
.row_t_test <- function(in.data, region, min_bins=50, ...){
    if (!is.factor(region)) {
        region <- factor(region)
    }
    if (nlevels(region) != 2) {
        stop("'region' must have exactly two levels (groups).")
    }
    group_levels <- levels(region)
    int <- which(as.character(group_levels)=="Interacting")
    interacting <- group_levels[int]
    patname <- group_levels[-int]
    idx_pat <- which(region == patname)
    idx_int <- which(region == interacting)
    t_scores <- sapply(rownames(in.data), function(r) {
        pat <- in.data[r, idx_pat]
        inter <- in.data[r, idx_int]
        if (length(pat) < min_bins || length(inter) < min_bins) {
            return(c(statistic=NA, p.value=NA, n1=0, n2=0, effect_size=NA))
        }
        tmp <- t.test(x=inter, y=pat, 
                      alternative = "two.sided", var.equal = FALSE, 
                      na.action = na.omit)
        effect_size <- effsize::cohen.d(inter, pat)$estimate
        return(c(statistic=tmp$statistic, p.value=tmp$p.value, n1=length(inter), n2=length(pat), effect_size=effect_size))
    })

    t_scores <- t(t_scores)
    colnames(t_scores) <- c("statistic", "p.value", "n1", "n2", "effect_size")
    return(t_scores)
}
#' @title Calculate interaction scores for a specific pattern pair
#' @description This function calculates interaction scores for a specific pattern pair
#' using the ` classify_spots` function to determine the region of each spot.
#' @param data A numeric matrix with genes as rows and barcodes as columns.
#' @param pat_hotspots A data frame with pattern hotspots, containing columns for x, y, and barcode.
#' @param influence_hotspots A data frame with influence hotspots, containing columns for x, y, and barcode.
#' @param patternpair A character vector of length 2 specifying the pattern pair to analyze.
#' @param avoid_confounders Logical (default=FALSE) indicating whether to avoid confounding effects due to colocalization.
#' @param ... Additional parameters to pass to lower level functions.
#' @return A data frame with interaction scores for the specified pattern pair.
.calc_IM_scores <- function(data, pat_hotspots, influence_hotspots, patternpair, avoid_confounders=FALSE,...) {
    spotClass <-  classify_spots(pat_hotspots, influence_hotspots, patternpair = patternpair)
    pat1 <- patternpair[1]
    pat2 <- patternpair[2]
    region <- spotClass[,1]
    if (avoid_confounders==TRUE)
        region[which(!is.na(spotClass[,2]))] <- NA
    if (sum(!is.na(unique(region)))==2){
        t1table <- .row_t_test(data, region=region)

    } else {
       t1table <- matrix(NA, nrow=nrow(data), ncol=5)
       colnames(t1table) <- c("statistic", "p.value", "n1", "n2", "effect_size")
    }
    t1table <- as.data.frame(t1table)
    t1table$gene <- rownames(t1table)
    t1table$cell_interaction <- paste0(pat1, "_near_", pat2)

    region <- spotClass[,2]
    if (avoid_confounders==TRUE)
        region[which(!is.na(spotClass[,1]))] <- NA
    if (sum(!is.na(unique(region)))==2){
        t2table <- .row_t_test(data, region=region)
    } else {
        t2table <- matrix(NA, nrow=nrow(data), ncol=5)
        colnames(t2table) <- c("statistic", "p.value", "n1", "n2", "effect_size")
    }
    t2table <- as.data.frame(t2table)
    t2table$gene <- rownames(t2table)
    t2table$cell_interaction <- paste0(pat2, "_near_", pat1)

    tscores <- rbind(t1table, t2table)
    tscores$p.adj <- p.adjust(tscores$p.value, method = "BH")
    return(tscores)
}

#' @title Calculate interaction scores for all pattern pairs
#' @description This function calculates interaction scores for all pattern pairs
#' using the `.calc_IM_scores` function. It can run in parallel if BiocParallel is available.
#' @param data A numeric matrix with genes as rows and barcodes as columns.
#' @param pat_hotspots A data frame with pattern hotspots, containing columns for x, y, and barcode.
#' @param influence_hotspots A data frame with influence hotspots, containing columns for x, y, and barcode.
#' @param pattern_pairs A data frame with pattern pairs to calculate interaction scores for. If NULL, 
#' all combinations of patterns in `pat_hotspots` will be used.
#' If provided, it should have two columns with pattern names. 
#' Each row should represent a pair of patterns for which interaction scores will be calculated.
#' @param ... Additional parameters to pass to lower level functions.
#' @return A data frame with interaction scores for all pattern pairs.
#' @export
calculate_gene_scores_directed <- function(data, pat_hotspots, influence_hotspots, pattern_pairs=NULL,...) {
    if (is.null(pattern_pairs)) {
        pattern_pairs <- utils::combn(setdiff(colnames(pat_hotspots), c("x", "y", "barcode")), 2, simplify = FALSE)
    }
    use_biocparallel <- (requireNamespace("BiocParallel", quietly = TRUE)) && 
                                        (BiocParallel::bpparam()$workers > 1) &&
                                        (!is.null(nrow(pattern_pairs)) && (nrow(pattern_pairs) > 1))
    if (use_biocparallel) {
        bpp <- BiocParallel::bpparam()
        IMscores_list <- BiocParallel::bplapply(
            seq_len(nrow(pattern_pairs)),
            function(i) {
                patternpair <- pattern_pairs[i,]
                IMscores.pair <- .calc_IM_scores(data, pat_hotspots, influence_hotspots, patternpair,...)
                return(IMscores.pair)
            },
            BPPARAM = bpp
        )
        IMscores <- do.call(rbind, IMscores_list)
        message("Processed all pattern pairs in parallel.\n")
    } else {
        IMscores <- c()
        for (i in 1:nrow(pattern_pairs)) {
            patternpair <- pattern_pairs[i,]
            IMscores.pair <- .calc_IM_scores(data, pat_hotspots, influence_hotspots, patternpair,...)
            IMscores <- rbind(IMscores, IMscores.pair)
            message("Processed pattern pair: ", patternpair[1], " and ", patternpair[2], "\n")
            if (i %% 10 == 0) {
                message("Processed ", i, "pattern pairs so far.\n")
            }
        }
    }
    return(IMscores)
}

#' Create a long-format ligand–receptor dataframe
#'
#' Convert ligand, receptor, and ligand–receptor (LR) score matrices into a
#' tidy long-format dataframe, with one row per ligand–receptor interaction and
#' source–target cell-type pair.
#'
#' @param ligand_scores A numeric matrix or data.frame of ligand scores with
#'   ligands in rows and source/target context in columns
#'   (e.g. \code{"B_near_Epi"}).
#' @param receptor_scores A numeric matrix or data.frame of receptor scores
#'   with receptor complexes in rows and target cell types in columns
#'   (e.g. \code{"Epi"}, \code{"Mast"}).
#' @param lrscores A numeric matrix of ligand–receptor scores with
#'   interactions in rows (e.g. \code{"TGFB1_TGFBR1_TGFBR2"}) and
#'   source-to-target directions in columns
#'   (e.g. \code{"B_to_Epi"}).
#' @param complex_sep A character string specifying the separator used between
#'   receptor components in \code{receptor_scores} rownames
#'   (default is \code{", "}).
#'
#' @return A data.frame with columns:
#'   \item{source_celltype}{Source cell type.}
#'   \item{target_celltype}{Target cell type.}
#'   \item{ligand}{Ligand gene/symbol.}
#'   \item{receptor}{Receptor complex string.}
#'   \item{source_to_target}{Source-to-target label (e.g. \code{"B_to_Epi"}).}
#'   \item{interaction}{Ligand–receptor interaction identifier.}
#'   \item{ligand_score}{Ligand score for this interaction/direction.}
#'   \item{receptor_score}{Receptor score for this interaction/target.}
#'   \item{lr_score}{Combined ligand–receptor score.}
#'
#' @importFrom reshape2 melt
#' @importFrom dplyr select distinct
#'
#' @examples
#' \dontrun{
#' df <- create_lr_dataframe(ligand_scores, receptor_scores, lrscores)
#' head(df)
#' }
create_lr_dataframe <- function(ligand_scores, receptor_scores, lrscores,complex_sep = ", ") {
  
  # --- 1. Melt lrscores to long format ---
  # 'interaction' is in rownames, 'source_to_target' is in colnames
  lr_long <- reshape2::melt(as.matrix(lrscores), varnames = c("interaction", "source_to_target"), value.name = "lr_score")
  
  # Remove rows where lr_score is NA or NaN (using is.finite is a good way to check for valid numbers)
  #lr_long <- lr_long[is.finite(lr_long$lr_score) & !is.na(lr_long$lr_score), ]
  
  # --- 2. Extract components from 'source_to_target' ---
  lr_long$source_celltype <- gsub("_to_.*", "", lr_long$source_to_target)
  lr_long$target_celltype <- gsub(".*_to_", "", lr_long$source_to_target)
  
  # --- 3. Extract components from 'interaction' ---
  # Split the interaction string (Ligand_Receptor1_Receptor2...) by underscore
  # This part assumes the ligand is the first component, and receptors are the rest.
  # The first component before the first underscore is the ligand.
  lr_long$ligand <- sub("_.*", "", lr_long$interaction)
  
  # The receptor is the rest of the string after the first underscore.
  # This is a simplification; in a multi-receptor complex, this column will contain
  # all receptors concatenated by underscores (e.g., TGFBR1_TGFBR2).
  lr_long$receptor <- sub("^[^_]*_", "", lr_long$interaction)
  
  # --- 4. Prepare and merge ligand_scores ---
  
  # Melt ligand_scores to long format
  ligand_long <- reshape2::melt(as.matrix(ligand_scores), varnames = c("interaction", "ligand_col"), value.name = "ligand_score")
  
  # Remove NA/NaN scores
  #ligand_long <- ligand_long[is.finite(ligand_long$ligand_score) & !is.na(ligand_long$ligand_score), ]
  
  
  # Match the column names in ligand_scores (e.g., B_near_Epi) to the lr_long format (e.g., B_to_Epi)
  ligand_long$source_to_target <- gsub("near_", "to_", ligand_long$ligand_col)
  
  # Select and rename columns for merging
  ligand_scores_to_merge <- ligand_long %>%
    dplyr::select(interaction, source_to_target, ligand_score)
  
  # Merge with lr_long
  # lr_df <- merge(lr_long, ligand_scores_to_merge, by = c("interaction", "source_to_target"), all.x = TRUE)
  #lr_long$source_to_target <- NULL
  # Merge with lr_long
  lr_df <- merge(lr_long, ligand_scores_to_merge, by = c("interaction", "source_to_target"), all.x = TRUE)
  
  
  # --- 5. Prepare and merge receptor_scores ---
  
  # The column names in receptor_scores are the target cell types.
  # The row names are the receptor components (comma-separated, e.g., TGFBR2, TGFBR1)
  
  # Convert rownames to a proper column before melting
  rownames(receptor_scores) <- gsub(pattern = complex_sep,replacement = "_",rownames(receptor_scores))
  receptor_scores <- data.frame(receptor_scores)
  receptor_scores$receptor_components <- rownames(receptor_scores)
  # Remove any suffixes after the last dot in receptor_components (if present) 
  # e.g., TGFBR2_TGFBR1.1 -> TGFBR2_TGFBR1
  receptor_scores$receptor_components <- sub("\\.[^.]*$", "", receptor_scores$receptor_components)
  
  # Melt receptor_scores
  receptor_long <- reshape2::melt(receptor_scores, id.vars = "receptor_components", 
                                  variable.name = "target_celltype_col", 
                                  value.name = "receptor_score")
  
  
  # Get everything after the first underscore since ligands typically don't have underscores
  lr_df$receptor_components_match <- sub("^[^_]*_", "", lr_df$interaction)
  
  # The target_celltype_col is the target_celltype from lr_df
  receptor_long$target_celltype <- as.character(receptor_long$target_celltype_col)
  
  # Select and rename columns for merging
  receptor_scores_to_merge <- receptor_long %>%
    dplyr::select(receptor_components, target_celltype, receptor_score)
  
  # Merge by 'receptor_components_match' and 'target_celltype'
  lr_df <- merge(lr_df, receptor_scores_to_merge, 
                 by.x = c("receptor_components_match", "target_celltype"), 
                 by.y = c("receptor_components", "target_celltype"), 
                 all.x = TRUE)
  
  # --- 6. Final cleanup and column selection ---
  
  # Re-order and select the specified columns
  final_df <- lr_df %>%
    dplyr::select(source_celltype, target_celltype, ligand, receptor, 
                  source_to_target, interaction, ligand_score, 
                  receptor_score, lr_score) %>%
    # Remove the intermediate 'receptor_components_match' column used for merging
    dplyr::distinct() # Use distinct to remove any duplicate rows introduced by the merging logic
  
  # Convert target_celltype/source_celltype to character if they were factors
  final_df$source_celltype <- as.character(final_df$source_celltype)
  final_df$target_celltype <- as.character(final_df$target_celltype)
  final_df$source_to_target <- as.character(final_df$source_to_target)
  final_df$interaction <- as.character(final_df$interaction)
  
  final_df$lr_score[is.na(final_df$lr_score)] <- 0
  final_df$lr_score[which(final_df$lr_score == "NaN")] <- 0
  return(final_df)
}
