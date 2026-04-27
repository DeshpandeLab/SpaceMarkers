#' @importFrom sparseMatrixStats rowSums2
#' @importFrom matrixStats rowRanks
#' @importFrom matrixTests row_kruskalwallis
#' @importFrom stats pnorm
#' @importFrom qvalue qvalue
## author: Atul Deshpande
## email: adeshpande@jhu.edu

# =============================================================================
# .row_dunn_test
# =============================================================================
# Receives a dense matrix (expressed genes x pair spots) and a region vector
# that may contain NAs. The NA-subsetting is retained for interface
# compatibility with get_interacting_genes; when called from
# .find_genes_of_interest the input is already pair-subsetted so
# keep_cols is all-TRUE, which is a no-op.
#


.row_dunn_test <- function(in.data, region, pattern1, pattern2) {

  keep_cols <- !is.na(region)
  rsub      <- region[keep_cols]

  # matrixStats::rowRanks handles dense matrices correctly.
  # in.data is always dense at this point (converted in .find_genes_of_interest).
  in.ranks <- matrixStats::rowRanks(
    in.data[, keep_cols, drop = FALSE],
    ties.method = "average"
  )

  N  <- length(rsub)
  NI <- sum(rsub == "Interacting")
  N1 <- sum(rsub == pattern1)
  N2 <- sum(rsub == pattern2)

  idx_int <- which(rsub == "Interacting")
  idx_p1  <- which(rsub == pattern1)
  idx_p2  <- which(rsub == pattern2)

  # ── Exact tie correction ───────────────────────────────────────────────────
  # For each gene, compute sum over tie groups of t*(t^2-1) where t is the
  # size of each tied group. This matches kruskal.test() and the original
  # matrixStats-based implementation exactly.
  tiesStat <- apply(in.ranks, 1, function(rr) {
    rr_table <- table(rr)
    sum(rr_table * (rr_table^2 - 1))
  })
  tiesStat2 <- sqrt(1 - tiesStat / N / (N - 1) / (N + 1))

  # Guard: constant genes (all tied) produce tiesStat2 = 0.
  # Set to 1 so z-scores are 0 rather than NaN; these genes are zeroed out
  # downstream by the zero_genes correction anyway.
  tiesStat2[tiesStat2 == 0] <- 1

  # ── Mean ranks per region ──────────────────────────────────────────────────
  mInt <- rowMeans(in.ranks[, idx_int, drop = FALSE])
  mP1  <- rowMeans(in.ranks[, idx_p1,  drop = FALSE])
  mP2  <- rowMeans(in.ranks[, idx_p2,  drop = FALSE])

  SE_const <- N * (N + 1) / 12
  SEI1 <- sqrt(SE_const * (1/NI + 1/N1))
  SEI2 <- sqrt(SE_const * (1/NI + 1/N2))
  SE12 <- sqrt(SE_const * (1/N1  + 1/N2))

  zP1_Int <- (mP1 - mInt) / SEI1 / tiesStat2
  zP2_Int <- (mP2 - mInt) / SEI2 / tiesStat2
  zP2_P1  <- (mP2 - mP1)  / SE12  / tiesStat2

  zVals <- cbind(zP1_Int, zP2_Int, zP2_P1)
  pvals <- cbind(
    pmin(1, 2 * pnorm(zP1_Int, lower.tail = TRUE)),
    pmin(1, 2 * pnorm(zP2_Int, lower.tail = TRUE)),
    2  * pnorm(abs(zP2_P1), lower.tail = FALSE)
  )

  colnames(zVals) <- c("zP1_Int", "zP2_Int", "zP2_P1")
  colnames(pvals) <- c("pval_1_Int", "pval_2_Int", "pval_2_1")
  cbind(zVals, pvals)
}

# =============================================================================
# .find_genes_of_interest
# =============================================================================
#' .find_genes_of_interest
#' Identify genes associated with pattern interaction.
#' This function identifies genes exhibiting significantly higher values of
#' testMat in the Interaction region of the two patterns compared to regions
#' with exclusive influence from either pattern. It uses a Kruskal-Wallis
#' test followed by posthoc analysis using Dunn's Test to identify the genes.
#'
#' @usage
#' .find_genes_of_interest(testMat, goodGenes, region, fdr.level=0.05,
#'        analysis=c("enrichment","overlap"),...)
#' @param    testMat A matrix of counts with genes as rows and spots as columns.
#'           Sparse (dgCMatrix) or dense. For residual mode, the reconstruction
#'           has already been subtracted before this function is called.
#' @param    goodGenes A vector of user specified genes expected to interact
#'           a priori. The default is NULL.
#' @param    region A character/factor vector of region labels for each spot.
#'           Values are the two pattern names, "Interacting", or NA.
#' @param    fdr.level False Discovery Rate threshold. Default 0.05.
#' @param    analysis "enrichment" (default) or "overlap".
#' @param    ... Additional arguments passed to lower level functions.
#' @return A list containing one data frame of per-gene statistics and
#'   SpaceMarkersMetric, sorted by the metric.
#'

.find_genes_of_interest <- function(testMat, goodGenes = NULL, region,
                                    fdr.level = 0.05,
                                    analysis  = c("enrichment", "overlap"),
                                    ...) {

  analysis <- match.arg(analysis)
  region   <- factor(region)

  patnames <- levels(region)[levels(region) != "Interacting"]
  if (length(patnames) < 2)
    stop("Region factor must have at least 'Interacting' and two distinct ",
         "pattern levels.")
  pattern1 <- patnames[1]
  pattern2 <- patnames[2]

  # ── Gene filter: sparse-aware, no dense conversion yet ────────────────────
  if (!is.null(goodGenes)) {
    subset_goodGenes <- intersect(rownames(testMat), goodGenes)
  } else {
    # rowSums2 dispatches efficiently on dgCMatrix
    subset_goodGenes <- names(which(sparseMatrixStats::rowSums2(testMat) > 0))
  }

  # ── Convert only the pair-relevant slice to dense ─────────────────────────
  # Dimensions: expressed_genes x pair_spots (spots with non-NA region only).
  # This is the minimum materializaton needed for exact rank-based tests.
  # The full G x S matrix is never converted.
  keep_cells   <- !is.na(region)
  testMat_pair <- as.matrix(
    testMat[subset_goodGenes, keep_cells, drop = FALSE]
  )
  region_pair  <- droplevels(region[keep_cells])

  # ── Exact Kruskal-Wallis via matrixTests ───────────────────────────────────
  # Produces identical results to the original implementation and exactly
  # matches kruskal.test() output including the tie-correction denominator.
  res_kruskal <- matrixTests::row_kruskalwallis(x = testMat_pair,
                                                g = region_pair)

  qq <- qvalue::qvalue(res_kruskal$pvalue, fdr.level = fdr.level,
                       pfdr = FALSE, pi0 = 1)
  res_kruskal <- cbind(res_kruskal, p.adj = qq$qvalues)

  # Force zero-variance genes to p = 1
  zero_genes <- names(which(rowSums(testMat_pair) == 0))
  res_kruskal[zero_genes, c("df", "statistic")] <- 0
  res_kruskal[zero_genes, c("pvalue", "p.adj")]  <- 1

  # ── Exact Dunn's test on the same dense slice ─────────────────────────────
  # region_pair has no NAs so keep_cols inside .row_dunn_test is all-TRUE;
  # the subsetting is a no-op, preserving interface compatibility.
  res_dunn_test <- .row_dunn_test(
    in.data  = testMat_pair,
    region   = region_pair,
    pattern1 = pattern1,
    pattern2 = pattern2
  )
  rownames(res_dunn_test) <- rownames(res_kruskal)

  # ── FDR correction for Dunn ───────────────────────────────────────────────
  qDunn <- qvalue::qvalue(res_dunn_test[, 4:6], fdr.level = fdr.level,
                          pfdr = FALSE, pi0 = 1)

  ind <- rownames(res_kruskal[res_kruskal$p.adj < fdr.level, ])

  if (length(ind) > 0) {
    if (length(ind) == 1) {
      qDunn$qvalues[ind, ] <- res_dunn_test[ind, 4:6]
    } else {
      qq2 <- qvalue::qvalue(res_dunn_test[ind, 4:6],
                            fdr.level = fdr.level, pfdr = FALSE, pi0 = 1)
      qDunn$qvalues[ind, ] <- qq2$qvalues
    }
  }

  res_dunn_test <- cbind(res_dunn_test, qDunn$qvalues)
  colnames(res_dunn_test)[7:9] <- paste0(colnames(res_dunn_test)[7:9], ".adj")

  res_dunn_test[zero_genes, 1:3] <- 0
  res_dunn_test[zero_genes, 4:9] <- 1

  .build_interact_genes_df(res_kruskal, res_dunn_test, ind,
                           fdr.level, pattern1, pattern2, analysis)
}

# =============================================================================
# .build_interact_genes_df  — unchanged
# =============================================================================
.build_interact_genes_df <- function(res_kruskal, res_dunn_test, ind,
                                     fdr.level = 0.05, pattern1, pattern2,
                                     analysis) {
    interact_patt1 <- res_dunn_test[ind, "pval_1_Int.adj"] < fdr.level
    interact_patt2 <- res_dunn_test[ind, "pval_2_Int.adj"] < fdr.level
    interacting_over_both_patterns <- interact_patt1 & interact_patt2
    not_pattern1_diff_pattern2 <- res_dunn_test[ind, "pval_2_1.adj"] >= fdr.level
    exc_interact_patt1 <- interact_patt1 & not_pattern1_diff_pattern2
    exc_interact_patt2 <- interact_patt2 & not_pattern1_diff_pattern2
    names(exc_interact_patt2)              <- ind
    names(exc_interact_patt1)              <- names(exc_interact_patt2)
    names(interacting_over_both_patterns)  <- names(exc_interact_patt1)
    interact_genes <- matrix(
        FALSE, nrow = nrow(res_dunn_test), ncol = 2,
        dimnames = list(rownames(res_dunn_test),
                        c("Gene", paste0(pattern1, " x ", pattern2)))
    )
    interact_genes[, 1] <- rownames(interact_genes)
    interact_genes[ind[which(exc_interact_patt1)], 2] <- paste0("vs", pattern1)
    interact_genes[ind[which(exc_interact_patt2)], 2] <- paste0("vs", pattern2)
    interact_genes[ind[which(interacting_over_both_patterns)], 2] <- "vsBoth"
    colnames(res_kruskal)   <- paste0("KW.",   colnames(res_kruskal))
    colnames(res_dunn_test) <- paste0("Dunn.", colnames(res_dunn_test))
    interact_genes <- cbind(interact_genes, res_kruskal, res_dunn_test)
    if (analysis == "overlap")
        interact_genes <- interact_genes[interact_genes[, 2] != "FALSE", ]
    return(list(interact_genes))
}
