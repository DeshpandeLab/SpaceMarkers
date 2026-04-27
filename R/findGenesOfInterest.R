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
# Receives a dense matrix and a region vector. region may contain NAs;
# keep_cols subsetting is retained for interface compatibility. When called
# from .find_genes_of_interest the slice is already pair-subsetted so
# keep_cols is all-TRUE (no-op).
#
# Tie correction: exact formula from the original paper implementation,
# written verbatim: sum(table(rr)^3 - table(rr)).

.row_dunn_test <- function(in.data, region, pattern1, pattern2) {

  keep_cols <- !is.na(region)
  rsub      <- region[keep_cols]

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

  # Exact tie correction — verbatim from original
  tiesStat  <- apply(in.ranks, 1, function(rr) sum(table(rr)^3 - table(rr)))
  tiesStat2 <- sqrt(1 - tiesStat / N / (N - 1) / (N + 1))
  tiesStat2[tiesStat2 == 0] <- 1   # guard: constant genes → z = 0

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

  # P-values — verbatim ifelse logic from original
  pval_1_Int <- ifelse(
    2 * pnorm(zP1_Int, lower.tail = TRUE) <= 1,
    2 * pnorm(zP1_Int, lower.tail = TRUE), 1
  )
  pval_2_Int <- ifelse(
    2 * pnorm(zP2_Int, lower.tail = TRUE) <= 1,
    2 * pnorm(zP2_Int, lower.tail = TRUE), 1
  )
  pval_2_1 <- ifelse(
    zP2_P1 > 0,
    2 * pnorm(zP2_P1, lower.tail = FALSE),
    2 * pnorm(zP2_P1, lower.tail = TRUE)
  )

  pvals <- cbind(pval_1_Int, pval_2_Int, pval_2_1)
  colnames(zVals) <- c("zP1_Int", "zP2_Int", "zP2_P1")
  colnames(pvals) <- c("pval_1_Int", "pval_2_Int", "pval_2_1")
  cbind(zVals, pvals)
}

# =============================================================================
# .find_genes_of_interest
# =============================================================================
#' .find_genes_of_interest
#' Identify genes associated with pattern interaction.
#'
#' @param testMat  genes x spots matrix (sparse dgCMatrix or dense).
#' @param goodGenes optional vector of genes to restrict analysis to.
#' @param region   region label vector per spot ("Interacting", p1, p2, NA).
#' @param fdr.level FDR threshold (default 0.05).
#' @param analysis "enrichment" (default) or "overlap".
#' @param ...      additional arguments.
#'
#' Bug fixed vs previous sparse implementation
#' --------------------------------------------
#' Bug 1 (KW qvalue denominator) — the only confirmed bug:
#'   The previous version pre-filtered to expressed genes before running
#'   qvalue on KW p-values. With pi0 = 1, qvalue is BH-equivalent; the
#'   threshold for rank-k gene among m tests is p <= (k/m)*alpha. Halving m
#'   by dropping zero genes doubles the effective threshold, producing
#'   too many genes passing KW FDR. The original passes ALL genes to
#'   row_kruskalwallis; zero-expression genes produce p = 1, contributing
#'   to a properly calibrated denominator.
#'   Fix: convert the full pair slice (all genes x pair spots) to dense
#'   before KW, with no gene pre-filtering.
#'
#' Mode-specific Dunn behaviour (by design, not bugs)
#' ---------------------------------------------------
#' enrichment: Dunn on ALL expressed genes, two-step qvalue, .adj
#'   classification. Returns all genes sorted by SpaceMarkersMetric —
#'   needed for GSEA-style ranking of the full gene list.
#'
#' overlap: Dunn on ind genes ONLY (matching original paper), single
#'   qvalue on ind p-values, .adj classification. Returns only non-FALSE
#'   classified genes — matches original paper output counts.

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

  # goodGenes filter — only when explicitly provided (matching original)
  if (!is.null(goodGenes)) {
    testMat <- testMat[intersect(rownames(testMat), goodGenes), , drop = FALSE]
  }

  # ── Bug 1 fix: materialize ALL genes x pair spots ─────────────────────────
  # No expressed-gene pre-filtering. Zero genes produce KW p = 1, which
  # calibrates the qvalue denominator correctly — matching the original's
  # row_kruskalwallis(x = as.matrix(testMat), g = region) call where
  # matrixTests excludes NA-group observations but uses all genes.
  keep_cells   <- !is.na(region)
  testMat_pair <- as.matrix(testMat[, keep_cells, drop = FALSE])
  region_pair  <- droplevels(region[keep_cells])

  # ── Exact KW on all genes ──────────────────────────────────────────────────
  res_kruskal <- matrixTests::row_kruskalwallis(x = testMat_pair,
                                                g = region_pair)
  qq_kw <- qvalue::qvalue(res_kruskal$pvalue, fdr.level = fdr.level,
                          pfdr = FALSE, pi0 = 1)
  res_kruskal <- cbind(res_kruskal, p.adj = qq_kw$qvalues)

  ind <- rownames(res_kruskal[res_kruskal$p.adj < fdr.level, ])

  if (length(ind) == 0) {
    empty <- matrix(
      FALSE, nrow = 0L, ncol = 2L,
      dimnames = list(NULL, c("Gene", paste0(pattern1, " x ", pattern2)))
    )
    return(list(empty))
  }

  # ── Mode-specific Dunn path ────────────────────────────────────────────────

  if (analysis == "enrichment") {

    # Dunn on all expressed genes — enables full-gene SpaceMarkersMetric
    # ranking needed for GSEA. Two-step qvalue: global correction first,
    # then re-corrected within ind for genes that passed KW.
    expressed <- names(which(sparseMatrixStats::rowSums2(
      testMat[, keep_cells, drop = FALSE]) > 0))

    res_dunn_test <- .row_dunn_test(
      in.data  = testMat_pair[expressed, , drop = FALSE],
      region   = region_pair,
      pattern1 = pattern1,
      pattern2 = pattern2
    )
    rownames(res_dunn_test) <- expressed

    qDunn <- qvalue::qvalue(res_dunn_test[, 4:6], fdr.level = fdr.level,
                            pfdr = FALSE, pi0 = 1)

    # Re-adjust within ind (expressed subset)
    ind_expressed <- intersect(ind, expressed)
    if (length(ind_expressed) > 1) {
      qq2 <- qvalue::qvalue(res_dunn_test[ind_expressed, 4:6],
                            fdr.level = fdr.level, pfdr = FALSE, pi0 = 1)
      qDunn$qvalues[ind_expressed, ] <- qq2$qvalues
    } else if (length(ind_expressed) == 1L) {
      qDunn$qvalues[ind_expressed, ] <- res_dunn_test[ind_expressed, 4:6]
    }

    res_dunn_test <- cbind(res_dunn_test, qDunn$qvalues)
    colnames(res_dunn_test)[7:9] <- paste0(colnames(res_dunn_test)[7:9], ".adj")

    # ind for classification is restricted to expressed genes
    .build_interact_genes_df(res_kruskal, res_dunn_test, ind_expressed,
                             fdr.level, pattern1, pattern2, analysis)

  } else {

    # overlap: Dunn on ind only — matches original paper behaviour.
    # Single qvalue on ind p-values only.
    res_dunn_test <- .row_dunn_test(
      in.data  = testMat_pair[ind, , drop = FALSE],
      region   = region_pair,
      pattern1 = pattern1,
      pattern2 = pattern2
    )
    rownames(res_dunn_test) <- ind

    qq_dunn <- qvalue::qvalue(res_dunn_test[, 4:6], fdr.level = fdr.level,
                              pfdr = FALSE, pi0 = 1)
    res_dunn_test <- cbind(res_dunn_test, qq_dunn$qvalues)
    colnames(res_dunn_test)[7:9] <- paste0(colnames(res_dunn_test)[7:9], ".adj")

    .build_interact_genes_df(res_kruskal, res_dunn_test, ind,
                             fdr.level, pattern1, pattern2, analysis)
  }
}

# =============================================================================
# .build_interact_genes_df
# =============================================================================
# Used by both modes. Key design:
#   - res_dunn_test rows define the gene set included in output.
#     enrichment: all expressed genes; overlap: ind genes only.
#   - ind defines which rows receive a non-FALSE classification label.
#     enrichment: ind_expressed; overlap: ind (same as res_dunn_test rows).
#   - res_kruskal is subsetted to rownames(res_dunn_test) for the cbind,
#     ensuring row counts always match.
#   - Classification uses .adj Dunn p-values (user-confirmed valid choice).
#   - overlap mode filters FALSE-classified rows at the end.

.build_interact_genes_df <- function(res_kruskal, res_dunn_test, ind,
                                     fdr.level = 0.05, pattern1, pattern2,
                                     analysis) {

  genes_for_output <- rownames(res_dunn_test)

  # Classification using .adj p-values
  interact_patt1 <- res_dunn_test[ind, "pval_1_Int.adj"] < fdr.level
  interact_patt2 <- res_dunn_test[ind, "pval_2_Int.adj"] < fdr.level
  interacting_over_both_patterns <- interact_patt1 & interact_patt2
  not_pattern1_diff_pattern2     <- res_dunn_test[ind, "pval_2_1.adj"] >= fdr.level
  exc_interact_patt1 <- interact_patt1 & not_pattern1_diff_pattern2
  exc_interact_patt2 <- interact_patt2 & not_pattern1_diff_pattern2

  names(exc_interact_patt2)             <- ind
  names(exc_interact_patt1)             <- names(exc_interact_patt2)
  names(interacting_over_both_patterns) <- names(exc_interact_patt1)

  interact_genes <- matrix(
    FALSE, nrow = length(genes_for_output), ncol = 2,
    dimnames = list(genes_for_output,
                    c("Gene", paste0(pattern1, " x ", pattern2)))
  )
  interact_genes[, 1] <- rownames(interact_genes)
  interact_genes[ind[which(exc_interact_patt1)], 2] <- paste0("vs", pattern1)
  interact_genes[ind[which(exc_interact_patt2)], 2] <- paste0("vs", pattern2)
  interact_genes[ind[which(interacting_over_both_patterns)], 2] <- "vsBoth"

  colnames(res_kruskal)   <- paste0("KW.",   colnames(res_kruskal))
  colnames(res_dunn_test) <- paste0("Dunn.", colnames(res_dunn_test))

  # res_kruskal subsetted to genes_for_output so row counts match
  interact_genes <- cbind(
    interact_genes,
    res_kruskal[genes_for_output, ],
    res_dunn_test
  )

  if (analysis == "overlap")
    interact_genes <- interact_genes[interact_genes[, 2] != "FALSE", ]

  return(list(interact_genes))
}
