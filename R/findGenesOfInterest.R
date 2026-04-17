#' @importFrom sparseMatrixStats rowRanks rowSums2 rowMeans2 rowVars
#' @importFrom Matrix rowSums
#' @importFrom stats pchisq pnorm
#' @importFrom qvalue qvalue
## author: Atul Deshpande
## email: adeshpande@jhu.edu

.row_dunn_test <- function(in.data, region, pattern1, pattern2) {
  keep_cols <- !is.na(region)
  rsub <- region[keep_cols]
  
  # ranking on sparse subset
  in.ranks <- sparseMatrixStats::rowRanks(in.data[, keep_cols], ties.method = "average")
  
  N <- length(rsub)
  idx_int <- which(rsub == "Interacting")
  idx_p1  <- which(rsub == pattern1)
  idx_p2  <- which(rsub == pattern2)
  NI <- length(idx_int); N1 <- length(idx_p1); N2 <- length(idx_p2)
  
  # Tie correction using variance (Sparse)
  row_vars <- sparseMatrixStats::rowVars(in.ranks)
  tiesStat2 <- sqrt(row_vars / ((N * (N + 1)) / 12))
  tiesStat2[row_vars == 0] <- 1 # Avoid division by zero
  
  # Mean Ranks
  mInt <- sparseMatrixStats::rowMeans2(in.ranks[, idx_int, drop=FALSE])
  mP1  <- sparseMatrixStats::rowMeans2(in.ranks[, idx_p1, drop=FALSE])
  mP2  <- sparseMatrixStats::rowMeans2(in.ranks[, idx_p2, drop=FALSE])
  
  SE_const <- N * (N + 1) / 12
  SEI1 <- sqrt(SE_const * (1/NI + 1/N1))
  SEI2 <- sqrt(SE_const * (1/NI + 1/N2))
  SE12 <- sqrt(SE_const * (1/N1 + 1/N2))
  
  zP1_Int <- (mP1-mInt)/SEI1/tiesStat2
  zP2_Int <- (mP2-mInt)/SEI2/tiesStat2
  zP2_P1 <- (mP2-mP1)/SE12/tiesStat2
  
  zVals <- cbind(zP1_Int, zP2_Int, zP2_P1)
  pvals <- cbind(
    pmin(1, 2 * pnorm(zP1_Int, lower.tail = TRUE)),
    pmin(1, 2 * pnorm(zP2_Int, lower.tail = TRUE)),
    2 * pnorm(abs(zP2_P1), lower.tail = FALSE)
  )
  
  colnames(zVals) <- c("zP1_Int", "zP2_Int", "zP2_P1")
  colnames(pvals) <- c("pval_1_Int", "pval_2_Int", "pval_2_1")
  return(cbind(zVals, pvals))
}

#===================
#' .find_genes_of_interest
#' Identify genes associated with pattern interaction.
#' This function identifies genes exhibiting significantly higher values of 
#' testMat in the Interaction region of the two 
#' patterns compared to regions with exclusive influence from either 
#' pattern. It uses Kruskal-Wallis test followed by
#' posthoc analysis using Dunn's Test to identify the genes.
#'
#' @usage
#' .find_genes_of_interest(testMat, goodGenes, region, fdr.level=0.05,
#'        analysis=c("enrichment","overlap"),...)
#' @param    testMat A matrix of counts with cells as columns and genes as rows
#' @param    goodGenes A vector of user specified genes expected to interact 
#' a priori. The default for this is NULL as the function can find these genes 
#' itself
#' @param    region A data frame of the reference pattern regions that overlap 
#' with the other patterns
#' @param    fdr.level False Discovery Rate. The default value is 0.05.
#' @param    analysis a character string that specifies the type of analysis to 
#' carry out, whether overlap or enrichment.
#' @param    ... Additional arguments to be passed to lower level functions
#' @return a list of genes exhibiting significantly higher values of testMat in 
#' the Interaction region of the two #' patterns compared to regions with 
#' exclusive influence from either pattern.

.find_genes_of_interest <- function(testMat, goodGenes=NULL, region, fdr.level=0.05,
                                    analysis=c("enrichment", "overlap"), ...) {
  
  # If it's a dense Matrix class (like dgeMatrix from residuals), 
  # convert to standard base matrix so sparseMatrixStats can handle it.
  if (is(testMat, "dgeMatrix")) {
    testMat <- as.matrix(testMat)
  }
  # --------------------------------------
  
  # 1. Validation and Setup
  analysis <- match.arg(analysis)
  region <- factor(region)
  
  # Identify the two exclusive patterns (excluding the Interaction level)
  patnames <- levels(region)[which(levels(region) != "Interacting")]
  if (length(patnames) < 2) {
    stop("Region factor must have at least 'Interacting' and two distinct pattern levels.")
  }
  pattern1 <- patnames[1]
  pattern2 <- patnames[2]
  
  # 2. Gene Filtering (Using Sparse-aware rowSums)
  if (!is.null(goodGenes)) {
    subset_goodGenes <- intersect(rownames(testMat), goodGenes)
  } else {
    # rowSums2 handles dgCMatrix efficiently
    subset_goodGenes <- names(which(sparseMatrixStats::rowSums2(testMat) > 0))
  }
  testMat <- testMat[subset_goodGenes, ]
  
  # 3. Sparse Manual Kruskal-Wallis with Tie Correction
  # Filter to non-NA spots for the test
  keep_cells <- !is.na(region)
  r_sub <- region[keep_cells]
  
  # Calculate Ranks across the subset of cells
  # Note: returns a dense matrix of ranks, but only for active genes/cells
  ranks <- sparseMatrixStats::rowRanks(testMat[, keep_cells], ties.method = "average")
  
  N <- ncol(ranks)
  global_mean_rank <- (N + 1) / 2
  
  # Calculate Tie Correction (C) based on rank variance
  # This ensures exact match with matrixTests::row_kruskalwallis
  row_vars <- sparseMatrixStats::rowVars(ranks)
  tie_correction <- row_vars / ((N * (N + 1)) / 12)
  
  # Calculate H-Statistic sum across groups
  H_sum <- numeric(nrow(testMat))
  for (lvl in levels(r_sub)) {
    idx <- which(r_sub == lvl)
    n_i <- length(idx)
    m_i <- sparseMatrixStats::rowMeans2(ranks[, idx, drop = FALSE])
    H_sum <- H_sum + n_i * (m_i - global_mean_rank)^2
  }
  
  # Final H statistic (Adjusted for Ties)
  h_stat <- ((12 / (N * (N + 1))) * H_sum) / tie_correction
  h_stat[tie_correction == 0] <- 0 # Handle zero-variance genes
  
  p_vals_kw <- pchisq(h_stat, df = length(levels(r_sub)) - 1, lower.tail = FALSE)
  
  res_kruskal <- data.frame(
    statistic = h_stat,
    pvalue = p_vals_kw,
    df = length(levels(r_sub)) - 1,
    row.names = rownames(testMat)
  )
  
  # 4. Correct for FDR (Kruskal)
  qq <- qvalue::qvalue(res_kruskal$pvalue, fdr.level = fdr.level, pfdr = FALSE, pi0 = 1)
  res_kruskal$p.adj <- qq$qvalues
  
  # Handle zero variance/expression genes explicitly
  zero_genes <- names(which(sparseMatrixStats::rowSums2(testMat[, keep_cells]) == 0))
  res_kruskal[zero_genes, c("df", "statistic")] <- 0
  res_kruskal[zero_genes, c("pvalue", "p.adj")] <- 1
  
  # 5. Sparse Post-hoc Dunn's Test
  res_dunn_test <- .row_dunn_test(in.data = testMat, region = region,
                                  pattern1 = pattern1, pattern2 = pattern2)
  rownames(res_dunn_test) <- rownames(res_kruskal)
  
  # 6. Adjust Dunn P-values
  qDunn <- qvalue::qvalue(res_dunn_test[, 4:6], fdr.level = fdr.level, pfdr = FALSE, pi0 = 1)
  
  # Only readjust for genes that passed Kruskal threshold
  ind <- rownames(res_kruskal[which(res_kruskal$p.adj < fdr.level), ])
  if (length(ind) > 0) {
    if (length(ind) == 1) {
      qDunn$qvalues[ind, ] <- res_dunn_test[ind, 4:6]
    } else {
      qq_dunn_sub <- qvalue::qvalue(res_dunn_test[ind, 4:6], fdr.level = fdr.level, pfdr = FALSE, pi0 = 1)
      qDunn$qvalues[ind, ] <- qq_dunn_sub$qvalues
    }
  }
  
  res_dunn_test <- cbind(res_dunn_test, qDunn$qvalues)
  colnames(res_dunn_test)[7:9] <- paste0(colnames(res_dunn_test)[7:9], ".adj")
  
  # Force zero gene results to 1
  res_dunn_test[zero_genes, 1:3] <- 0
  res_dunn_test[zero_genes, 4:9] <- 1
  
  # 7. Final Output Construction
  interactGenes <- .build_interact_genes_df(res_kruskal, res_dunn_test, ind,
                                            fdr.level, pattern1, pattern2,
                                            analysis)
  return(interactGenes)
}

.build_interact_genes_df <- function(res_kruskal,res_dunn_test,ind,
                                fdr.level=0.05,pattern1,pattern2,analysis) {
    interact_patt1 <- res_dunn_test[ind,"pval_1_Int.adj"]<fdr.level
    interact_patt2 <- res_dunn_test[ind,"pval_2_Int.adj"]<fdr.level
    interacting_over_both_patterns <- interact_patt1 & interact_patt2
    not_pattern1_diff_pattern2 <- res_dunn_test[ind,"pval_2_1.adj"]>=fdr.level
    exc_interact_patt1 <- interact_patt1 & not_pattern1_diff_pattern2
    exc_interact_patt2 <- interact_patt2 & not_pattern1_diff_pattern2
    names(exc_interact_patt2)<-ind
    names(exc_interact_patt1)<-names(exc_interact_patt2)
    names(interacting_over_both_patterns)<-names(exc_interact_patt1)
    interact_genes<-matrix(FALSE,nrow=nrow(res_dunn_test),ncol=2,dimnames=list(
        rownames(res_dunn_test),c("Gene",paste0(pattern1,' x ',pattern2))))
    interact_genes[,1] <- rownames(interact_genes)
    interact_genes[ind[which(exc_interact_patt1)],2]<-paste0("vs",pattern1)
    interact_genes[ind[which(exc_interact_patt2)],2]<-paste0("vs",pattern2)
    interact_genes[ind[which(interacting_over_both_patterns)],2]<-"vsBoth"
    colnames(res_kruskal) <- paste0("KW.",colnames(res_kruskal))
    colnames(res_dunn_test) <- paste0("Dunn.",colnames(res_dunn_test))
    interact_genes <- cbind(interact_genes,res_kruskal,res_dunn_test)
    if (analysis=="overlap"){
        interact_genes <-interact_genes[interact_genes[,2]!="FALSE",] 
    }
    return(list(interact_genes))
}
