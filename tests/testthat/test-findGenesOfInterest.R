test_that("findGenesOfInterest handles no significant genes without error", {
    # Mock input data
    fdr <- 0.05
    testMat <- matrix(rnorm(100), nrow=10, ncol=10)
    rownames(testMat) <- paste0("Gene", 1:10)
    region <- factor(rep(c("Pattern1", "Pattern2", "Interacting"), length.out=10))
    
    # Mock result of row_kruskalwallis to have no significant genes
    res_kruskal <- data.frame(
        pvalue = runif(10, 0.5, 1),  # p-values all above typical significance threshold
        row.names = rownames(testMat)
    )

    # Mock qvalue adjustment to have no significant adjusted p-values
    qq <- qvalue::qvalue(res_kruskal$pvalue, fdr.level = fdr, pfdr = FALSE, pi0 = 1)
    res_kruskal$p.adj <- qq$qvalues

    # Call the function and expect no errors
    expect_no_error(findGenesOfInterest(testMat, region = region, fdr.level = fdr))
})

test_that("findGenesOfInterest returns statistics for all non-zero genes", {
    # Mock input data
    fdr <- 0.05
    testMat <- matrix(runif(80), nrow=8, ncol=20)
    testMat <- rbind(testMat, matrix(0, nrow=4, ncol=20))  # Add 4 more rows
    testMat[10, 6] <- .2; # Add 3 non-zero genes
    testMat[11, 2] <- .4;
    testMat[12,9] <- .3;
    rownames(testMat) <- paste0("Gene", 1:12)
    region <- factor(rep(c("Pattern1", "Pattern2", "Interacting"), length.out=20))
    
    region[c(2,6,9)] <- NA  # Set some spots to NA
    region <- factor(region)
    suppressWarnings(
        sm_test <- findGenesOfInterest(testMat, region = region, fdr.level = fdr)
    )
    zero_genes <- c("Gene10", "Gene11", "Gene12")
    coi <- c("KW.df","KW.statistic","KW.pvalue","KW.p.adj","Dunn.zP1_Int",
                "Dunn.zP2_Int","Dunn.zP2_P1","Dunn.pval_1_Int",
                "Dunn.pval_2_Int","Dunn.pval_2_1","Dunn.pval_1_Int.adj",
                "Dunn.pval_2_Int.adj","Dunn.pval_2_1.adj")
    # Check that all non-zero genes are included in the output
    expect_equal(nrow(sm_test[[1]]), 11)
    expect_equal(any(is.na(sm_test[[1]])), FALSE)
    expect_equal(sum(abs(t(as.matrix(sm_test[[1]][zero_genes,coi])) 
        - c(0,0,1,1,0,0,0,1,1,1,1,1,1))), 0)
})
