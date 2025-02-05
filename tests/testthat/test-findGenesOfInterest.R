test_that("find_genes_of_interest handles empty significant genes without error", {
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
    expect_no_error(find_genes_of_interest(testMat, region = region, fdr.level = fdr))
})