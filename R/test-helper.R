# helper file for tests

# create sample data for getPairwiseInteractingGenes
createSampleData <- function(...) {
    data <- matrix(stats::rnorm(100), nrow = 10)
    rownames(data) <- paste0("gene_", seq_len(10))
    colnames(data) <- paste0("spot_", seq_len(10))
    spPatterns <- as.data.frame(matrix(stats::rnorm(40), nrow = 10))
    optParams <- data.matrix(data.frame("pattern1" = c(.2,3),
                                           "pattern2" = c(.2,3), 
                                            "pattern3" = c(.2,3), 
                                            "pattern4" = c(.2,3)))
    rownames(optParams) <- c("sigmaOpt","threshOpt")
    patnames <- colnames(spPatterns) <- 
        c("pattern1", "pattern2", "pattern3", "pattern4")
    coords <- data.frame(x = stats::runif(10), y = stats::runif(10), 
        barcode = paste0("spot_",seq_len(10)))
    spPatterns <- cbind(coords,spPatterns)
    return(list(data = data, spPatterns = spPatterns, 
                patnames = patnames, optParams = optParams))
}
