# 1. Test with no 'patternPairs' provided
test_that("No patternPairs - computes all pairs", {
    # Create sample data (adjust to match your actual data structure)
    args <- createSampleData()
    # Call the function
    
    suppressMessages(
        result <- getPairwiseInteractingGenes(data=args$data,
        spPatterns=args$spPatterns,
        optParams=args$optParams, minOverlap = 0)
    )
    patnames <- setdiff(colnames(args$spPatterns), c("x", "y", "barcode"))
    # Check the number of results (should be 1 since there are 2 patterns)
    expect_equal(length(result), choose(length(patnames), 2))
    patcombs <- t(utils::combn(patnames, 2))
    patcombs <- apply(patcombs, 1, function(x) paste(x, collapse = "_"))
    # Check the names of the results (should be "pattern1_pattern2")
    expect_equal(names(result), patcombs)
})

# 2. Test with a matrix of 'patternPairs'
test_that("Matrix patternPairs - specific pairs", {
    # Create sample data
    args <- createSampleData()
    # Define specific pattern pairs
    patternPairs <- rbind(c("pattern1", "pattern2"), c("pattern1", "pattern3"))
    
    # Call the function
    suppressMessages(
        result <- getPairwiseInteractingGenes(data=args$data, 
        spPatterns=args$spPatterns, 
        optParams=args$optParams, 
        patternPairs = patternPairs, minOverlap = 0)
    )
    # Check the number of results (should be 1)
    expect_equal(length(result), 2)
    
    # Check the names of the results
    expect_equal(names(result), c("pattern1_pattern2", "pattern1_pattern3"))
})

# 3. Test with a list of 'patternPairs'
test_that("List patternPairs - specific pairs", {
    # Create sample data
    args <- createSampleData()
    # Define specific pattern pairs as a list
    patternPairs <- list(c("pattern1", "pattern2"))
    
    # Call the function
    suppressMessages(
    result <- getPairwiseInteractingGenes(data=args$data, 
        spPatterns=args$spPatterns, 
        optParams=args$optParams, 
        patternPairs = patternPairs, minOverlap = 0)
    )
    # Check the number of results (should be 1)
    expect_equal(length(result), 1)
    
    # Check the names of the results
    expect_equal(names(result), "pattern1_pattern2")
})

# 4. Test with invalid 'patternPairs'
test_that("Invalid patternPairs - throws error", {
    # Create sample data
    args <- createSampleData()
    
    # Define invalid pattern pairs
    patternPairs <- matrix(c("pattern1", "invalid_pattern"), ncol = 2)
    
    # Expect an error
    expect_error(result <- getPairwiseInteractingGenes(data=args$data, 
        spPatterns=args$spPatterns, 
        optParams=args$optParams, 
        patternPairs = patternPairs, minOverlap = 0), "not pattern names")
})

# 5. Test without BiocParallel
test_that("No BiocParallel - sequential execution", {
    # Mock 'requireNamespace' to always return FALSE
    myRequire <- function(...) FALSE
    local_mocked_bindings(requireNamespace = myRequire, .package = "base")
        # Create sample data
        args <- createSampleData()
        # Define specific pattern pairs
        patternPairs <- rbind(c("pattern1", "pattern2"), 
                                    c("pattern1", "pattern3"))
        
        # Call the function
        suppressMessages(
            result <- getPairwiseInteractingGenes(data=args$data, 
                spPatterns=args$spPatterns, 
                optParams=args$optParams, 
                patternPairs = patternPairs, minOverlap = 0)
        )
        # Check the number of results (should be 1)
        expect_equal(length(result), 2)
        
        # Check the names of the results
        expect_equal(names(result), c("pattern1_pattern2", 
            "pattern1_pattern3"))
})

#6. Test for single pair of patterns
test_that("Single pair of patterns", {
    # Create sample data
    args <- createSampleData()
    # Define specific pattern pairs
    patternPairs <- c("pattern1", "pattern2")
    
    # Call the function
    suppressMessages(
        result <- getPairwiseInteractingGenes(data=args$data, 
            spPatterns=args$spPatterns, 
            optParams=args$optParams, 
            patternPairs = patternPairs, minOverlap = 0)    
    )
    
    # Check the number of results (should be 1)
    expect_equal(length(result), 1)
    
    # Check the names of the results
    expect_equal(names(result), "pattern1_pattern2")
})

test_that("getPairwiseInteractingGenes works with specified hotspots", {
    args <- createSampleData()
    spHotspots <- findAllHotspots(spPatterns = spPatterns, params = optParams,
        outlier = "positive", nullSamples = 1000, includeSelf = TRUE)

    suppressMessages({
        interactingGenes <- getPairwiseInteractingGenes(
            data=args$data,
            spPatterns=args$spPatterns,
            optParams=args$optParams,
            minOverlap = 0,
            hotspots = spHotspots)})
    expect_true(is.list(interactingGenes))

})