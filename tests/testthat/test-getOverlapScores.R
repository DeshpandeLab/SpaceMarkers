test_that("getOverlapScores works correctly", {
    # Create a sample hotspots data frame
    hotspots <- data.frame(
        x = c(1, 2, 3, 4, 5),
        y = c(1, 2, 3, 4, 5),
        barcode = c("A", "B", "C", "D", "E"),
        pattern1 = c("pattern1", NA, "pattern1", NA, "pattern1"),
        pattern2 = c("pattern2", "pattern2", NA, "pattern2", "pattern2")
    )
    
    # Test default method (Szymkiewicz–Simpson)
    result <- getOverlapScores(hotspots)
    expect_equal(nrow(result), 1)
    expect_equal(as.character(result$pattern1), "pattern1")
    expect_equal(as.character(result$pattern2), "pattern2")
    expect_equal(result$overlapScore, 0.6666667,
                 tolerance = 1e-7) 
    # Test Jaccard method
    result_jaccard <- getOverlapScores(hotspots, method = "Jaccard")
    expect_equal(nrow(result_jaccard), 1)
    expect_equal(result_jaccard$overlapScore, 0.4)
    # Test Sørensen–Dice method
    result_sorensen <- getOverlapScores(hotspots, method = "Sorensen–Dice")
    expect_equal(nrow(result_sorensen), 1)
    expect_equal(result_sorensen$overlapScore, 0.5714286,
                 tolerance = 1e-7)
    # Test Ochiai method
    result_ochiai <- getOverlapScores(hotspots, method = "Ochiai")
    expect_equal(nrow(result_ochiai), 1)
    expect_equal(result_ochiai$overlapScore, 0.5773503,
                 tolerance = 1e-7)
    # Test absolute method
    result_absolute <- getOverlapScores(hotspots, method = "absolute")
    expect_equal(nrow(result_absolute), 1)
    expect_equal(result_absolute$overlapScore, 2)

    # Test with patternList
    result_patternList <- getOverlapScores(hotspots, patternList = c("pattern1", "pattern2"))
    expect_equal(nrow(result_patternList), 1)
    expect_equal(as.character(result_patternList$pattern1), "pattern1")
    expect_equal(as.character(result_patternList$pattern2), "pattern2")
    expect_equal(result_patternList$overlapScore, 0.6666667,
                 tolerance = 1e-7)

    # Test invalid pattern names
    expect_error(getOverlapScores(hotspots, patternList = c("invalidPattern")))
    
    # Test with multiple patterns
    hotspots$pattern3 <- c(NA, "pattern3", "pattern3", NA, NA)
    result_multiple <- getOverlapScores(hotspots, patternList = c("pattern1", "pattern2", "pattern3"))
    expect_true(nrow(result_multiple) > 1)

    # Test with multiple methods
    expect_warning(getOverlapScores(hotspots, method = c("Jaccard", "absolute")))

})