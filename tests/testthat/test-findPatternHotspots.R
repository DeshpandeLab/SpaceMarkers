# Sample data for testing
spPatterns <- data.frame(
  x = runif(100, 0, 10),
  y = runif(100, 0, 10),
  Pattern_1 = rnorm(100)
)

# Test with includeSelf = FALSE
test_that("findPatternHotspots works with includeSelf = FALSE", {
  result <- find_pattern_hotspots(spPatterns, patternName = "Pattern_1",
                                includeSelf = FALSE)
  expect_type(result, "character")
})

# Test behavior with different outlier types
test_that("findPatternHotspots handles outlier parameter correctly", {
  result_positive <- find_pattern_hotspots(spPatterns, patternName = "Pattern_1",
                                         outlier = "positive")
  result_two_sided <- find_pattern_hotspots(spPatterns, patternName = "Pattern_1",
                                          outlier = "two.sided")
  
  expect_true(all(is.na(result_positive) | result_positive == "Pattern_1"))
  expect_true(all(is.na(result_two_sided) | result_two_sided == "Pattern_1"))
})