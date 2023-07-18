test_that("getSpatialParameters returns optimal parameters", {
  # Create test data
  cells <- c()
  test_num <- 500
  for(i in 1:test_num){
    cells[length(cells)+1] <- paste0("cell_",i)
  }
  spatialPatterns <- readRDS(test_path("testdata", "spatialPatterns.rds"))
  # Call the getSpatialParameters function with the test data
  optParams <- getSpatialParameters(spatialPatterns)
  optParams_test <- readRDS(test_path("testdata", "optParams_test.rds"))
  # Perform assertions to check if the optimal parameters are returned correctly
  expect_equal(optParams["sigmaOpt",] ,optParams_test["sigmaOpt",], tolerance = 3,
               info = "Parameters are outside the reasonable range")
  expect_equal(optParams["threshOpt",] ,optParams_test["threshOpt",], tolerance = 3,
               info = "Parameters are outside the reasonable range")
  expect_equal(ncol(optParams), 2, info = "Incorrect number of pattern parameters")
  expect_equal(colnames(optParams), c("Pattern_1", "Pattern_2"), info = "Incorrect pattern names")

})
