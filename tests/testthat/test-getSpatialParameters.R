test_that("getSpatialParameters returns optimal parameters", {
  # Create test data
  cells <- c()
  test_num <- 500
  for(i in 1:test_num){
    cells[length(cells)+1] <- paste0("cell_",i)
    }
  spatialPatterns <- data.frame(barcode = cells,
                                y = runif(test_num, min=0, max=test_num),
                                x = runif(test_num, min=0, max=500),
                                Pattern_1 = runif(test_num, min=0, max=1),
                                Pattern_2 = runif(test_num, min=0, max=1) )
  
  # Call the getSpatialParameters function with the test data
  optParams <- getSpatialParameters(spatialPatterns)
  
  # Perform assertions to check if the optimal parameters are returned correctly
  expect_equal(ncol(optParams), 2, info = "Incorrect number of pattern parameters")
  expect_equal(colnames(optParams), c("Pattern_1", "Pattern_2"), info = "Incorrect pattern names")

})
