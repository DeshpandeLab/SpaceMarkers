
# Sample data to use in tests
spatialPatterns <- data.frame(
  barcode = c("A", "B", "C"),
  x = c(1, 2, 3),
  y = c(4, 5, 6),
  pattern1 = c(0.5, 0.6, 0.7),
  pattern2 = c(0.8, 0.9, 1.0)
)
temp <- "mock_visiumDir"
# Mock JSON data for testing
mock_json <- function(temp_dir = temp) {
  spatial_dir <- file.path(temp_dir, "spatial")
  
  dir.create(spatial_dir, showWarnings = FALSE, recursive = TRUE)
  
  json_data <- list(
    scale_factors = 1,
    tissue_hires_scalef = 0.75,
    tissue_lowres_scalef = 0.25,
    spot_diameter_fullres = 4
  )
  jsonlite::write_json(json_data,file.path(spatial_dir,
                                           "scalefactors_json.json"), 
                       simplifyVector =TRUE)
  return(list(visiumDir = temp_dir, spatialDir = spatial_dir,
              pattern = "scalefactors_json.json"))
}

# Test cases for getSpatialParamsExternal
test_that("getSpatialParameters works with provided sigma", {
  result <- getSpatialParameters(spatialPatterns, sigma = 5)
  
  expect_equal(nrow(result), 2)
  expect_equal(ncol(result), 2)
  expect_equal(result[["sigmaOpt", "pattern1"]], 5)
  expect_equal(result[["threshOpt", "pattern1"]], 4)
  expect_equal(result[["sigmaOpt", "pattern2"]], 5)
  expect_equal(result[["threshOpt", "pattern2"]], 4)
})

test_that("getSpatialParameters works by reading from JSON file", {
  paths <- mock_json()
  
result <- getSpatialParameters(spatialPatterns,
                                         visiumDir = paths$visiumDir,
                                         spatialDir = "spatial")
  
  expect_equal(nrow(result), 2)
  expect_equal(ncol(result), 2)
  expect_equal(result[["sigmaOpt", "pattern1"]], 4)
  expect_equal(result[["threshOpt", "pattern1"]], 4)
  expect_equal(result[["sigmaOpt", "pattern2"]], 4)
  expect_equal(result[["threshOpt", "pattern2"]], 4)
})

test_that("getSpatialParameters works with threshold", {
  result <- getSpatialParameters(spatialPatterns, sigma = 6,
                                         threshold = 10)
  
  expect_equal(result[["sigmaOpt", "pattern1"]], 6)
  expect_equal(result[["threshOpt", "pattern1"]], 10)
  expect_equal(result[["sigmaOpt", "pattern2"]], 6)
  expect_equal(result[["threshOpt", "pattern2"]], 10)
})
unlink(temp, recursive = TRUE)
