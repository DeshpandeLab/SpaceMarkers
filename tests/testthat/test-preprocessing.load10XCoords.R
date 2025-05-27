#a mock function for the header/non header cases, can be extended
#to other cases (above) as well
create_visium_mock <- function(version = "1.0", probe_set = TRUE) {
  stopifnot(version %in% c("1.0", "2.0","HD"))
  if (version == "1.0") {
    header <- FALSE
    tissue_pos_name <- "tissue_positions_list.csv"
  } else if (version == "2.0") {
    header <- TRUE
    tissue_pos_name <- "tissue_positions.csv"
  } else if (version == "HD"){
    tissue_pos_name <- "tissue_positions.parquet"
  }
  
  dir.create("mock_visiumDir/spatial", recursive = TRUE)
  
  #Create probe_set file
  if (probe_set) {
    writeLines(c(sprintf("#probe_set_file_version=%s", version),
                 "#panel_name=TEST",
                 "gene_id",
                 "ENSMUSG00000000001",
                 "ENSMUSG00000000003",
                 "ENSMUSG00000000028"),
               "mock_visiumDir/probe_set.csv")
  }
  
  
  # Create mock scalefactors_json.json
  scalefactors_json <- list(
    spot_diameter_fullres = 1,
    lowres = 0.5
  )
  jsonlite::write_json(scalefactors_json,
                       "mock_visiumDir/spatial/scalefactors_json.json")
  
  # Create mock tissue_positions_list.csv
  tissue_positions <- data.frame(
    barcode = c("ACGCCTGACACGCGCT-1",
                "TACCGATCCAACACTT-1",
                "CGGAGGCTCTCGTCTG-1"),
    in_tissue = c(1, 1, 1),
    array_row = c(0, 1, 0),
    array_col = 0:2,
    pxl_row_in_fullres = c(1, 2, 3),
    pxl_col_in_fullres = c(4, 5, 6)
  )
  
  if (version == "HD"){
    nanoparquet::write_parquet(tissue_positions,
                               file = sprintf(
                                 "mock_visiumDir/spatial/%s.parquet",
                                 tissue_pos_name))
  } else {
    write.table(tissue_positions,
                sprintf("mock_visiumDir/spatial/%s.csv", tissue_pos_name),
                row.names = FALSE,
                col.names = header,
                sep = ",")
    
  }
  
  
  # Create expected output
  expected_output <- data.frame(
    barcode = c("ACGCCTGACACGCGCT-1",
                "TACCGATCCAACACTT-1",
                "CGGAGGCTCTCGTCTG-1"),
    y = c(0.5, 1, 1.5),
    x = c(2, 2.5, 3)
  )
  return(expected_output)
}

test_that("load10XCoords returns the expected output for spaceranger v1.0", {
  expected_output <- create_visium_mock(version = "1.0")
  output <- load10XCoords("mock_visiumDir", 
            resolution="lowres",
            version = "1.0")
  expect_equal(output, expected_output)
  unlink("mock_visiumDir", recursive = TRUE)
})

test_that("load10XCoords returns the expected output for spaceranger v2.0", {
  expected_output <- create_visium_mock(version = "2.0")
  output <- load10XCoords("mock_visiumDir", resolution="lowres",
  version = "2.0")
  expect_equal(output, expected_output)
  unlink("mock_visiumDir", recursive = TRUE)
})

test_that("load10XCoords returns the expected output for visium HD", {
  expected_output <- create_visium_mock(version = "HD")
  output <- load10XCoords("mock_visiumDir", resolution="lowres",
                          version = "HD")
  expect_equal(as.data.frame(output), expected_output)
  unlink("mock_visiumDir", recursive = TRUE)
})

test_that("load10XCoords returns the expected output without probe_set file", {
  expected_output <- create_visium_mock(version = "1.0", probe_set = FALSE)
  output <- load10XCoords("mock_visiumDir", resolution="lowres")
  expect_equal(output, expected_output)
  unlink("mock_visiumDir", recursive = TRUE)
})
