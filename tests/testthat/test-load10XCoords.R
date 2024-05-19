
create_mocks <- function(header = TRUE) {
  dir.create("mock_visiumDir/spatial", recursive = TRUE)

  # Create mock scalefactors_json.json
  scalefactors_json <- list(
    spot_diameter_fullres = 1,
    lowres = 0.5
  )
  jsonlite::write_json(scalefactors_json,
                       "mock_visiumDir/spatial/scalefactors_json.json")

  # Create mock tissue_positions_list.csv
  tissue_positions_list <- data.frame(
    barcode = c("ACGCCTGACACGCGCT-1", "TACCGATCCAACACTT-1",
                "CGGAGGCTCTCGTCTG-1"),
    in_tissue = c(1, 1, 1),
    array_row = c(0, 1, 0),
    array_col = 0:2,
    pxl_row_in_fullres = c(1, 2, 3),
    pxl_col_in_fullres = c(4, 5, 6)
  )

  write.table(tissue_positions_list,
              "mock_visiumDir/spatial/tissue_positions_list.csv",
              row.names = FALSE,
              col.names = header,
              sep = ",")

  # Create expected output
  expected_output <- data.frame(
    barcode = c("ACGCCTGACACGCGCT-1", "TACCGATCCAACACTT-1",
                "CGGAGGCTCTCGTCTG-1"),
    y = c(0.5, 1, 1.5),
    x = c(2, 2.5, 3)
  )
  return(expected_output)
}

test_that("load10XCoords returns the expected output with header = FALSE", {
  expected_output <- create_mocks(header = FALSE)
  output <- load10XCoords("mock_visiumDir")
  expect_equal(output, expected_output)
  unlink("mock_visiumDir", recursive = TRUE)
})

test_that("load10XCoords returns the expected output with header = TRUE", {
  expected_output <- create_mocks(header = TRUE)
  output <- load10XCoords("mock_visiumDir")
  expect_equal(output, expected_output)
  unlink("mock_visiumDir", recursive = TRUE)
})