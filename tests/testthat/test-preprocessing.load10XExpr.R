
test_that("load10XExpr loads data correctly", {
  # Create a temporary directory for testing
  temp_dir <- "."
  # Generate a test h5 file
  file.create("test.h5")
  h5_file <- file.path("test.h5")
  hf <- hdf5r::H5File$new(h5_file , mode = "w")

  # Create test data
  counts <- c(1, 2, 3, 4, 5, 6)
  indices <- c(0, 1, 2, 0, 1, 2)
  indptr <- c(0, 3, 6)
  shp <- c(3, 2)
  features <- c("gene1", "gene2", "gene3")
  barcodes <- c("cell1", "cell2")

  # Write test data to h5 file
  hf$create_group("matrix")
  hf$create_group("matrix/features")

  hf[["matrix/data"]] <- counts
  hf[["matrix/indices"]] <- indices
  hf[["matrix/indptr"]] <- indptr
  hf[["matrix/shape"]] <- shp
  hf[["matrix/features/name"]] <- features 
  hf[["matrix/barcodes"]] <- barcodes
  
  # Close the h5 file
  hf$close_all()

  # Call the load10XExpr function with the test h5 file
  spMat <- load10XExpr(visiumDir = temp_dir, h5filename = "test.h5")
  
  # Perform assertions to check if the data is loaded correctly
  expect_equal(length(grepl(h5_file,list.files())[grepl(h5_file,list.files())]),1,
               info = paste0("More than one file with ", h5_file))
  expect_equal(dim(spMat), c(3, 2), info = "Incorrect matrix dimensions")
  expect_equal(rownames(spMat), make.unique(features), info = "Incorrect row names")
  expect_equal(colnames(spMat), barcodes, info = "Incorrect column names")
  expect_equal(spMat[1, 1], log2(1 + counts[1]), info = "Incorrect matrix value")
  expect_equal(spMat[3, 2], log2(1 + counts[6]), info = "Incorrect matrix value")
  expect_true(is(spMat, "dgCMatrix"), info = "Incorrect matrix class")
  
  # Clean up the temporary directory
  unlink("test.h5")
})
