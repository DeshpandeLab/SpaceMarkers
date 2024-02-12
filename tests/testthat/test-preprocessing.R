
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
  
  # Clean up the temporary directory
  unlink("test.h5")
})

# Define the test cases
test_that("load10XCoords loads data correctly", {
  # Create a temporary directory for testing
  temp_dir <- "spatial"
  # Create test data
  dir.create(temp_dir)
  scale_json <- file.path(temp_dir, "scalefactors_json.json")
  coord_file <- file.path(temp_dir, "tissue_positions_list.csv")
  # Write test data to files
  write('{"spot_diameter_fullres": 10, "scalefactors_lowres": 0.2}', scale_json)
  writeLines("spot_1,0.5,0.5,0,200,600\nspot_2,1.5,1.5,0,300,700", coord_file)
  
  #Check for expected directories and file names
  expect_equal(dir.exists('spatial'),TRUE,info = "No directory callled 'spatial'")
  expect_equal(file.exists('spatial/scalefactors_json.json'),
               TRUE,info ="No file called scalefactors_json.json")
  expect_equal(file.exists('spatial/tissue_positions_list.csv'),
               TRUE,info ="No file called tissue_positions_list.csv")
  
  # Call the load10XCoords function with the test directory
  coord_values <- load10XCoords(visiumDir = ".", resolution = "lowres")
  
  # Perform assertions to check if the data is loaded correctly
  expect_equal(coord_values$barcode, c("spot_1", "spot_2"), info = "Incorrect barcode values")
  expect_equal(coord_values$y, c(40, 60), info = "Incorrect y coordinates")
  expect_equal(coord_values$x, c(120, 140), info = "Incorrect x coordinates")
  
  # Clean up the temporary directory
  unlink(temp_dir, recursive = TRUE)
})

