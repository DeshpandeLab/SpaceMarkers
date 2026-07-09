# .get_seurat_assay_data() and load_seurat() rely on the Seurat package,
# which is an optional (Suggests) dependency. All tests below skip cleanly
# when it is not installed.

test_that("load_seurat errors for input that is neither a path nor a Seurat object", {
  skip_if_not_installed("Seurat")
  expect_error(load_seurat(list(a = 1)),
               "must be a path")
  expect_error(load_seurat(42),
               "must be a path")
})

test_that("load_seurat errors when the .rds file does not contain a Seurat object", {
  skip_if_not_installed("Seurat")
  temp <- tempfile(fileext = ".rds")
  on.exit(unlink(temp))
  saveRDS(data.frame(x = 1), temp)
  expect_error(load_seurat(temp), "did not contain a Seurat object")
})

.make_fixture_seurat <- function(n_spots = 12, n_genes = 15) {
  counts <- matrix(rpois(n_genes * n_spots, lambda = 3), nrow = n_genes,
                   dimnames = list(paste0("G", seq_len(n_genes)),
                                   paste0("spot", seq_len(n_spots))))
  # Use a sparse matrix and populate the "data"/"scale.data" layers so
  # Seurat::as.SingleCellExperiment() (called inside load_seurat()) has
  # something to read for every layer instead of warning about empty
  # ones -- purely fixture hygiene, unrelated to load_seurat() itself.
  counts <- Matrix::Matrix(counts, sparse = TRUE)
  seurat_object <- Seurat::CreateSeuratObject(counts = counts, assay = "Spatial")
  seurat_object <- Seurat::NormalizeData(seurat_object, verbose = FALSE)
  seurat_object <- Seurat::ScaleData(seurat_object,
                                     features = rownames(seurat_object),
                                     verbose = FALSE)
  seurat_object
}

# Creating an assay from data= alone (no counts=) is exactly how
# deconvolution/cell-type-fraction assays normally look (e.g. RCTD output),
# but Seurat's own `[[<-` assignment warns "Layer counts isn't present..."
# every time -- expected Seurat-internal noise, not something load_seurat()
# can or should silence for real callers, so it's suppressed only here.
.add_deconv_assay <- function(seurat_object, deconv) {
  suppressWarnings(
    seurat_object[["deconv"]] <- Seurat::CreateAssayObject(data = deconv)
  )
  seurat_object
}

test_that("load_seurat(deconv_assay = NULL) adds no patterns and points to get_spatial_features", {
  skip_if_not_installed("Seurat")
  seurat_object <- .make_fixture_seurat()
  expect_message(sme <- load_seurat(seurat_object, deconv_assay = NULL),
                 "get_spatial_features")
  expect_s4_class(sme, "SpaceMarkersExperiment")
  expect_null(spatial_patterns(sme))
})

test_that("load_seurat errors when deconv_assay is not found on the object", {
  skip_if_not_installed("Seurat")
  seurat_object <- .make_fixture_seurat()
  expect_error(load_seurat(seurat_object, deconv_assay = "nope"),
               "not found in the Seurat object")
})

test_that("load_seurat adds spatial patterns from the named deconv assay", {
  skip_if_not_installed("Seurat")
  seurat_object <- .make_fixture_seurat()
  n_spots <- ncol(seurat_object)
  
  deconv <- matrix(
    runif(3 * n_spots),
    nrow = 3,
    dimnames = list(c("CT1", "CT2", "CT3"), colnames(seurat_object))
  )
  
  seurat_object <- .add_deconv_assay(seurat_object, deconv)
  
  expect_warning(
    sme <- load_seurat(seurat_object, deconv_assay = "deconv"),
    "no images"
  )
  
  expect_s4_class(sme, "SpaceMarkersExperiment")
  expect_equal(ncol(sme), n_spots)
  expect_setequal(colnames(spatial_patterns(sme)), c("CT1", "CT2", "CT3"))
  expect_equal(
    as.matrix(spatial_patterns(sme))[, c("CT1", "CT2", "CT3")],
    t(deconv)[colnames(sme), ],
    ignore_attr = TRUE
  )
})
test_that("load_seurat warns (not errors) when the object has no images", {
  skip_if_not_installed("Seurat")
  seurat_object <- .make_fixture_seurat()
  n_spots <- ncol(seurat_object)
  deconv <- matrix(runif(n_spots), nrow = 1,
                   dimnames = list("CT1", colnames(seurat_object)))
  seurat_object <- .add_deconv_assay(seurat_object, deconv)
  
  expect_warning(sme <- load_seurat(seurat_object, deconv_assay = "deconv"),
                 "no images")
  expect_s4_class(sme, "SpaceMarkersExperiment")
  expect_null(spatial_params(sme))
})