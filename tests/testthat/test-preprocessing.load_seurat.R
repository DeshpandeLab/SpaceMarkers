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

test_that("load_seurat requires the Seurat package", {
    with_mocked_bindings(
        requireNamespace = function(package, ...) {
            if (identical(package, "Seurat")) FALSE else base::requireNamespace(package, ...)
        },
        .package = "SpaceMarkers",
        {
            expect_error(load_seurat("some/path.rds"), "requires the 'Seurat' package")
        }
    )
})

.make_fixture_seurat <- function(n_spots = 12, n_genes = 15) {
    counts <- matrix(rpois(n_genes * n_spots, lambda = 3), nrow = n_genes,
                     dimnames = list(paste0("G", seq_len(n_genes)),
                                     paste0("spot", seq_len(n_spots))))
    Seurat::CreateSeuratObject(counts = counts, assay = "Spatial")
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
    deconv <- matrix(runif(3 * n_spots), nrow = 3,
                     dimnames = list(c("CT1", "CT2", "CT3"),
                                     colnames(seurat_object)))
    seurat_object[["deconv"]] <- Seurat::CreateAssayObject(data = deconv)

    sme <- load_seurat(seurat_object, deconv_assay = "deconv")
    expect_s4_class(sme, "SpaceMarkersExperiment")
    expect_equal(ncol(sme), n_spots)
    expect_setequal(colnames(spatial_patterns(sme)), c("CT1", "CT2", "CT3"))
    expect_equal(
        as.matrix(spatial_patterns(sme))[, c("CT1", "CT2", "CT3")],
        t(deconv)[colnames(sme), ],
        ignore_attr = TRUE
    )
})

test_that("load_seurat derives spatial_params from the image scale factors", {
    skip_if_not_installed("Seurat")
    seurat_object <- .make_fixture_seurat()
    n_spots <- ncol(seurat_object)
    deconv <- matrix(runif(2 * n_spots), nrow = 2,
                     dimnames = list(c("CT1", "CT2"), colnames(seurat_object)))
    seurat_object[["deconv"]] <- Seurat::CreateAssayObject(data = deconv)

    # load_seurat() only relies on the image object having a `scale.factors`
    # slot (read generically via methods::slot()); it doesn't depend on the
    # real Seurat/SeuratObject image S4 class, so a minimal stand-in keeps
    # this test stable across Seurat versions.
    methods::setClass("SpaceMarkersFakeSeuratImage",
                      representation(scale.factors = "list"),
                      where = globalenv())
    on.exit(suppressWarnings(methods::removeClass(
        "SpaceMarkersFakeSeuratImage", where = globalenv())), add = TRUE)
    fake_image <- methods::new("SpaceMarkersFakeSeuratImage",
                               scale.factors = list(spot = 10, lowres = 0.05))
    seurat_object@images[["slice1"]] <- fake_image

    sme <- load_seurat(seurat_object, deconv_assay = "deconv", threshold = 4)
    op <- spatial_params(sme)
    expect_false(is.null(op))
    expect_equal(unname(op["sigmaOpt", ]), rep(10 * 0.05, 2))
    expect_equal(unname(op["threshOpt", ]), rep(4, 2))
})

test_that("load_seurat warns (not errors) when the object has no images", {
    skip_if_not_installed("Seurat")
    seurat_object <- .make_fixture_seurat()
    n_spots <- ncol(seurat_object)
    deconv <- matrix(runif(n_spots), nrow = 1,
                     dimnames = list("CT1", colnames(seurat_object)))
    seurat_object[["deconv"]] <- Seurat::CreateAssayObject(data = deconv)

    expect_warning(sme <- load_seurat(seurat_object, deconv_assay = "deconv"),
                   "no images")
    expect_s4_class(sme, "SpaceMarkersExperiment")
    expect_null(spatial_params(sme))
})
