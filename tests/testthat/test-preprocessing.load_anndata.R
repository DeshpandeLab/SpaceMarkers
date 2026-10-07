test_that("load_anndata loads staple.h5ad into a SpaceMarkersExperiment", {
    skip_if_not_installed("anndataR")

    file <- testthat::test_path("assets", "staple.h5ad")
    original_read_h5ad <- get("read_h5ad", envir = asNamespace("anndataR"))
    sce <- suppressWarnings(original_read_h5ad(file, as = "SingleCellExperiment"))

    expected_spots <- rownames(S4Vectors::metadata(sce)$cell_type_composition)
    expected_patterns <- as.matrix(
        S4Vectors::metadata(sce)$cell_type_composition[expected_spots, , drop = FALSE]
    )
    expected_coords <- SingleCellExperiment::reducedDim(sce, "spatial")[
        expected_spots, , drop = FALSE
    ]
    colnames(expected_coords) <- c("y", "x")
    sigma <- S4Vectors::metadata(sce)$spatial[[1]][["scalefactors"]][["spot_diameter_fullres"]]
    expected_params <- matrix(
        c(sigma, 4),
        nrow = 2,
        ncol = ncol(expected_patterns),
        dimnames = list(c("sigmaOpt", "threshOpt"), colnames(expected_patterns))
    )
    colnames(expected_patterns) <- make.names(colnames(expected_patterns))
    colnames(expected_params) <- make.names(colnames(expected_params))

    with_mocked_bindings(
        read_h5ad = function(path, ...) {
            list(
                as_SingleCellExperiment = function() {
                    suppressWarnings(original_read_h5ad(path, as = "SingleCellExperiment"))
                }
            )
        },
        .package = "anndataR",
        {
            sme <- suppressWarnings(load_anndata(file, reader = "anndataR"))

            expect_s4_class(sme, "SpaceMarkersExperiment")
            expect_equal(dim(sme), c(nrow(sce), length(expected_spots)))
            expect_equal(colnames(sme), expected_spots)
            expect_equal(SpatialExperiment::spatialCoords(sme), expected_coords)
            expect_equal(as.matrix(spatial_patterns(sme)), expected_patterns)
            expect_equal(spatial_params(sme), expected_params)
        }
    )
})
