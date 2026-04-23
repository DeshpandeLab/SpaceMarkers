test_that("SpaceMarkersExperiment constructor works with minimal input", {
    mat <- matrix(rnorm(50), nrow = 5, ncol = 10)
    rownames(mat) <- paste0("gene_", seq_len(5))
    colnames(mat) <- paste0("spot_", seq_len(10))
    coords <- matrix(runif(20), ncol = 2)
    colnames(coords) <- c("y", "x")
    rownames(coords) <- colnames(mat)

    sme <- SpaceMarkersExperiment(
        assays = list(logcounts = mat),
        spatialCoords = coords
    )

    expect_s4_class(sme, "SpaceMarkersExperiment")
    expect_s4_class(sme, "SpatialExperiment")
    expect_equal(nrow(sme), 5L)
    expect_equal(ncol(sme), 10L)
})

test_that("SpaceMarkersExperiment constructor accepts list for spaceMarkers", {
    mat <- matrix(0, nrow = 2, ncol = 3)
    rownames(mat) <- c("g1", "g2")
    colnames(mat) <- c("s1", "s2", "s3")
    coords <- matrix(1:6, ncol = 2)
    colnames(coords) <- c("y", "x")
    rownames(coords) <- colnames(mat)

    sme <- SpaceMarkersExperiment(
        assays = list(logcounts = mat),
        spatialCoords = coords,
        spaceMarkers = list(analysis = "undirected")
    )

    expect_s4_class(sme, "SpaceMarkersExperiment")
    expect_equal(sme@spacemarkers$analysis, "undirected")
})

test_that("validity rejects invalid analysis type", {
    mat <- matrix(0, nrow = 2, ncol = 3)
    rownames(mat) <- c("g1", "g2")
    colnames(mat) <- c("s1", "s2", "s3")
    coords <- matrix(1:6, ncol = 2)
    colnames(coords) <- c("y", "x")
    rownames(coords) <- colnames(mat)

    sme <- SpaceMarkersExperiment(
        assays = list(logcounts = mat),
        spatialCoords = coords
    )
    sme@spacemarkers$analysis <- "invalid_type"
    expect_error(validObject(sme))
})

test_that("validity rejects unknown spaceMarkers elements", {
    mat <- matrix(0, nrow = 2, ncol = 3)
    rownames(mat) <- c("g1", "g2")
    colnames(mat) <- c("s1", "s2", "s3")
    coords <- matrix(1:6, ncol = 2)
    colnames(coords) <- c("y", "x")
    rownames(coords) <- colnames(mat)

    sme <- SpaceMarkersExperiment(
        assays = list(logcounts = mat),
        spatialCoords = coords
    )
    sme@spacemarkers$unknown_slot <- "bad"
    expect_error(validObject(sme))
})

test_that("spatial_patterns accessor works", {
    mat <- matrix(rnorm(20), nrow = 2, ncol = 10)
    rownames(mat) <- c("g1", "g2")
    colnames(mat) <- paste0("s", seq_len(10))
    coords <- matrix(runif(20), ncol = 2)
    colnames(coords) <- c("y", "x")
    rownames(coords) <- colnames(mat)

    cd <- S4Vectors::DataFrame(
        Pattern_1 = runif(10),
        Pattern_5 = runif(10),
        row.names = colnames(mat)
    )

    sme <- SpaceMarkersExperiment(
        assays = list(logcounts = mat),
        colData = cd,
        spatialCoords = coords,
        spaceMarkers = list(
            params = list(pattern_names = c("Pattern_1", "Pattern_5"))
        )
    )

    sp <- spatial_patterns(sme)
    expect_s4_class(sp, "DataFrame")
    expect_equal(ncol(sp), 2L)
    expect_equal(colnames(sp), c("Pattern_1", "Pattern_5"))
})

test_that("spatial_patterns returns NULL when no patterns set", {
    mat <- matrix(0, nrow = 2, ncol = 3)
    rownames(mat) <- c("g1", "g2")
    colnames(mat) <- c("s1", "s2", "s3")
    coords <- matrix(1:6, ncol = 2)
    colnames(coords) <- c("y", "x")
    rownames(coords) <- colnames(mat)

    sme <- SpaceMarkersExperiment(
        assays = list(logcounts = mat),
        spatialCoords = coords
    )
    expect_null(spatial_patterns(sme))
})

test_that("spatial_patterns<- setter works", {
    mat <- matrix(0, nrow = 2, ncol = 3)
    rownames(mat) <- c("g1", "g2")
    colnames(mat) <- c("s1", "s2", "s3")
    coords <- matrix(1:6, ncol = 2)
    colnames(coords) <- c("y", "x")
    rownames(coords) <- colnames(mat)

    sme <- SpaceMarkersExperiment(
        assays = list(logcounts = mat),
        spatialCoords = coords
    )

    pats <- data.frame(P1 = c(0.1, 0.2, 0.3), P2 = c(0.4, 0.5, 0.6))
    spatial_patterns(sme) <- pats

    sp <- spatial_patterns(sme)
    expect_equal(ncol(sp), 2L)
    expect_equal(colnames(sp), c("P1", "P2"))
    expect_equal(sme@spacemarkers$params$pattern_names, c("P1", "P2"))
})

test_that("spatial_patterns<- reorders by rownames and rejects mismatches", {
    mat <- matrix(0, nrow = 2, ncol = 3)
    rownames(mat) <- c("g1", "g2")
    colnames(mat) <- c("s1", "s2", "s3")
    coords <- matrix(1:6, ncol = 2); colnames(coords) <- c("y", "x")
    rownames(coords) <- colnames(mat)
    sme <- SpaceMarkersExperiment(
        assays = list(logcounts = mat), spatialCoords = coords)

    # Rownames in reverse order should be reordered to match colnames(x)
    pats <- data.frame(P1 = c(0.3, 0.2, 0.1), row.names = c("s3", "s2", "s1"))
    spatial_patterns(sme) <- pats
    expect_equal(spatial_patterns(sme)$P1, c(0.1, 0.2, 0.3))

    # Rowname set mismatch -> error
    bad <- data.frame(P1 = c(0.1, 0.2, 0.3),
                      row.names = c("s1", "s2", "sX"))
    expect_error(spatial_patterns(sme) <- bad, "must match colnames")

    # No rownames, wrong length -> error
    wrong_len <- data.frame(P1 = c(0.1, 0.2))
    expect_error(spatial_patterns(sme) <- wrong_len, "nrow\\(value\\)")
})

test_that("spatial_params accessor works", {
    mat <- matrix(0, nrow = 2, ncol = 3)
    rownames(mat) <- c("g1", "g2")
    colnames(mat) <- c("s1", "s2", "s3")
    coords <- matrix(1:6, ncol = 2)
    colnames(coords) <- c("y", "x")
    rownames(coords) <- colnames(mat)

    optP <- matrix(c(6, 2, 6, 2), nrow = 2)
    rownames(optP) <- c("sigmaOpt", "threshOpt")
    colnames(optP) <- c("P1", "P2")

    sme <- SpaceMarkersExperiment(
        assays = list(logcounts = mat),
        spatialCoords = coords,
        spaceMarkers = list(params = list(spatial_params = optP))
    )

    expect_equal(spatial_params(sme), optP)
})

test_that("spatial_params returns NULL when not set", {
    mat <- matrix(0, nrow = 2, ncol = 3)
    rownames(mat) <- c("g1", "g2")
    colnames(mat) <- c("s1", "s2", "s3")
    coords <- matrix(1:6, ncol = 2)
    colnames(coords) <- c("y", "x")
    rownames(coords) <- colnames(mat)

    sme <- SpaceMarkersExperiment(
        assays = list(logcounts = mat),
        spatialCoords = coords
    )
    expect_null(spatial_params(sme))
})

test_that("spatial_params<- setter works", {
    mat <- matrix(0, nrow = 2, ncol = 3)
    rownames(mat) <- c("g1", "g2")
    colnames(mat) <- c("s1", "s2", "s3")
    coords <- matrix(1:6, ncol = 2)
    colnames(coords) <- c("y", "x")
    rownames(coords) <- colnames(mat)

    sme <- SpaceMarkersExperiment(
        assays = list(logcounts = mat),
        spatialCoords = coords
    )

    optP <- matrix(c(6, 2, 6, 2), nrow = 2)
    rownames(optP) <- c("sigmaOpt", "threshOpt")
    colnames(optP) <- c("P1", "P2")

    spatial_params(sme) <- optP
    expect_equal(spatial_params(sme), optP)
})

test_that("params returns full hyperparameter list", {
    mat <- matrix(0, nrow = 2, ncol = 3)
    rownames(mat) <- c("g1", "g2")
    colnames(mat) <- c("s1", "s2", "s3")
    coords <- matrix(1:6, ncol = 2)
    colnames(coords) <- c("y", "x")
    rownames(coords) <- colnames(mat)

    optP <- matrix(c(6, 2, 6, 2), nrow = 2)
    rownames(optP) <- c("sigmaOpt", "threshOpt")
    colnames(optP) <- c("P1", "P2")

    sme <- SpaceMarkersExperiment(
        assays = list(logcounts = mat),
        spatialCoords = coords,
        spaceMarkers = list(params = list(
            pattern_names = c("P1", "P2"),
            spatial_params = optP,
            min_gene_expr = 10,
            directed = FALSE
        ))
    )

    p <- params(sme)
    expect_true(is.list(p))
    expect_equal(p$spatial_params, optP)
    expect_equal(p$pattern_names, c("P1", "P2"))
    expect_equal(p$min_gene_expr, 10)
    expect_equal(p$directed, FALSE)
})

test_that("params returns NULL when not set", {
    mat <- matrix(0, nrow = 2, ncol = 3)
    rownames(mat) <- c("g1", "g2")
    colnames(mat) <- c("s1", "s2", "s3")
    coords <- matrix(1:6, ncol = 2)
    colnames(coords) <- c("y", "x")
    rownames(coords) <- colnames(mat)

    sme <- SpaceMarkersExperiment(
        assays = list(logcounts = mat),
        spatialCoords = coords
    )
    expect_null(params(sme))
})

test_that("hotspots accessor works", {
    mat <- matrix(0, nrow = 2, ncol = 3)
    rownames(mat) <- c("g1", "g2")
    colnames(mat) <- c("s1", "s2", "s3")
    coords <- matrix(1:6, ncol = 2)
    colnames(coords) <- c("y", "x")
    rownames(coords) <- colnames(mat)

    hs <- data.frame(
        barcode = c("s1", "s2", "s3"),
        P1 = c("P1", NA, "P1"),
        P2 = c(NA, "P2", "P2")
    )

    sme <- SpaceMarkersExperiment(
        assays = list(logcounts = mat),
        spatialCoords = coords
    )
    S4Vectors::metadata(sme)$hotspots <- list(undirected = hs)

    expect_equal(hotspots(sme, type = "undirected"), hs)
    expect_null(hotspots(sme, type = "pattern"))
})

test_that("undirected_scores and directed_scores accessors work", {
    mat <- matrix(0, nrow = 2, ncol = 3)
    rownames(mat) <- c("g1", "g2")
    colnames(mat) <- c("s1", "s2", "s3")
    coords <- matrix(1:6, ncol = 2)
    colnames(coords) <- c("y", "x")
    rownames(coords) <- colnames(mat)

    ims <- data.frame(Gene = c("g1", "g2"), P1_P2 = c(1.5, 0.8))
    ds <- data.frame(
        gene = c("g1", "g2"),
        cell_interaction = c("P1_near_P2", "P1_near_P2"),
        statistic = c(2.1, -0.5)
    )

    sme <- SpaceMarkersExperiment(
        assays = list(logcounts = mat),
        spatialCoords = coords,
        spaceMarkers = list(
            results = list(undirected_scores = ims, directed_scores = ds),
            analysis = "both"
        )
    )

    expect_equal(undirected_scores(sme), ims)
    expect_equal(directed_scores(sme), ds)
    expect_equal(analysis_type(sme), "both")
})

test_that("lr_scores accessor works", {
    mat <- matrix(0, nrow = 2, ncol = 3)
    rownames(mat) <- c("g1", "g2")
    colnames(mat) <- c("s1", "s2", "s3")
    coords <- matrix(1:6, ncol = 2)
    colnames(coords) <- c("y", "x")
    rownames(coords) <- colnames(mat)

    sme <- SpaceMarkersExperiment(
        assays = list(logcounts = mat),
        spatialCoords = coords
    )
    expect_null(lr_scores(sme))

    lr <- matrix(c(0.5, 0.3), nrow = 2, ncol = 1)
    rownames(lr) <- c("LIG1_REC1", "LIG2_REC2")
    colnames(lr) <- "to_P1"
    sme <- SpaceMarkersExperiment(
        assays = list(logcounts = mat),
        spatialCoords = coords,
        spaceMarkers = list(results = list(lr_scores = lr))
    )
    expect_equal(lr_scores(sme), lr)
})

test_that("overlap_scores accessor works", {
    mat <- matrix(0, nrow = 2, ncol = 3)
    rownames(mat) <- c("g1", "g2")
    colnames(mat) <- c("s1", "s2", "s3")
    coords <- matrix(1:6, ncol = 2)
    colnames(coords) <- c("y", "x")
    rownames(coords) <- colnames(mat)

    sme <- SpaceMarkersExperiment(
        assays = list(logcounts = mat),
        spatialCoords = coords
    )
    expect_null(overlap_scores(sme))

    ov <- data.frame(
        pattern1 = "P1", pattern2 = "P2", overlapScore = 0.67
    )
    sme <- SpaceMarkersExperiment(
        assays = list(logcounts = mat),
        spatialCoords = coords,
        spaceMarkers = list(results = list(overlap_scores = ov))
    )
    expect_equal(overlap_scores(sme), ov)
})

test_that("interactions accessor works", {
    mat <- matrix(0, nrow = 2, ncol = 3)
    rownames(mat) <- c("g1", "g2")
    colnames(mat) <- c("s1", "s2", "s3")
    coords <- matrix(1:6, ncol = 2)
    colnames(coords) <- c("y", "x")
    rownames(coords) <- colnames(mat)

    pair_result <- list(
        interacting_genes = list(data.frame(Gene = "g1", KW.p = 0.01)),
        hotspots = data.frame(barcode = c("s1", "s2", "s3")),
        patterns = c("P1", "P2")
    )

    sme <- SpaceMarkersExperiment(
        assays = list(logcounts = mat),
        spatialCoords = coords
    )
    S4Vectors::metadata(sme)$interactions <- list(P1_P2 = pair_result)

    expect_equal(interactions(sme, pair = "P1_P2"), pair_result)
    expect_equal(length(interactions(sme)), 1L)
    expect_null(interactions(sme, pair = "nonexistent"))
})

test_that("influence_map accessor works", {
    mat <- matrix(0, nrow = 2, ncol = 3)
    rownames(mat) <- c("g1", "g2")
    colnames(mat) <- c("s1", "s2", "s3")
    coords <- matrix(1:6, ncol = 2)
    colnames(coords) <- c("y", "x")
    rownames(coords) <- colnames(mat)

    inf_df <- data.frame(
        barcode = c("s1", "s2", "s3"),
        x = 1:3, y = 4:6,
        P1 = c(0.1, 0.5, 0.9),
        P2 = c(0.3, 0.2, 0.8)
    )

    sme <- SpaceMarkersExperiment(
        assays = list(logcounts = mat),
        spatialCoords = coords
    )
    S4Vectors::metadata(sme)$influence <- inf_df

    expect_equal(influence_map(sme), inf_df)
})

test_that("show method runs without error", {
    mat <- matrix(0, nrow = 2, ncol = 3)
    rownames(mat) <- c("g1", "g2")
    colnames(mat) <- c("s1", "s2", "s3")
    coords <- matrix(1:6, ncol = 2)
    colnames(coords) <- c("y", "x")
    rownames(coords) <- colnames(mat)

    sme <- SpaceMarkersExperiment(
        assays = list(logcounts = mat),
        spatialCoords = coords
    )

    expect_output(show(sme), "spacemarkers analysis: none")
})

test_that("show method displays analysis info", {
    mat <- matrix(0, nrow = 2, ncol = 3)
    rownames(mat) <- c("g1", "g2")
    colnames(mat) <- c("s1", "s2", "s3")
    coords <- matrix(1:6, ncol = 2)
    colnames(coords) <- c("y", "x")
    rownames(coords) <- colnames(mat)

    ims <- data.frame(Gene = c("g1", "g2"), P1_P2 = c(1.5, 0.8))

    sme <- SpaceMarkersExperiment(
        assays = list(logcounts = mat),
        spatialCoords = coords,
        spaceMarkers = list(
            results = list(undirected_scores = ims),
            analysis = "undirected",
            params = list(pattern_names = c("P1", "P2"))
        )
    )

    expect_output(show(sme), "undirected")
    expect_output(show(sme), "P1, P2")
})

test_that("SpaceMarkersExperiment accepts SpatialExperiment input", {
    mat <- matrix(rnorm(30), nrow = 3, ncol = 10)
    rownames(mat) <- paste0("gene_", 1:3)
    colnames(mat) <- paste0("spot_", 1:10)
    coords <- matrix(runif(20), ncol = 2)
    colnames(coords) <- c("y", "x")
    rownames(coords) <- colnames(mat)

    spe <- SpatialExperiment::SpatialExperiment(
        assays = list(logcounts = mat),
        spatialCoords = coords
    )

    sme <- SpaceMarkersExperiment(spe)
    expect_s4_class(sme, "SpaceMarkersExperiment")
    expect_s4_class(sme, "SpatialExperiment")
    expect_equal(nrow(sme), 3L)
    expect_equal(ncol(sme), 10L)
    expect_null(params(sme))
})

test_that("SpaceMarkersExperiment coerces SPE with features", {
    mat <- matrix(rnorm(30), nrow = 3, ncol = 10)
    rownames(mat) <- paste0("gene_", 1:3)
    colnames(mat) <- paste0("spot_", 1:10)
    coords <- matrix(runif(20), ncol = 2)
    colnames(coords) <- c("y", "x")
    rownames(coords) <- colnames(mat)

    spe <- SpatialExperiment::SpatialExperiment(
        assays = list(logcounts = mat),
        spatialCoords = coords
    )

    feats <- data.frame(
        P1 = runif(10), P2 = runif(10),
        row.names = colnames(mat)
    )

    sme <- SpaceMarkersExperiment(spe, features = feats)
    expect_s4_class(sme, "SpaceMarkersExperiment")
    sp <- spatial_patterns(sme)
    expect_equal(ncol(sp), 2L)
    expect_equal(colnames(sp), c("P1", "P2"))
})

test_that("as() coercion from SPE to SME works", {
    mat <- matrix(rnorm(30), nrow = 3, ncol = 10)
    rownames(mat) <- paste0("gene_", 1:3)
    colnames(mat) <- paste0("spot_", 1:10)
    coords <- matrix(runif(20), ncol = 2)
    colnames(coords) <- c("y", "x")
    rownames(coords) <- colnames(mat)

    spe <- SpatialExperiment::SpatialExperiment(
        assays = list(logcounts = mat),
        spatialCoords = coords
    )

    sme <- as(spe, "SpaceMarkersExperiment")
    expect_s4_class(sme, "SpaceMarkersExperiment")
    expect_equal(nrow(sme), 3L)
    expect_equal(ncol(sme), 10L)
})

test_that("SCE -> SME promotes reducedDim('spatial') to spatialCoords", {
    mat <- matrix(runif(30), nrow = 3, ncol = 10)
    rownames(mat) <- paste0("g", 1:3)
    colnames(mat) <- paste0("c", 1:10)
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(logcounts = mat)
    )
    SingleCellExperiment::reducedDim(sce, "spatial") <-
        matrix(runif(20), ncol = 2,
               dimnames = list(colnames(sce), c("row", "col")))
    sme <- as(sce, "SpaceMarkersExperiment")
    expect_s4_class(sme, "SpaceMarkersExperiment")
    expect_equal(nrow(SpatialExperiment::spatialCoords(sme)), 10L)
    expect_equal(ncol(SpatialExperiment::spatialCoords(sme)), 2L)
})

test_that("SCE -> SME leaves spatialCoords empty when no spatial reducedDim", {
    mat <- matrix(runif(30), nrow = 3, ncol = 10)
    rownames(mat) <- paste0("g", 1:3)
    colnames(mat) <- paste0("c", 1:10)
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(logcounts = mat)
    )
    sme <- as(sce, "SpaceMarkersExperiment")
    expect_s4_class(sme, "SpaceMarkersExperiment")
    # SPE default from SCE yields a matrix with ncells rows and 0 columns.
    expect_equal(ncol(SpatialExperiment::spatialCoords(sme)), 0L)
})

test_that("load_anndata errors clearly when zellkonverter is unavailable", {
    # Simulate zellkonverter missing by shadowing requireNamespace via the
    # function's own environment. Faster and more portable than actually
    # unloading the package.
    with_mocked_bindings(
        requireNamespace = function(package, ...) {
            if (identical(package, "zellkonverter")) FALSE
            else base::requireNamespace(package, ...)
        },
        .package = "base",
        expect_error(load_anndata("ignored.h5ad"),
                     "requires the 'zellkonverter' package")
    )
})
