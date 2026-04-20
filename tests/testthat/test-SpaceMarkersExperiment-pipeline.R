# Integration tests for SpaceMarkersExperiment pipeline

# Helper to create a minimal SME for testing
.create_test_sme <- function(n_genes = 10, n_spots = 20) {
    mat <- matrix(abs(rnorm(n_genes * n_spots)), nrow = n_genes, ncol = n_spots)
    rownames(mat) <- paste0("gene_", seq_len(n_genes))
    colnames(mat) <- paste0("spot_", seq_len(n_spots))

    coords <- matrix(runif(n_spots * 2) * 100, ncol = 2)
    colnames(coords) <- c("y", "x")
    rownames(coords) <- colnames(mat)

    # Create two patterns with some spatial structure
    pat1 <- runif(n_spots)
    pat2 <- runif(n_spots)
    cd <- S4Vectors::DataFrame(
        Pattern_1 = pat1,
        Pattern_2 = pat2,
        row.names = colnames(mat)
    )

    SpaceMarkersExperiment(
        assays = list(logcounts = mat),
        colData = cd,
        spatialCoords = coords,
        spaceMarkers = list(
            params = list(pattern_names = c("Pattern_1", "Pattern_2"))
        )
    )
}

test_that("SME constructed from test helper is valid", {
    sme <- .create_test_sme()
    expect_s4_class(sme, "SpaceMarkersExperiment")
    expect_equal(nrow(sme), 10L)
    expect_equal(ncol(sme), 20L)
    expect_equal(spatial_patterns(sme) |> ncol(), 2L)
    expect_null(spatial_params(sme))
    expect_null(undirected_scores(sme))
    expect_null(directed_scores(sme))
    expect_null(analysis_type(sme))
})

test_that(".sme_to_spPatterns produces correct legacy format", {
    sme <- .create_test_sme()
    sp <- SpaceMarkers:::.sme_to_spPatterns(sme)

    expect_true(is.data.frame(sp))
    expect_true(all(c("barcode", "y", "x", "Pattern_1", "Pattern_2") %in%
                    colnames(sp)))
    expect_equal(nrow(sp), 20L)
    expect_equal(sp$barcode, colnames(sme))
})

test_that(".sme_get_expr extracts expression matrix", {
    sme <- .create_test_sme()
    expr <- SpaceMarkers:::.sme_get_expr(sme)
    expect_equal(dim(expr), c(10L, 20L))
    expect_equal(rownames(expr), paste0("gene_", seq_len(10)))
})

test_that("storing and retrieving results round-trips correctly", {
    sme <- .create_test_sme()

    # Simulate storing undirected results
    ims <- data.frame(
        Gene = paste0("gene_", 1:5),
        Pattern_1_Pattern_2 = runif(5)
    )
    sm <- sme@spacemarkers
    sm$results <- list(undirected_scores = ims)
    sm$analysis <- "undirected"
    sme@spacemarkers <- sm

    expect_equal(undirected_scores(sme), ims)
    expect_equal(analysis_type(sme), "undirected")

    # Store hotspots in metadata
    hs <- data.frame(
        barcode = colnames(sme)[1:5],
        y = 1:5, x = 6:10,
        Pattern_1 = c("Pattern_1", NA, "Pattern_1", NA, "Pattern_1"),
        Pattern_2 = c(NA, "Pattern_2", "Pattern_2", NA, NA)
    )
    S4Vectors::metadata(sme)$hotspots$undirected <- hs
    expect_equal(hotspots(sme, type = "undirected"), hs)
})

test_that("storing and retrieving directed results", {
    sme <- .create_test_sme()

    ds <- data.frame(
        gene = paste0("gene_", 1:4),
        cell_interaction = rep(c("P1_near_P2", "P2_near_P1"), each = 2),
        statistic = rnorm(4),
        p.value = runif(4),
        effect_size = rnorm(4)
    )

    sm <- sme@spacemarkers
    sm$results <- list(directed_scores = ds)
    sm$analysis <- "directed"
    sme@spacemarkers <- sm

    expect_equal(directed_scores(sme), ds)
    expect_equal(analysis_type(sme), "directed")
})

test_that("interactions stored in metadata are accessible", {
    sme <- .create_test_sme()

    pair_result <- list(
        interacting_genes = list(data.frame(
            Gene = paste0("gene_", 1:3),
            SpaceMarkersMetric = c(2.1, 1.5, 0.3)
        )),
        hotspots = data.frame(
            barcode = colnames(sme),
            Pattern_1 = rep(c("Pattern_1", NA), length.out = ncol(sme)),
            Pattern_2 = rep(c(NA, "Pattern_2"), length.out = ncol(sme))
        ),
        patterns = c("Pattern_1", "Pattern_2")
    )

    S4Vectors::metadata(sme)$interactions <- list(
        Pattern_1_Pattern_2 = pair_result
    )

    result <- interactions(sme, pair = "Pattern_1_Pattern_2")
    expect_equal(result$patterns, c("Pattern_1", "Pattern_2"))
    expect_equal(length(interactions(sme)), 1L)
})

test_that("influence map stored in metadata is accessible", {
    sme <- .create_test_sme()

    inf_df <- data.frame(
        barcode = colnames(sme),
        y = runif(ncol(sme)),
        x = runif(ncol(sme)),
        Pattern_1 = runif(ncol(sme)),
        Pattern_2 = runif(ncol(sme))
    )

    S4Vectors::metadata(sme)$influence <- inf_df
    expect_equal(influence_map(sme), inf_df)
})

test_that("SME can be subset and accessors still work", {
    sme <- .create_test_sme()

    # Subset to first 5 spots
    sub <- sme[, 1:5]
    expect_s4_class(sub, "SpaceMarkersExperiment")
    expect_equal(ncol(sub), 5L)
    # Patterns should still be accessible
    sp <- spatial_patterns(sub)
    expect_equal(nrow(sp), 5L)
})

test_that("get_spatial_features accepts SpatialExperiment objects", {
    mat <- matrix(abs(rnorm(20)), nrow = 2, ncol = 10)
    rownames(mat) <- c("g1", "g2")
    colnames(mat) <- paste0("s", seq_len(10))
    coords <- matrix(runif(20), ncol = 2)
    colnames(coords) <- c("y", "x")
    rownames(coords) <- colnames(mat)

    cd <- S4Vectors::DataFrame(
        PatternA = runif(10),
        PatternB = runif(10),
        row.names = colnames(mat)
    )

    spe <- SpatialExperiment::SpatialExperiment(
        assays = list(logcounts = mat),
        colData = cd,
        spatialCoords = coords
    )

    features <- get_spatial_features(spe)
    expect_true(is.data.frame(features))
    expect_true(all(c("PatternA", "PatternB") %in% colnames(features)))
    expect_equal(nrow(features), 10L)
})

# ---- add_features tests ----

test_that("add_features with aligned spots adds features without warning", {
    mat <- matrix(abs(rnorm(50)), nrow = 5, ncol = 10)
    rownames(mat) <- paste0("gene_", seq_len(5))
    colnames(mat) <- paste0("spot_", seq_len(10))
    coords <- matrix(runif(20), ncol = 2)
    colnames(coords) <- c("y", "x")
    rownames(coords) <- colnames(mat)

    sme <- SpaceMarkersExperiment(
        assays = list(logcounts = mat),
        spatialCoords = coords
    )

    feats <- data.frame(
        P1 = runif(10), P2 = runif(10),
        row.names = colnames(mat)
    )

    expect_silent(sme2 <- add_features(sme, feats))
    expect_s4_class(sme2, "SpaceMarkersExperiment")
    expect_equal(ncol(sme2), 10L)
    sp <- spatial_patterns(sme2)
    expect_equal(ncol(sp), 2L)
    expect_equal(colnames(sp), c("P1", "P2"))
    expect_equal(sme2@spacemarkers$params$pattern_names, c("P1", "P2"))
})

test_that("add_features with mismatched spots warns and intersects", {
    mat <- matrix(abs(rnorm(50)), nrow = 5, ncol = 10)
    rownames(mat) <- paste0("gene_", seq_len(5))
    colnames(mat) <- paste0("spot_", seq_len(10))
    coords <- matrix(runif(20), ncol = 2)
    colnames(coords) <- c("y", "x")
    rownames(coords) <- colnames(mat)

    sme <- SpaceMarkersExperiment(
        assays = list(logcounts = mat),
        spatialCoords = coords
    )

    # Features with only 7 matching spots + 3 extras
    feat_spots <- c(paste0("spot_", 1:7), "extra_1", "extra_2", "extra_3")
    feats <- data.frame(
        P1 = runif(10), P2 = runif(10),
        row.names = feat_spots
    )

    expect_warning(
        sme2 <- add_features(sme, feats),
        "Spot mismatch"
    )
    expect_equal(ncol(sme2), 7L)
    expect_equal(sme2@spacemarkers$params$pattern_names, c("P1", "P2"))
})

test_that("add_features errors with no common spots", {
    mat <- matrix(abs(rnorm(20)), nrow = 2, ncol = 10)
    rownames(mat) <- c("g1", "g2")
    colnames(mat) <- paste0("spot_", seq_len(10))
    coords <- matrix(runif(20), ncol = 2)
    colnames(coords) <- c("y", "x")
    rownames(coords) <- colnames(mat)

    sme <- SpaceMarkersExperiment(
        assays = list(logcounts = mat),
        spatialCoords = coords
    )

    feats <- data.frame(P1 = 1:3, row.names = c("a", "b", "c"))
    expect_error(add_features(sme, feats), "No common spots")
})

test_that("add_features accepts matrix input", {
    mat <- matrix(abs(rnorm(30)), nrow = 3, ncol = 10)
    rownames(mat) <- paste0("gene_", seq_len(3))
    colnames(mat) <- paste0("spot_", seq_len(10))
    coords <- matrix(runif(20), ncol = 2)
    colnames(coords) <- c("y", "x")
    rownames(coords) <- colnames(mat)

    sme <- SpaceMarkersExperiment(
        assays = list(logcounts = mat),
        spatialCoords = coords
    )

    feat_mat <- matrix(runif(20), ncol = 2)
    colnames(feat_mat) <- c("Pattern_A", "Pattern_B")
    rownames(feat_mat) <- colnames(mat)

    sme2 <- add_features(sme, feat_mat)
    expect_equal(sme2@spacemarkers$params$pattern_names,
                 c("Pattern_A", "Pattern_B"))
})
