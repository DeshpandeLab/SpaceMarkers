test_that("SpaceMarkers() on a manually-built SME without visiumDir uses default params with a message", {
    set.seed(2)
    nb <- 30; ng <- 8
    counts <- matrix(rpois(ng * nb, 3), nrow = ng,
                     dimnames = list(paste0("G", seq_len(ng)),
                                     paste0("s", seq_len(nb))))
    coords <- matrix(runif(2 * nb, 0, 100), ncol = 2,
                     dimnames = list(paste0("s", seq_len(nb)),
                                     c("y", "x")))
    pat <- data.frame(P1 = runif(nb), P2 = runif(nb),
                      row.names = paste0("s", seq_len(nb)))

    sme <- SpaceMarkersExperiment(
        assays = list(logcounts = counts),
        colData = pat,
        spatialCoords = coords)
    # SpaceMarkersExperiment ctor doesn't infer pattern names from colData
    # on its own; declare them explicitly.
    spatial_patterns(sme) <- pat

    # No visiumDir, no spatial_params -> we expect a *message*, not a hard
    # stop, and the SME slot should get filled with defaults.
    expect_message(
        sme_filt <- SpaceMarkers:::.apply_sme_filters(sme,
            min.gene.expr = 0,
            sigma = NULL, threshold = 4,
            resolution = "lowres", spatialDir = "spatial",
            pattern = "scalefactors_json.json"),
        "spatial_params\\(sme\\) is NULL and no visiumDir")
    sp <- spatial_params(sme_filt)
    expect_false(is.null(sp))
    expect_equal(rownames(sp), c("sigmaOpt", "threshOpt"))
    expect_equal(unname(sp["sigmaOpt", ]), c(10, 10))
    expect_equal(unname(sp["threshOpt", ]), c(4, 4))
})
