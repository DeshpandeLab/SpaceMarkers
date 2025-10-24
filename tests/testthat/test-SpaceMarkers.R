hdpath <- "./assets/hdvisium/binned_outputs/square_016um/"
hdfeatures <- "./assets/hdvisium/hdvisium_cellular_compositions.csv"
sdpath <- "./assets/visium/"
sdfeatures <- "./assets/visium/visium_cellular_compositions.csv"

test_that("Directed Visium HD throws no error", {
    expect_no_error(
        sm <- SpaceMarkers(features = hdfeatures, data = hdpath, cpus = 4)
    )
})

test_that("UnDirected Visium HD throws no error", {
    expect_no_error(
        sm <- SpaceMarkers(features = hdfeatures, data = hdpath,
                           cpus = 4, directed = FALSE)
    )
})

test_that("UnDirected Visium SD throws no error", {
    expect_no_error(
        suppressWarnings(
        #1: In Smooth.ppp(X, at = "points", sigma = sigma[1], ...):
        #Bandwidth is close to zero: nearest-neighbour interpolation performed
            sm <- SpaceMarkers(features = sdfeatures, data = sdpath,
                               cpus = 4, directed = FALSE,
                               avoid_confounders = FALSE)
        )
    )
})

test_that("Directed Visium SD throws no error", {
    expect_no_error(
        sm <- SpaceMarkers(features = sdfeatures, data = sdpath,
                           cpus = 4, directed = TRUE,
                           avoid_confounders = FALSE)
    )
})