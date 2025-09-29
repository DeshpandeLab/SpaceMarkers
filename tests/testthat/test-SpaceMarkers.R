
test_that("Directed Visium HD throws no error", {
    expect_no_error(
        sm <- SpaceMarkers(features='./assets/hdvisium/hdvisium_cellular_compositions.csv',
                    data='./assets/hdvisium/binned_outputs/square_016um/', cpus=4)
    )
})

test_that("UnDirected Visium HD throws no error", {
    expect_no_error(
        sm <- SpaceMarkers(features='./assets/hdvisium/hdvisium_cellular_compositions.csv',
                    data='./assets/hdvisium/binned_outputs/square_016um/', cpus=4,
                    directed = FALSE)
    )
})

test_that(" Visium SD throws no error", {
    expect_no_error(
        suppressWarnings(
        #1: In Smooth.ppp(X, at = "points", sigma = sigma[1], ...):
        #Bandwidth is close to zero: nearest-neighbour interpolation performed
            sm <- SpaceMarkers(features='./assets/visium/visium_cellular_compositions.csv',
                         avoid_confounders=FALSE, data='./assets/visium', cpus=4)
        )
    )
})

test_that("UnDirected Visium SD throws no error", {
    expect_no_error(
        sm <- SpaceMarkers(features='./assets/hdvisium/hdvisium_cellular_compositions.csv',
                    data='./assets/hdvisium/binned_outputs/square_016um/', cpus=4,
                    directed = FALSE)
    )
})