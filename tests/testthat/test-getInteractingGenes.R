test_that("test no interacting genes in nmf output", {
    expect_no_error(getSpaceMarkersMetric(list()))
})