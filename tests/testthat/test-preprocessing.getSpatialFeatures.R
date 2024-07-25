test_that("getCogaps can read a CoGAPS object",{
    sf <- getCogaps("assets/cogaps.rds")
    expect_equal(ncol(sf), 3)
})

test_that("getBtme can read a BayesTME object",{
    sf <- getBtme("assets/btme.h5ad")
    expect_equal(ncol(sf), 3)
})

test_that("getSpatialFeatures fails with unsupported method",{
    expect_error(getSpatialFeatures("assets/cogaps.rds", method = "unsupported"))
})

test_that("getSpatialFeatures fails with no matching feature names",{
    expect_error(getSpatialFeatures("assets/cogaps.rds", method = "CoGAPS", featureNames = "no_match"),
                 "Regex no_match does not match any feature.")
})

test_that("getSpatialFeatures work with custom regex", {
    sf <- getSpatialFeatures("assets/cogaps.rds", method = "CoGAPS", featureNames = "_1")
    expect_equal(ncol(sf), 1)
})

test_that("getSpatialFeatures works with a feature name set", {
    sf <- getSpatialFeatures("assets/cogaps.rds", method = "CoGAPS", featureNames = c("Pattern_1", "Pattern_2"))
    expect_equal(ncol(sf), 2)
})

test_that("getSpatialFeatures warns if some of the features are not found", {
    expect_error(getSpatialFeatures("assets/cogaps.rds", method = "CoGAPS", 
                                    featureNames = c("Pattern_1", "Pattern_2", "Pattern_4")),
                                    "Some of the features were not found: Pattern_4")
})

test_that("getSpatialFeatures warns if some of the features are not found", {
    expect_error(getSpatialFeatures("assets/cogaps.rds", method = "CoGAPS",
                                    featureNames = c("Pattern_1", "Pattern_2", "Pattern_3", "qq", "ww")),
                                    "Some of the features were not found: qq ww")
})