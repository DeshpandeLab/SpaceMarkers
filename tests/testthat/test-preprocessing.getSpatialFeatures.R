#reading input files and extracting features

test_that(".read_format can read cogaps and btme files", {
    expect_no_error(.read_format("assets/cogaps.rds"))
    expect_no_error(.read_format("assets/btme.h5ad"))
})

test_that(".get_cogaps_features can read features from a CoGAPS object",{
    cg <- .read_format("assets/cogaps.rds")
    sf <- .get_cogaps_features(cg)
    expect_equal(ncol(sf), 3)
})

test_that(".get_BTME_features can read features from Anndata object",{
    bt <- .read_format("assets/btme.h5ad")
    sf <- .get_BTME_features(bt)
    expect_equal(ncol(sf), 3)
})

test_that(".infer_method can infer the method used to obtain spatial features", {
    expect_equal(.infer_method(.read_format("assets/cogaps.rds"), NULL), "CoGAPS")
    expect_equal(.infer_method(.read_format("assets/btme.h5ad"), NULL), "BayesTME")
})

#main function tests
test_that("get_spatial_features fails with unsupported method",{
    expect_error(get_spatial_features("assets/cogaps.rds", method = "unsupported"))
})

test_that("get_spatial_features fails with no matching feature names",{
    expect_error(get_spatial_features("assets/cogaps.rds", method = "CoGAPS", featureNames = "no_match"),
                 "Regex no_match does not match any feature.")
})

test_that("get_spatial_features work with custom regex", {
    sf <- get_spatial_features("assets/cogaps.rds", method = "CoGAPS", featureNames = "_1")
    expect_equal(ncol(sf), 1)
})

test_that("get_spatial_features works with a feature name set", {
    sf <- get_spatial_features("assets/cogaps.rds", method = "CoGAPS", featureNames = c("Pattern_1", "Pattern_2"))
    expect_equal(ncol(sf), 2)
})

test_that("get_spatial_features warns if some of the features are not found", {
    expect_error(get_spatial_features("assets/cogaps.rds", method = "CoGAPS", 
                                    featureNames = c("Pattern_1", "Pattern_2", "Pattern_4")),
                                    "Some of the features were not found: Pattern_4")
})

test_that("get_spatial_features warns if some of the features are not found", {
    expect_error(get_spatial_features("assets/cogaps.rds", method = "CoGAPS",
                                    featureNames = c("Pattern_1", "Pattern_2", "Pattern_3", "qq", "ww")),
                                    "Some of the features were not found: qq ww")
})

test_that("get_spatial_features works in infer mode with Cogaps object", {
    sf <- get_spatial_features("assets/cogaps.rds")
    expect_equal(ncol(sf), 3)
})

test_that("get_spatial_features works in infer mode with Cogaps object", {
    sf <- get_spatial_features("assets/btme.h5ad")
    expect_equal(ncol(sf), 3)
})

test_that(".complete_pattern_rows flags and reports incomplete rows", {
    df <- data.frame(P1 = c(1, NA, 3), P2 = c(1, 2, NA))
    expect_message(keep <- .complete_pattern_rows(df),
                   "2 NA rows out of 3 total rows removed")
    expect_equal(keep, c(TRUE, FALSE, FALSE))
})

test_that(".complete_pattern_rows is silent when there are no NAs", {
    df <- data.frame(P1 = c(1, 2, 3), P2 = c(1, 2, 3))
    expect_no_message(keep <- .complete_pattern_rows(df))
    expect_true(all(keep))
})

test_that("get_spatial_features drops NA rows (e.g. RCTD-unestimated barcodes)
    and messages the count", {
    temp <- tempfile(fileext = ".csv")
    on.exit(unlink(temp))
    write.csv(data.frame(barcode = c("A1", "B2", "C3"),
                         P1 = c(1, NA, 3), P2 = c(1, 2, 3)),
             temp, row.names = FALSE)
    expect_message(sf <- get_spatial_features(temp, method = "CSV"),
                   "1 NA row out of 3 total rows removed")
    expect_equal(rownames(sf), c("A1", "C3"))
})