test_that("getSpatialFeatures can read a CoGAPS object",{
    sf <- getSpatialFeatures("assets/cogaps.rds")
    expect_equal(ncol(sf), 3)
})

test_that("getSpatialFeatures can read a BayesTME object",{
    sf <- getSpatialFeatures("assets/btme.h5ad", method = "BayesTME")
    expect_equal(ncol(sf), 3)
})

test_that("getSpatialFeatures can read a Seurat object",{
    #to be defined

})