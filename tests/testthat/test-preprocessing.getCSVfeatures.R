test_that("getCSVfeatures works with a barcode column", {
    obj <- data.frame(barcode = c("A1", "B2"), x = c(1, 2), y = c(3, 4))
    sf <- .getCSVfeatures(obj)
    expect_equal(rownames(sf), c("A1", "B2"))
})

test_that("getCSVfeatures works without a barcode column", {
    #simulate empty first colname
    temp <- tempfile(fileext = ".csv")
    write.csv(data.frame(x = c(1, 2), y = c(3, 4)), temp, row.names = TRUE)
    obj <- read.csv(temp)
    sf <- .getCSVfeatures(obj)
    expect_equal(rownames(sf), c("1", "2"))
})

test_that("getCSVfeatures fails without a barcode column or rownames", {
    obj <- data.frame(x = c(1, 2), y = c(3, 4))
    expect_error(.getCSVfeatures(obj), "No barcode column")
})
