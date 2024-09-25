test_that('row dunn test works with both sparse and dense matrices',{
    dense <- matrix(rnorm(100), ncol = 10)
    sparse <- Matrix::Matrix(dense, sparse = TRUE)

    region <- matrix(sample(c(NA, "Interacting"), 10, replace = TRUE), ncol = 1)
    region <- factor(region)

    expect_true("dgCMatrix" %in% class(sparse))
    expect_error(row.dunn.test(sparse, region),
        "Argument 'x' must be a matrix or a vector")

    expect_true("matrix" %in% class(dense))
    expect_no_error(row.dunn.test(dense, region))
})