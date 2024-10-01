test_that('row dunn test works with dense matrix but not sparse',{
    dense <- matrix(rnorm(100), ncol = 10)
    sparse <- Matrix::Matrix(dense, sparse = TRUE)

    region <- matrix(sample(c(NA, "Interacting"), 10, replace = TRUE), ncol = 1)
    region <- factor(region)

    patnames <- levels(region)[which(levels(region)!="Interacting")]
    
    pattern1 <- patnames[1]
    pattern2 <- patnames[2]

    expect_true("dgCMatrix" %in% class(sparse))
    expect_error(row.dunn.test(sparse, region, pattern1, pattern2),
        "Argument 'x' must be a matrix or a vector")

    expect_true("matrix" %in% class(dense))
    expect_no_error(row.dunn.test(dense, region, pattern1, pattern2))
})