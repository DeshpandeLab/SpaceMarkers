test_that("create_lr_dataframe round-trips a small fixture", {
    lrpairs <- data.frame(
        ligand   = c("L1", "L2"),
        receptor = c("R1, R2", "R3"),
        row.names = c("L1_R1_R2", "L2_R3"),
        stringsAsFactors = FALSE)
    lr <- matrix(c(0.4, 0.7, 0.1, 0.2),
                 nrow = 2,
                 dimnames = list(rownames(lrpairs), c("A_to_B", "B_to_A")))
    lig <- matrix(c(0.5, 0.6, 0.2, 0.3),
                  nrow = 2,
                  dimnames = list(rownames(lrpairs), c("A_near_B", "B_near_A")))
    rec <- matrix(c(0.7, 0.8, 0.9, 0.6),
                  nrow = 2,
                  dimnames = list(rownames(lrpairs), c("A", "B")))

    out <- create_lr_dataframe(lr, lig, rec, lrpairs)

    expect_named(out, c("source_cell_type", "ligand", "receptor",
                        "target_cell_type", "ligand_score",
                        "receptor_score", "score", "interaction",
                        "source_to_target"))
    expect_equal(nrow(out), 4L)
    # complexes flattened with underscore
    expect_true("R1_R2" %in% out$receptor)
    # source_to_target == "<source>_to_<target>"
    expect_equal(out$source_to_target,
                 paste0(out$source_cell_type, "_to_", out$target_cell_type))
})

test_that("create_lr_dataframe preserves NA in score (no silent NA -> 0)", {
    lrpairs <- data.frame(ligand = "L1", receptor = "R1",
                          row.names = "L1_R1",
                          stringsAsFactors = FALSE)
    lr  <- matrix(c(0.5, NA), nrow = 1,
                  dimnames = list("L1_R1", c("A_to_B", "B_to_A")))
    lig <- matrix(c(0.4, 0.4), nrow = 1,
                  dimnames = list("L1_R1", c("A_near_B", "B_near_A")))
    rec <- matrix(c(0.7, 0.7), nrow = 1,
                  dimnames = list("L1_R1", c("A", "B")))

    out <- create_lr_dataframe(lr, lig, rec, lrpairs)
    expect_equal(sum(is.na(out$score)), 1L)
    expect_true(any(out$score > 0, na.rm = TRUE))
})
