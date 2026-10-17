test_that(".calc_IM_scores labels fallback (single-group) rows with gene names", {
    set.seed(1)
    n <- 300
    bc <- paste0("s", seq_len(n))
    genes <- paste0("GENE", 1:6)
    data <- matrix(rnorm(length(genes) * n), nrow = length(genes),
                   dimnames = list(genes, bc))

    mk <- function(a_idx, b_idx) {
        df <- data.frame(barcode = bc, y = seq_len(n), x = seq_len(n),
                         A = NA_character_, B = NA_character_)
        df$A[a_idx] <- "A"
        df$B[b_idx] <- "B"
        df
    }
    pat <- mk(1:100, 101:250)
    # A's influence reaches into B's hotspot; B's influence reaches nothing.
    # So direction B_near_A has two groups, but A_near_B has only one and
    # takes the all-NA fallback table.
    inf <- mk(101:175, integer(0))

    res <- .calc_IM_scores(data, pat, inf, c("A", "B"))

    expect_equal(nrow(res), 2L * length(genes))
    expect_equal(sort(unique(res$cell_interaction)),
                 c("A_near_B", "B_near_A"))

    # Every row, including the fallback direction, is labelled with a real
    # gene name -- never a row number like "1", "2", ...
    expect_type(res$gene, "character")
    expect_false(any(grepl("^[0-9]+$", res$gene)))
    expect_setequal(res$gene, genes)
    expect_equal(sort(res$gene[res$cell_interaction == "A_near_B"]),
                 sort(genes))

    # The fallback direction is all NA by design; the other has real values.
    expect_true(all(is.na(res$statistic[res$cell_interaction == "A_near_B"])))
    expect_false(any(is.na(res$statistic[res$cell_interaction == "B_near_A"])))
})
