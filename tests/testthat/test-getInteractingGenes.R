test_that("getInteracting genes return empty interacting_genes object when no interacting genes found", {
  library(testthat)
  # Create some sample data and patterns for testing
  set.seed(123)
  cells <- c()
  cell_num <- 500
  gene_num <- 100
  for(i in 1:cell_num){
    cells[length(cells)+1] <- paste0("cell_",i)
  }
  data <- matrix(rnorm(5), ncol = gene_num,nrow = cell_num)
  rownames(data) <- cells
  colnames(data) <- paste0("gene",1:gene_num)
  
  spPatterns <- data.frame(barcode = cells,
                           y = runif(cell_num, min=0, max=cell_num),
                           x = runif(cell_num , min=0, max=cell_num ),
                           "Pattern_1" = runif(cell_num , min=0, max=1),
                           "Pattern_2" = runif(cell_num , min=0, max=1))
  
  optParams <- NULL
  output <- getInteractingGenes(data=data, spPatterns, reconstruction=NULL,
                                optParams = optParams,
                                spPatterns = spPatterns,
                                refPattern = "Pattern_1",
                                mode="DE",analysis="overlap")
  # Checkin result if no interacting genes found
  res<-evaluate_promise(output$interacting_genes)
  expect_equal(res$result, NULL)
  
})

test_that("test no interacting genes in nmf output", {
    expect_no_error(getSpaceMarkersMetric(list()))
})
