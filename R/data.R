#' Optimal paramters of 5 patterns from CoGAPS.
#'
#' A dataset with the optimal width of the gaussian distribution (sigmaOpt) and 
#' the outlier threshold around the set of spots (thresOpt) for each pattern
#' obtained from CoGAPS. CoGAPS was ran on spatial transcriptomic data from a 
#' breast cancer sample.
#' 
#' 
#' @name optParams
#' @format A data frame with 2 rows and 5 columns:
#' \describe{
#'   \item{Pattern_1}{immune cell pattern paramters }
#'   \item{Pattern_2}{Disp.1 parameters}
#'   \item{Pattern_3}{intraductal carcinoma (DCIS) parameters }
#'   \item{Pattern_2}{Disp.2 parameters }
#'   \item{Pattern_5}{invasive carcinoma lesion pattern paramters }
#'
#' }
#' @return A matrix of optimal parameters for patterns identified by CoGAPS
NULL

#' Latent Feature Space for each pattern
#' 
#' A CoGAPS object where the major requirements for SpaceMarkers are the
#' matrices of genes, barcodes and patterns learned of the latent-feature space 
#' 
#' @name cogaps_result
#' @format CogapsResult object with 24228 features and 6 samples:
#' \describe{
#'   \item{featureLoadings}{Data frame of Gene for each pattern}
#'   \item{sampleFactors}{Data frame of cell barcodes and the 5 patterns}
#' }
#' @return A matrix of statistics for each pattern across each barcode
NULL

#' Curated Genes for example purposes
#' 
#' A vector with genes selected based on previous runs of SpaceMarkers on the
#' Visium 10x breast ductal carcinoma spatial transcriptomics dataset
#'
#' @name curated_genes
#' @format A vector with 114 pre-selected genes
#' @return a vector of genes
NULL
