#' Interacting Genes including the KW-test genes
#'
#' A dataframe with:
#' Interactions between Pattern_1 and Pattern_5 from CoGAPS and
#' statistics about the interacting genes associated with them.
#'  
#' @name residual_p1_p5_enrichment_df
#' @format A dataframe of 21331 rows and 18 columns:
#' \describe{
#'   \item{Gene}{interacting genes}
#'   \item{Pattern_1 x Pattern_5}{vsBoth,vsPattern_1,vsPattern5,FALSE}
#'   \item{KW stats}{...}
#'   \item{Dunns stats}{...}
#'   \item{SpaceMarkersMetric}{|Dunn.zP1|x|Dunn.zP2|x2^-min(Dunn.zP1,Dunn.zP2)}
#' }
#' @return Dataframe summarizing latent feature interactions
NULL

#' Interacting Genes including the KW-test genes
#'
#' A dataframe with:
#' Interactions between Pattern_1 and Pattern_5 from CoGAPS and
#' statistics about the interacting genes associated with them.
#'  
#' @name DE_p1_p5_enrichment_df
#' @format A dataframe of 21145 rows and 18 columns:
#' \describe{
#'   \item{Gene}{interacting genes}
#'   \item{Pattern_1 x Pattern_5}{vsBoth,vsPattern_1,vsPattern5,FALSE}
#'   \item{KW stats}{...}
#'   \item{Dunns stats}{...}
#'   \item{SpaceMarkersMetric}{|Dunn.zP1|x|Dunn.zP2|x2^-min(Dunn.zP1,Dunn.zP2)}
#' }
#' @return Dataframe summarizing latent feature interactions
NULL

#' Interacting Genes
#'
#' A dataframe with:
#' Interactions between Pattern_1 and Pattern_5 from CoGAPS and
#' statistics about the interacting genes associated with them.
#'  
#' @name DE_p1_p5_overlap_df
#' @format A dataframe of 2799 rows and 18 columns:
#' \describe{
#'   \item{Gene}{interacting genes}
#'   \item{Pattern_1 x Pattern_5}{vsBoth,vsPattern_1,vsPattern5,FALSE}
#'   \item{KW stats}{...}
#'   \item{Dunns stats}{...}
#'   \item{SpaceMarkersMetric}{|Dunn.zP1|x|Dunn.zP2|x2^-min(Dunn.zP1,Dunn.zP2)}
#' }
#' @return Dataframe summarizing latent feature interactions
NULL

#' Interacting Genes for each Pattern.
#'
#' A list with:
#' 1. A list of 4 each containing a data frame of Pattern interactions from 
#' CoGAPS and statistics about the interacting genes associated with them.
#' 
#' 2. A matrix of the hotspot regions from CoGAPS for each spot
#'  
#' @name SpaceMarkers_Residualmode_overlap
#' @format A list of a list of 4 data frames and 1 matrix:
#' \describe{
#'   \item{interacting genes list of 4 dfs}{genes,CoGAPS patterns and stats}
#'   \item{hotspot Regions matrix}{CoGAPS hotspots for each spot}
#' }
#' @return List of data frames summarizing latent feature interactions
NULL

#' Optimal paramters of 5 patterns from CoGAPS.
#'
#' A dataset with the optimal width of the gaussian distribution (sigmaOpt) and 
#' the outlier threshold around the set of spots (thresOpt) for each pattern
#' obtained from CoGAPS. CoGAPS was ran on spatial transcriptomic data from a 
#' breast cancer sample.
#' 
#' 
#' @name optParams_breast_cancer
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



