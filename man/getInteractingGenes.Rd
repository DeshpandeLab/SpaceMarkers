% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getInteractingGenes.R
\name{getInteractingGenes}
\alias{getInteractingGenes}
\title{Calculate Interaction Regions and Associated Genes}
\usage{
getInteractingGenes(
  data,
  spPatterns,
  refPattern = "Pattern_1",
  mode = c("DE", "residual"),
  optParams = NULL,
  reconstruction = NULL,
  hotspots = NULL,
  analysis = c("enrichment", "overlap"),
  minOverlap = 50,
  ...
)
}
\arguments{
\item{data}{original spatial data matrix.}

\item{spPatterns}{A data frame that contains the spatial coordinates
and metrics for spatial features (cell types/cell processes). The column
names must include 'x' and 'y' as well as the spatially varying features.}

\item{refPattern}{a character string that specifies the pattern whose
"interaction" with every other pattern we want
to study. The default value is "Pattern_1".}

\item{mode}{SpaceMarkers mode of operation. Possible values are
"DE" (the default) or "residual".}

\item{optParams}{a matrix with dimensions 2 X N, where N is the number
of spatial patterns with optimal parameters. The first row contains the
kernel width 'sigmaOpt' for each pattern, and the second row is the
threshOpt (outlier threshold) for each pattern. Users can also input their
preferred param values. The default value is NULL.}

\item{reconstruction}{reconstruction of the data matrix from latent
spaces. Required for "residual" mode.}

\item{hotspots}{a vector that specifies the patterns to compare
to the 'refPattern'. The default is NULL which indicates that all patterns
would be compared to the 'refPattern'.}

\item{analysis}{a character string that specifies the type of downstream
analysis to be performed. Possible values are "enrichment" (default)
and "overlap". In enrichment mode, all genes are returned, ranked by
the SpaceMarkers metric. In overlap mode, only the genes which are
significantly overexpressed in the interaction region are returned.}

\item{minOverlap}{a number that specifies the minimum overlap between
genes in two patterns to be considered for the statistical tests.
The default is 50.}

\item{...}{Arguments passed to methods}
}
\value{
a list of data frames with information about the interacting genes
of the refPattern and each latent feature pattern matrix
(interacting_genes object). There is also a data frame with all of the
regions of influence for any two of patterns (the hotspots object).
}
\description{
This function calculates statistically significant genes using a
non-parametric Kruskal-Wallis test for genes in any one region
of influence and a post hoc Dunn's test is used for analysis of
genes between regions.
}
\examples{
library(SpaceMarkers)
#Visium data links
urls <- read.csv(system.file("extdata","visium_data.txt",
package="SpaceMarkers",mustWork = TRUE))
counts_url <- urls[["visium_url"]][1]
sp_url <- urls[["visium_url"]][2]
#Remove present Directories if any
unlink(basename(sp_url))
unlink("spatial", recursive = TRUE)
files <- list.files(".")[grepl(basename(counts_url),list.files("."))]
unlink(files)
download.file(counts_url,basename(counts_url), mode = "wb")
counts_matrix<-load10XExpr(visiumDir=".",h5filename = basename(counts_url))
#Obtaining CoGAPS Patterns
cogaps_result <- readRDS(system.file("extdata","CoGAPS_result.rds",
package="SpaceMarkers",mustWork = TRUE))
features <- intersect(rownames(counts_matrix),rownames(
    slot(cogaps_result,"featureLoadings")))
barcodes <- intersect(colnames(counts_matrix),rownames(
    slot(cogaps_result,"sampleFactors")))
counts_matrix <- counts_matrix[features,barcodes]
cogaps_matrix <- slot(cogaps_result,"featureLoadings")[features,]\%*\%
    t(slot(cogaps_result,"sampleFactors")[barcodes,])
#Obtaining Spatial Coordinates
download.file(sp_url, basename(sp_url), mode = "wb")
untar(basename(sp_url))
spCoords <- load10XCoords(visiumDir = ".", version = "1.0")
rownames(spCoords) <- spCoords$barcode
spCoords <- spCoords[barcodes,]
spPatterns <- cbind(spCoords,slot(cogaps_result,
"sampleFactors")[barcodes,])
data("curated_genes")
spPatterns<-spPatterns[c("barcode","y","x","Pattern_1","Pattern_5")]
counts_matrix <- counts_matrix[curated_genes,]
cogaps_matrix <- cogaps_matrix[curated_genes, ]
data("optParams")
SpaceMarkersMode <- "DE"
ref_Pattern <- "Pattern_1"
SpaceMarkers_test <- getInteractingGenes(
    data=counts_matrix,reconstruction=NULL,
    optParams = optParams,
    spPatterns = spPatterns,
    refPattern = "Pattern_1",
    mode="DE",analysis="overlap")
#Remove present Directories if any
unlink(basename(sp_url))
unlink("spatial", recursive = TRUE)
files <- list.files(".")[grepl(basename(counts_url),list.files("."))]
unlink(files)

}
\seealso{
Other getIntGenes: 
\code{\link{findPatternHotspots}()},
\code{\link{getPairwiseInteractingGenes}()}
}
\concept{getIntGenes}
