% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getSpatialParameters.R
\name{getSpatialParamsMoransI}
\alias{getSpatialParamsMoransI}
\title{Calculate the optimal parameters from spatial kernel density for
cell-cell interactions}
\usage{
getSpatialParamsMoransI(spatialPatterns, ...)
}
\arguments{
\item{spatialPatterns}{A data frame that contains the spatial coordinates
for each cell type. The column names must include 'x' and 'y' as well as a
set of numbered columns named  'Pattern_1.....N'.}

\item{...}{Arguments passed to methods}
}
\value{
a numeric matrix of sigmaOpts - the optimal width of the gaussian
distribution, and the threshOpt - outlier threshold around the set of spots
for each pattern
}
\description{
This function uses Morans.I to calculate the optimal width of
the kernel density (sigmaOpt) as well as the outlier threshold around the set
of spots (threshOpt) for a null distribution.
}
\examples{
library(SpaceMarkers)
# Create test data
cells <- c()
test_num <- 500
for(i in 1:test_num){
    cells[length(cells)+1] <- paste0("cell_",i)
}
spPatterns <- data.frame(barcode = cells,
y = runif(test_num, min=0, max=test_num),
x = runif(test_num, min=0, max=test_num),
Pattern_1 = runif(test_num, min=0, max=1),
Pattern_2 = runif(test_num, min=0, max=1) )
# Call the getSpatialParamsMoransI function with the test data
optParams <- getSpatialParamsMoransI(spPatterns)

}
