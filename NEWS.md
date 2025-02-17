# SpaceMarkers - unreleased

* Added `findAllHotspots` to use with `getPairwiseInteractingGenes`
* Return empty interacting genes list instead of failing with error
when no genes pass fdr threshold.

# SpaceMarkers 1.2.0

* Updated SpaceMarkersMetric by fixing signage and log transformed to scale
magnitude
* Added getSpatialParamsExternal which enables getting spatial parameters
from file or from the user.
* Deprecated getSpatialParameters
* Enabled includeSelf = TRUE in getInteractingGenes.R to improve hotspot 
detection
* Enabled load10XCoords to read coordinates from VisiumHD directory
* Optimized the long running row.dunn.test() function
* Corrected sparse -> dense conversions
* Added getPairwiseInteractingGenes which enables pairwise analysis of 
interacting patterns 
* `getSpatialFeatures`: add default method to infer the object passed to it. 

# SpaceMarkers 0.1.0

* Added a `NEWS.md` file to track changes to the package.
