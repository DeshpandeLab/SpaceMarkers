# SpaceMarkers - unreleased

# SpaceMarkers 2.0.0
* Added support for directed cell-cell interaction (see calculate_gene_scores_directed, calculate_influence)
* Functions for supporting ligand-receptor interactions based on 
directed cell-cell interactions
* Added ligand–receptor and gene-set scoring utilities: calculate_lr_scores,
calculate_gene_set_score, and calculate_gene_set_specificity.
* Added support for handling Visium HD
* plotting functions using circlize for showing ligand-receiver interactions.
* API and naming standardization (backwards-incompatible) Major public functions  renamed to snake_case for consistency (examples: calcInfluence → calculate_influence, getInteractingGenes → get_interacting_genes, getPairwiseInteractingGenes → get_pairwise_interacting_genes, getSpatialFeatures → get_spatial_features, getSpatialParameters → get_spatial_parameters)

# SpaceMarkers 1.5.0
* Added `findAllHotspots` to use with `getPairwiseInteractingGenes`
* Return empty interacting genes list instead of failing with error
when no genes pass fdr threshold.

# SpaceMarkers 1.2.0

* Updated SpaceMarkersMetric by fixing signage and log transformed to scale
magnitude
* Added get_spatial_paramsExternal which enables getting spatial parameters
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
