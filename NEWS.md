# SpaceMarkers 2.3.3

## Breaking changes
* Directed overlap scores now use `target` / `source` column names in place of
`pattern` / `influence` (`relAbundance` is unchanged). `target` is the pattern
whose own GMM hotspot is being examined; `source` (prefixed `near.`) is the
pattern whose influence zone overlaps into that hotspot. Code that reads the
old `pattern` / `influence` columns from directed results must be updated.

## Bug fixes
* `SpaceMarkers()` on a hand-built `SpaceMarkersExperiment` no longer errors
with "spatial_params(sme) is NULL and no visiumDir is stored" — an SME
assembled without `load10X()` now falls back to computing default spatial
parameters instead of failing.
* `add_features()`: fixed a positional misalignment when the SME and the
features share the same barcodes in a different order. The NA mask and the
subsequent `colData` assignment were applied against the wrong spot, so values
were attached to the wrong barcodes and the wrong spot was dropped. Features
are now reordered to the SME's column order first.
* `load10X()` / `add_features()`: closed an NA-padding gap left by the earlier
spatial-pattern NA handling.
* Directed scores no longer mix numeric and character gene labels. When one
direction of a pattern pair did not have two groups, the all-NA fallback table
in `.calc_IM_scores()` was built without row names, so its `gene` column
became row numbers ("1", "2", ...) while the other direction used real gene
names. The fallback tables now carry gene names.
* `calculate_fc_score()`: reworked binning logic so `low_bins` and `high_bins`
can no longer come back empty.
* The parallel backend is no longer registered globally from inside a closure,
leaving a user's own registered `BiocParallel` backend intact.

## Internal
* `getSpatialParameters` and the image-creation path refactored; AnnData loader
improvements; missing documentation added for several functions.
* Tests skip gracefully when optional AnnData reader dependencies are absent.
* No changes to the exported API in this release (`NAMESPACE` is unchanged
since 2.3.1).

# SpaceMarkers 2.3.1

## New SpaceMarkersExperiment class
* Added `SpaceMarkersExperiment` (SME), an S4 class extending `SpatialExperiment` that carries hotspots, influence maps, and pattern metadata alongside expression and spatial coordinates.
* Added SME-aware methods and dispatch for `get_spatial_features`, `get_interacting_genes`, `get_pairwise_interacting_genes`, and related accessors so users can stay in a single Bioconductor object end-to-end.
* Re-exported `SpatialExperiment` / `SummarizedExperiment` accessors so SME workflows do not require attaching upstream packages.

## AnnData I/O
* Added `load_anndata()` to read `.h5ad` files directly into a `SpaceMarkersExperiment`.
* Added `save_anndata()` to write SME objects back to `.h5ad` with round-trip fidelity.
* Added direct `SingleCellExperiment` → `SpaceMarkersExperiment` coercion.

## Directed interaction analysis and visualization
* Added `overlap_map()` for per-spot directed interaction classification (3-level factor per direction).
* `plot_spatial()` gained a `source` argument (`colData`/`assay`/`hotspots`/`influence_map`/`interaction`) and directed interaction labels for three-way hotspot overlays.
* Hotspot plots in vignettes now use distinct colors per pattern.

## Vignettes
* New `SpaceMarkersExperiment_vignette` walks through the SME-first workflow.
* `SpaceMarkersStepByStep_vignette` restructured as a tutorial with per-step plots interleaved with narrative.
* Vignette data downloads now route through `BiocFileCache` instead of ad hoc temp files.
* Compact JPEG raster output keeps rendered HTML vignettes under Bioconductor's per-file size guideline.

## Bioconductor readiness
* Runnable examples and `\value` sections added across exported functions for `BiocCheck`.
* Tightened `.Rbuildignore` and `.gitignore` for a clean source tarball.
* Renamed `R/OneSpaceMarkers.R` → `R/SpaceMarkers.R` and regenerated `Collate` order.

## Performance
* Avoided dense coercion in the gene filter and skipped materializing the `spPatterns` data frame, reducing memory pressure on large samples.

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
