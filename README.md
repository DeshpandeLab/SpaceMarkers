<p align="center">
  <img src="SpaceMarkersHexWhite.png" width="200" title="SpaceMarkers hex logo">
</p>

# SpaceMarkers
An R/Bioconductor software tool to identify genes associated with latent space interactions in spatial transcriptomics.

## Citation
If you use the SpaceMarkers software please cite:

Atul Deshpande, Melanie Loth, et al.,
[Uncovering the spatial landscape of molecular interactions within the tumor microenvironment through latent spaces](https://doi.org/10.1016/j.cels.2023.03.004).
*Cell Systems,* April 2023. https://doi.org/10.1016/j.cels.2023.03.004

## Installation
You can install the latest Bioconductor release from SpaceMarkers using the code below.
```
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("SpaceMarkers")
```
You can install the development version of SpaceMarkers directly from the Github source.
```
install.packages("remotes")
remotes::install_github("DeshpandeLab/SpaceMarkers", dependencies = TRUE, build_vignettes = TRUE)
```
