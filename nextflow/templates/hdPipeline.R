#!/usr/bin/env Rscript
# @dimalvovs: template to be able to run SpaceMarkers as part of
# btc-spatial-pipelines while it is in dev and until it is added
# to the SpaceMarkers package
# run example:
# nextflow run nextflow/visiumhd.nf --input nextflow/hd-samplesheet.csv -profile docker -c nextflow/nextflow.config -resume


# script start
#' @title HD Pipeline for Cell-Cell Interactions
#' @description This script processes spatial transcriptomics data to identify and visualize cell-cell interactions using the HD method.
## @author Atul Deshpande
#' @date 2025-06-07
#' 
#' 
# Load necessary libraries

library("dplyr")
library(SpaceMarkers)
library(effsize)

data_dir <- "${data}"           # example: "sample1/binned_outputs/"
patternpath <- "${features}"    # example: "rctd_cell_types.csv" # 
output_dir <- "${prefix}"       # example: "hd_pipeline_output" # 
figure_dir <- file.path(output_dir, "figures")
set.seed(${params.seed})
useLigandReceptorGenes <- as.logical("${params.use_ligand_receptor_genes}") # limit to ligand-receptor genes
goodgeneThreshold <- ${params.good_gene_threshold} # limit to genes with high expression


dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, showWarnings = FALSE)

coords <- load10XCoords(data_dir)
rownames(coords) <- coords[["barcode"]]

spPatterns <- get_spatial_features(patternpath)
barcodes <- intersect(rownames(coords), rownames(spPatterns))
spPatterns <- cbind(coords[barcodes,],spPatterns[barcodes,])
optParams <- get_spatial_parameters(spPatterns,visiumDir=data_dir)
sigmaPair <- optParams[,1]
patnames <- setdiff(colnames(spPatterns),c("x","y","barcode"))

# Calculate hotspots and influence for each pattern
print("Calculating hotspots and influence for each pattern...")
patthresholds <- calculate_thresholds(spPatterns, minvals=0.1, maxvals=0.8)
pat_hotspots <- find_all_hotspots(spPatterns, threshold=patthresholds)

spInfluence <-  calculate_influence(spPatterns,optParams)
infthresholds <- calculate_thresholds(spInfluence, minvals=0.01, maxvals=0.5)
influence_hotspots <- find_all_hotspots(spInfluence, threshold=infthresholds)

#create a table of pattern pairs
pattern_pairs <- t(combn(patnames,2))

data <- load10XExpr(data_dir)
data <- data[,barcodes]

if (useLigandReceptorGenes) {
  print("Limiting data to ligand-receptor genes...")
    # Filter data to include only ligand-receptor genes
    library(CellChat)
    lrdf <- CellChat::CellChatDB.human
    lrpairs <- lrdf[["interaction"]][,c("ligand.symbol","receptor.symbol")]
    ligands <- sapply(lrpairs[["ligand.symbol"]],function(i) strsplit(i,split=", "))
    receptors <- sapply(lrpairs[["receptor.symbol"]],function(i) strsplit(i,split=", "))
    names(ligands) <- names(receptors) <- rownames(lrpairs)

    lrgenes <- union(unlist(ligands),unlist(receptors))
    data <- data[rownames(data) %in% lrgenes,]

    # Calculate interaction scores for all pattern pairs
    IMscores <- calculate_gene_scores_directed(data=data,
                                pat_hotspots=pat_hotspots,
                                influence_hotspots=influence_hotspots,
                                pattern_pairs=pattern_pairs,
                                avoid_confounders=TRUE)

    saveRDS(IMscores, file = sprintf("%s/IMscores.rds", output_dir))

    ligand_scores <- calculate_gene_set_score(IMscores,gene_sets = ligands, weighted = TRUE, method = "arithmetic_mean")
    receptor_scores <- calculate_gene_set_specificity(data, spPatterns, gene_sets=receptors, weighted = TRUE, method = "arithmetic_mean")
    lr_scores <- calculate_lr_scores(ligand_scores,receptor_scores,lr_pairs=lrpairs, method = "geometric_mean", weighted = TRUE)

    saveRDS(ligand_scores, file = sprintf("%s/ligand_scores.rds", output_dir))
    saveRDS(receptor_scores, file = sprintf("%s/receptor_scores.rds", output_dir))
    saveRDS(lr_scores, file = sprintf("%s/LRscores.rds", output_dir))
} else {
    # Use top "good" genes
    goodgenes <- apply(data,1,sum) |> sort(decreasing=TRUE) |> head(goodgeneThreshold) |> names()
    data <- data[goodgenes,]
    IMscores <- calculate_gene_scores_directed(data=data,
                                pat_hotspots=pat_hotspots,
                                influence_hotspots=influence_hotspots,
                                pattern_pairs=pattern_pairs,
                                avoid_confounders=TRUE)
    saveRDS(IMscores, file = sprintf("%s/IMscores.rds", output_dir))
}


#output versions
#versions
message("reading session info")
sinfo <- sessionInfo()
versions <- lapply(sinfo[["otherPkgs"]], function(x) {sprintf("  %s: %s",x[["Package"]],x[["Version"]])})
versions[['R']] <- sprintf("  R: %s
",packageVersion("base"))
cat(paste0("process",":
"), file="versions.yml")
cat(unlist(versions), file="versions.yml", append=TRUE, sep="
")
