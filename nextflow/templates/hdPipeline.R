#!/usr/bin/env Rscript
# @dimalvovs: template to be able to run SpaceMarkers as part of
# btc-spatial-pipelines while it is in dev and until it is added
# to she SpaceMarkers package
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
library("SpaceMarkers")

# read data_dir and patternpath form arguments
data_dir <- "${data}"           # example: /HDsample/binned_outputs/square_016um/
patternpath <- "${features}"    # example: /rctd_cell_types-2.csv
output_dir <- "${prefix}"       # example: /HDsample/hd_pipeline_output/
figure_dir <- file.path(output_dir, "figures")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, showWarnings = FALSE)

# limit the analysis to specific cell types
include_cells <- c("FIBROBLASTS","CYCLING.DUCTAL", "DUCTAL","MYELOID", "TNK", "B.CELLS", "CYCLING.TNK", "MAST", "CYCLING.MYELOID")
reference_cells <- c("FIBROBLASTS", "CYCLING.DUCTAL")

# limit the analysis to ligand-receptor genes
use.ligand.receptor.genes <- TRUE

# limit the analysis to genes with high expression
# goodgeneThreshold <- 1000

# identify the top n interactions to report for each cell pair (n=10)
top_n_table <- 25
top_n_plots <- 10
coords <- load10XCoords(data_dir)
rownames(coords) <- coords[["barcode"]]

spPatterns <- getSpatialFeatures(patternpath)
barcodes <- intersect(rownames(coords), rownames(spPatterns))
spPatterns <- cbind(coords[barcodes,],spPatterns[barcodes,])
optParams <- getSpatialParameters(spPatterns,visiumDir=data_dir)
sigmaPair <- optParams[,1]
patnames <- setdiff(colnames(spPatterns),c("x","y","barcode"))

# Calculate hotspots and influence for each pattern
print("Calculating hotspots and influence for each pattern...")
patthresholds <- calcAllThresholds(spPatterns)
patHotspots <- findAllHotspots.value(spPatterns, threshold=patthresholds)

spInfluence <- calcInfluence(spPatterns,optParams)
infthresholds <- calcAllThresholds(spInfluence)
infHotspots <- findAllHotspots.value(spInfluence, threshold=infthresholds)

#create a table of pattern pairs
patternPairs <- t(combn(patnames,2))

data <- load10XExpr(data_dir)
data <- data[,barcodes]

data(lrdf)
lrpairs <- lrdf[["interaction"]][,c("ligand.symbol","receptor.symbol")]
ligands <- sapply(lrpairs[["ligand.symbol"]],function(i) strsplit(i,split=", "))
receptors <- sapply(lrpairs[["receptor.symbol"]],function(i) strsplit(i,split=", "))

lrgenes <- union(unlist(ligands),unlist(receptors))

if (use.ligand.receptor.genes) {
  print("Limiting data to ligand-receptor genes...")
    # Filter data to include only ligand-receptor genes
    data <- data[rownames(data) %in% lrgenes,]
} else{
    # Use "good" genes
    goodgenes <- which(apply(data,1,sum)>goodgeneThreshold)
    data <- data[goodgenes,]
}

# Calculate interaction scores for all pattern pairs
IMscores <- calcAllIMscores.HD(data=data, 
                                patHotspots=patHotspots, 
                                infHotspots=infHotspots, 
                                patternPairs=patternPairs,
                                min_bins=100)

saveRDS(IMscores, file = sprintf("%s/IMscores.rds", output_dir))

LScores <- getGeneSetScore(IMscores,genes = ligands)
RScores <- getGeneSetScore(IMscores,genes = receptors)

ind1 <- seq(1,ncol(LScores),by=2)
ind2 <- seq(2,ncol(LScores),by=2)
LR1 <- LScores[,ind1] + RScores[,ind2]
rownames(LR1) <- rownames(lrpairs)
LR2 <- LScores[,ind2] + RScores[,ind1]
rownames(LR2) <- rownames(lrpairs)

LRcomb <- cbind(LR1,LR2)

saveRDS(LRcomb, file = sprintf("%s/LRscores.rds", output_dir))
# Select top entries from each column
top_entries <- lapply(colnames(LRcomb), function(cc) { x<- LRcomb[,cc];
    if (length(x) < top_n_table) {
        return(x[!is.na(x)])
    } else {
        cc <- gsub("t_","",cc)
        cc <- strsplit(cc, split = "_near_")
        source_cell <- cc[[1]][1]
        target_cell <- cc[[1]][2]
        x <- x[!is.na(x)]
        ind <- order(x,decreasing = TRUE)
        ind <- names(x[ind])
        lls <- lrpairs[ind[1:top_n_table],'ligand.symbol']
        lls <- gsub(", ","|", lls)
        rrs <- lrpairs[ind[1:top_n_table],'receptor.symbol']
        rrs <- gsub(", ","|", rrs)
        x <- x[ind]
        return(data.frame(score=x[1:top_n_table], interaction=ind[1:top_n_table],source_cell_type=source_cell,target_cell_type=target_cell,ligand=lls,receptor=rrs))
    }
})
topinteractions <- Reduce(rbind,top_entries)

# order the interactions by score

topinteractions <- topinteractions %>%
    dplyr::arrange(desc(score))

# Save the top interactions to a CSV file
write.csv(topinteractions, 
          file = sprintf("%s/top_%s_pairwise_LR_interactions.csv", output_dir, top_n_table),
          row.names = FALSE)

score_pal <- circlize::colorRamp2(
  breaks = c(min(topinteractions[["score"]]), 
            median(topinteractions[["score"]]), 
            max(topinteractions[["score"]])), # Vector of length 3
  colors = c("grey90", "orange", "red"),                 # Vector of length 3
  transparency = 0.4
)

topNinteractions <- topinteractions %>% 
filter(source_cell_type %in% include_cells | 
target_cell_type %in% include_cells) %>%
    group_by(source_cell_type,target_cell_type) %>% #limit to top 5 interactions per pair
    arrange(desc(score)) %>%
    slice_head(n = top_n_plots) %>%
    ungroup()


summarized_interactions <- topNinteractions %>% 
    summarise(score = sum(score)) %>%
    arrange(desc(score))
png_filename <- paste0("circos_top_",top_n_plots,".png")
png_path <- file.path(figure_dir, png_filename)
png(png_path, width = 2000, height = 2000, res = 300)
plot.new() # If not in RStudio or if you get "R newBella page" error
  plotCellInteractionCircos(topNinteractions,
                            link_buffer_fraction = 0.05,
                            link_connection_rou = 0.85,
                            split_segments_for_links = T,
                            scale_link_width_by_score = TRUE,
                            score_transform_for_width = function(s) log1p(s) + 0.2, # Log transform for width
                            score_color_palette_fun = score_pal,
                            inter_molecule_segment_gap = 0.1,
                            link_arrowhead_length=0.1,
                            link_arrowhead_width=0.2)
dev.off()

for(source_cell in reference_cells) {
  print(paste("Processing source cell type:", source_cell))
  target_cells <- setdiff(include_cells, source_cell)
  # Filter interactions for the current source cell type
  filtered_interactions <- topNinteractions %>% 
                            filter(source_cell_type == source_cell) %>% 
                            filter(target_cell_type %in% target_cells)
  
  # If there are no interactions, skip to the next iteration
  if (nrow(filtered_interactions) == 0) {
    next
  }

  score_pal <- circlize::colorRamp2(
    breaks = c(min(filtered_interactions[["score"]]), 
    median(filtered_interactions[["score"]]), 
    max(filtered_interactions[["score"]])), # Vector of length 3
  colors = c("grey90", "orange", "red"),                 # Vector of length 3
  transparency = 0.4
  )
  all_plot_cells <- c(source_cell, target_cells)
  custom_gaps <- setNames(rep(2, length(all_plot_cells)), all_plot_cells)
  custom_gaps[length(custom_gaps)] <- 20 
  # Create a plot for the current source cell type
  png_filename <- paste0("circos_interactions_", source_cell, "_source.png")
  png_path <- file.path(figure_dir, png_filename)
  png(png_path, width = 2000, height = 2000, res = 300)
  par(mar = c(3, 3, 3, 3)) # Provides 2 lines of margin on all sides
  plot.new() # If not in RStudio or if you get "R newBella page" error
  plotSourceToTargetCircos(
  filtered_interactions,
  source_cell_name = source_cell,
  target_cell_names = target_cells,
  # Balance sector sizes by transforming score for width
  scale_link_width_by_score = TRUE,
  score_transform_for_width = function(s) log1p(s) + 0.2, # Key change!
  
  # Improve label readability
  molecule_label_cex = 0.5,
  inter_molecule_segment_gap = 0.15,
  track_height_molecules = 0.12,
  
  # Enhance link and layout aesthetics
  score_color_palette_fun = score_pal,
  gap_degree_after_sector = custom_gaps,
  link_connection_rou = 0.8, # Slightly more curve
  link_transparency = 0.4, # Set to 0 if transparency is in the colorRamp2 function

  # Keep other successful parameters
  split_segments_for_links = TRUE,
  link_buffer_fraction = 0.05,
  link_arrowhead_length = 0.05,
  link_arrowhead_width = 0.1
)
  dev.off()
}


for (target_cell in reference_cells) {
  print(paste("Processing target cell type:", target_cell))
  source_cells <- setdiff(include_cells, target_cell)
  
  # Filter interactions for the current target cell type
  filtered_interactions <- topNinteractions %>% 
  filter(target_cell_type == target_cell) %>%
  filter(source_cell_type %in% source_cells)
  
  # If there are no interactions, skip to the next iteration
  if (nrow(filtered_interactions) == 0) {
    next
  }
  
  score_pal <- circlize::colorRamp2(
    breaks = c(min(filtered_interactions[["score"]]), 
              median(filtered_interactions[["score"]]), 
              max(filtered_interactions[["score"]])), # Vector of length 3
    colors = c("grey90", "orange", "red"),                 # Vector of length 3
    transparency = 0.4
  )
 
  png_filename <- paste0("circos_interactions_", target_cell, "_target.png")
  png_path <- file.path(figure_dir, png_filename)
  png(png_path, width = 2000, height = 2000, res = 300)
  par(mar = c(1, 1, 1, 1)) # Provides 2 lines of margin on all sides

  all_plot_cells <- c(source_cells, target_cell)
  custom_gaps <- setNames(rep(2, length(all_plot_cells)), all_plot_cells)
  custom_gaps[length(custom_gaps)] <- 20 # Increase gap after the last sector

  plotTargetFromSourcesCircos(
    topinteractions,
    target_cell_name = target_cell,
    source_cell_names = source_cells,
    # Balance sector sizes by transforming score for width
    scale_link_width_by_score = TRUE,
    score_transform_for_width = function(s) log1p(s) + 0.2, # Key change!
    
    # Improve label readability
    molecule_label_cex = 0.5,
    inter_molecule_segment_gap = 0.15,
    track_height_molecules = 0.12,
    
    # Enhance link and layout aesthetics
    score_color_palette_fun = score_pal,
    gap_degree_after_sector = custom_gaps,
    link_connection_rou = 0.9, # Slightly more curve
    link_transparency = 0.4, # Set to 0 if transparency is in the colorRamp2 function

    # Keep other successful parameters
    split_segments_for_links = TRUE,
    link_buffer_fraction = 0.05,
    link_arrowhead_length = 0.05,
    link_arrowhead_width = 0.1
  )
  dev.off()
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