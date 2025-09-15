#!/usr/bin/env Rscript
library("dplyr")
library("SpaceMarkers")

# limit the analysis to specific cell types
include_cells <- c("FIBROBLASTS", "PDAC", "MYELOID", "TNK", "B.CELLS", "CYCLING.TNK", "MAST", "CYCLING.MYELOID")
reference_cells <- c("FIBROBLASTS", "PDAC")

# identify the top n interactions to report for each cell pair (n=10)
top_n_table <- 25
top_n_plots <- 10

LRcomb <- readRDS("${lrscores}")  # example: "sample1/binned_outputs/hdLR_combined.rds"

if (all(dim(LRcomb)==0)){
    message("No ligand-receptor interactions found. Exiting script.")
    quit(status = 0)
}

output_dir <- "${prefix}"         # example: "hd_pipeline_output" #
figure_dir <- file.path(output_dir, "figures")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, showWarnings = FALSE)

data(lrdf)
lrpairs <- lrdf[["interaction"]][,c("ligand.symbol","receptor.symbol")]

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
