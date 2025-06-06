
data_dir <- "~/0_Projects/02_BTC/DPT/HDsample/binned_outputs/square_016um/"
coords <- load10XCoords(data_dir)
rownames(coords) <- coords$barcode

spPatterns <- getSpatialFeatures("~/Downloads/rctd_cell_types-2.csv")
barcodes <- intersect(rownames(coords), rownames(spPatterns))
spPatterns <- cbind(coords[barcodes,],spPatterns[barcodes,])
optParams <- getSpatialParameters(spPatterns,visiumDir=data_dir)
sigmaPair <- optParams[,1]
patnames <- setdiff(colnames(spPatterns),c("x","y","barcode"))


patthresholds <- calcAllThresholds(spPatterns)
patHotspots <- findAllHotspots.value(spPatterns, threshold=patthresholds)

spInfluence <- calcInfluence(spPatterns,optParams)
infthresholds <- calcAllThresholds(spInfluence)
infHotspots <- findAllHotspots.value(spInfluence, threshold=infthresholds)

#create a table of pattern pairs
patternPairs <- t(combn(patnames,2))

data <- load10XExpr(data_dir)
data <- data[,barcodes]

lrlist <- data(lrlist)
lrgenes <- Reduce(union,lrlist)
lrgenes <- intersect(lrgenes,rownames(data))
data <- data[lrgenes,]

# Calculate interaction scores for all pattern pairs
IMscores <- calcAllIMscores.HD(data, patHotspots, infHotspots, patternpairs)

data(lrdf)
lrpairs <- lrdf$interaction[,c("ligand.symbol","receptor.symbol")]
ligands <- sapply(lrpairs$ligand.symbol,function(i) strsplit(i,split=", "))
receptors <- sapply(lrpairs$receptor.symbol,function(i) strsplit(i,split=", "))

LScores <- getGeneSetScore(IMscores,genes = ligands)
RScores <- getGeneSetScore(IMscores,genes = receptors)

LRScores <- LScores + RScores
rownames(LRScores) <- rownames(lrpairs)
ind1 <- seq(1,ncol(LScores),by=2)
ind2 <- seq(2,ncol(LScores),by=2)
LR1 <- LScores[,ind1] + RScores[,ind2]
rownames(LR1) <- rownames(lrpairs)
LR2 <- LScores[,ind2] + RScores[,ind1]
rownames(LR2) <- rownames(lrpairs)

LRcomb <- cbind(LR1,LR2)
apply(LRcomb, 2, function(x) x[!is.na(x)])

# Select top entries from each column
top_n <- 10
top_entries <- lapply(colnames(LRcomb), function(cc) { x<- LRcomb[,cc];
    if (length(x) < top_n) {
        return(x[!is.na(x)])
    } else {
        cc <- gsub("t_","",cc)
        cc <- strsplit(cc, split = "_near_")
        source_cell <- cc[[1]][1]
        target_cell <- cc[[1]][2]
        x <- x[!is.na(x)]
        ind <- order(x,decreasing = TRUE)
        ind <- names(x[ind])
        lls <- lrpairs[ind[1:top_n],'ligand.symbol']
        lls <- gsub(", ","|", lls)
        rrs <- lrpairs[ind[1:top_n],'receptor.symbol']
        rrs <- gsub(", ","|", rrs)
        x <- x[ind]
        return(data.frame(score=x[1:top_n], interaction=ind[1:top_n],source_cell_type=source_cell,target_cell_type=target_cell,ligand=lls,receptor=rrs))
    }
})
topinteractions <- Reduce(rbind,top_entries)
gsub("G..M","G.M", topinteractions$source_cell_type) -> topinteractions$source_cell_type
gsub("G..M","G.M", topinteractions$target_cell_type) -> topinteractions$target_cell_type
target_cells <- c("FIBROBLASTS","CYCLING.DUCTAL", "DUCTAL", "TNK", "B.CELLS", "CYCLING.TNK", "MAST")
source_cell <- "CYCLING.MYELOID"

topinteractions %>% filter(source_cell_type == source_cell) %>%
    #filter(target_cell_type %in% target_cells) %>%
    group_by(target_cell_type) %>%
    summarise(score = sum(score)) %>%
    arrange(desc(score))
topinteractions %>% filter(source_cell_type == "FIBROBLASTS") %>%
    arrange(desc(score)) -> plotinteractions

plot.new() # If not in RStudio or if you get "R newBella page" error
plotCellInteractionCircos(topinteractions,
                            link_buffer_fraction = 0.05,
                            link_connection_rou = 0.85,
                            split_segments_for_links = T,
                            scale_link_width_by_score = TRUE,
                            score_transform_for_width = function(s) log1p(s) + 0.2, # Log transform for width
                            score_color_palette_fun = score_pal_func_test3,
                            inter_molecule_segment_gap = 0.1,
                            link_arrowhead_length=0.1,
                            link_arrowhead_width=0.2)
# Correct way to create a three-color gradient
score_pal <- circlize::colorRamp2(
  breaks = c(min(topinteractions$score), 
            median(topinteractions$score), 
            max(topinteractions$score)), # Vector of length 3
  colors = c("grey90", "orange", "red"),                 # Vector of length 3
  transparency = 0.4
)

target_cells <- c("CYCLING.DUCTAL", "DUCTAL", "TNK", "B.CELLS", "CYCLING.TNK", "MAST")
all_plot_cells <- c("FIBROBLASTS", target_cells)

custom_gaps <- setNames(rep(2, length(all_plot_cells)), all_plot_cells)
custom_gaps["MAST"] <- 20 

plotSourceToTargetCircos(
  topinteractions,
  source_cell_name = "FIBROBLASTS",
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

plotTargetFromSourcesCircos(
  topinteractions,
  target_cell_name = "FIBROBLASTS",
  source_cell_names = target_cells,
  
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

lgd_links = ComplexHeatmap::Legend(
  col_fun = score_pal, # Use the *exact same* color function
  title = "Interaction Score",
  direction = "vertical" # Can also be "horizontal"
)

ComplexHeatmap::draw(
  lgd_links,
  x = unit(1, "npc") - unit(2, "mm"), # 2mm from the right edge
  y = unit(2, "mm"),                  # 2mm from the bottom edge
  just = c("right", "bottom")         # Justification of the legend box
)
plotSourceToTargetCircos(topinteractions,
                        source_cell_name = "FIBROBLASTS",
                        target_cell_names = c("CYCLING.DUCTAL", "DUCTAL", "TNK", "B.CELLS", "CYCLING.TNK","MAST"), # Include autocrin
                        scale_link_width_by_score = TRUE,
                        #score_transform_for_width = function(s) log1p(s) + 0.2, # Log transform for width
                        score_color_palette_fun = score_pal_func_test3,
                        inter_molecule_segment_gap = 0.1,
                        link_connection_rou = 0.85,
                        split_segments_for_links = T,
                        link_buffer_fraction = 0.05,
                        link_arrowhead_length=0.05,
                        link_arrowhead_width=0.1)

plotTargetFromSourcesCircos(topinteractions,
                        target_cell_name = "FIBROBLASTS",
                        source_cell_names = c("CYCLING.DUCTAL", "DUCTAL", "TNK", "B.CELLS", "CYCLING.TNK","MAST"), # Include autocrin
                        scale_link_width_by_score = TRUE,
                        #score_transform_for_width = function(s) log1p(s) + 0.2, # Log transform for width
                        score_color_palette_fun = score_pal_func_test3,
                        inter_molecule_segment_gap = 0.1,
                        link_connection_rou = 0.85,
                        split_segments_for_links = T,
                        link_buffer_fraction = 0.05,
                        link_arrowhead_length=0.05,
                        link_arrowhead_width=0.1)
score_palette <- colorRamp2::colorRamp2(c(min(lr_sample_data$score), max(lr_sample_data$score)), c("lightyellow", "orangered"))
plotLRCircos(lr_sample_data, score_color_palette = score_palette, ligand_color = "lightgreen", receptor_color="lightblue")
plotLRCircos_v2(lr_sample_data, score_color_palette = score_palette)
