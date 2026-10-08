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
lr_reference_path <- "${params.lr_reference}" # optional CSV overriding the CellChat ligand-receptor reference below
workers <- min(as.numeric("${task.cpus}"), parallel::detectCores())

BiocParallel::register(BiocParallel::MulticoreParam(workers = workers)) #register backend

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, showWarnings = FALSE)

# ---- Setup: build the SpaceMarkersExperiment in one call ----
# resolution = "fullres" matches this script's previous behavior: the old
# load10XCoords(data_dir) call never specified a resolution, and "fullres"
# is load10XCoords()'s own default.
sme <- load10X(
  visiumDir  = data_dir,
  features   = patternpath,
  method     = "CSV",
  resolution = "fullres",
  version    = "HD"
)

# ---- Gene restriction, and the ligand-receptor reference if requested ----
# Same useLigandReceptorGenes toggle and behavior as before, just computed
# up front so it can be passed straight into the one-line directed call
# below via genes=/lr_pairs= instead of manually subsetting a data matrix.
# If params.lr_reference points to an existing CSV, it overrides the
# CellChatDB.human-derived lrpairs below; otherwise CellChatDB.human is
# used exactly as before.
lrpairs <- NULL
if (useLigandReceptorGenes) {
  print("Limiting data to ligand-receptor genes...")
  if (nzchar(lr_reference_path) && file.exists(lr_reference_path)) {
    message(
      "Using user-supplied lr_reference: ", lr_reference_path, ". Expected ",
      "format: a CSV with ligand-receptor pair IDs as row names (first ",
      "column) and columns 'ligand.symbol'/'receptor.symbol' (or ",
      "'ligand'/'receptor', which will be renamed). Each entry may be a ",
      "single gene symbol or multiple symbols separated by ', ' (e.g. a ",
      "receptor complex) -- the same convention CellChatDB.human uses below."
    )
    lrpairs <- read.csv(lr_reference_path, row.names = 1)
    if (!all(c("ligand.symbol", "receptor.symbol") %in% colnames(lrpairs))) {
      if (all(c("ligand", "receptor") %in% colnames(lrpairs))) {
        lrpairs[["ligand.symbol"]]   <- lrpairs[["ligand"]]
        lrpairs[["receptor.symbol"]] <- lrpairs[["receptor"]]
      } else {
        stop(
          "lr_reference CSV at '", lr_reference_path, "' must have columns ",
          "'ligand.symbol'/'receptor.symbol' (or 'ligand'/'receptor'). Found: ",
          paste(colnames(lrpairs), collapse = ", ")
        )
      }
    }
  } else if (nzchar(lr_reference_path)) {
    warning(
      "params.lr_reference was set to '", lr_reference_path, "' but that ",
      "file does not exist; falling back to CellChatDB.human."
    )
  }
  if (is.null(lrpairs)) {
    library(CellChat)
    lrdf <- CellChat::CellChatDB.human
    lrpairs <- lrdf[["interaction"]][, c("ligand.symbol", "receptor.symbol")]
  }
  ligands   <- sapply(lrpairs[["ligand.symbol"]],   function(i) strsplit(i, split = ", "))
  receptors <- sapply(lrpairs[["receptor.symbol"]], function(i) strsplit(i, split = ", "))
  names(ligands) <- names(receptors) <- rownames(lrpairs)
  gene_list <- union(unlist(ligands), unlist(receptors))
} else {
  # Use the top "good" genes by total expression, same ranking as before.
  # SpaceMarkers(genes = ...) only intersects a supplied list against what's
  # present in the data -- it doesn't rank by expression -- so that ranking
  # still has to happen against the raw counts here.
  expr_for_ranking <- load10XExpr(data_dir)
  expr_for_ranking <- expr_for_ranking[, colnames(sme)]
  gene_list <- apply(expr_for_ranking, 1, sum) |>
    sort(decreasing = TRUE) |>
    head(goodgeneThreshold) |>
    names()
}

# ---- Directed SpaceMarkers (one line) ----
sme <- SpaceMarkers(
  sme,
  directed          = TRUE,
  genes             = gene_list,
  lr_pairs          = lrpairs,
  avoid_confounders = TRUE
)

# ---- New output: the full SpaceMarkersExperiment ----
saveRDS(sme, file = sprintf("%s/sme.rds", output_dir))

# ---- Backward-compatible legacy outputs (same file names/shapes as before) ----
IMscores <- directed_scores(sme)
saveRDS(IMscores, file = sprintf("%s/IMscores.rds", output_dir))

if (useLigandReceptorGenes) {
  ligand_scores   <- S4Vectors::metadata(sme)[["ligand_scores"]]
  receptor_scores <- S4Vectors::metadata(sme)[["receptor_scores"]]
  lr_scores_df    <- lr_scores(sme)
  
  saveRDS(ligand_scores,   file = sprintf("%s/ligand_scores.rds", output_dir))
  saveRDS(receptor_scores, file = sprintf("%s/receptor_scores.rds", output_dir))
  saveRDS(lr_scores_df,    file = sprintf("%s/LRscores.rds", output_dir))
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