#!/usr/bin/env Rscript
# SpaceMarkers >= 2.2 on AnnData via SpaceMarkersExperiment.
# NB: avoid dollar signs and backslashes here, Nextflow templates interpolate them.

suppressPackageStartupMessages({
  devtools::load_all("/spacemarkers")
  library(SummarizedExperiment)
})

adata_path <- "${adata}"
output_dir <- "${prefix}"
patterns_key <- "${params.sm_patterns_uns}"
directed_param <- "${params.sm_directed}"
spot_diameter_param <- "${params.sm_spot_diameter}"
use_lr_genes <- as.logical("${params.sm_use_lr_genes}")
good_gene_threshold <- ${params.sm_good_gene_thd}
n_spots_directed <- 10000
set.seed(${params.seed})

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

message('setting up parallel param..')
workers <- min(as.numeric("${task.cpus}"), parallel::detectCores())
BiocParallel::register(BiocParallel::MulticoreParam(workers = workers)) #register backend

message('creating sme object..')
sme <- load_anndata(adata_path, patterns_meta_table = patterns_key)

message('extracting spatial patterns..')
pats <- spatial_patterns(sme)
if (is.null(pats) || ncol(pats) < 2) {
  stop(sprintf("Need at least 2 latent features in uns[['%s']] of %s", patterns_key, adata_path))
}

# Mode: explicit param wins, otherwise by total spot count
n_spots <- ncol(sme)
if (nzchar(directed_param)) {
  directed <- as.logical(directed_param)
} else {
  directed <- n_spots >= n_spots_directed
}
message(sprintf("spots: %d, mode: %s", n_spots, ifelse(directed, "directed", "undirected")))

# Spatial kernel width: param override, then uns scalefactor (set by load_anndata), then estimate
estimate_spot_diameter <- function(coords, n_sample = 2000) {
  coords <- as.matrix(coords[, c("x", "y")])
  idx <- sample.int(nrow(coords), min(n_sample, nrow(coords)))
  nn <- vapply(idx, function(i) {
    d <- sqrt((coords[, 1] - coords[i, 1])^2 + (coords[, 2] - coords[i, 2])^2)
    min(d[d > 0])
  }, numeric(1))
  # 55 um spot diameter over 100 um centre-to-centre pitch on Visium
  0.55 * median(nn)
}

sigma <- NULL
sigma_method <- NULL
if (nzchar(spot_diameter_param)) {
  sigma <- as.numeric(spot_diameter_param)
  sigma_method <- "param"
} else if (is.null(spatial_params(sme))) {
  sigma <- estimate_spot_diameter(spatialCoords(sme))
  sigma_method <- "estimated"
}
if (!is.null(sigma)) {
  if (!is.finite(sigma) || sigma <= 0) stop("Could not determine a valid spot diameter")
  message(sprintf("spot diameter (%s): %.2f", sigma_method, sigma))
  sp_patterns <- cbind(as.data.frame(pats),
                       barcode = colnames(sme),
                       x = spatialCoords(sme)[, "x"],
                       y = spatialCoords(sme)[, "y"])
  spatial_params(sme) <- get_spatial_parameters(spatialPatterns = sp_patterns, sigma = sigma)
} else {
  message("spot diameter (scalefactor): from AnnData uns")
}

expr <- assay(sme, if ("logcounts" %in% assayNames(sme)) "logcounts" else 1L)

if (directed) {
  if (use_lr_genes) {
    library(CellChat)
    lrdf <- CellChat::CellChatDB.human
    lrpairs <- lrdf[["interaction"]][, c("ligand.symbol", "receptor.symbol")]
    ligands <- sapply(lrpairs[["ligand.symbol"]], function(i) strsplit(i, split = ", "))
    receptors <- sapply(lrpairs[["receptor.symbol"]], function(i) strsplit(i, split = ", "))
    names(ligands) <- names(receptors) <- rownames(lrpairs)
    lrgenes <- union(unlist(ligands), unlist(receptors))
    sme <- sme[rownames(sme) %in% lrgenes, ]
  } else {
    goodgenes <- names(head(sort(rowSums(expr), decreasing = TRUE), good_gene_threshold))
    sme <- sme[rownames(sme) %in% goodgenes, ]
  }
  message('calculating directed influence scores...')
  sme <- calculate_influence(sme)
  message('finding pattern hotspots using GMM...')
  sme <- find_hotspots_gmm(sme, type = "pattern")
  message('finding influence hotspots using GMM...')
  sme <- find_hotspots_gmm(sme, type = "influence")
  message('calculating gene scores...')
  sme <- calculate_gene_scores_directed(sme, avoid_confounders = TRUE)
  saveRDS(directed_scores(sme), file = file.path(output_dir, "IMscores.rds"))

  if (use_lr_genes) {
    message('calculating ligand-receptor scores...')
    sme <- calculate_gene_set_score(sme, gene_sets = ligands, weighted = TRUE, method = "arithmetic_mean")
    sme <- calculate_gene_set_specificity(sme, gene_sets = receptors, weighted = TRUE, method = "arithmetic_mean")
    sme <- calculate_lr_scores(sme, lr_pairs = lrpairs, method = "geometric_mean", weighted = TRUE)
    saveRDS(lr_scores(sme), file = file.path(output_dir, "LRscores.rds"))
  }
} else {
  # drop barely expressed genes
  sme <- sme[rowSums(expr) > 10, ]

  sme <- find_all_hotspots(sme)
  sme <- get_pairwise_interacting_genes(sme, mode = "DE", analysis = "enrichment",
                                        minOverlap = 10, workers = ${task.cpus})
  sme <- get_im_scores(sme)
  imscores <- undirected_scores(sme)
  if ("Gene" %in% colnames(imscores)) {
    rownames(imscores) <- imscores[["Gene"]]
    imscores[["Gene"]] <- NULL
  }
  saveRDS(imscores, file = file.path(output_dir, "IMscores.rds"))
}

sinfo <- sessionInfo()
versions <- lapply(sinfo[["otherPkgs"]], function(x) sprintf("  %s: %s", x[["Package"]], x[["Version"]]))
versions[["R"]] <- sprintf("  R: %s", as.character(packageVersion("base")))
writeLines(c('"${task.process}":', unlist(versions)), "versions.yml")
