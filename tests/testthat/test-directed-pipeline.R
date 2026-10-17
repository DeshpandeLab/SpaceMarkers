test_that("there is no numeric genes in the output", {
  spatial_path <- "assets/staple.h5ad"
  patterns_key <- "cell_type_composition"
  directed_param <- TRUE
  spot_diameter_param <- ""
  use_lr_genes <- TRUE
  good_gene_threshold <- 1000
  n_spots_directed <- 10000
  set.seed(123)
  features_path <- ""
  lr_reference_path <- "assets/ligrec_db.csv"
  output_dir <- "."


  message("creating sme object..")
  if (grepl(".h5ad", spatial_path)) {
    message("loading from AnnData..")
    sme <- load_anndata(spatial_path, patterns_meta_table = patterns_key)
  } else {
    message("loading from visium dir..")
    sme <- load10X(spatial_path, features_path)
  }

  message("extracting spatial patterns..")
  pats <- spatial_patterns(sme)
  if (is.null(pats) || ncol(pats) < 2) {
    stop(sprintf("Need at least 2 latent features in uns[['%s']] of %s", patterns_key, spatial_path))
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
  if (nzchar(spot_diameter_param) && !is.na(as.numeric(spot_diameter_param))) {
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

  expr <- assay(sme, if ("logcounts" %in% SummarizedExperiment::assayNames(sme)) "logcounts" else 1L)

  if (directed) {
    if (use_lr_genes) {
      lrpairs <- NULL
      if (nzchar(lr_reference_path)) {
        message("using user-supplied lr_reference: ", lr_reference_path)
        lrpairs <- read.csv(lr_reference_path)
        # CSV with pair IDs as row names (first column) and columns
        # ligand.symbol/receptor.symbol (or ligand/receptor); entries may hold
        # several symbols separated by ', ' (same convention as CellChatDB.human)
        if (!all(c("ligand.symbol", "receptor.symbol") %in% colnames(lrpairs))) {
          if (all(c("ligand", "receptor") %in% colnames(lrpairs))) {
            message("Using 'ligand' and 'receptor' columns as ligand.symbol and receptor.symbol")
            lrpairs[["ligand.symbol"]] <- lrpairs[["ligand"]]
            lrpairs[["receptor.symbol"]] <- lrpairs[["receptor"]]
          } else if (all(c("source","target") %in% colnames(lrpairs))) {
            message("Omnipath format detected: using 'source' and 'target' columns")
            lrpairs[["ligand.symbol"]] <- lrpairs[["source"]]
            lrpairs[["receptor.symbol"]] <- lrpairs[["target"]]
          }
            else {
            stop("lr_reference CSV must have columns 'ligand.symbol'/'receptor.symbol' ",
                "(or 'ligand'/'receptor'). Found: ", paste(colnames(lrpairs), collapse = ", "))
          }
        }
        lrpairs <- lrpairs[, c("ligand.symbol", "receptor.symbol")]
      } else {
        library(CellChat)
        lrdf <- CellChat::CellChatDB.human
        lrpairs <- lrdf[["interaction"]][, c("ligand.symbol", "receptor.symbol")]
      }
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
                                          minOverlap = 10, workers = workers)
    sme <- get_im_scores(sme)
    imscores <- undirected_scores(sme)
    if ("Gene" %in% colnames(imscores)) {
      rownames(imscores) <- imscores[["Gene"]]
      imscores[["Gene"]] <- NULL
    }
    saveRDS(imscores, file = file.path(output_dir, "IMscores.rds"))
  }

  # full SpaceMarkersExperiment, saved in both directed and undirected modes
  saveRDS(sme, file = file.path(output_dir, "sme.rds"))

  expect_false(directed_scores(sme)[['gene']] |> as.numeric() |> sum(na.rm=TRUE) > 0)
})