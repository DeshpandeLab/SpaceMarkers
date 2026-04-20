#' @title Main dispatcher for SpaceMarkers
#' @description Main dispatcher for SpaceMarkers. Accepts either a
#'   \code{\link{SpaceMarkersExperiment}} object or file paths to features
#'   and a Visium directory, and returns interaction scores.
#'
#' @param x An optional \code{\link{SpaceMarkersExperiment}} object. When
#'   provided, \code{features} and \code{data} are ignored.
#' @param features A path to a csv features file. Ignored when \code{x} is
#'   provided.
#' @param data A path to a 10X Visium directory. Ignored when \code{x} is
#'   provided.
#' @param directed Logical; run directed analysis (`TRUE`) or undirected
#'   analysis (`FALSE`).
#' @param genes Optional character vector of genes to retain. If `NULL`, genes
#'   are filtered by `min.gene.expr`.
#' @param min.gene.expr Minimum summed expression threshold for retaining genes
#'   when `genes` is `NULL`.
#' @param resolution Resolution passed to [load10XCoords()]. One of `"fullres"`,
#'   `"lowres"`, or `"hires"`.
#' @param version Optional Spaceranger version passed to [load10XCoords()].
#' @param h5filename Name of the 10X H5 expression file passed to
#'   [load10XExpr()].
#' @param spatialDir Name of the spatial subdirectory passed to
#'   [get_spatial_parameters()].
#' @param pattern Name of the JSON scale-factor file passed to
#'   [get_spatial_parameters()].
#' @param sigma Optional numeric sigma passed to [get_spatial_parameters()].
#' @param threshold Numeric threshold passed to [get_spatial_parameters()].
#' @param cpus Number of workers used in the undirected workflow passed to
#'   [get_pairwise_interacting_genes()].
#' @param returnSME Logical; if `TRUE` (default), return a
#'   \code{SpaceMarkersExperiment} with results stored in the object. If
#'   `FALSE`, return the legacy data frame of scores.
#' @param ... Additional arguments passed only to
#'   [get_pairwise_interacting_genes()] or [calculate_gene_scores_directed()].
#'
#' @return If \code{returnSME = TRUE}, a \code{SpaceMarkersExperiment} with
#'   analysis results. If \code{returnSME = FALSE}, a matrix or data frame of
#'   interaction scores (legacy format).
#' @export
SpaceMarkers <- function(x = NULL,
                         features = NULL,
                         data = NULL,
                         directed = FALSE,
                         genes = NULL,
                         min.gene.expr = 10,
                         resolution = c("fullres", "lowres", "hires"),
                         version = NULL,
                         h5filename = "filtered_feature_bc_matrix.h5",
                         spatialDir = "spatial",
                         pattern = "scalefactors_json.json",
                         sigma = NULL,
                         threshold = 4,
                         cpus = 1,
                         lr_pairs = NULL,
                         returnSME = TRUE,
                         ...) {
  resolution <- match.arg(resolution)

  # Determine input mode
  use_sme <- is(x, "SpaceMarkersExperiment")

  if (!use_sme && is.null(features) && is.null(data)) {
    stop("Provide either a SpaceMarkersExperiment 'x', or both 'features' and 'data'.")
  }

  if (isTRUE(directed)) {
    if (use_sme) {
      return(.directed_SpaceMarkers_sme(
        sme = x, genes = genes, min.gene.expr = min.gene.expr,
        spatialDir = spatialDir, pattern = pattern, sigma = sigma,
        threshold = threshold, resolution = resolution,
        lr_pairs = lr_pairs, returnSME = returnSME, ...))
    }
    result <- .directed_SpaceMarkers(
      features = features, data = data, genes = genes,
      min.gene.expr = min.gene.expr, resolution = resolution,
      version = version, h5filename = h5filename,
      spatialDir = spatialDir, pattern = pattern,
      sigma = sigma, threshold = threshold, ...)
    if (returnSME) {
      return(.wrap_directed_result(result, features, data, resolution,
                                   version, h5filename))
    }
    return(result)
  } else if (isFALSE(directed)) {
    if (use_sme) {
      return(.undirected_SpaceMarkers_sme(
        sme = x, cpus = cpus, genes = genes,
        min.gene.expr = min.gene.expr, spatialDir = spatialDir,
        pattern = pattern, sigma = sigma, threshold = threshold,
        resolution = resolution, returnSME = returnSME, ...))
    }
    result <- .undirected_SpaceMarkers(
      features = features, data = data, cpus = cpus, genes = genes,
      min.gene.expr = min.gene.expr, resolution = resolution,
      version = version, h5filename = h5filename,
      spatialDir = spatialDir, pattern = pattern,
      sigma = sigma, threshold = threshold, ...)
    if (returnSME) {
      return(.wrap_undirected_result(result, features, data, resolution,
                                    version, h5filename))
    }
    return(result)
  } else {
    stop("Check input: `directed` must be TRUE or FALSE.")
  }
}

# ---- Legacy file-path workflow (unchanged) ----

.process_input_data <- function(features,
                                data,
                                genes = NULL,
                                min.gene.expr = 10,
                                resolution = c("fullres", "lowres", "hires"),
                                version = NULL,
                                h5filename = "filtered_feature_bc_matrix.h5",
                                spatialDir = "spatial",
                                pattern = "scalefactors_json.json",
                                sigma = NULL,
                                threshold = 4) {
  resolution <- match.arg(resolution)

  message("Loading data...")

  visium_dir <- data
  coords <- load10XCoords(
    visiumDir = visium_dir,
    resolution = resolution,
    version = version
  )
  expr <- load10XExpr(
    visiumDir = visium_dir,
    h5filename = h5filename
  )
  spPatterns <- get_spatial_features(features)

  message("Preparing data...")
  rownames(coords) <- coords[["barcode"]]

  barcodes <- intersect(rownames(coords), rownames(spPatterns))
  coords <- coords[barcodes, , drop = FALSE]
  spPatterns <- spPatterns[barcodes, , drop = FALSE]

  spPatterns <- cbind(coords, spPatterns)

  message("Filtering data...")
  keepBarcodes <- intersect(colnames(expr), barcodes)
  expr <- expr[, keepBarcodes, drop = FALSE]
  spPatterns <- spPatterns[keepBarcodes, , drop = FALSE]

  if (!is.null(genes)) {
    keepGenes <- intersect(genes, rownames(expr))
  } else {
    keepGenes <- rownames(expr)[which(apply(expr, 1, sum) > min.gene.expr)]
  }

  expr <- expr[keepGenes, , drop = FALSE]

  message("Computing optimal parameters...")
  optParams <- get_spatial_parameters(
    spatialPatterns = spPatterns,
    visiumDir = visium_dir,
    spatialDir = spatialDir,
    pattern = pattern,
    sigma = sigma,
    threshold = threshold,
    resolution = resolution
  )

  list(
    data = expr,
    spPatterns = spPatterns,
    optParams = optParams
  )
}


#' @title Undirected SpaceMarkers workflow
#' @description Internal function to run the undirected SpaceMarkers analysis.
#'
#' @inheritParams SpaceMarkers
#'
#' @return Interaction scores for the undirected workflow.
.undirected_SpaceMarkers <- function(features,
                                     data,
                                     cpus = 1,
                                     genes = NULL,
                                     min.gene.expr = 10,
                                     resolution = c("fullres", "lowres", "hires"),
                                     version = NULL,
                                     h5filename = "filtered_feature_bc_matrix.h5",
                                     spatialDir = "spatial",
                                     pattern = "scalefactors_json.json",
                                     sigma = NULL,
                                     threshold = 4,
                                     ...) {
  inputs <- .process_input_data(
    features = features,
    data = data,
    genes = genes,
    min.gene.expr = min.gene.expr,
    resolution = resolution,
    version = version,
    h5filename = h5filename,
    spatialDir = spatialDir,
    pattern = pattern,
    sigma = sigma,
    threshold = threshold
  )

  expr <- inputs[["data"]]
  spPatterns <- inputs[["spPatterns"]]
  optParams <- inputs[["optParams"]]

  message("Finding hotspots...")
  hotspots <- find_all_hotspots(spPatterns)

  message("Finding interacting genes...")
  spaceMarkers <- get_pairwise_interacting_genes(
    data = expr,
    optParams = optParams,
    spPatterns = spPatterns,
    hotspots = hotspots,
    mode = "DE",
    analysis = "enrichment",
    minOverlap = 10,
    workers = cpus,
    ...
  )

  message("Calculating Interaction Scores...")
  IMScores <- get_im_scores(spaceMarkers)

  return(IMScores)
}


#' @title Directed SpaceMarkers workflow
#' @description Internal function to run the directed SpaceMarkers analysis.
#'
#' @inheritParams SpaceMarkers
#'
#' @return Interaction scores for the directed workflow.
.directed_SpaceMarkers <- function(features,
                                   data,
                                   genes = NULL,
                                   min.gene.expr = 10,
                                   resolution = c("fullres", "lowres", "hires"),
                                   version = NULL,
                                   h5filename = "filtered_feature_bc_matrix.h5",
                                   spatialDir = "spatial",
                                   pattern = "scalefactors_json.json",
                                   sigma = NULL,
                                   threshold = 4,
                                   ...) {
  inputs <- .process_input_data(
    features = features,
    data = data,
    genes = genes,
    min.gene.expr = min.gene.expr,
    resolution = resolution,
    version = version,
    h5filename = h5filename,
    spatialDir = spatialDir,
    pattern = pattern,
    sigma = sigma,
    threshold = threshold
  )

  expr <- inputs[["data"]]
  spPatterns <- inputs[["spPatterns"]]
  optParams <- inputs[["optParams"]]

  patnames <- setdiff(colnames(spPatterns), c("x", "y", "barcode"))

  message("Calculating hotspots and influence for each pattern...")
  patthresholds <- calculate_thresholds(
    spPatterns,
    minvals = 0.1,
    maxvals = 0.8
  )
  patHotspots <- find_hotspots_gmm(
    spPatterns,
    threshold = patthresholds
  )

  spInfluence <- calculate_influence(spPatterns, optParams)
  infthresholds <- calculate_thresholds(
    spInfluence,
    minvals = 0.01,
    maxvals = 0.5
  )
  infHotspots <- find_hotspots_gmm(
    spInfluence,
    threshold = infthresholds
  )

  if (length(patnames) < 2) {
    stop("At least two spatial patterns are required to compute directed pattern pairs.")
  }

  patternPairs <- t(combn(patnames, 2))

  message("Calculating directed interaction scores...")
  IMscores <- calculate_gene_scores_directed(
    data = expr,
    pat_hotspots = patHotspots,
    influence_hotspots = infHotspots,
    pattern_pairs = patternPairs,
    ...
  )

  return(IMscores)
}

# ---- SME-based workflows ----

#' Process input from a SpaceMarkersExperiment
#' @keywords internal
.process_sme_input <- function(sme,
                               genes = NULL,
                               min.gene.expr = 10,
                               spatialDir = "spatial",
                               pattern = "scalefactors_json.json",
                               sigma = NULL,
                               threshold = 4,
                               resolution = "fullres") {
  message("Preparing data from SpaceMarkersExperiment...")
  expr <- .sme_get_expr(sme)
  spPatterns <- .sme_to_spPatterns(sme)

  message("Filtering data...")
  all_expr <- colnames(expr)
  all_pats <- rownames(spPatterns)
  keepBarcodes <- intersect(all_expr, all_pats)
  dropped <- length(all_expr) + length(all_pats) - 2L * length(keepBarcodes)
  if (dropped > 0L) {
    warning(sprintf(
      "Dropped %d spots not present in both features and expression. Keeping %d common spots.",
      dropped, length(keepBarcodes)
    ))
  }
  expr <- expr[, keepBarcodes, drop = FALSE]
  spPatterns <- spPatterns[keepBarcodes, , drop = FALSE]

  if (!is.null(genes)) {
    keepGenes <- intersect(genes, rownames(expr))
  } else {
    keepGenes <- rownames(expr)[which(apply(expr, 1, sum) > min.gene.expr)]
  }
  expr <- expr[keepGenes, , drop = FALSE]

  message("Computing optimal parameters...")
  optParams <- spatial_params(sme)
  if (is.null(optParams)) {
    # Fall back to visiumDir stored at load time, or sigma if provided
    stored_visiumDir <- sme@spacemarkers$params$visiumDir
    stored_resolution <- if (!is.null(sme@spacemarkers$params$resolution))
        sme@spacemarkers$params$resolution else "fullres"
    optParams <- get_spatial_parameters(
      spatialPatterns = spPatterns,
      visiumDir = if (!is.null(stored_visiumDir)) stored_visiumDir else ".",
      spatialDir = spatialDir,
      pattern = pattern,
      sigma = sigma,
      threshold = threshold,
      resolution = stored_resolution
    )
  }

  list(data = expr, spPatterns = spPatterns, optParams = optParams)
}

#' Undirected workflow on a SpaceMarkersExperiment
#' @keywords internal
.undirected_SpaceMarkers_sme <- function(sme, cpus = 1, genes = NULL,
                                         min.gene.expr = 10,
                                         spatialDir = "spatial",
                                         pattern = "scalefactors_json.json",
                                         sigma = NULL, threshold = 4,
                                         resolution = "fullres",
                                         returnSME = TRUE, ...) {
  inputs <- .process_sme_input(sme, genes, min.gene.expr, spatialDir,
                               pattern, sigma, threshold, resolution)
  expr <- inputs$data
  spPatterns <- inputs$spPatterns
  optParams <- inputs$optParams

  message("Finding hotspots...")
  hs <- find_all_hotspots(spPatterns)

  message("Finding interacting genes...")
  spaceMarkersResult <- get_pairwise_interacting_genes(
    data = expr, optParams = optParams, spPatterns = spPatterns,
    hotspots = hs, mode = "DE", analysis = "enrichment",
    minOverlap = 10, workers = cpus, ...
  )

  message("Calculating Interaction Scores...")
  IMScores <- get_im_scores(spaceMarkersResult)

  if (!returnSME) return(IMScores)

  message("Calculating overlap scores...")
  overlap <- calculate_overlap_undirected(hs)

  # Store results in SME
  sm <- sme@spacemarkers
  if (is.null(sm$results)) sm$results <- list()
  sm$results$undirected_scores <- IMScores
  sm$results$overlap_scores <- overlap
  sm$analysis <- if (is.null(sm$analysis) || sm$analysis == "directed") {
    if (!is.null(sm$results$directed_scores)) "both" else "undirected"
  } else { sm$analysis }
  sm$params$spatial_params <- optParams
  sm$params$min_gene_expr <- min.gene.expr
  sm$params$mode <- "DE"
  sm$params$analysis_method <- "enrichment"
  sm$params$min_overlap <- 10
  sm$params$directed <- FALSE
  sm$params$genes <- genes
  sm$params$cpus <- cpus
  sme@spacemarkers <- sm

  # Store detailed intermediates in metadata
  S4Vectors::metadata(sme)$hotspots$undirected <- hs
  S4Vectors::metadata(sme)$interactions <- spaceMarkersResult

  sme
}

#' Directed workflow on a SpaceMarkersExperiment
#' @keywords internal
.directed_SpaceMarkers_sme <- function(sme, genes = NULL,
                                       min.gene.expr = 10,
                                       spatialDir = "spatial",
                                       pattern = "scalefactors_json.json",
                                       sigma = NULL, threshold = 4,
                                       resolution = "fullres",
                                       lr_pairs = NULL,
                                       returnSME = TRUE, ...) {
  inputs <- .process_sme_input(sme, genes, min.gene.expr, spatialDir,
                               pattern, sigma, threshold, resolution)
  expr <- inputs$data
  spPatterns <- inputs$spPatterns
  optParams <- inputs$optParams

  patnames <- setdiff(colnames(spPatterns), c("x", "y", "barcode"))

  message("Calculating hotspots and influence for each pattern...")
  patthresholds <- calculate_thresholds(spPatterns, minvals = 0.1, maxvals = 0.8)
  patHotspots <- find_hotspots_gmm(spPatterns, threshold = patthresholds)

  spInfluence <- calculate_influence(spPatterns, optParams)
  infthresholds <- calculate_thresholds(spInfluence, minvals = 0.01, maxvals = 0.5)
  infHotspots <- find_hotspots_gmm(spInfluence, threshold = infthresholds)

  if (length(patnames) < 2) {
    stop("At least two spatial patterns are required to compute directed pattern pairs.")
  }

  patternPairs <- t(combn(patnames, 2))

  message("Calculating directed interaction scores...")
  IMscores <- calculate_gene_scores_directed(
    data = expr, pat_hotspots = patHotspots,
    influence_hotspots = infHotspots,
    pattern_pairs = patternPairs, ...
  )

  if (!returnSME) return(IMscores)

  message("Calculating overlap scores...")
  overlap <- calculate_overlap_directed(patHotspots, infHotspots)

  # Store results in SME
  sm <- sme@spacemarkers
  if (is.null(sm$results)) sm$results <- list()
  sm$results$directed_scores <- IMscores
  sm$results$overlap_scores <- overlap
  sm$analysis <- if (is.null(sm$analysis) || sm$analysis == "undirected") {
    if (!is.null(sm$results$undirected_scores)) "both" else "directed"
  } else { sm$analysis }
  sm$params$spatial_params <- optParams
  sm$params$min_gene_expr <- min.gene.expr
  sm$params$directed <- TRUE
  sm$params$genes <- genes

  # LR scores (if lr_pairs provided)
  if (!is.null(lr_pairs)) {
    message("Calculating LR scores...")
    ligand_scores <- calculate_gene_set_score(IMscores, lr_pairs$ligand.symbol)
    receptor_scores <- calculate_gene_set_specificity(
      IMscores, lr_pairs$receptor.symbol)
    lr <- calculate_lr_scores(ligand_scores, receptor_scores, lr_pairs)
    sm$results$lr_scores <- lr
    sm$params$lr_pairs <- lr_pairs
    S4Vectors::metadata(sme)$ligand_scores <- ligand_scores
    S4Vectors::metadata(sme)$receptor_scores <- receptor_scores
  }

  sme@spacemarkers <- sm

  # Store detailed intermediates in metadata
  S4Vectors::metadata(sme)$hotspots$pattern <- patHotspots
  S4Vectors::metadata(sme)$hotspots$influence <- infHotspots
  S4Vectors::metadata(sme)$influence <- spInfluence

  sme
}

# ---- Legacy result wrapping helpers ----

#' Wrap legacy undirected result as SME (when file-path input used with returnSME=TRUE)
#' @keywords internal
.wrap_undirected_result <- function(IMScores, features, data, resolution,
                                   version, h5filename) {
  sme <- tryCatch(
    load10X(visiumDir = data, features = features,
                 h5filename = h5filename, resolution = resolution,
                 version = version),
    error = function(e) NULL
  )
  if (is.null(sme)) return(IMScores)

  sm <- sme@spacemarkers
  sm$results <- list(undirected_scores = IMScores)
  sm$analysis <- "undirected"
  sme@spacemarkers <- sm
  sme
}

#' Wrap legacy directed result as SME (when file-path input used with returnSME=TRUE)
#' @keywords internal
.wrap_directed_result <- function(IMscores, features, data, resolution,
                                 version, h5filename) {
  sme <- tryCatch(
    load10X(visiumDir = data, features = features,
                 h5filename = h5filename, resolution = resolution,
                 version = version),
    error = function(e) NULL
  )
  if (is.null(sme)) return(IMscores)

  sm <- sme@spacemarkers
  sm$results <- list(directed_scores = IMscores)
  sm$analysis <- "directed"
  sme@spacemarkers <- sm
  sme
}
