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
#' @param lr_pairs A data.frame with \code{ligand.symbol}, \code{receptor.symbol},
#'   and \code{pair} columns, or \code{NULL}. Used by the directed SME workflow.
#' @param returnSME Logical; if `TRUE` (default), return a
#'   \code{SpaceMarkersExperiment} with results stored in the object. If
#'   `FALSE`, return the legacy data frame of scores.
#' @param ... Additional arguments passed only to
#'   [get_pairwise_interacting_genes()] or [calculate_gene_scores_directed()].
#'
#' @return If \code{returnSME = TRUE}, a \code{SpaceMarkersExperiment} with
#'   analysis results. If \code{returnSME = FALSE}, a matrix or data frame of
#'   interaction scores (legacy format).
#' @examples
#' \donttest{
#' # End-to-end on an SME (small synthetic example):
#' set.seed(1)
#' nb <- 30
#' bc <- paste0("s", seq_len(nb))
#' sme <- SpaceMarkersExperiment(
#'     assays = list(logcounts = matrix(rpois(20 * nb, 3), 20, nb,
#'         dimnames = list(paste0("G", 1:20), bc))),
#'     spatialCoords = matrix(runif(2 * nb), nb, 2,
#'         dimnames = list(bc, c("y", "x"))))
#' spatial_patterns(sme) <- data.frame(Pattern_1 = runif(nb),
#'                                     Pattern_2 = runif(nb),
#'                                     row.names = bc)
#' spatial_params(sme) <- matrix(c(0.1, 4, 0.1, 4), nrow = 2,
#'     dimnames = list(c("sigmaOpt","threshOpt"), c("Pattern_1","Pattern_2")))
#' result <- SpaceMarkers(sme, directed = FALSE)
#' }
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
    # Matrix::rowSums stays sparse-aware on dgCMatrix; apply() would densify
    # the whole expression matrix (~1 GB transient on a typical Visium load).
    keepGenes <- rownames(expr)[Matrix::rowSums(expr) > min.gene.expr]
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

#' Filter an SME by gene expression / barcode match and (re)compute spatial params
#' @return The filtered \code{SpaceMarkersExperiment} (genes/spots restricted,
#'   \code{spatial_params(sme)} populated when not already set).
#' @keywords internal
.apply_sme_filters <- function(sme, genes = NULL, min.gene.expr = 10,
                               spatialDir = "spatial",
                               pattern = "scalefactors_json.json",
                               sigma = NULL, threshold = 4,
                               resolution = "fullres") {
    message("Preparing data from SpaceMarkersExperiment...")

    # Barcode consistency: intersect expression colnames with pattern rownames
    expr <- .sme_expr(sme)
    pat_barcodes <- rownames(spatial_patterns(sme))
    if (!is.null(pat_barcodes)) {
        keep <- intersect(colnames(expr), pat_barcodes)
        if (length(keep) < ncol(sme)) {
            sme <- sme[, keep]
        }
    }

    # Gene filtering — Matrix::rowSums avoids densifying a dgCMatrix
    # (apply(e, 1, sum) would force a ~1 GB transient on a typical Visium load).
    e <- .sme_expr(sme)
    keep_genes <- if (!is.null(genes))
        intersect(genes, rownames(e))
    else
        rownames(e)[Matrix::rowSums(e) > min.gene.expr]
    sme <- sme[keep_genes, ]

    # Spatial params fallback
    if (is.null(spatial_params(sme))) {
        stored_visiumDir <- sme@spacemarkers$params$visiumDir
        if (is.null(stored_visiumDir)) {
            # No visiumDir to read scalefactors from. Fall back to the
            # same defaults find_pattern_hotspots() uses when params is
            # NULL: sigmaOpt = sigma (or 10), threshOpt = threshold.
            patternList <- setdiff(colnames(.sme_spPatterns(sme)),
                                   c("x", "y", "barcode"))
            sigma_default <- if (is.null(sigma)) 10 else as.numeric(sigma[1])
            op <- matrix(c(sigma_default, threshold),
                         nrow = 2, ncol = length(patternList),
                         dimnames = list(c("sigmaOpt", "threshOpt"),
                                         patternList))
            message("spatial_params(sme) is NULL and no visiumDir is ",
                    "stored. Using sigma = ", sigma_default,
                    ", threshold = ", threshold,
                    " for all patterns. Set spatial_params(sme) ",
                    "explicitly before running SpaceMarkers() if you ",
                    "want tuned values per pattern.")
            spatial_params(sme) <- op
        } else {
            message("Computing optimal parameters...")
            stored_res <- sme@spacemarkers$params$resolution %||% resolution
            op <- get_spatial_parameters(
                spatialPatterns = .sme_spPatterns(sme),
                visiumDir = stored_visiumDir,
                spatialDir = spatialDir, pattern = pattern,
                sigma = sigma, threshold = threshold, resolution = stored_res)
            spatial_params(sme) <- op
        }
    }

    sme
}

#' Undirected workflow on a SpaceMarkersExperiment
#' @return If \code{returnSME = TRUE}, the input SME with hotspot,
#'   overlap, and interaction results attached; otherwise the legacy
#'   data frame produced by the underlying pipeline.
#' @keywords internal
.undirected_SpaceMarkers_sme <- function(sme, cpus = 1, genes = NULL,
                                         min.gene.expr = 10,
                                         spatialDir = "spatial",
                                         pattern = "scalefactors_json.json",
                                         sigma = NULL, threshold = 4,
                                         resolution = "fullres",
                                         returnSME = TRUE, ...) {
    sme <- .apply_sme_filters(sme, genes = genes,
                              min.gene.expr = min.gene.expr,
                              spatialDir = spatialDir, pattern = pattern,
                              sigma = sigma, threshold = threshold,
                              resolution = resolution)
    message("Finding hotspots...")
    sme <- find_all_hotspots(sme)
    message("Calculating overlap scores...")
    sme <- calculate_overlap_undirected(sme)
    message("Finding interacting genes...")
    sme <- get_pairwise_interacting_genes(
        sme, mode = "DE", analysis = "enrichment",
        minOverlap = 10, workers = cpus, ...)
    message("Calculating Interaction Scores...")
    sme <- get_im_scores(sme)
    if (!returnSME) return(undirected_scores(sme))
    sm <- sme@spacemarkers
    if (is.null(sm$params)) sm$params <- list()
    sm$params$min_gene_expr <- min.gene.expr
    sm$params$directed <- FALSE
    sm$params$genes <- genes
    sm$params$cpus <- cpus
    sme@spacemarkers <- sm
    sme
}

#' Directed workflow on a SpaceMarkersExperiment
#' @return If \code{returnSME = TRUE}, the input SME with directed
#'   hotspots, gene scores, and LR scores attached; otherwise the
#'   legacy data frame produced by the underlying pipeline.
#' @keywords internal
.directed_SpaceMarkers_sme <- function(sme, genes = NULL,
                                       min.gene.expr = 10,
                                       spatialDir = "spatial",
                                       pattern = "scalefactors_json.json",
                                       sigma = NULL, threshold = 4,
                                       resolution = "fullres",
                                       lr_pairs = NULL,
                                       returnSME = TRUE, ...) {
    sme <- .apply_sme_filters(sme, genes = genes,
                              min.gene.expr = min.gene.expr,
                              spatialDir = spatialDir, pattern = pattern,
                              sigma = sigma, threshold = threshold,
                              resolution = resolution)
    message("Calculating influence and pattern hotspots...")
    sme <- calculate_influence(sme)
    sme <- find_hotspots_gmm(sme, type = "pattern")
    sme <- find_hotspots_gmm(sme, type = "influence")
    message("Calculating overlap scores...")
    sme <- calculate_overlap_directed(sme)
    message("Calculating directed interaction scores...")
    sme <- calculate_gene_scores_directed(sme, ...)
    if (!returnSME) return(directed_scores(sme))
    if (!is.null(lr_pairs)) {
        message("Calculating LR scores...")
        sm <- sme@spacemarkers
        if (is.null(sm$params)) sm$params <- list()
        sm$params$lr_pairs <- lr_pairs
        sme@spacemarkers <- sm
        sme <- calculate_gene_set_score(sme)
        sme <- calculate_gene_set_specificity(sme)
        sme <- calculate_lr_scores(sme)
    }
    sm <- sme@spacemarkers
    if (is.null(sm$params)) sm$params <- list()
    sm$params$min_gene_expr <- min.gene.expr
    sm$params$directed <- TRUE
    sm$params$genes <- genes
    sme@spacemarkers <- sm
    sme
}

# ---- Legacy result wrapping helpers ----

#' Wrap legacy undirected result as SME (when file-path input used with returnSME=TRUE)
#' @return A \code{SpaceMarkersExperiment} carrying the legacy
#'   undirected result (\code{IMScores}, hotspots, params) in its
#'   \code{spacemarkers} slot.
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
#' @return A \code{SpaceMarkersExperiment} carrying the legacy
#'   directed result (\code{IMScores}, hotspots, LR scores) in its
#'   \code{spacemarkers} slot.
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
