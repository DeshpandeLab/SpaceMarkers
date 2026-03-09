#' @title Main dispatcher for SpaceMarkers
#' @description Main dispatcher for SpaceMarkers. It takes a features object/path
#'   and a Visium directory and returns interaction scores.
#'
#' @param features A path to a csv features file.
#' @param data A path to a 10X Visium directory.
#' @param directed Logical; run directed analysis (`TRUE`) or undirected analysis (`FALSE`).
#' @param genes Optional character vector of genes to retain. If `NULL`, genes are
#'   filtered by `min.gene.expr`.
#' @param min.gene.expr Minimum summed expression threshold for retaining genes
#'   when `genes` is `NULL`.
#' @param resolution Resolution passed to [load10XCoords()]. One of `"fullres"`,
#'   `"lowres"`, or `"hires"`.
#' @param version Optional Spaceranger version passed to [load10XCoords()].
#' @param h5filename Name of the 10X H5 expression file passed to [load10XExpr()].
#' @param spatialDir Name of the spatial subdirectory passed to
#'   [get_spatial_parameters()].
#' @param pattern Name of the JSON scale-factor file passed to
#'   [get_spatial_parameters()].
#' @param sigma Optional numeric sigma passed to [get_spatial_parameters()].
#' @param threshold Numeric threshold passed to [get_spatial_parameters()].
#' @param cpus Number of workers used in the undirected workflow passed to
#'   [get_pairwise_interacting_genes()].
#' @param ... Additional arguments passed only to
#'   [get_pairwise_interacting_genes()] or [calculate_gene_scores_directed()].
#'
#' @return A matrix or data frame of interaction scores, depending on the selected
#'   workflow.
#' @export
SpaceMarkers <- function(features,
                         data,
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
                         ...) {
  resolution <- match.arg(resolution)
  
  if (isTRUE(directed)) {
    return(
      .directed_SpaceMarkers(
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
        threshold = threshold,
        ...
      )
    )
  } else if (isFALSE(directed)) {
    return(
      .undirected_SpaceMarkers(
        features = features,
        data = data,
        cpus = cpus,
        genes = genes,
        min.gene.expr = min.gene.expr,
        resolution = resolution,
        version = version,
        h5filename = h5filename,
        spatialDir = spatialDir,
        pattern = pattern,
        sigma = sigma,
        threshold = threshold,
        ...
      )
    )
  } else {
    stop("Check input: `directed` must be TRUE or FALSE.")
  }
}

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
  spPatterns <- get_spatial_features(features )
  
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