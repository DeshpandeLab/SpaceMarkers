#' @title SpaceMarkersExperiment class
#' @description An S4 class extending SpatialExperiment to hold SpaceMarkers
#' analysis results alongside spatial transcriptomics data.
#'
#' @slot spacemarkers A \code{SimpleList} containing primary analysis results:
#'   \describe{
#'     \item{params}{A list of all hyperparameters used in analysis, including:
#'       \code{pattern_names} (character vector),
#'       \code{spatial_params} (2 x N matrix of sigmaOpt/threshOpt per pattern),
#'       \code{min_gene_expr} (numeric), \code{mode} (character),
#'       \code{analysis_method} (character), \code{min_overlap} (numeric),
#'       \code{directed} (logical), \code{genes} (character or NULL), etc.
#'       Access via \code{params(sme)} for the full list, or
#'       \code{spatial_params(sme)} for just the kernel parameter matrix.}
#'     \item{results}{A list with \code{undirected_scores} (data.frame),
#'       \code{directed_scores} (data.frame), \code{lr_scores} (matrix of
#'       ligand-receptor pair scores), and \code{overlap_scores} (data.frame
#'       of pattern overlap scores).}
#'     \item{analysis}{A character string: "undirected", "directed", or "both".}
#'   }
#'
#'   Detailed intermediate results are stored in \code{metadata()}:
#'   \describe{
#'     \item{hotspots}{A list with elements \code{undirected}, \code{pattern},
#'       and \code{influence} (each a data.frame).}
#'     \item{interactions}{A named list of per-pair results from
#'       \code{get_pairwise_interacting_genes()}.}
#'     \item{influence}{A data.frame of per-spot influence values.}
#'     \item{ligand_scores}{Output of \code{calculate_gene_set_score()}.}
#'     \item{receptor_scores}{Output of \code{calculate_gene_set_specificity()}.}
#'   }
#'
#' @importClassesFrom SpatialExperiment SpatialExperiment
#' @importFrom S4Vectors SimpleList
#' @importFrom methods setClass setValidity setMethod new is validObject
#'   callNextMethod show
#' @exportClass SpaceMarkersExperiment
setClass(
    "SpaceMarkersExperiment",
    contains = "SpatialExperiment",
    slots = c(spacemarkers = "SimpleList")
)

setValidity("SpaceMarkersExperiment", function(object) {
    msg <- character()
    sm <- object@spacemarkers
    if (!is(sm, "SimpleList")) {
        msg <- c(msg, "'spacemarkers' slot must be a SimpleList")
    }
    valid_names <- c("params", "results", "analysis")
    unknown <- setdiff(names(sm), valid_names)
    if (length(unknown) > 0) {
        msg <- c(msg, paste("Unknown spacemarkers elements:",
                            paste(unknown, collapse = ", ")))
    }
    if (!is.null(sm$analysis)) {
        if (!sm$analysis %in% c("undirected", "directed", "both")) {
            msg <- c(msg, "'analysis' must be 'undirected', 'directed', or 'both'")
        }
    }
    if (length(msg)) msg else TRUE
})

#' @title Create a SpaceMarkersExperiment object
#' @description Constructor for the SpaceMarkersExperiment class. Accepts
#'   either a \code{SpatialExperiment} object to coerce or raw components
#'   (assays, coordinates, etc.) to build from scratch. When an SPE is
#'   provided, features can optionally be added in the same call via
#'   \code{\link{add_features}}.
#'
#' @param assays A \code{SpatialExperiment} object to coerce, or a list of
#'   assay matrices (e.g., logcounts) when building from scratch.
#' @param colData A DataFrame of per-spot metadata. Ignored when \code{assays}
#'   is a SpatialExperiment.
#' @param spatialCoords A numeric matrix of spatial coordinates (spots x 2).
#'   Ignored when \code{assays} is a SpatialExperiment.
#' @param rowData Optional DataFrame of per-gene metadata. Ignored when
#'   \code{assays} is a SpatialExperiment.
#' @param spaceMarkers A SimpleList or list for SpaceMarkers results.
#'   Defaults to an empty SimpleList.
#' @param features Optional: a file path, matrix/data.frame, or R object
#'   (e.g., CoGAPS result) of spatial features to add via
#'   \code{\link{add_features}}.
#' @param ... Additional arguments passed to \code{SpatialExperiment()} when
#'   building from scratch, or to \code{\link{get_spatial_features}} when
#'   \code{features} is provided.
#' @return A \code{SpaceMarkersExperiment} object.
#' @examples
#' set.seed(1)
#' nb <- 20
#' sme <- SpaceMarkersExperiment(
#'     assays = list(logcounts = matrix(rpois(10 * nb, 2), 10, nb,
#'         dimnames = list(paste0("G", 1:10), paste0("s", seq_len(nb))))),
#'     spatialCoords = matrix(runif(2 * nb), nb, 2,
#'         dimnames = list(paste0("s", seq_len(nb)), c("y", "x"))))
#' sme
#' @export
#' @importFrom SpatialExperiment SpatialExperiment spatialCoords
#' @importFrom SummarizedExperiment assays colData rowData
#' @importFrom S4Vectors SimpleList DataFrame
SpaceMarkersExperiment <- function(
    assays,
    colData = NULL,
    spatialCoords = NULL,
    rowData = NULL,
    spaceMarkers = S4Vectors::SimpleList(),
    features = NULL,
    ...
) {
    if (is.list(spaceMarkers) && !is(spaceMarkers, "SimpleList")) {
        spaceMarkers <- as(spaceMarkers, "SimpleList")
    }

    # Coerce from SpatialExperiment
    if (is(assays, "SpatialExperiment")) {
        sme <- new("SpaceMarkersExperiment", assays,
                    spacemarkers = spaceMarkers)
        if (!is.null(features)) {
            sme <- add_features(sme, features, ...)
        }
        return(sme)
    }

    # Build from raw components
    spe_args <- list(assays = assays, spatialCoords = spatialCoords, ...)
    if (!is.null(colData)) spe_args$colData <- colData
    if (!is.null(rowData)) spe_args$rowData <- rowData
    spe <- do.call(SpatialExperiment::SpatialExperiment, spe_args)
    sme <- new("SpaceMarkersExperiment", spe, spacemarkers = spaceMarkers)
    if (!is.null(features)) {
        sme <- add_features(sme, features, ...)
    }
    sme
}

#' Coercion from SpatialExperiment to SpaceMarkersExperiment
#'
#' @return A \code{SpaceMarkersExperiment} wrapping \code{from}, with an
#'   empty \code{spacemarkers} slot ready to be populated by the
#'   pipeline.
#' @name as-SpatialExperiment-SpaceMarkersExperiment
#' @importFrom methods setAs
setAs("SpatialExperiment", "SpaceMarkersExperiment", function(from) {
    new("SpaceMarkersExperiment", from,
        spacemarkers = S4Vectors::SimpleList())
})

#' Coercion from SingleCellExperiment to SpaceMarkersExperiment
#'
#' AnnData-sourced objects (e.g. via `zellkonverter::readH5AD()` or
#' `anndataR::as_SingleCellExperiment()`) typically arrive as
#' `SingleCellExperiment` with spatial coordinates stored as a
#' reducedDim, commonly named "spatial". This coerce method promotes
#' that reducedDim into `spatialCoords()` and wraps the result as an
#' SME, so a full AnnData -> SME path reads:
#'
#'     sme <- as(zellkonverter::readH5AD(path), "SpaceMarkersExperiment")
#'
#' If no suitable spatial reducedDim is found, the SME is built with
#' empty `spatialCoords`; the caller can populate them afterwards.
#'
#' @return A \code{SpaceMarkersExperiment} carrying the same assays /
#'   colData / rowData as \code{from}, with the spatial reducedDim
#'   promoted to \code{spatialCoords()} when present.
#' @name as-SingleCellExperiment-SpaceMarkersExperiment
setAs("SingleCellExperiment", "SpaceMarkersExperiment", function(from) {
    # Look for a spatial reducedDim under common AnnData conventions.
    rd_names <- SingleCellExperiment::reducedDimNames(from)
    spatial_rd <- intersect(c("spatial", "X_spatial", "SPATIAL"), rd_names)
    spatial_coords <- NULL
    if (length(spatial_rd) > 0L) {
        spatial_coords <- as.matrix(
            SingleCellExperiment::reducedDim(from, spatial_rd[1L])
        )
        if (ncol(spatial_coords) >= 2L) {
            colnames(spatial_coords)[seq_len(2L)] <- c("x", "y")
        }
        rownames(spatial_coords) <- colnames(from)
    }
    spe <- methods::as(from, "SpatialExperiment")
    if (!is.null(spatial_coords) &&
        ncol(SpatialExperiment::spatialCoords(spe)) == 0L) {
        SpatialExperiment::spatialCoords(spe) <- spatial_coords
    }
    methods::as(spe, "SpaceMarkersExperiment")
})

#' @rdname SpaceMarkersExperiment-class
#' @aliases show,SpaceMarkersExperiment-method
#' @param object A \code{SpaceMarkersExperiment}.
setMethod("show", "SpaceMarkersExperiment", function(object) {
    callNextMethod()
    sm <- object@spacemarkers
    analysis <- sm$analysis
    cat("spacemarkers analysis:",
        if (is.null(analysis)) "none" else analysis, "\n")
    if (!is.null(sm$params$pattern_names)) {
        cat("  patterns:", paste(sm$params$pattern_names, collapse = ", "), "\n")
    }
    if (!is.null(sm$results$undirected_scores)) {
        n_pairs <- ncol(sm$results$undirected_scores) - 1L
        cat("  undirected scores:", n_pairs, "pattern pairs\n")
    }
    if (!is.null(sm$results$directed_scores)) {
        n_int <- length(unique(sm$results$directed_scores$cell_interaction))
        cat("  directed scores:", n_int, "cell interactions\n")
    }
    if (!is.null(sm$results$overlap_scores)) {
        cat("  overlap scores:", nrow(sm$results$overlap_scores), "pairs\n")
    }
    if (!is.null(sm$results$lr_scores)) {
        cat("  LR scores:", nrow(sm$results$lr_scores), "pairs\n")
    }
    md <- S4Vectors::metadata(object)
    if (!is.null(md$interactions)) {
        cat("  detailed interactions:", length(md$interactions), "pairs\n")
    }
})
