#' Internal imports for SpaceMarkersExperiment methods
#'
#' @importFrom S4Vectors metadata metadata<-
#' @importFrom SummarizedExperiment assayNames colData<-
#' @importFrom methods setGeneric validObject
#' @name SpaceMarkersExperiment-imports
#' @keywords internal
NULL

#' @include AllClasses.R AllGenerics.R

#' @rdname spatial_patterns
#' @aliases spatial_patterns,SpaceMarkersExperiment-method
#' @param x A \code{SpaceMarkersExperiment} object.
#' @return A \code{DataFrame} of spatial pattern values per spot, or NULL if
#'   pattern names have not been set.
#' @export
setMethod("spatial_patterns", "SpaceMarkersExperiment", function(x) {
    pn <- x@spacemarkers$params$pattern_names
    if (is.null(pn)) return(NULL)
    pn <- intersect(pn, colnames(SummarizedExperiment::colData(x)))
    if (length(pn) == 0L) return(NULL)
    SummarizedExperiment::colData(x)[, pn, drop = FALSE]
})

#' @rdname spatial_patterns
#' @aliases spatial_patterns<-,SpaceMarkersExperiment-method
#' @param value A data.frame or DataFrame of spatial pattern values. Row names
#'   must match \code{colnames(x)}.
#' @export
setMethod("spatial_patterns<-", "SpaceMarkersExperiment", function(x, value) {
    value <- S4Vectors::DataFrame(value)
    cd <- SummarizedExperiment::colData(x)
    for (col in colnames(value)) {
        cd[[col]] <- value[[col]]
    }
    SummarizedExperiment::colData(x) <- cd
    sm <- x@spacemarkers
    if (is.null(sm$params)) sm$params <- list()
    sm$params$pattern_names <- colnames(value)
    x@spacemarkers <- sm
    x
})

#' @rdname spatial_params
#' @aliases spatial_params,SpaceMarkersExperiment-method
#' @param x A \code{SpaceMarkersExperiment} object.
#' @return The optimal kernel parameters matrix (2 x N with rows sigmaOpt,
#'   threshOpt), or NULL if not computed.
#' @export
setMethod("spatial_params", "SpaceMarkersExperiment", function(x) {
    x@spacemarkers$params$spatial_params
})

#' @rdname spatial_params
#' @aliases spatial_params<-,SpaceMarkersExperiment-method
#' @param value A matrix of spatial parameters (2 x N, rows sigmaOpt and
#'   threshOpt).
#' @export
setMethod("spatial_params<-", "SpaceMarkersExperiment", function(x, value) {
    sm <- x@spacemarkers
    if (is.null(sm$params)) sm$params <- list()
    sm$params$spatial_params <- value
    x@spacemarkers <- sm
    x
})

#' @rdname params
#' @aliases params,SpaceMarkersExperiment-method
#' @param x A \code{SpaceMarkersExperiment} object.
#' @return The full list of hyperparameters stored during analysis, including
#'   \code{spatial_params}, \code{pattern_names}, \code{min_gene_expr},
#'   \code{mode}, \code{analysis}, \code{min_overlap}, \code{directed},
#'   \code{outlier}, \code{null_samples}, \code{genes}, etc. Returns NULL if
#'   no parameters have been set.
#' @export
setMethod("params", "SpaceMarkersExperiment", function(x) {
    x@spacemarkers$params
})

#' @rdname hotspots
#' @aliases hotspots,SpaceMarkersExperiment-method
#' @param x A \code{SpaceMarkersExperiment} object.
#' @param type Character: one of "undirected", "pattern", or "influence".
#' @return A data.frame of hotspot assignments, or NULL.
#' @export
setMethod("hotspots", "SpaceMarkersExperiment",
    function(x, type = c("undirected", "pattern", "influence")) {
        type <- match.arg(type)
        S4Vectors::metadata(x)$hotspots[[type]]
    }
)

#' @rdname influence_map
#' @aliases influence_map,SpaceMarkersExperiment-method
#' @param x A \code{SpaceMarkersExperiment} object.
#' @return A data.frame of per-spot influence values, or NULL.
#' @export
setMethod("influence_map", "SpaceMarkersExperiment", function(x) {
    S4Vectors::metadata(x)$influence
})

#' @rdname interactions
#' @aliases interactions,SpaceMarkersExperiment-method
#' @param x A \code{SpaceMarkersExperiment} object.
#' @param pair Optional character string specifying a pattern pair name
#'   (e.g., "Pattern_1_Pattern_5"). If NULL, returns all pairs.
#' @return A named list of per-pair interaction results, or a single pair's
#'   results if \code{pair} is specified.
#' @export
setMethod("interactions", "SpaceMarkersExperiment",
    function(x, pair = NULL) {
        ints <- S4Vectors::metadata(x)$interactions
        if (is.null(pair)) return(ints)
        ints[[pair]]
    }
)

#' @rdname undirected_scores
#' @aliases undirected_scores,SpaceMarkersExperiment-method
#' @param x A \code{SpaceMarkersExperiment} object.
#' @return A data.frame of undirected interaction scores (genes x pattern
#'   pairs), or NULL.
#' @export
setMethod("undirected_scores", "SpaceMarkersExperiment", function(x) {
    x@spacemarkers$results$undirected_scores
})

#' @rdname directed_scores
#' @aliases directed_scores,SpaceMarkersExperiment-method
#' @param x A \code{SpaceMarkersExperiment} object.
#' @return A data.frame of directed interaction scores, or NULL.
#' @export
setMethod("directed_scores", "SpaceMarkersExperiment", function(x) {
    x@spacemarkers$results$directed_scores
})

#' @rdname lr_scores
#' @aliases lr_scores,SpaceMarkersExperiment-method
#' @param x A \code{SpaceMarkersExperiment} object.
#' @return A matrix of ligand-receptor pair scores, or NULL.
#' @export
setMethod("lr_scores", "SpaceMarkersExperiment", function(x) {
    x@spacemarkers$results$lr_scores
})

#' @rdname overlap_scores
#' @aliases overlap_scores,SpaceMarkersExperiment-method
#' @param x A \code{SpaceMarkersExperiment} object.
#' @return A data.frame of pattern overlap scores, or NULL.
#' @export
setMethod("overlap_scores", "SpaceMarkersExperiment", function(x) {
    x@spacemarkers$results$overlap_scores
})

#' @rdname analysis_type
#' @aliases analysis_type,SpaceMarkersExperiment-method
#' @param x A \code{SpaceMarkersExperiment} object.
#' @return A character string: "undirected", "directed", "both", or NULL.
#' @export
setMethod("analysis_type", "SpaceMarkersExperiment", function(x) {
    x@spacemarkers$analysis
})

# ---- Internal bridge helper ----

#' Reconstruct legacy spPatterns data.frame from a SpaceMarkersExperiment
#'
#' @param sme A \code{SpaceMarkersExperiment} object.
#' @return A data.frame with columns barcode, y, x, and pattern columns.
#' @keywords internal
.sme_to_spPatterns <- function(sme) {
    coords <- as.data.frame(SpatialExperiment::spatialCoords(sme))
    # spatialCoords returns a matrix with columns matching the original names
    # Ensure we have x and y columns
    if (!all(c("x", "y") %in% colnames(coords))) {
        # Try standard column renaming
        if (ncol(coords) >= 2) {
            colnames(coords) <- c("y", "x")
        }
    }
    coords$barcode <- colnames(sme)

    pn <- sme@spacemarkers$params$pattern_names
    if (!is.null(pn)) {
        pats <- as.data.frame(
            SummarizedExperiment::colData(sme)[, pn, drop = FALSE]
        )
        result <- cbind(coords[, c("barcode", "y", "x"), drop = FALSE], pats)
    } else {
        result <- coords[, c("barcode", "y", "x"), drop = FALSE]
    }
    rownames(result) <- result$barcode
    result
}

#' Extract expression matrix from a SpaceMarkersExperiment
#'
#' @param sme A \code{SpaceMarkersExperiment} object.
#' @param assay_name Character name of the assay to extract. Default "logcounts".
#' @return A matrix or sparse matrix of expression values.
#' @keywords internal
.sme_get_expr <- function(sme, assay_name = "logcounts") {
    if (assay_name %in% SummarizedExperiment::assayNames(sme)) {
        SummarizedExperiment::assay(sme, assay_name)
    } else {
        SummarizedExperiment::assay(sme, 1L)
    }
}

# ---- add_features ----

#' @title Add spatial features to a SpaceMarkersExperiment
#' @description Adds spatial features (patterns) to a
#'   \code{\link{SpaceMarkersExperiment}} object's colData after object
#'   initialization.  When spots differ between the SME and the features,
#'   only the common spots are kept and a warning is emitted.
#'
#' @param sme A \code{\link{SpaceMarkersExperiment}} object.
#' @param features A matrix or data.frame of features (rows = spots,
#'   cols = patterns), a file path, or an R object accepted by
#'   \code{\link{get_spatial_features}}.
#' @param ... Additional arguments passed to \code{get_spatial_features()} when
#'   \code{features} is a file path or object.
#' @return The \code{SpaceMarkersExperiment} with features added to
#'   \code{colData} and \code{pattern_names} updated in the
#'   \code{spacemarkers} slot.
#' @export
add_features <- function(sme, features, ...) {
    if (!is(sme, "SpaceMarkersExperiment")) {
        stop("'sme' must be a SpaceMarkersExperiment object.")
    }

    # Load features if needed
    if (is.character(features) ||
        (!is.data.frame(features) && !is.matrix(features))) {
        features <- get_spatial_features(features, ...)
    }
    features <- as.data.frame(features)

    # Spot alignment: intersect + warn
    sme_spots <- colnames(sme)
    feat_spots <- rownames(features)
    common <- intersect(sme_spots, feat_spots)

    if (length(common) == 0L) {
        stop("No common spots between SME and features.")
    }

    dropped_sme <- length(sme_spots) - length(common)
    dropped_feat <- length(feat_spots) - length(common)
    if (dropped_sme > 0L || dropped_feat > 0L) {
        warning(sprintf(
            paste0("Spot mismatch: %d spots in SME and %d spots in features ",
                   "not shared. Keeping %d common spots."),
            dropped_sme, dropped_feat, length(common)
        ))
        sme <- sme[, common]
        features <- features[common, , drop = FALSE]
    }

    # Add to colData
    for (col in colnames(features)) {
        SummarizedExperiment::colData(sme)[[col]] <- features[[col]]
    }

    # Update pattern_names in spacemarkers slot
    sm <- sme@spacemarkers
    if (is.null(sm$params)) sm$params <- list()
    sm$params$pattern_names <- colnames(features)
    sme@spacemarkers <- sm

    sme
}

# NULL-coalescing helper (used by SME pipeline methods for optional-arg fallbacks)
`%||%` <- function(a, b) if (!is.null(a)) a else b

# Unexported helpers used by SME pipeline methods
.sme_expr <- function(sme) {
    if ("logcounts" %in% SummarizedExperiment::assayNames(sme))
        SummarizedExperiment::assay(sme, "logcounts")
    else
        SummarizedExperiment::assay(sme, 1L)
}

.sme_spPatterns <- function(sme) {
    coords <- as.data.frame(SpatialExperiment::spatialCoords(sme))
    pats <- as.data.frame(spatial_patterns(sme))
    df <- data.frame(barcode = colnames(sme),
                     coords, pats,
                     check.names = FALSE,
                     stringsAsFactors = FALSE)
    rownames(df) <- colnames(sme)
    df
}

#' @rdname hotspots
#' @aliases hotspots<-,SpaceMarkersExperiment-method
#' @param value A data.frame of hotspot assignments.
#' @export
setMethod("hotspots<-", "SpaceMarkersExperiment",
    function(x, type = c("undirected", "pattern", "influence"), value) {
        type <- match.arg(type)
        md <- S4Vectors::metadata(x)
        if (is.null(md$hotspots)) md$hotspots <- list()
        md$hotspots[[type]] <- value
        S4Vectors::metadata(x) <- md
        x
    }
)

#' @rdname interactions
#' @aliases interactions<-,SpaceMarkersExperiment-method
#' @param value A named list of per-pair interaction results.
#' @export
setMethod("interactions<-", "SpaceMarkersExperiment", function(x, value) {
    md <- S4Vectors::metadata(x)
    md$interactions <- value
    S4Vectors::metadata(x) <- md
    x
})

#' @rdname influence_map
#' @aliases influence_map<-,SpaceMarkersExperiment-method
#' @param value A data.frame of per-spot influence values.
#' @export
setMethod("influence_map<-", "SpaceMarkersExperiment", function(x, value) {
    md <- S4Vectors::metadata(x)
    md$influence <- value
    S4Vectors::metadata(x) <- md
    x
})

#' @rdname undirected_scores
#' @aliases undirected_scores<-,SpaceMarkersExperiment-method
#' @param value A data.frame of undirected interaction scores.
#' @export
setMethod("undirected_scores<-", "SpaceMarkersExperiment", function(x, value) {
    sm <- x@spacemarkers
    if (is.null(sm$results)) sm$results <- list()
    sm$results$undirected_scores <- value
    x@spacemarkers <- sm
    x
})

#' @rdname directed_scores
#' @aliases directed_scores<-,SpaceMarkersExperiment-method
#' @param value A data.frame of directed interaction scores.
#' @export
setMethod("directed_scores<-", "SpaceMarkersExperiment", function(x, value) {
    sm <- x@spacemarkers
    if (is.null(sm$results)) sm$results <- list()
    sm$results$directed_scores <- value
    x@spacemarkers <- sm
    x
})

#' @rdname lr_scores
#' @aliases lr_scores<-,SpaceMarkersExperiment-method
#' @param value A matrix of ligand-receptor pair scores.
#' @export
setMethod("lr_scores<-", "SpaceMarkersExperiment", function(x, value) {
    sm <- x@spacemarkers
    if (is.null(sm$results)) sm$results <- list()
    sm$results$lr_scores <- value
    x@spacemarkers <- sm
    x
})

#' @rdname overlap_scores
#' @aliases overlap_scores<-,SpaceMarkersExperiment-method
#' @param value A data.frame of pattern overlap scores.
#' @export
setMethod("overlap_scores<-", "SpaceMarkersExperiment", function(x, value) {
    sm <- x@spacemarkers
    if (is.null(sm$results)) sm$results <- list()
    sm$results$overlap_scores <- value
    x@spacemarkers <- sm
    x
})

#' @rdname analysis_type
#' @aliases analysis_type<-,SpaceMarkersExperiment-method
#' @param value Character: "undirected", "directed", or "both".
#' @export
setMethod("analysis_type<-", "SpaceMarkersExperiment", function(x, value) {
    sm <- x@spacemarkers
    sm$analysis <- value
    x@spacemarkers <- sm
    methods::validObject(x)
    x
})

#' @rdname find_all_hotspots
#' @aliases find_all_hotspots,SpaceMarkersExperiment-method
#' @export
setMethod("find_all_hotspots", "SpaceMarkersExperiment",
    function(spPatterns, params = NULL, outlier = "positive",
             nullSamples = 1000, includeSelf = TRUE, ...) {
        sme <- spPatterns
        if (is.null(params)) params <- spatial_params(sme)
        hs <- find_all_hotspots(.sme_spPatterns(sme),
                                params = params, outlier = outlier,
                                nullSamples = nullSamples,
                                includeSelf = includeSelf, ...)
        hotspots(sme, type = "undirected") <- hs
        sme
    }
)

#' @rdname get_pairwise_interacting_genes
#' @aliases get_pairwise_interacting_genes,SpaceMarkersExperiment-method
#' @export
setMethod("get_pairwise_interacting_genes", "SpaceMarkersExperiment",
    function(data, spPatterns = NULL, mode = "residual", optParams = NULL,
             reconstruction = NULL, hotspots = NULL, minOverlap = 50,
             analysis = c("overlap", "enrichment"),
             pattern_pairs = NULL, ..., workers = 1) {
        sme <- data
        hs <- hotspots %||% hotspots(sme, "undirected")
        if (is.null(hs)) {
            stop("Run find_all_hotspots(x) before get_pairwise_interacting_genes().")
        }
        analysis <- match.arg(analysis)
        res <- get_pairwise_interacting_genes(
            data = .sme_expr(sme),
            spPatterns = spPatterns %||% .sme_spPatterns(sme),
            mode = mode,
            optParams = optParams %||% spatial_params(sme),
            reconstruction = reconstruction,
            hotspots = hs,
            minOverlap = minOverlap,
            analysis = analysis,
            pattern_pairs = pattern_pairs,
            workers = workers,
            ...
        )
        interactions(sme) <- res
        sm <- sme@spacemarkers
        if (is.null(sm$params)) sm$params <- list()
        sm$params$mode <- mode
        sm$params$analysis_method <- analysis
        sm$params$min_overlap <- minOverlap
        sme@spacemarkers <- sm
        sme
    }
)

#' @rdname get_im_scores
#' @aliases get_im_scores,SpaceMarkersExperiment-method
#' @export
setMethod("get_im_scores", "SpaceMarkersExperiment", function(SpaceMarkers) {
    sme <- SpaceMarkers
    ints <- interactions(sme)
    if (is.null(ints)) {
        stop("Run get_pairwise_interacting_genes(x) before get_im_scores().")
    }
    undirected_scores(sme) <- get_im_scores(ints)
    sme
})

#' @rdname calculate_overlap_undirected
#' @aliases calculate_overlap_undirected,SpaceMarkersExperiment-method
#' @export
setMethod("calculate_overlap_undirected", "SpaceMarkersExperiment",
    function(hotspots, patternList = NULL,
             method = c("Szymkiewicz-Simpson", "Jaccard", "Sorensen-Dice",
                        "Ochiai", "absolute")) {
        sme <- hotspots
        hs <- hotspots(sme, "undirected")
        if (is.null(hs)) {
            stop("Run find_all_hotspots(x) before calculate_overlap_undirected().")
        }
        ov <- calculate_overlap_undirected(hs, patternList = patternList,
                                           method = method)
        overlap_scores(sme) <- ov
        cur <- analysis_type(sme)
        new_type <- if (identical(cur, "directed") || identical(cur, "both"))
            "both" else "undirected"
        analysis_type(sme) <- new_type
        sme
    }
)

#' @rdname calculate_influence
#' @aliases calculate_influence,SpaceMarkersExperiment-method
#' @export
setMethod("calculate_influence", "SpaceMarkersExperiment",
    function(spPatterns, optParams = NULL, ...) {
        sme <- spPatterns
        op <- optParams %||% spatial_params(sme)
        inf <- calculate_influence(.sme_spPatterns(sme), optParams = op, ...)
        influence_map(sme) <- inf
        sme
    }
)

#' @rdname find_hotspots_gmm
#' @aliases find_hotspots_gmm,SpaceMarkersExperiment-method
#' @export
setMethod("find_hotspots_gmm", "SpaceMarkersExperiment",
    function(df, threshold = 0.1, ..., type = c("pattern", "influence"),
             minvals = NULL, maxvals = NULL) {
        sme <- df
        type <- match.arg(type)
        src <- if (type == "pattern") .sme_spPatterns(sme) else influence_map(sme)
        if (is.null(src)) {
            stop("Run calculate_influence(x) before find_hotspots_gmm(x, 'influence').")
        }
        # Choose sensible default thresholds per type when user didn't pass one
        thr <- if (missing(threshold) || identical(threshold, 0.1)) {
            default_min <- if (type == "pattern") 0.1 else 0.01
            default_max <- if (type == "pattern") 0.8 else 0.5
            calculate_thresholds(src,
                minvals = minvals %||% default_min,
                maxvals = maxvals %||% default_max)
        } else threshold
        hs <- find_hotspots_gmm(src, threshold = thr, ...)
        hotspots(sme, type = type) <- hs
        sme
    }
)

#' @rdname calculate_gene_scores_directed
#' @aliases calculate_gene_scores_directed,SpaceMarkersExperiment-method
#' @export
setMethod("calculate_gene_scores_directed", "SpaceMarkersExperiment",
    function(data, pat_hotspots = NULL, influence_hotspots = NULL,
             pattern_pairs = NULL, ...) {
        sme <- data
        pat <- pat_hotspots %||% hotspots(sme, "pattern")
        inf <- influence_hotspots %||% hotspots(sme, "influence")
        if (is.null(pat) || is.null(inf)) {
            stop("Run find_hotspots_gmm(x, 'pattern') and find_hotspots_gmm(x, 'influence') before calculate_gene_scores_directed().")
        }
        if (is.null(pattern_pairs)) {
            patnames <- setdiff(colnames(.sme_spPatterns(sme)),
                                c("barcode", "x", "y"))
            pattern_pairs <- t(utils::combn(patnames, 2))
        }
        scores <- calculate_gene_scores_directed(
            data = .sme_expr(sme),
            pat_hotspots = pat,
            influence_hotspots = inf,
            pattern_pairs = pattern_pairs,
            ...
        )
        directed_scores(sme) <- scores
        sme
    }
)
