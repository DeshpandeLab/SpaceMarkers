#' Internal imports for SpaceMarkersExperiment methods
#'
#' Documentation stub whose only job is to declare package imports that
#' roxygen2 turns into \code{importFrom} entries in \code{NAMESPACE}.
#'
#' @importFrom S4Vectors metadata metadata<-
#' @importFrom SummarizedExperiment assayNames colData<-
#' @importFrom methods setGeneric validObject
#' @name SpaceMarkersExperiment-imports
#' @return No return value; used only for its roxygen side effects.
#' @keywords internal
NULL

#' @include AllClasses.R AllGenerics.R
NULL

#' Access hotspot assignments on a SpaceMarkersExperiment
#'
#' Read (\code{hotspots(x, type)}) or write
#' (\code{hotspots(x, type) <- value}) the per-type hotspot data.frame stored
#' in \code{metadata(x)$hotspots[[type]]}. \code{type} selects one of
#' \code{"undirected"}, \code{"pattern"}, or \code{"influence"}.
#'
#' @name hotspots
#' @aliases hotspots hotspots<-
#' @param x A \code{SpaceMarkersExperiment}.
#' @param type Character; one of \code{"undirected"}, \code{"pattern"},
#'   \code{"influence"}.
#' @param value A data.frame of hotspot assignments.
#' @param ... Not used.
#' @return For the getter, a data.frame or \code{NULL}. For the setter, the
#'   modified \code{SpaceMarkersExperiment}.
#' @examples
#' sme <- SpaceMarkersExperiment(
#'     assays = list(logcounts = matrix(rpois(200, 2), 10, 20,
#'         dimnames = list(paste0("G", 1:10), paste0("s", 1:20)))),
#'     spatialCoords = matrix(runif(40), 20, 2,
#'         dimnames = list(NULL, c("y", "x"))))
#' colnames(sme) <- paste0("s", seq_len(20))
#' hs <- data.frame(barcode = colnames(sme), x = 1:20, y = 1:20,
#'                  Pattern_1 = c(rep("hot", 10), rep(NA, 10)))
#' hotspots(sme, type = "undirected") <- hs
#' head(hotspots(sme, "undirected"))
NULL

#' Access pairwise interaction results on a SpaceMarkersExperiment
#'
#' Read (\code{interactions(x, pair)}) or write
#' (\code{interactions(x) <- value}) the per-pattern-pair interaction results
#' stored in \code{metadata(x)$interactions}. When \code{pair} is supplied to
#' the getter, only that pair's results are returned; otherwise the full
#' named list is returned.
#'
#' @name interactions
#' @aliases interactions interactions<-
#' @param x A \code{SpaceMarkersExperiment}.
#' @param pair Optional character string specifying a pattern-pair name
#'   (e.g. \code{"Pattern_1_Pattern_5"}). If \code{NULL}, all pairs are
#'   returned.
#' @param value A named list of per-pair interaction results.
#' @param ... Not used.
#' @return For the getter, a named list of per-pair results (or a single
#'   pair's results when \code{pair} is specified). For the setter, the
#'   modified \code{SpaceMarkersExperiment}.
#' @examples
#' sme <- SpaceMarkersExperiment(
#'     assays = list(logcounts = matrix(rpois(200, 2), 10, 20,
#'         dimnames = list(paste0("G", 1:10), paste0("s", 1:20)))),
#'     spatialCoords = matrix(runif(40), 20, 2,
#'         dimnames = list(NULL, c("y", "x"))))
#' interactions(sme) <- list(
#'     Pattern_1_Pattern_2 = list(interacting_genes = list()))
#' names(interactions(sme))
NULL

#' Access the per-spot influence map on a SpaceMarkersExperiment
#'
#' Read (\code{influence_map(x)}) or write
#' (\code{influence_map(x) <- value}) the per-spot influence values stored in
#' \code{metadata(x)$influence}.
#'
#' @name influence_map
#' @aliases influence_map influence_map<-
#' @param x A \code{SpaceMarkersExperiment}.
#' @param value A data.frame of per-spot influence values.
#' @param ... Not used.
#' @return For the getter, a data.frame of per-spot influence values or
#'   \code{NULL}. For the setter, the modified \code{SpaceMarkersExperiment}.
#' @examples
#' sme <- SpaceMarkersExperiment(
#'     assays = list(logcounts = matrix(rpois(200, 2), 10, 20,
#'         dimnames = list(paste0("G", 1:10), paste0("s", 1:20)))),
#'     spatialCoords = matrix(runif(40), 20, 2,
#'         dimnames = list(NULL, c("y", "x"))))
#' colnames(sme) <- paste0("s", seq_len(20))
#' influence_map(sme) <- data.frame(barcode = colnames(sme),
#'     x = 1:20, y = 1:20, Pattern_1 = runif(20))
#' head(influence_map(sme))
NULL

#' Access undirected interaction scores on a SpaceMarkersExperiment
#'
#' Read (\code{undirected_scores(x)}) or write
#' (\code{undirected_scores(x) <- value}) the data.frame of undirected
#' interaction scores stored in \code{x@spacemarkers$results$undirected_scores}.
#'
#' @name undirected_scores
#' @aliases undirected_scores undirected_scores<-
#' @param x A \code{SpaceMarkersExperiment}.
#' @param value A data.frame of undirected interaction scores
#'   (genes x pattern pairs).
#' @param ... Not used.
#' @return For the getter, a data.frame of undirected scores or \code{NULL}.
#'   For the setter, the modified \code{SpaceMarkersExperiment}.
#' @examples
#' sme <- SpaceMarkersExperiment(
#'     assays = list(logcounts = matrix(rpois(200, 2), 10, 20,
#'         dimnames = list(paste0("G", 1:10), paste0("s", 1:20)))),
#'     spatialCoords = matrix(runif(40), 20, 2,
#'         dimnames = list(NULL, c("y", "x"))))
#' undirected_scores(sme) <- data.frame(
#'     Gene = paste0("G", 1:5), Pattern_1_Pattern_2 = runif(5))
#' head(undirected_scores(sme))
NULL

#' Access directed interaction scores on a SpaceMarkersExperiment
#'
#' Read (\code{directed_scores(x)}) or write
#' (\code{directed_scores(x) <- value}) the data.frame of directed
#' interaction scores stored in \code{x@spacemarkers$results$directed_scores}.
#'
#' @name directed_scores
#' @aliases directed_scores directed_scores<-
#' @param x A \code{SpaceMarkersExperiment}.
#' @param value A data.frame of directed interaction scores.
#' @param ... Not used.
#' @return For the getter, a data.frame of directed scores or \code{NULL}.
#'   For the setter, the modified \code{SpaceMarkersExperiment}.
#' @examples
#' sme <- SpaceMarkersExperiment(
#'     assays = list(logcounts = matrix(rpois(200, 2), 10, 20,
#'         dimnames = list(paste0("G", 1:10), paste0("s", 1:20)))),
#'     spatialCoords = matrix(runif(40), 20, 2,
#'         dimnames = list(NULL, c("y", "x"))))
#' directed_scores(sme) <- data.frame(
#'     gene = paste0("G", 1:5),
#'     cell_interaction = "Pattern_1->Pattern_5",
#'     effect_size = runif(5))
#' head(directed_scores(sme))
NULL

#' Access ligand-receptor scores on a SpaceMarkersExperiment
#'
#' Read (\code{lr_scores(x)}) or write (\code{lr_scores(x) <- value}) the
#' matrix of ligand-receptor pair scores stored in
#' \code{x@spacemarkers$results$lr_scores}.
#'
#' @name lr_scores
#' @aliases lr_scores lr_scores<-
#' @param x A \code{SpaceMarkersExperiment}.
#' @param value A matrix of ligand-receptor pair scores.
#' @param ... Not used.
#' @return For the getter, a matrix of LR pair scores or \code{NULL}. For the
#'   setter, the modified \code{SpaceMarkersExperiment}.
#' @examples
#' sme <- SpaceMarkersExperiment(
#'     assays = list(logcounts = matrix(rpois(200, 2), 10, 20,
#'         dimnames = list(paste0("G", 1:10), paste0("s", 1:20)))),
#'     spatialCoords = matrix(runif(40), 20, 2,
#'         dimnames = list(NULL, c("y", "x"))))
#' lr_scores(sme) <- matrix(runif(4), 2, 2,
#'     dimnames = list(c("LR1", "LR2"), c("P1_P2", "P2_P1")))
#' lr_scores(sme)
NULL

#' Access pattern overlap scores on a SpaceMarkersExperiment
#'
#' Read (\code{overlap_scores(x)}) or write
#' (\code{overlap_scores(x) <- value}) the data.frame of pattern overlap
#' scores stored in \code{x@spacemarkers$results$overlap_scores}.
#'
#' @name overlap_scores
#' @aliases overlap_scores overlap_scores<-
#' @param x A \code{SpaceMarkersExperiment}.
#' @param value A data.frame of pattern overlap scores.
#' @param ... Not used.
#' @return For the getter, a data.frame of pattern overlap scores or
#'   \code{NULL}. For the setter, the modified \code{SpaceMarkersExperiment}.
#' @examples
#' sme <- SpaceMarkersExperiment(
#'     assays = list(logcounts = matrix(rpois(200, 2), 10, 20,
#'         dimnames = list(paste0("G", 1:10), paste0("s", 1:20)))),
#'     spatialCoords = matrix(runif(40), 20, 2,
#'         dimnames = list(NULL, c("y", "x"))))
#' overlap_scores(sme) <- data.frame(pattern1 = "Pattern_1",
#'     pattern2 = "Pattern_2", overlapScore = 0.5)
#' overlap_scores(sme)
NULL

#' Access or set the analysis type on a SpaceMarkersExperiment
#'
#' Read (\code{analysis_type(x)}) or write
#' (\code{analysis_type(x) <- value}) the analysis-type tag stored in
#' \code{x@spacemarkers$analysis}. Valid values are \code{"undirected"},
#' \code{"directed"}, or \code{"both"}. The setter delegates validation to
#' \code{setValidity("SpaceMarkersExperiment", ...)}.
#'
#' @name analysis_type
#' @aliases analysis_type analysis_type<-
#' @param x A \code{SpaceMarkersExperiment}.
#' @param value Character; one of \code{"undirected"}, \code{"directed"},
#'   or \code{"both"}.
#' @param ... Not used.
#' @return For the getter, a character string or \code{NULL}. For the setter,
#'   the modified \code{SpaceMarkersExperiment}.
#' @examples
#' sme <- SpaceMarkersExperiment(
#'     assays = list(logcounts = matrix(rpois(200, 2), 10, 20,
#'         dimnames = list(paste0("G", 1:10), paste0("s", 1:20)))),
#'     spatialCoords = matrix(runif(40), 20, 2,
#'         dimnames = list(NULL, c("y", "x"))))
#' analysis_type(sme) <- "undirected"
#' analysis_type(sme)
NULL

#' Access analysis hyperparameters on a SpaceMarkersExperiment
#'
#' Read-only accessor that returns the full list stored in
#' \code{x@spacemarkers$params}, including \code{spatial_params},
#' \code{pattern_names}, \code{min_gene_expr}, \code{mode}, \code{analysis},
#' \code{min_overlap}, \code{directed}, \code{outlier}, \code{null_samples},
#' \code{genes}, and any other hyperparameters set during analysis.
#'
#' @name params
#' @aliases params
#' @param x A \code{SpaceMarkersExperiment}.
#' @param ... Not used.
#' @return The full list of hyperparameters, or \code{NULL} if none has been
#'   set.
#' @examples
#' sme <- SpaceMarkersExperiment(
#'     assays = list(logcounts = matrix(rpois(200, 2), 10, 20,
#'         dimnames = list(paste0("G", 1:10), paste0("s", 1:20)))),
#'     spatialCoords = matrix(runif(40), 20, 2,
#'         dimnames = list(NULL, c("y", "x"))))
#' params(sme)
NULL

#' Access or set the spatial-kernel parameters on a SpaceMarkersExperiment
#'
#' Read (\code{spatial_params(x)}) or write
#' (\code{spatial_params(x) <- value}) the 2 x N kernel-parameter matrix
#' (rows \code{sigmaOpt}, \code{threshOpt}; one column per pattern) stored in
#' \code{x@spacemarkers$params$spatial_params}.
#'
#' @name spatial_params
#' @aliases spatial_params spatial_params<-
#' @param x A \code{SpaceMarkersExperiment}.
#' @param value A 2 x N numeric matrix of spatial kernel parameters (rows
#'   \code{sigmaOpt}, \code{threshOpt}).
#' @param ... Not used.
#' @return For the getter, a 2 x N matrix of spatial parameters or
#'   \code{NULL}. For the setter, the modified \code{SpaceMarkersExperiment}.
#' @examples
#' sme <- SpaceMarkersExperiment(
#'     assays = list(logcounts = matrix(rpois(200, 2), 10, 20,
#'         dimnames = list(paste0("G", 1:10), paste0("s", 1:20)))),
#'     spatialCoords = matrix(runif(40), 20, 2,
#'         dimnames = list(NULL, c("y", "x"))))
#' spatial_params(sme) <- matrix(c(1, 2, 1, 2), 2, 2,
#'     dimnames = list(c("sigmaOpt", "threshOpt"),
#'                     c("Pattern_1", "Pattern_2")))
#' spatial_params(sme)
NULL

#' Access or set spatial pattern columns on a SpaceMarkersExperiment
#'
#' Read (\code{spatial_patterns(x)}) or write
#' (\code{spatial_patterns(x) <- value}) the per-spot spatial-pattern values
#' stored in the columns of \code{colData(x)} named by
#' \code{x@spacemarkers$params$pattern_names}. The setter also updates
#' \code{pattern_names} to the column names of \code{value}.
#'
#' @name spatial_patterns
#' @aliases spatial_patterns spatial_patterns<-
#' @param x A \code{SpaceMarkersExperiment}.
#' @param value A data.frame or \code{DataFrame} of spatial pattern values
#'   (rows matching \code{colnames(x)}).
#' @param ... Not used.
#' @return For the getter, a \code{DataFrame} of spatial patterns or
#'   \code{NULL}. For the setter, the modified \code{SpaceMarkersExperiment}.
#' @examples
#' sme <- SpaceMarkersExperiment(
#'     assays = list(logcounts = matrix(rpois(200, 2), 10, 20,
#'         dimnames = list(paste0("G", 1:10), paste0("s", 1:20)))),
#'     spatialCoords = matrix(runif(40), 20, 2,
#'         dimnames = list(NULL, c("y", "x"))))
#' colnames(sme) <- paste0("s", seq_len(20))
#' spatial_patterns(sme) <- S4Vectors::DataFrame(
#'     Pattern_1 = runif(20), Pattern_2 = runif(20),
#'     row.names = colnames(sme))
#' head(spatial_patterns(sme))
NULL

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

# ---- overlap_map ----

#' Per-spot overlap classification map for a pattern pair
#'
#' The spot-level counterpart to \code{\link{overlap_scores}}. Returns a
#' factor with one label per spot, summarizing how each spot relates to a
#' pattern pair given the hotspot information stored on the SME.
#'
#' \strong{Undirected} (when \code{hotspots(x, "undirected")} is present and
#'   \code{directed = FALSE} or auto-detected). Three levels: \code{<pat1>},
#'   \code{"interacting"} (override via \code{interaction_label}),
#'   \code{<pat2>}.
#'
#' \strong{Directed} (when \code{hotspots(x, "pattern")} AND
#'   \code{hotspots(x, "influence")} are both present). Because a directed
#'   analysis runs a separate t-test per direction (\code{pat1 near pat2}
#'   vs \code{pat1}, and \code{pat2 near pat1} vs \code{pat2}), the helper
#'   returns one direction at a time — the three-level factor that
#'   directly matches that direction's t-test:
#'   \itemize{
#'     \item \code{direction = "forward"} (default; source = \code{pat1},
#'       target = \code{pat2}). Levels: \code{<pat1>} (the "control" in
#'       the t-test), \code{"<pat1> near <pat2>"} (the interaction region
#'       - the "treatment"), \code{"<pat2> influence"} (context: where
#'       \code{pat2}'s influence extends beyond \code{pat1}'s hotspot).
#'     \item \code{direction = "reverse"} (source = \code{pat2}, target =
#'       \code{pat1}). Levels are the mirror image: \code{<pat2>},
#'       \code{"<pat2> near <pat1>"}, \code{"<pat1> influence"}.
#'   }
#'
#' Downstream callers (including \code{\link{plot_spatial}(..., source =
#' "interaction")}) use this helper so the plot and any user-side
#' analysis stay consistent.
#'
#' @param x A \code{SpaceMarkersExperiment}.
#' @param interaction_patterns Length-2 character vector of pattern names
#'   \code{c(pat1, pat2)}. For directed mode these are interpreted as
#'   \code{c(source, target)} under \code{direction = "forward"}.
#' @param direction One of \code{"forward"} or \code{"reverse"} — controls
#'   which direction's classification is returned in directed mode. Ignored
#'   for undirected.
#' @param directed Logical or \code{NULL}. When \code{NULL} (default),
#'   picks directed mode iff pattern and influence hotspots are both
#'   stored on the SME.
#' @param interaction_label Optional override for the middle label in
#'   undirected mode (default \code{"interacting"}). Ignored in directed
#'   mode where labels are derived from \code{interaction_patterns}.
#' @return A factor with one entry per spot, \code{NA} for spots outside
#'   any hotspot. Names are the spot barcodes.
#' @export
#' @examples
#' sme <- SpaceMarkersExperiment(
#'     assays = list(logcounts = matrix(rpois(200, 2), 10, 20,
#'         dimnames = list(paste0("G", 1:10), paste0("s", 1:20)))),
#'     spatialCoords = matrix(runif(40), 20, 2,
#'         dimnames = list(NULL, c("y", "x"))))
#' colnames(sme) <- paste0("s", seq_len(20))
#' hs <- data.frame(barcode = colnames(sme), x = 1:20, y = 1:20,
#'                  Pattern_1 = c(rep("hot", 10), rep(NA, 10)),
#'                  Pattern_2 = c(rep(NA, 5), rep("hot", 10), rep(NA, 5)))
#' hotspots(sme, type = "undirected") <- hs
#' table(overlap_map(sme, c("Pattern_1", "Pattern_2")))
overlap_map <- function(x, interaction_patterns,
                        direction = c("forward", "reverse"),
                        directed = NULL,
                        interaction_label = NULL) {
    stopifnot(is(x, "SpaceMarkersExperiment"))
    if (length(interaction_patterns) != 2L) {
        stop("interaction_patterns must be length-2.")
    }
    direction <- match.arg(direction)
    pat1 <- interaction_patterns[1]; pat2 <- interaction_patterns[2]

    pat_hs <- hotspots(x, "pattern")
    inf_hs <- hotspots(x, "influence")
    und_hs <- hotspots(x, "undirected")

    if (is.null(directed)) {
        directed <- !is.null(pat_hs) && !is.null(inf_hs)
    }

    if (directed) {
        if (is.null(pat_hs) || is.null(inf_hs)) {
            stop("Directed classification needs both pattern and influence hotspots; run calculate_influence() + find_hotspots_gmm(type='pattern') and find_hotspots_gmm(type='influence').")
        }
        for (p in interaction_patterns) {
            if (!(p %in% colnames(pat_hs)))
                stop(sprintf("Pattern '%s' not in pattern hotspots.", p))
            if (!(p %in% colnames(inf_hs)))
                stop(sprintf("Pattern '%s' not in influence hotspots.", p))
        }
        # Source / target depend on direction
        if (direction == "forward") {
            src <- pat1; tgt <- pat2
        } else {
            src <- pat2; tgt <- pat1
        }
        src_pat <- !is.na(pat_hs[[src]])
        tgt_inf <- !is.na(inf_hs[[tgt]])

        src_only          <- src
        src_near_tgt      <- sprintf("%s near %s", src, tgt)
        tgt_influence_lbl <- sprintf("%s influence", tgt)

        labs <- rep(NA_character_, length(src_pat))
        labs[src_pat & !tgt_inf] <- src_only
        labs[src_pat &  tgt_inf] <- src_near_tgt
        labs[!src_pat &  tgt_inf] <- tgt_influence_lbl
        lvls <- c(src_only, src_near_tgt, tgt_influence_lbl)
        out <- factor(labs, levels = lvls)
        nm <- pat_hs$barcode %||% rownames(pat_hs) %||% colnames(x)
        names(out) <- nm
        return(out)
    }

    # Undirected
    if (is.null(und_hs)) {
        stop("Undirected classification needs undirected hotspots; run find_all_hotspots() first.")
    }
    for (p in interaction_patterns) {
        if (!(p %in% colnames(und_hs)))
            stop(sprintf("Pattern '%s' not in undirected hotspots.", p))
    }
    mid <- interaction_label %||% "interacting"
    in1 <- !is.na(und_hs[[pat1]])
    in2 <- !is.na(und_hs[[pat2]])
    labs <- rep(NA_character_, length(in1))
    labs[in1 & !in2] <- pat1
    labs[!in1 & in2] <- pat2
    labs[in1 &  in2] <- mid
    out <- factor(labs, levels = c(pat1, mid, pat2))
    names(out) <- und_hs$barcode %||% rownames(und_hs) %||% colnames(x)
    out
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
    function(data, spPatterns = NULL, mode = c("DE", "residual"),
             optParams = NULL, reconstruction = NULL, hotspots = NULL,
             minOverlap = 50,
             analysis = c("enrichment", "overlap"),
             pattern_pairs = NULL, ..., workers = NULL) {
        sme <- data
        hs <- hotspots %||% hotspots(sme, "undirected")
        if (is.null(hs)) {
            stop("Run find_all_hotspots(x) before get_pairwise_interacting_genes().")
        }
        mode <- match.arg(mode)
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
#' @param type Character; one of \code{"pattern"} or \code{"influence"}.
#'   Selects which SME slot to operate on.
#' @param minvals,maxvals Numeric thresholds passed to
#'   \code{calculate_thresholds} when the default \code{threshold = 0.1} is
#'   used.
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
            # Read pattern names from params() instead of materializing the full
            # spPatterns data.frame just to take colnames().
            patnames <- params(sme)$pattern_names %||%
                setdiff(colnames(SummarizedExperiment::colData(sme)),
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

#' @rdname calculate_overlap_directed
#' @aliases calculate_overlap_directed,SpaceMarkersExperiment-method
#' @export
setMethod("calculate_overlap_directed", "SpaceMarkersExperiment",
    function(pat_hotspots, influence_hotspots = NULL,
             patternList = NULL,
             method = c("relative-abundance", "differential-abundance", "absolute")) {
        sme <- pat_hotspots
        pat <- hotspots(sme, "pattern")
        inf <- influence_hotspots %||% hotspots(sme, "influence")
        if (is.null(pat) || is.null(inf)) {
            stop("Run find_hotspots_gmm(x, 'pattern') and find_hotspots_gmm(x, 'influence') before calculate_overlap_directed().")
        }
        ov <- calculate_overlap_directed(pat, inf,
                                         patternList = patternList,
                                         method = method)
        overlap_scores(sme) <- ov
        cur <- analysis_type(sme)
        new_type <- if (identical(cur, "undirected") || identical(cur, "both"))
            "both" else "directed"
        analysis_type(sme) <- new_type
        sme
    }
)

#' @rdname calculate_gene_set_score
#' @aliases calculate_gene_set_score,SpaceMarkersExperiment-method
#' @export
setMethod("calculate_gene_set_score", "SpaceMarkersExperiment",
    function(IMscores, gene_sets = NULL, weighted = TRUE,
             method = c("geometric_mean", "arithmetic_mean")) {
        sme <- IMscores
        method <- match.arg(method)
        scores <- directed_scores(sme)
        if (is.null(scores)) {
            stop("Run calculate_gene_scores_directed(x) before calculate_gene_set_score().")
        }
        gs <- gene_sets %||% params(sme)$lr_pairs$ligand.symbol
        if (is.null(gs)) {
            stop("No gene_sets supplied and no params(x)$lr_pairs$ligand.symbol available.")
        }
        ls <- calculate_gene_set_score(scores, gene_sets = gs,
                                       weighted = weighted, method = method)
        md <- S4Vectors::metadata(sme)
        md$ligand_scores <- ls
        S4Vectors::metadata(sme) <- md
        sme
    }
)

#' @rdname calculate_gene_set_specificity
#' @aliases calculate_gene_set_specificity,SpaceMarkersExperiment-method
#' @export
setMethod("calculate_gene_set_specificity", "SpaceMarkersExperiment",
    function(data, spPatterns = NULL, gene_sets = NULL, weighted = TRUE,
             method = c("geometric_mean", "arithmetic_mean")) {
        sme <- data
        method <- match.arg(method)
        gs <- gene_sets %||% params(sme)$lr_pairs$receptor.symbol
        if (is.null(gs)) {
            stop("No gene_sets supplied and no params(x)$lr_pairs$receptor.symbol available.")
        }
        rs <- calculate_gene_set_specificity(
            data = .sme_expr(sme),
            spPatterns = spPatterns %||% .sme_spPatterns(sme),
            gene_sets = gs,
            weighted = weighted,
            method = method
        )
        md <- S4Vectors::metadata(sme)
        md$receptor_scores <- rs
        S4Vectors::metadata(sme) <- md
        sme
    }
)

#' @rdname calculate_lr_scores
#' @aliases calculate_lr_scores,SpaceMarkersExperiment-method
#' @export
setMethod("calculate_lr_scores", "SpaceMarkersExperiment",
    function(ligand_scores, receptor_scores = NULL, lr_pairs = NULL,
             ligand_test = NULL, method = "geometric_mean", weighted = TRUE) {
        sme <- ligand_scores
        md <- S4Vectors::metadata(sme)
        ls <- md$ligand_scores
        rs <- receptor_scores %||% md$receptor_scores
        lr <- lr_pairs %||% params(sme)$lr_pairs
        if (is.null(ls) || is.null(rs) || is.null(lr)) {
            stop("Run calculate_gene_set_score(x) and calculate_gene_set_specificity(x), and ensure params(x)$lr_pairs is set, before calculate_lr_scores().")
        }
        out <- calculate_lr_scores(ls, rs, lr, ligand_test = ligand_test,
                                   method = method, weighted = weighted)
        lr_scores(sme) <- out
        sme
    }
)
