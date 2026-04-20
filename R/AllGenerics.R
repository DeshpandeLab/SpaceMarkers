#' @include AllClasses.R

#' @title Accessor generics for SpaceMarkersExperiment
#' @name SpaceMarkersExperiment-accessors
#' @description Generics for accessing SpaceMarkersExperiment results.

#' @rdname spatial_patterns
#' @export
setGeneric("spatial_patterns", function(x, ...) {
    standardGeneric("spatial_patterns")
})

#' @rdname spatial_patterns
#' @export
setGeneric("spatial_patterns<-", function(x, value) {
    standardGeneric("spatial_patterns<-")
})

#' @rdname spatial_params
#' @export
setGeneric("spatial_params", function(x, ...) {
    standardGeneric("spatial_params")
})

#' @rdname spatial_params
#' @export
setGeneric("spatial_params<-", function(x, value) {
    standardGeneric("spatial_params<-")
})

#' @rdname params
#' @export
setGeneric("params", function(x, ...) {
    standardGeneric("params")
})

#' @rdname hotspots
#' @export
setGeneric("hotspots", function(x, ...) {
    standardGeneric("hotspots")
})

#' @rdname influence_map
#' @export
setGeneric("influence_map", function(x, ...) {
    standardGeneric("influence_map")
})

#' @rdname interactions
#' @export
setGeneric("interactions", function(x, ...) {
    standardGeneric("interactions")
})

#' @rdname undirected_scores
#' @export
setGeneric("undirected_scores", function(x, ...) {
    standardGeneric("undirected_scores")
})

#' @rdname directed_scores
#' @export
setGeneric("directed_scores", function(x, ...) {
    standardGeneric("directed_scores")
})

#' @rdname lr_scores
#' @export
setGeneric("lr_scores", function(x, ...) {
    standardGeneric("lr_scores")
})

#' @rdname overlap_scores
#' @export
setGeneric("overlap_scores", function(x, ...) {
    standardGeneric("overlap_scores")
})

#' @rdname analysis_type
#' @export
setGeneric("analysis_type", function(x, ...) {
    standardGeneric("analysis_type")
})

#' @rdname hotspots
#' @export
setGeneric("hotspots<-", function(x, type = "undirected", value) {
    standardGeneric("hotspots<-")
})

#' @rdname interactions
#' @export
setGeneric("interactions<-", function(x, value) {
    standardGeneric("interactions<-")
})

#' @rdname influence_map
#' @export
setGeneric("influence_map<-", function(x, value) {
    standardGeneric("influence_map<-")
})

#' @rdname undirected_scores
#' @export
setGeneric("undirected_scores<-", function(x, value) {
    standardGeneric("undirected_scores<-")
})

#' @rdname directed_scores
#' @export
setGeneric("directed_scores<-", function(x, value) {
    standardGeneric("directed_scores<-")
})

#' @rdname lr_scores
#' @export
setGeneric("lr_scores<-", function(x, value) {
    standardGeneric("lr_scores<-")
})

#' @rdname overlap_scores
#' @export
setGeneric("overlap_scores<-", function(x, value) {
    standardGeneric("overlap_scores<-")
})

#' @rdname analysis_type
#' @export
setGeneric("analysis_type<-", function(x, value) {
    standardGeneric("analysis_type<-")
})

#' @rdname find_all_hotspots
#' @export
setGeneric("find_all_hotspots", function(spPatterns, params = NULL,
                                         outlier = "positive",
                                         nullSamples = 1000,
                                         includeSelf = TRUE, ...) {
    standardGeneric("find_all_hotspots")
})

#' @rdname get_pairwise_interacting_genes
#' @export
setGeneric("get_pairwise_interacting_genes",
    function(data, spPatterns = NULL, mode = "residual", optParams = NULL,
             reconstruction = NULL, hotspots = NULL, minOverlap = 50,
             analysis = c("overlap", "enrichment"),
             pattern_pairs = NULL, ..., workers = 1) {
        standardGeneric("get_pairwise_interacting_genes")
    }
)

#' @rdname get_im_scores
#' @export
setGeneric("get_im_scores", function(SpaceMarkers) {
    standardGeneric("get_im_scores")
})

#' @rdname calculate_overlap_undirected
#' @export
setGeneric("calculate_overlap_undirected",
    function(hotspots, patternList = NULL,
             method = c("Szymkiewicz-Simpson", "Jaccard", "Sorensen-Dice",
                        "Ochiai", "absolute")) {
        standardGeneric("calculate_overlap_undirected")
    }
)

#' @rdname calculate_influence
#' @export
setGeneric("calculate_influence", function(spPatterns, optParams = NULL, ...) {
    standardGeneric("calculate_influence")
})

#' @rdname find_hotspots_gmm
#' @export
setGeneric("find_hotspots_gmm", function(df, threshold = 0.1, ...) {
    standardGeneric("find_hotspots_gmm")
})

#' @rdname calculate_gene_scores_directed
#' @export
setGeneric("calculate_gene_scores_directed",
    function(data, pat_hotspots = NULL, influence_hotspots = NULL,
             pattern_pairs = NULL, ...) {
        standardGeneric("calculate_gene_scores_directed")
    }
)

#' @rdname calculate_overlap_directed
#' @export
setGeneric("calculate_overlap_directed",
    function(pat_hotspots, influence_hotspots = NULL,
             patternList = NULL,
             method = c("relative-abundance", "differential-abundance", "absolute")) {
        standardGeneric("calculate_overlap_directed")
    }
)

#' @rdname calculate_gene_set_score
#' @export
setGeneric("calculate_gene_set_score",
    function(IMscores, gene_sets = NULL, weighted = TRUE,
             method = c("geometric_mean", "arithmetic_mean")) {
        standardGeneric("calculate_gene_set_score")
    }
)
