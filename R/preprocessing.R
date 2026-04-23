#' @importFrom hdf5r h5file
#' @importFrom jsonlite read_json
#' @importFrom nanoparquet read_parquet
#' @importFrom utils read.csv
#' @importFrom Matrix sparseMatrix
#' @importFrom methods slot
#import description end
0

## author: Atul Deshpande
## email: adeshpande@jhu.edu

#
#===================
#' @title Load 10X Visium Expression Data
#' @description This loads log-transformed 10X Visium expression data from
#'  standard 10X 
#' Visium folder.
#' @export
#' @param visiumDir  A string path to the h5 file with expression information.
#' @param h5filename  A string of the name of the h5 file in the directory.
#' @return A matrix of class dgeMatrix or Matrix that contains the expression 
#' info for each sample (cells) across multiple features (genes)
#' @examples
#' library(SpaceMarkers)
#' #Visium data links
#' urls <- read.csv(system.file("extdata","visium_data.txt",
#' package = "SpaceMarkers",mustWork = TRUE))
#' counts_url <- urls[["visium_url"]][1]
#' #Remove present Directories if any
#' files <- list.files(".")[grepl(basename(counts_url),list.files("."))]
#' unlink(files)
#' download.file(counts_url,basename(counts_url), mode = "wb")
#' counts_matrix<-load10XExpr(visiumDir=".",h5filename = basename(counts_url))
#' files <- list.files(".")[grepl(basename(counts_url),list.files("."))]
#' unlink(files)
#' 

load10XExpr<- function(visiumDir=NULL,
                       h5filename='filtered_feature_bc_matrix.h5'){
  
    h5FilePath <- dir(path = visiumDir,pattern = h5filename,full.names = TRUE)
    hf <- hdf5r::h5file(filename = h5FilePath, mode='r')
    mat <- names(hf)
    counts <- hf[[paste0(mat, '/data')]]
    indices <- hf[[paste0(mat, '/indices')]]
    indptr <- hf[[paste0(mat, '/indptr')]]
    shp <- hf[[paste0(mat, '/shape')]]
    features <- hf[[paste0(mat, '/features/name')]][]
    barcodes <- hf[[paste0(mat, '/barcodes')]][]
    spMat <- Matrix::sparseMatrix(
        i = indices[] + 1,
        p = indptr[],
        x = as.numeric(x = counts[]),
        dims = shp[],
        repr = "C"
    )
    spMat@x <- log2(1+spMat@x)
    features <- make.unique(names = features)
    rownames(spMat) <- features
    colnames(spMat) <- barcodes
    hf$close_all()
    return(spMat)
}

#===================
#' @title Load 10x Visium Spatial Coordinates
#' @description This function loads spatial coordinates for each cell from a
#' 10X Visium 
#' spatial folder.
#' @export
#' @param visiumDir A string path to the location of the folder containing the 
#' spatial coordinates. The folder in your visiumDir must be named 'spatial' 
#' and must contain files 'scalefactors_json.json' 
#' and 'tissue_positions_list.csv.'
#' @param resolution A string specifying which values to look for in the .json 
#' object. Can be either fullres (default), lowres or hires.
#' @param version A string specifying the version of the spaceranger data.
#' @return a data frame of the spatial coordinates 
#' ( x and y) for each spot/cell
#' @examples
#' library(SpaceMarkers)
#' #Visium data links
#' urls <- read.csv(system.file("extdata","visium_data.txt",
#' package = "SpaceMarkers",mustWork = TRUE))
#' sp_url <- urls[["visium_url"]][2]
#' # Spatial Coordinates
#' download.file(sp_url, basename(sp_url), mode = "wb")
#' untar(basename(sp_url))
#' spCoords <- load10XCoords(visiumDir = ".", version = "1.0")
#' unlink("spatial", recursive = TRUE)
#' unlink("Visium_Human_Breast_Cancer_spatial.tar.gz")
#' 

load10XCoords <- function(visiumDir, resolution = c("fullres","lowres","hires"), version = NULL){

    #resolve resolution parameter
    resolution <- match.arg(resolution, several.ok = FALSE)
    message("resolution: ", resolution)
    #determine spacerager version
    if(is.null(version)){
        message("Version not provided. Trying to infer.")
        if("probe_set.csv" %in% dir(visiumDir)){
            config_line <- readLines(paste0(visiumDir,"/probe_set.csv"), 1)
            version <- strsplit(config_line, "=")[[1]][2]
        } else if ("tissue_positions.parquet" %in% 
        dir(file.path(visiumDir,"spatial"))) {
          version <- "HD"
          visiumDir <- file.path(visiumDir)
          message(".parquet file found.
                  Detected VisiumHD.")
        } else {
            message(
              "probe_set.csv or .parquet not found.Assuming version 1.0.")
            version <- "1.0"
        }
    }
    spatial_dir <- paste0(visiumDir,"/spatial")
    #account for different versions of visium data
    if(version == "1.0"){
      has_header <- FALSE
      tissue_pos_name <- "tissue_positions_list.csv"
      coord_reader <- read.csv
    } else if (version == "2.0") {
      has_header <- TRUE
      tissue_pos_name <- "tissue_positions.csv"
      coord_reader <- read.csv
    } else if (version == "HD") {
      has_header <- TRUE
      tissue_pos_name <- "tissue_positions.parquet"
      coord_reader <- nanoparquet::read_parquet
    }
    coord_file <- dir(spatial_dir,pattern = tissue_pos_name, full.names = TRUE)
    if(!is.null(formals(coord_reader)["header"])) 
      formals(coord_reader)["header"] <- has_header
    coord_values <- coord_reader(coord_file)
    scale_json <- dir(spatial_dir,
                        pattern = "scalefactors_json.json",full.names = TRUE)
    scale_values <- jsonlite::read_json(scale_json)
    if (resolution == "fullres"){
        scale_factor <- 1
    } else {
        scale_factor <- scale_values[grepl(resolution, names(scale_values))][[1]]
    }
    
    coord_values <- coord_values[,c(1,5,6)]
    coord_values[,2:3] <- coord_values[,2:3]*scale_factor
    names(coord_values) <- c("barcode","y","x")
    return(coord_values)
}

#===================
#' @title Load spatial features
#' @description This function loads spatial features from a file containing
#' spatial features
#' @export
#' @param filePath A string path to the location of the file containing the 
#' spatial features. 
#' @param method A string specifying the type of object to obtain spatial
#' feature from. Default NULL, where the method is inferred based on object
#' type. Other methods are: "CoGAPS", "Seurat", or "BayesTME".
#' @param featureNames An array of strings specifying the column names 
#' corresponding to the feature names or a regex string. In the case of Seurat,
#' all metadata columns with "_Feature" suffix are selected.
#' @return a matrix of spatial features with barcodes associated 
#' with individual coordinates
#' @examples
#' library(SpaceMarkers)
#' #CoGAPS data filePath
#' filePath <- system.file("extdata","CoGAPS_result.rds", 
#' package = "SpaceMarkers",mustWork = TRUE)
#' spFeatures <- get_spatial_features(filePath, method = "CoGAPS")
#' head(spFeatures)
#' 

get_spatial_features <- function(filePath, method = NULL, featureNames = "."){

    # If already an R object (e.g., SpatialExperiment), use directly
    if (is.object(filePath) && !is.character(filePath)) {
        spObject <- filePath
    } else {
        #read the features object based on the format
        spObject <- .read_format(filePath)
    }

    #determine the method to use for feature extractioin
    method <- .infer_method(spObject, method)

    spFun <- c("CoGAPS"=.get_cogaps_features,
                "BayesTME"=.get_BTME_features,
                "Seurat"=.get_seurat_features,
                "CSV"=.get_csv_features,
                "SpatialExperiment"=.get_spe_features)

    spFeatures <- spFun[[method]](spObject)

    dataNames <- colnames(spFeatures)

    #subset the features based on the featureNames
    if(length(featureNames) == 1) {
        #assume regex is provided
        namePattern <- featureNames
        featureNames <- dataNames[grepl(pattern = namePattern,
                                        dataNames, ignore.case = TRUE)]
        if(length(featureNames) == 0) {
           stop(sprintf("Regex %s does not match any feature.", namePattern))
        }
    } else if(!all(featureNames %in% dataNames)) {
    stop("Some of the features were not found:",
            sprintf(" %s", setdiff(featureNames, dataNames)))
    }

    featureNames <- intersect(featureNames, dataNames)
    spFeatures <- spFeatures[,featureNames, drop = FALSE]

    return(spFeatures)
}

#' readFormat
#' Reads a format into an R object
#' @keywords internal
#' 
.read_format <- function(path){
    if(grepl(".rds",path)){
        obj <- readRDS(path)
    } else if (grepl(".h5ad",path)){
        obj <- hdf5r::h5file(filename = path, mode='r')
    } else if (grepl(".csv",path)){
        obj <- read.csv(path)
    } else
    {
        stop("File format not supported.")
    }
    return(obj)
}

#' inferMethod
#' Infer the method used to obtain spatial features
#' @keywords internal
.infer_method <- function(spObject, method){
    if(is.null(method)){
        if(inherits(spObject, "SpatialExperiment")){
            method <- "SpatialExperiment"
        } else if(inherits(spObject, "H5File")){
            method <- "BayesTME"
        } else if(inherits(spObject, "CogapsResult")){
            method <- "CoGAPS"
        } else if(inherits(spObject, "Seurat")){
            method <- "Seurat"
        } else if(inherits(spObject, "data.frame")){
            method <- "CSV"
        } else {
            stop("Method not supported.")
        }
    }
    return(method)
}

#' .get_cogaps_features
#' Load features CoGAPS object
#' @keywords internal
#' 
.get_cogaps_features <- function(obj){
    spFeatures <- slot(obj, "sampleFactors")
    return(spFeatures)
}

#' .get_BTME_features
#' Load features BayesTME object
#' 
#' @keywords internal
#' 
.get_BTME_features <- function(hf){
    feat_loc <- "obsm/bayestme_cell_type_counts"
    barc_loc <- "obs/_index"
    spFeatures <- t(hdf5r::readDataSet(hf[[feat_loc]]))
    barcodes <- hdf5r::readDataSet(hf[[barc_loc]])
    rownames(spFeatures) <- barcodes
    if (is.null(colnames(spFeatures))) {
        colnames(spFeatures)<-paste0("BayesTME_",seq(1,ncol(spFeatures)))
    }
    
    return(spFeatures)
}

#' .get_seurat_features
#' Load features Seurat object
#' @keywords internal
#' 
.get_seurat_features <- function(obj){
    spFeatures <- slot(obj, "meta.data")
    selection <- grepl("_Feature",colnames(spFeatures), ignore.case = TRUE)
    if (!any(selection)){
        stop("No _Feature columns found in Seurat object.")
    }
    spFeatures <- spFeatures[,selection, drop = FALSE]
    return(spFeatures)
}

#' .get_csv_features
#' Load features from dataframe
#' @keywords internal
#'
.get_csv_features <- function(obj){
    spFeatures <- obj
    if ("barcode" %in% colnames(spFeatures)){
        rownames(spFeatures) <- spFeatures$barcode
    } else {
        if(!colnames(spFeatures)[1]=="X"){
            stop("No barcode column found and first colname is not blank.
                    Stopping.")
        } else {
            rownames(spFeatures) <- spFeatures[,"X"]
        }
    }
    removeCols <- c("NA","barcode","in_tissue","array_row","array_col","pxl_col_in_fullres","pxl_row_in_fullres")
    spFeatures <- spFeatures[,-which(startsWith(colnames(spFeatures),"X") | colnames(spFeatures) %in% removeCols)]
    return(spFeatures)
}

#' .get_spe_features
#' Extract spatial features from a SpatialExperiment object's colData
#' @keywords internal
.get_spe_features <- function(obj) {
    cd <- as.data.frame(SummarizedExperiment::colData(obj))
    # Remove standard SE/SPE columns
    remove_cols <- c("sample_id", "in_tissue", "array_row", "array_col",
                     "pxl_col_in_fullres", "pxl_row_in_fullres",
                     "sizeFactor")
    keep_cols <- setdiff(colnames(cd), remove_cols)
    # Keep only numeric columns as features
    is_numeric <- vapply(cd[, keep_cols, drop = FALSE], is.numeric,
                         logical(1))
    spFeatures <- cd[, keep_cols[is_numeric], drop = FALSE]
    if (ncol(spFeatures) == 0) {
        stop("No numeric feature columns found in SpatialExperiment colData.")
    }
    return(as.data.frame(spFeatures))
}

#===================
#' @title Load 10X Visium data as a SpaceMarkersExperiment
#' @description Convenience function that loads expression, coordinates, and
#'   optionally spatial features from a 10X Visium directory and assembles them
#'   into a \code{\link{SpaceMarkersExperiment}} object.
#' @param visiumDir A string path to the 10X Visium directory.
#' @param features Optional: a file path or object for spatial features, passed
#'   to \code{\link{get_spatial_features}}.
#' @param h5filename Name of the H5 expression file. Default
#'   "filtered_feature_bc_matrix.h5".
#' @param resolution Resolution for coordinates. One of "fullres", "lowres",
#'   "hires".
#' @param version Optional Spaceranger version.
#' @param ... Additional arguments passed to \code{get_spatial_features}.
#' @return A \code{\link{SpaceMarkersExperiment}} object.
#' @export
load10X <- function(visiumDir,
                         features = NULL,
                         h5filename = "filtered_feature_bc_matrix.h5",
                         resolution = c("fullres", "lowres", "hires"),
                         version = NULL,
                         ...) {
    resolution <- match.arg(resolution)
    expr <- load10XExpr(visiumDir = visiumDir, h5filename = h5filename)
    coords <- load10XCoords(visiumDir = visiumDir,
                            resolution = resolution, version = version)
    rownames(coords) <- coords$barcode

    all_expr <- colnames(expr)
    all_coords <- coords$barcode
    common <- intersect(all_expr, all_coords)
    dropped <- length(all_expr) + length(all_coords) - 2L * length(common)
    if (dropped > 0L) {
        warning(sprintf(
            "Dropped %d spots not present in both expression and coordinates. Keeping %d common spots.",
            dropped, length(common)
        ))
    }
    expr <- expr[, common, drop = FALSE]
    coords <- coords[common, , drop = FALSE]

    coord_mat <- as.matrix(coords[, c("y", "x")])
    rownames(coord_mat) <- common

    cd <- S4Vectors::DataFrame(row.names = common)
    pattern_names <- NULL

    if (!is.null(features)) {
        spFeatures <- get_spatial_features(features, ...)
        shared <- intersect(common, rownames(spFeatures))
        dropped_feat <- length(rownames(spFeatures)) - length(shared)
        dropped_sme <- length(common) - length(shared)
        if (dropped_feat > 0L || dropped_sme > 0L) {
            warning(sprintf(
                paste0("Spot mismatch: %d spots in expression/coords and ",
                       "%d spots in features not shared. ",
                       "Keeping %d common spots for features (NA-padded)."),
                dropped_sme, dropped_feat, length(shared)
            ))
        }
        feat_df <- S4Vectors::DataFrame(
            spFeatures[shared, , drop = FALSE],
            row.names = shared
        )
        # Pad unmatched spots with NA
        for (col in colnames(feat_df)) {
            cd[[col]] <- NA_real_
            cd[shared, col] <- feat_df[shared, col]
        }
        pattern_names <- colnames(spFeatures)
    }

    # Pre-compute spatial parameters from the scalefactors JSON if available
    spatial_par <- NULL
    if (!is.null(pattern_names) && length(shared) > 0L) {
        spCoords <- coords[shared, c("barcode", "y", "x"), drop = FALSE]
        spPats <- as.data.frame(spFeatures[shared, , drop = FALSE])
        spPatterns_combined <- cbind(spCoords, spPats)
        tryCatch({
            spatial_par <- get_spatial_parameters(
                spatialPatterns = spPatterns_combined,
                visiumDir = visiumDir,
                resolution = resolution
            )
        }, error = function(e) {
            warning(sprintf(
                "Pre-computing spatial_params failed: %s", conditionMessage(e)
            ))
            NULL
        })
    }

    sm <- as(list(
        params = list(
            pattern_names = pattern_names,
            spatial_params = spatial_par,
            visiumDir = visiumDir,
            resolution = resolution
        )
    ), "SimpleList")

    SpaceMarkersExperiment(
        assays = list(logcounts = expr),
        colData = cd,
        spatialCoords = coord_mat,
        spaceMarkers = sm
    )
}

#===================
#' @title Load an AnnData file as a SpaceMarkersExperiment
#' @description Convenience wrapper that reads an `.h5ad` file and
#'   coerces the resulting \code{SingleCellExperiment} directly into a
#'   \code{\link{SpaceMarkersExperiment}}. If the AnnData stores spatial
#'   coordinates under \code{obsm["spatial"]}, they are promoted to
#'   \code{spatialCoords()}.
#'
#'   Two readers are supported:
#'   \itemize{
#'     \item \pkg{anndataR} — pure-R, no Python/conda dependency, faster
#'       on first call. Recommended when available.
#'     \item \pkg{zellkonverter} — Bioconductor-standard, uses a
#'       \pkg{basilisk}-managed Python environment under the hood.
#'   }
#'   Both are \code{Suggests} dependencies. With \code{reader = "auto"}
#'   (default), anndataR is preferred when installed; otherwise
#'   zellkonverter is used. Pass \code{reader = "zellkonverter"} or
#'   \code{reader = "anndataR"} to force a specific backend.
#' @param file Path to an \code{.h5ad} file.
#' @param reader One of \code{"auto"}, \code{"anndataR"},
#'   \code{"zellkonverter"}.
#' @param ... Additional arguments forwarded to the chosen reader's
#'   read function.
#' @return A \code{\link{SpaceMarkersExperiment}} object.
#' @export
load_anndata <- function(file,
                         reader = c("auto", "anndataR", "zellkonverter"),
                         ...) {
    reader <- match.arg(reader)
    if (reader == "auto") {
        reader <- if (requireNamespace("anndataR", quietly = TRUE)) {
            "anndataR"
        } else if (requireNamespace("zellkonverter", quietly = TRUE)) {
            "zellkonverter"
        } else {
            stop(
                "load_anndata() requires either 'anndataR' (pure R) or ",
                "'zellkonverter' (Bioconductor). Install one of: ",
                "install.packages('anndataR')  or  ",
                "BiocManager::install('zellkonverter')."
            )
        }
    } else if (!requireNamespace(reader, quietly = TRUE)) {
        stop(sprintf("load_anndata(reader = '%s') requires the '%s' package.",
                     reader, reader))
    }
    sce <- switch(reader,
        anndataR      = anndataR::as_SingleCellExperiment(
                            anndataR::read_h5ad(file, ...)),
        zellkonverter = zellkonverter::readH5AD(file, ...)
    )
    sme <- methods::as(sce, "SpaceMarkersExperiment")
    .unpack_spacemarkers_state(sme)
}

#===================
#' @title Save a SpaceMarkersExperiment as an AnnData (.h5ad) file
#' @description Writes a \code{\link{SpaceMarkersExperiment}} out to
#'   an \code{.h5ad} file, preserving the full SpaceMarkers analysis
#'   state (hotspots, overlap scores, interactions, influence map,
#'   kernel parameters) under \code{uns["spacemarkers"]}. Reloading
#'   with \code{load_anndata()} restores that state.
#'
#'   The reader selection logic is identical to
#'   \code{\link{load_anndata}}: \pkg{anndataR} is preferred when
#'   installed, \pkg{zellkonverter} is the fallback.
#' @param sme A \code{\link{SpaceMarkersExperiment}} object.
#' @param file Path to the \code{.h5ad} file to write.
#' @param reader One of \code{"auto"}, \code{"anndataR"},
#'   \code{"zellkonverter"}.
#' @param ... Additional arguments forwarded to the chosen reader's
#'   write function.
#' @return The path \code{file}, returned invisibly.
#' @export
save_anndata <- function(sme, file,
                         reader = c("auto", "anndataR", "zellkonverter"),
                         ...) {
    reader <- match.arg(reader)
    if (reader == "auto") {
        reader <- if (requireNamespace("anndataR", quietly = TRUE)) {
            "anndataR"
        } else if (requireNamespace("zellkonverter", quietly = TRUE)) {
            "zellkonverter"
        } else {
            stop(
                "save_anndata() requires either 'anndataR' (pure R) or ",
                "'zellkonverter' (Bioconductor). Install one of: ",
                "install.packages('anndataR')  or  ",
                "BiocManager::install('zellkonverter')."
            )
        }
    } else if (!requireNamespace(reader, quietly = TRUE)) {
        stop(sprintf("save_anndata(reader = '%s') requires the '%s' package.",
                     reader, reader))
    }
    state <- .pack_spacemarkers_state(sme)
    # Coerce SME down to SCE, but carry the packed state in metadata so it
    # round-trips as uns["spacemarkers"] when the reader writes it out.
    sce <- methods::as(sme, "SingleCellExperiment")
    md <- S4Vectors::metadata(sce)
    md$spacemarkers <- state
    S4Vectors::metadata(sce) <- md
    switch(reader,
        anndataR      = anndataR::write_h5ad(
                            anndataR::as_AnnData(sce), file, ...),
        zellkonverter = zellkonverter::writeH5AD(sce, file, ...)
    )
    invisible(file)
}

#===================
# Pack / unpack SpaceMarkers analysis state as a plain nested list so it
# survives AnnData serialization (uns accepts lists of atomic / matrix /
# data.frame values but not arbitrary S4 objects).
.pack_spacemarkers_state <- function(sme) {
    sm <- sme@spacemarkers
    md <- S4Vectors::metadata(sme)
    as_plain_df <- function(x) {
        if (is.null(x)) NULL
        else if (methods::is(x, "DataFrame")) as.data.frame(x)
        else x
    }
    list(
        params = list(
            pattern_names  = sm$params$pattern_names,
            spatial_params = sm$params$spatial_params,
            visiumDir      = sm$params$visiumDir,
            resolution     = sm$params$resolution
        ),
        results = list(
            undirected_scores = as_plain_df(sm$results$undirected_scores),
            directed_scores   = as_plain_df(sm$results$directed_scores),
            overlap_scores    = as_plain_df(sm$results$overlap_scores),
            lr_scores         = as_plain_df(sm$results$lr_scores)
        ),
        analysis       = sm$analysis,
        hotspots       = lapply(md$hotspots, as_plain_df),
        interactions   = md$interactions,
        influence_map  = as_plain_df(md$influence)
    )
}

.unpack_spacemarkers_state <- function(sme) {
    md <- S4Vectors::metadata(sme)
    state <- md$spacemarkers
    if (is.null(state)) return(sme)
    sm <- sme@spacemarkers
    if (is.null(sm)) sm <- S4Vectors::SimpleList()
    if (!is.null(state$params))  sm$params  <- state$params
    if (!is.null(state$results)) sm$results <- state$results
    if (!is.null(state$analysis)) sm$analysis <- state$analysis
    sme@spacemarkers <- sm
    if (!is.null(state$hotspots))      md$hotspots      <- state$hotspots
    if (!is.null(state$interactions))  md$interactions  <- state$interactions
    if (!is.null(state$influence_map)) md$influence     <- state$influence_map
    md$spacemarkers <- NULL  # scrub the sidecar so we don't re-pack on next save
    S4Vectors::metadata(sme) <- md
    sme
}
