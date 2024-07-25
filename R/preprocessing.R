#' @importFrom hdf5r h5file
#' @importFrom jsonlite read_json
#' @importFrom utils read.csv
#' @importFrom Matrix sparseMatrix
#' @importFrom methods slot
#import description end
0

## author: Atul Deshpande
## email: adeshpande@jhu.edu

#
#===================
#' load10XExpr
#' load 10X Visium Expression Data
#'
#' This loads log-transformed 10X Visium expression data from standard 10X 
#' Visium folder.
#'
#' @export
#'
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
        repr = "T"
    )
    spMat <- log2(1+spMat)
    features <- make.unique(names = features)
    rownames(spMat) <- features
    colnames(spMat) <- barcodes
    hf$close_all()
    return(spMat)
}

#===================
#' load10XCoords
#' Load 10x Visium Spatial Coordinates
#'
#' This function loads spatial coordinates for each cell from a 10X Visium 
#' spatial folder.
#'
#' @export
#'
#' @param visiumDir A string path to the location of the folder containing the 
#' spatial coordinates. The folder in your visiumDir must be named 'spatial' 
#' and must contain files 'scalefactors_json.json' 
#' and 'tissue_positions_list.csv.'
#' @param resolution A string specifying which values to look for in the .json 
#' object. Can be either lowres or highres.
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

load10XCoords <- function(visiumDir, resolution = "lowres", version = NULL){
    #determine spacerager version
    if(is.null(version)){
        message("Version not provided. Trying to infer.")
        if("probe_set.csv" %in% dir(visiumDir)){
            config_line <- readLines(paste0(visiumDir,"/probe_set.csv"), 1)
            version <- strsplit(config_line, "=")[[1]][2]
        } else {
            message("probe_set.csv not found. Assuming version 1.0.")
            version <- "1.0"
        }
    }
    #account for different versions of visium data
    if(version == "1.0"){
        has_header <- FALSE
        tissue_pos_name <- "tissue_positions_list.csv"
    } else if (version == "2.0") {
        has_header <- TRUE
        tissue_pos_name <- "tissue_positions.csv"
  }
    spatial_dir <- paste0(visiumDir,'/spatial')
    scale_json <- dir(spatial_dir,
                        pattern = "scalefactors_json.json",full.names = TRUE)
    scale_values <- jsonlite::read_json(scale_json)
    scale_factor <- scale_values[grepl(resolution, names(scale_values))][[1]]
    coord_file <- dir(spatial_dir,
                        pattern = tissue_pos_name, full.names = TRUE)

    coord_values <- read.csv(coord_file, header = has_header)
    coord_values <- coord_values[,c(1,5,6)]
    coord_values[,2:3] <- coord_values[,2:3]*scale_factor
    names(coord_values) <- c("barcode","y","x")
    return(coord_values)
}

#===================
#' getSpatialFeatures
#' Load spatial features
#'
#' This function loads spatial features from a file containing spatial features
#'
#' @export
#'
#' @param filePath A string path to the location of the file containing the 
#' spatial features. 
#' @param method A string specifying the method used to obtain spatial 
#' features. e.g., "CoGAPS", "Seurat", or "BayesTME".
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
#' spFeatures <- getSpatialFeatures(filePath, method = "CoGAPS")
#' head(spFeatures)
#' 

getSpatialFeatures <- function(filePath,method = "CoGAPS",featureNames = "."){

    if(method == "CoGAPS") {
        spFeatures <- getCogaps(filePath)
    } else if(method == "BayesTME") {
        spFeatures <- getBtme(filePath)
    } else if(method == "Seurat") {
        spFeatures <- getSeurat(filePath)
    } else {
        stop("Method not supported.")
    }

    if (method == "Seurat") {
        namePattern <- "_Feature"
    } else {
        namePattern <- featureNames
    }

    if (length(featureNames) == 1) {
        dataNames <- colnames(spFeatures)
        featureNames <- dataNames[grepl(pattern = namePattern,
                                        dataNames, ignore.case = TRUE)]
        message(sprintf("%d features match expression '%s'.",
                        length(featureNames), namePattern))
    }

    featureNames <- intersect(featureNames,colnames(spFeatures))

    if(length(featureNames)>0) {
        spFeatures <- spFeatures[,featureNames,drop = FALSE]
    } else {stop("No features match provided feature names. Check your input.")}

    return(spFeatures)
}

#' getCogaps
#' Load features CoGAPS object
#' @keywords internal
#' 

getCogaps <- function(path){
    spFeatures <- readRDS(path)
    spFeatures <- slot(spFeatures,"sampleFactors")
    return(spFeatures)
}

#' getBtme
#' Load features BayesTME object
#' 
#' @keywords internal
#' 

getBtme <- function(path){
    hf <- hdf5r::h5file(filename = path, mode='r')
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

#' getSeurat
#' Load features Seurat object
#' @keywords internal
#' 

getSeurat <- function(path){
    spFeatures <- readRDS(path)
    spFeatures <- spFeatures[[]]
    return(spFeatures)
}