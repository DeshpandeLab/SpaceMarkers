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
#' download.file(counts_url,basename(counts_url))
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
#' @return a data frame of the spatial coordinates 
#' ( x and y) for each spot/cell
#' @examples
#' library(SpaceMarkers)
#' #Visium data links
#' urls <- read.csv(system.file("extdata","visium_data.txt",
#' package = "SpaceMarkers",mustWork = TRUE))
#' sp_url <- urls[["visium_url"]][2]
#' # Spatial Coordinates
#' download.file(sp_url, basename(sp_url))
#' untar(basename(sp_url))
#' spCoords <- load10XCoords(visiumDir = ".")
#' unlink("spatial", recursive = TRUE)
#' unlink("Visium_Human_Breast_Cancer_spatial.tar.gz")
#' 

load10XCoords <- function(visiumDir, resolution = "lowres"){
    scale_json <- dir(paste0(visiumDir,'/spatial'),
                        pattern = "scalefactors_json.json",full.names = TRUE)
    scale_values <- jsonlite::read_json(scale_json)
    scale_dia <- scale_values$spot_diameter_fullres
    scale_factor <- scale_values[grepl(resolution, names(scale_values))][[1]]
    coord_file <- dir(paste0(visiumDir,'/spatial'),
                        pattern="tissue_positions_list.csv",full.names = TRUE)
    coord_values <- read.csv(coord_file,header = FALSE)
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
#' corresponding to the feature names. If input is NULL, in the case of CoGAPS 
#' and BayesTME, all features are selected In the case of Seurat, all metadata 
#' columns with "_Feature" suffix are selected.
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

getSpatialFeatures <- function(filePath,method = "CoGAPS",featureNames = NULL){
    if(method=="CoGAPS"){
        spFeatures <- readRDS(filePath)
        spFeatures <- slot(spFeatures,"sampleFactors")
    } else if(method=="BayesTME"){
        hf <- hdf5r::h5file(filename = filePath, mode='r')
        spFeatures <- t(hdf5r::readDataSet(
            hf[["obsm/bayestme_cell_type_counts"]]))
        barcodes <- hdf5r::readDataSet(hf[["obs/_index"]])
        rownames(spFeatures) <- barcodes
        if (is.null(colnames(spFeatures)))
            colnames(spFeatures)<-paste0("BayesTME_",seq(1,ncol(spFeatures)))
    } else if(method=="Seurat"){
        spFeatures <- readRDS(filePath)
        spFeatures <- spFeatures[[]]
    } else {stop("Method not supported.")}
    if(is.null(featureNames)){
        featureNames <- colnames(spFeatures)
        message("No feature names provided. Using all available features.")
        if(method=="Seurat") {
            featureNames <- grepl(colnames(spFeatures),pattern = "_feature", 
                                    ignore.case = TRUE)
            message("Using all metadata columns with '_Feature' suffix.")
        }
    } else{
        if (length(featureNames) == 1){
            featureNames <- colnames(spFeatures)[grepl(pattern=featureNames,
                                                        colnames(spFeatures),
                                                        ignore.case = TRUE)]
            message("Only one featureName provided.
                    Assuming input is regular expression.")
            if (length(featureNames) == 0)
                stop("No features found with matching regular expression.
                        Please check your input.")
            else
                message("Found ",length(featureNames),
                        " features matching the regular expression.")
        }
        else
            featureNames <- intersect(featureNames,colnames(spFeatures))
        if(!is.null(featureNames))
            spFeatures <- spFeatures[,featureNames]
        else
            stop("No features found in the spatial
                    data with provided feature names.")
    }
    spFeatures <- spFeatures[,featureNames]
    return(spFeatures)
}
