#' @import hdf5r
#' @import jsonlite
#' @import dplyr
#import description end
0
 
## author: Atul Deshpande
## email: adeshpande@jhu.edu
## Need Documentation

# 
#===================
#' load10XExpr
#' load log-transformed 10X Visium expression data from standard 10X Visium folder
#'
#' This function loads ...
#'
#' @export
#'
#' @param visiumDir 	...
#' @param h5filename 	...
#'
#'
load10XExpr<- function(visiumDir = NULL,h5filename= 'filtered_feature_bc_matrix.h5'){
    h5FilePath <- dir(path = visiumDir,pattern = h5filename,full.names = T)
    
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
#' Load ...
#'
#' This function loads ...
#'
#' @export
#'
#' @param visiumDir 	...
#' @param resolution 	...
#'
#'
load10XCoords <- function(visiumDir, resolution = "lowres"){
    scale_json <- dir(paste0(visiumDir,'/spatial'),pattern = "scalefactors_json.json",full.names = T)
    scale_values <- jsonlite::read_json(scale_json)
    scale_dia <- scale_values$spot_diameter_fullres
    scale_factor <- scale_values[[contains(vars = names(scale_values),match=resolution)]]
    coord_file <- dir(paste0(visiumDir,'/spatial'),pattern = "tissue_positions_list.csv",full.names = T)
    coord_values <- read.csv(coord_file,header = F)
    coord_values <- coord_values[,c(1,5,6)]
    coord_values[,2:3] <- coord_values[,2:3]*scale_factor
    names(coord_values) <- c("barcode","y","x")
    return(coord_values)
}
  
