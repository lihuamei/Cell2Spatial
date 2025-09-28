#' @title recoverSpatialExpr
#'
#' @description Aggregate cell-level expression (e.g., output from Cell2Spatial) back to the spot level based on `SpotName`, and reconstruct a Seurat spatial object with either SCT or standard RNA normalization.
#'
#' @param sce  A Seurat object produced by Cell2Spatial (cell-level resolution).
#' @param sp.obj A Seurat object with spatial data.
#' @param assay Character, either `"SCT"` or `"RNA"`.
#' @return A Seurat object.
#' @export recoverSpatialExpr

recoverSpatialExpr <- function(sce, sp.obj, assay = c("SCT", "RNA")) {
    sc.mat <- GetAssayData(sce, slot = "counts")
    spots <- unique(sce$SpotName)
    mat <- future.apply::future_lapply(spots, function(sp) {
        idxes <- which(sce$SpotName == sp)
        sc.mat[, idxes, drop = FALSE] %>% Matrix::rowSums(.)
    }, future.seed = TRUE) %>%
        do.call(cbind, .) %>%
        as.data.frame() %>%
        `colnames<-`(spots)
    obj <- CreateSeuratObject(count = mat, assay = "Spatial")
    if (assay == "SCT") {
        obj <- obj %>% SCTransform(assay = "Spatial", verbose = FALSE)
    } else {
        obj <- obj %>%
            NormalizeData(verbose = FALSE) %>%
            FindVariableFeatures(., nfeatures = 3000, verbose = FALSE) %>%
            ScaleData(., verbose = FALSE)
    }
    obj@images <- sp.obj@images
    return(obj)
}
