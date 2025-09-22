recoverSpatialExpr <- function(sce, sp.obj) {
    sc.mat <- GetAssayData(sce, slot = "counts")
    spots <- unique(sce$SpotName)
    mat <- future.apply::future_lapply(spots, function(sp) {
        idxes <- which(sce$SpotName == sp)
        sc.mat[, idxes, drop = FALSE] %>% rowSums(.)
    }, future.seed = TRUE) %>%
        do.call(cbind, .) %>%
        as.data.frame() %>%
        `colnames<-`(spots)
    obj <- CreateSeuratObject(count = mat, assay = "Spatial") %>% SCTransform(verbose = FALSE, assay = "Spatial")
    obj@images <- sp.obj@images
    return(obj)
}
