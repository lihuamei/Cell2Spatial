#' @title getGsetScore
#'
#' @description Get signature scores in each cell or spot using a list of markers.
#' @param obj.seu Seurat object.
#' @param sc.markers A list of markers for cell types.
#' @param assay The assay used to calculate signature scores of cell types. Default: SCT.
#' @param method Method for scoring the signature of cell types in ST data: AddModuleScore, UCell, or AverageExpr. Default: AddModuleScore.
#' @return A data.frame of signature scores.
#' @export getGsetScore

getGsetScore <- function(obj.seu, sc.markers, assay = "SCT", method = c("AddModuleScore", "UCell", "AverageExpr")) {
    score.df <- switch(match.arg(method),
        AddModuleScore = {
            obj.seu <- tryCatch(
                {
                    obj.seu <- Seurat::AddModuleScore(obj = obj.seu, assay = assay, features = sc.markers, name = "CELL2SPATIAL")
                },
                error = function(e) {
                    obj.seu <- Seurat::AddModuleScore(obj = obj.seu, assay = assay, features = sc.markers, name = "CELL2SPATIAL", nbin = 10)
                }
            )
            score.df <- obj.seu@meta.data[, grep("CELL2SPATIAL", colnames(obj.seu@meta.data)), drop = FALSE]
            colnames(score.df) <- names(sc.markers)[1:ncol(score.df)]
            return(score.df)
        },
        UCell = {
            score.df <- UCell::ScoreSignatures_UCell(GetAssayData(obj.seu), features = sc.markers) %>% `colnames<-`(names(sc.markers))
        },
        AverageExpr = {
            expr <- GetAssayData(obj.seu, slot = "data", assay = assay)
            score.df <- future.apply::future_lapply(names(sc.markers), function(ctype) {
                gset <- sc.markers[[ctype]]
                expr.sub <- expr[gset, ] %>% as.data.frame()
                score <- colMeans(expr.sub, na.rm = TRUE)
            }, future.seed = TRUE) %>%
                do.call(rbind, .) %>%
                t() %>%
                `colnames<-`(names(sc.markers))
        }
    )
    return(score.df)
}

#' @title dectHotSpotsByGetisOrdGi
#'
#' @description Detect hot spots using Getis-Ord Gi* index
#' @param sp.obj Seurat object of spatial transcriptome (ST) data.
#' @param sp.score Signature score of cell types.
#' @param knn K nearest cells used. Default: 5.
#' @param p.value.threshold P-value threadhold for filtering out hot-spots. Default: 0.05.
#' @return A list of data frames containing hotspot information and P-values. A value of TRUE indicates a detected hotspot, while any other value indicates a non-hotspot.

dectHotSpotsByGetisOrdGi <- function(sp.obj, sp.score, knn = 5, p.value.threshold = 0.05) {
    image.coord <- GetTissueCoordinates(sp.obj)[colnames(sp.obj), c(1, 2)]
    nb <- spdep::knn2nb(spdep::knearneigh(image.coord, k = knn))

    listw <- spdep::nb2listw(nb, style = "W")
    gi.values <- apply(sp.score, 2, function(score) {
        spdep::localG(score, listw)
    })
    rownames(gi.values) <- rownames(sp.score)
    p.values <- apply(gi.values, 2, function(gi) {
        pnorm(gi, lower.tail = FALSE)
    })
    hotspot.spots <- p.values <= p.value.threshold
    return(list(x = hotspot.spots, y = p.values))
}
