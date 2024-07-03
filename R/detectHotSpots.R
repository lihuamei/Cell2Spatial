#' @title getGsetScore

#' @description Get signature scores in each cell or spot using a list of markers.
#' @param obj.seu Seurat object.
#' @param sc.markers A list of markers for cell types.
#' @param assay The assay used to calculate signature scores of cell types. Default: SCT.
#' @return A data.frame of signature scores.
#' @export getGsetScore

getGsetScore <- function(obj.seu, sc.markers, assay = "SCT") {
    expr <- GetAssayData(obj.seu, slot = "data", assay = assay)
    score.df <- future.apply::future_lapply(names(sc.markers), function(ctype) {
        gset <- sc.markers[[ctype]]
        expr.sub <- expr[gset, ] %>% as.data.frame()
        score <- colMeans(expr.sub, na.rm = TRUE)
    }, future.seed = TRUE) %>%
        do.call(rbind, .) %>%
        t() %>%
        `colnames<-`(names(sc.markers))
    return(score.df)
}

#' @title dectHotSpotsByTtest

#' @description Infer whether the signature scores in spots significant higher than background.
#' @param sp.obj Seurat object of spatial transcriptome (ST) data.
#' @param sp.score Signature score of cell types.
#' @param knn K nearest cells used. Default: 5.
#' @param p.value.threshold P-value threadhold for filtering out hot-spots. Default: 0.05.
#' @return A data frame containing hotspot information. A value of TRUE indicates a detected hotspot, while FALSE indicates a non-hotspot.

dectHotSpotsByTtest <- function(sp.obj, sp.score, knn = 5, p.value.threshold = 0.05) {
    dist.lst <- calcSpotsDist(sp.obj)
    avgs.vec <- colMeans(sp.score)
    sd.vec <- apply(sp.score, 2, var)
    nctypes <- ncol(sp.score)
    n1 <- knn + 1
    n2 <- nrow(sp.score)
    sd.vec.new <- (sd.vec^2 / (n2^2 * (n2 - 1)))
    fix.value.n1 <- n1^2 * (n1 - 1)

    spot.pvals <- future.apply::future_lapply(rownames(sp.score), function(spot) {
        spot.ids <- dist.lst[[Idents(sp.obj)[spot]]]
        knn.spots <- spot.ids[spot, ] %>%
            rownames(spot.ids)[.] %>%
            c(spot, .)

        pval.ctype <- rep(1, nctypes) %>% `names<-`(colnames(sp.score))
        vec1.mean <- colMeans(sp.score[knn.spots, ])
        bool.vec <- vec1.mean - avgs.vec
        idxes <- which(bool.vec > 0)
        if (length(idxes) == 0) {
            return(pval.ctype)
        }
        sig.spots <- lapply(names(idxes), function(ctype) {
            var1 <- var(sp.score[knn.spots, ctype])
            var2 <- sd.vec[ctype]
            df <- ((var1 / n1 + var2 / n2)^2) / ((var1^2 / fix.value.n1) + sd.vec.new[ctype])
            t.statistic <- bool.vec[ctype] / sqrt(var1 / n1 + var2 / n2)
            p.value <- pt(t.statistic, df, lower.tail = FALSE)
        }) %>%
            unlist() %>%
            as.vector() %>%
            `names<-`(names(idxes))
        pval.ctype[names(sig.spots)] <- sig.spots
        pval.ctype
    }, future.seed = TRUE) %>%
        do.call(rbind, .) %>%
        `rownames<-`(rownames(sp.score))
    hotspot.spots <- spot.pvals <= p.value.threshold
    return(spot.pvals)
}

#' @title dectHotSpotsByGetisOrdGi

#' @description Detect hot spots using Getis-Ord Gi* index
#' @param sp.obj Seurat object of spatial transcriptome (ST) data.
#' @param sp.score Signature score of cell types.
#' @param knn K nearest cells used. Default: 5.
#' @param p.value.threshold P-value threadhold for filtering out hot-spots. Default: 0.05.
#' @return A data frame containing hotspot information. A value of TRUE indicates a detected hotspot, while any other value indicates a non-hotspot.

dectHotSpotsByGetisOrdGi <- function(sp.obj, sp.score, knn = 5, p.value.threshold = 0.05) {
    image.coord <- GetTissueCoordinates(sp.obj)[colnames(sp.obj), ]
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
    return(hotspot.spots)
}
