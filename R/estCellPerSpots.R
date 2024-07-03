#' @title estCellPerSpots

#' @description Estimate the number of cells per spot using regression method.
#' @param sp.obj Seurat object of spatial transcriptomics (ST) data.
#' @param max.cells Maximum number of cells in Spots. Default: 10 (10x Genomics).
#' @param fix.cells.in.spot Fixed number of cells assigned to each spot. Default: NULL.
#' @param quantile.cut Spots with a number of UMIs greater than the quantile cutoff are set as a reference. Default: 0.99.
#' @param num.genes Number of stable genes used in the estimation. Default: 2000.
#' @return A vector of estimated cells per spot.
#' @export estCellPerSpots

estCellPerSpots <- function(sp.obj, max.cells = 10, fix.cells.in.spot = NULL, quantile.cut = 0.999, num.genes = 2000) {
    options(warn = -1)
    if (!is.null(fix.cells.in.spot)) {
        pred.cnt <- rep(fix.cells.in.spot, ncol(sp.obj)) %>% `names<-`(colnames(sp.obj))
    } else {
        if (max.cells == 1) {
            pred.cnt <- rep(1, ncol(sp.obj)) %>% `names<-`(colnames(sp.obj))
        } else {
            count.mat <- GetAssayData(sp.obj, slot = "count", assay = "Spatial") %>% as.data.frame()
            count.mat <- count.mat[!grepl("^MT|^RP[SL]", rownames(count.mat) %>% toupper(.)), ]
            umi.count <- rowMeans(count.mat)
            count.mat <- count.mat[umi.count > quantile(umi.count, 0.1), ]
            rank.mat <- apply(count.mat, 2, rank)
            rank.sd <- apply(rank.mat, 1, sd)
            genes <- rank.sd[order(rank.sd)][1:num.genes] %>% names()
            count.mat <- count.mat[genes, ]
            cnt.sum <- colSums(count.mat)
            ref.cut <- quantile(cnt.sum, quantile.cut)
            ref.df <- rowMeans(count.mat[, which(cnt.sum >= ref.cut), drop = FALSE])
            pred.cnt <- apply(count.mat, 2, function(gg) ceiling(coef(MASS::rlm(gg ~ ref.df - 1)) * max.cells))
            pred.cnt[pred.cnt <= 0] <- 1
            pred.cnt[pred.cnt > max.cells] <- max.cells
        }
    }
    return(pred.cnt)
}

#' @title adjustScObj

#' @description Adjust single-cell (SC) data by incorporating inferred proportions and cell counts for mapping to spatial transcriptomics (ST) spots
#' @param sc.obj Seurat object of single-cell (SC) data.
#' @param st.prop A data.frame containing proportions for ST data.
#' @param num.cells A vector of cell counts for ST spots.
#' @return Adjusted Seurat object of SC data.

adjustScObj <- function(sc.obj, st.prop, num.cells) {
    prop.vec <- colSums(st.prop) / sum(colSums(st.prop))
    ctype.cnt <-
        {
            sum(num.cells) * prop.vec
        } %>% ceiling()
    sc.idents <- Idents(sc.obj)
    ctype.cnt[which.max(ctype.cnt)] <- ctype.cnt[which.max(ctype.cnt)] + sum(num.cells) - sum(ctype.cnt)

    sys.mat <- future.apply::future_lapply(names(ctype.cnt), function(xx) {
        cells.tmp <- sc.idents[sc.idents == xx]
        rep.bool <- ifelse(length(cells.tmp) >= ctype.cnt[xx], FALSE, TRUE)
        if (rep.bool) {
            count.mat <- subset(sc.obj, cells = names(cells.tmp)) %>%
                GetAssayData(., slot = "count") %>%
                as.data.frame()
            diff.cnt <- ctype.cnt[xx] - length(cells.tmp)
            rand.mat <- count.mat[, sample(names(cells.tmp), diff.cnt, replace = TRUE)] %>% as.matrix()

            noise <- matrix(rnorm(n = length(rand.mat), mean = 0, sd = 1), nrow = nrow(rand.mat), ncol = ncol(rand.mat))
            rand.mat <- round(rand.mat + noise)
            rand.mat[rand.mat < 0] <- 0
            colnames(rand.mat) <- paste0(sprintf("Pseudo_%s", xx), "_", 1:diff.cnt)
            mat <- cbind.data.frame(count.mat, rand.mat)
        } else {
            cells.sample <- sample(cells.tmp, ctype.cnt[xx], replace = rep.bool) %>% names()
            mat <- GetAssayData(sc.obj, slot = "count")[, cells.sample] %>% as.data.frame()
        }
        return(list(EXPR = mat, CT = rep(xx, ncol(mat))))
    }, future.seed = TRUE)

    expr.mat <- lapply(sys.mat, function(xx) xx$EXPR) %>% do.call(cbind, .)
    cell.types <- lapply(sys.mat, function(xx) xx$CT) %>%
        unlist(.) %>%
        as.vector()
    meta.data <-
        {
            sc.obj@meta.data[colnames(expr.mat), ]
        } %>% `rownames<-`(colnames(expr.mat))
    sc.obj.new <- CreateSeuratObject(count = expr.mat, meta.data = meta.data, assay = "RNA") %>% SCTransform(verbose = FALSE, method = "glmGamPoi")
    Idents(sc.obj.new) <- cell.types
    return(sc.obj.new)
}
