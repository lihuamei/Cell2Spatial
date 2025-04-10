#' @title inferCellNumbers
#' @description Estimates the number of cells per spot in a spatial transcriptomics Seurat object using a corrected saturation model.
#'
#' @param sp.obj A Seurat object containing spatial transcriptomics data (Q matrix in the "counts" slot of the "Spatial" assay).
#' @param max.cells Maximum number of cells assumed in any spot (default = 10).
#' @param umi.thresh Minimum UMI count to consider a gene as expressed (default = 2).
#' @param w.L Weight for library size in the composite score (default = 0.5).
#' @param w.T Weight for gene count in the composite score (default = 0.5).
#' @param fix.cells.in.spot Fixed number of cells assigned to each spot. Default: NULL.
#' @return A named vector of estimated cell numbers per spot.
#' @export inferCellNumbers

inferCellNumbers <- function(sp.obj, max.cells = 10, umi.thresh = 2, w.L = 0.5, w.T = 0.5, fix.cells.in.spot = NULL) {
    if (!inherits(sp.obj, "Seurat")) {
        stop("Input 'sp.obj' must be a Seurat object.")
    }
    if (!("Spatial" %in% names(sp.obj@assays))) {
        stop("Seurat object must contain a 'Spatial' assay.")
    }
    if (!is.null(fix.cells.in.spot)) {
        pred.cnt <- rep(fix.cells.in.spot, ncol(sp.obj)) %>% `names<-`(colnames(sp.obj))
    } else {
        if (max.cells == 1) {
            pred.cnt <- rep(1, ncol(sp.obj)) %>% `names<-`(colnames(sp.obj))
        } else {
            if (w.L < 0 || w.T < 0 || (w.L + w.T) == 0) {
                stop("Weights 'w.L' and 'w.T' must be non-negative and sum to a positive value.")
            }
            Q <- GetAssayData(sp.obj, slot = "counts", assay = "Spatial") %>% as.data.frame()
            expr.genes <- Q >= umi.thresh
            T.j <- colSums(expr.genes)
            L.j <- colSums(Q)
            T.max <- max(T.j)

            # Compute composite score S_j
            L.max <- max(L.j)
            S_j <- w.L * (L.j / L.max) + w.T * (T.j / T.max)

            # Identify reference spot (j_ref) with highest composite score
            j_ref <- which.max(S_j)
            L.j_ref <- L.j[j_ref]
            T.j_ref <- T.j[j_ref]

            m <- L.j_ref / max.cells
            t <- T.j_ref / max.cells

            C.j <- numeric(length(T.j))
            C.j_lib <- numeric(length(T.j))
            for (j in seq_along(T.j)) {
                if (T.j[j] == 0) {
                    C.j[j] <- 1 # Avoid log(0); assume 1 cell for empty spots
                } else {
                    C.j[j] <- log(1 - T.j[j] / T.max) / log(1 - t / T.max)
                }
                C.j_lib[j] <- L.j[j] / m # Library size-based estimate
            }
            pred.cnt <- 2 / (1 / pmax(C.j, 1) + 1 / pmax(C.j_lib, 1))
            pred.cnt <- pmax(1, pmin(max.cells, round(pred.cnt)))
            names(pred.cnt) <- colnames(Q)
        }
    }
    return(pred.cnt)
}

#' @title adjustScObj
#'
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
            cells.rand <- sample(names(cells.tmp), diff.cnt, replace = TRUE)
            rand.mat <- count.mat[, cells.rand] %>% as.matrix()

            noise <- matrix(rnorm(n = length(rand.mat), mean = 0, sd = 1), nrow = nrow(rand.mat), ncol = ncol(rand.mat))
            rand.mat <- round(rand.mat + noise)
            rand.mat[rand.mat < 0] <- 0
            colnames(rand.mat) <- paste0(sprintf("Pseudo_%s", cells.rand), "_XYZ", 1:diff.cnt)
            mat <- cbind.data.frame(count.mat, rand.mat)
        } else {
            cells.sample <- sample(cells.tmp, ctype.cnt[xx], replace = rep.bool) %>% names()
            mat <- GetAssayData(sc.obj, slot = "count")[, cells.sample, drop = FALSE] %>% as.data.frame()
        }
        return(list(EXPR = mat, CT = rep(xx, ncol(mat))))
    }, future.seed = TRUE)

    expr.mat <- lapply(sys.mat, function(xx) xx$EXPR) %>% do.call(cbind, .)
    cell.types <- lapply(sys.mat, function(xx) xx$CT) %>%
        unlist(.) %>%
        as.vector()
    meta.data <-
        {
            sc.obj@meta.data[gsub("Pseudo_|_XYZ.*", "", colnames(expr.mat)), ]
        } %>% `rownames<-`(colnames(expr.mat))
    sc.obj.new <- CreateSeuratObject(count = expr.mat, meta.data = meta.data, assay = "RNA") %>% SCTransform(verbose = FALSE, method = "glmGamPoi")
    Idents(sc.obj.new) <- cell.types
    return(sc.obj.new)
}
