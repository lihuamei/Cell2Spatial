#' @title inferCellNumbers
#'
#' @description Estimates the number of cells per spot in a spatial transcriptomics Seurat object using a corrected saturation model.
#'
#' @param sp.obj A Seurat object containing spatial transcriptomics data (Q matrix in the "counts" slot of the "Spatial" assay).
#' @param max.cells Maximum number of cells assumed in any spot. Default: 10.
#' @param umi.thresh Minimum UMI count to consider a gene as expressed. Default: 1.
#' @param w.L Weight for library size in the composite score. Default: 0.5.
#' @param w.T Weight for gene count in the composite score: Default: 0.5.
#' @param fix.cells.in.spot Fixed number of cells assigned to each spot or not. Default: FALSE.
#' @param q.ref Quantile threshold used to select spots for reference. Default: 0.90.
#' @param winsor.p Quantile for winsorization. Default: 0.01.
#' @return A named vector of estimated cell numbers per spot.
#' @export inferCellNumbers

inferCellNumbers <- function(sp.obj, max.cells = 10, umi.thresh = 1, w.L = 0.5, w.T = 0.5, fix.cells.in.spot = FALSE, q.ref = 0.90, winsor.p = 0.01) {
    if (!inherits(sp.obj, "Seurat")) stop("Input 'sp.obj' must be a Seurat object.")
    if (!("Spatial" %in% names(sp.obj@assays))) stop("Seurat object must contain a 'Spatial' assay.")

    if (fix.cells.in.spot) {
        pred.cnt <- rep(max.cells.in.spot, ncol(sp.obj)) %>% `names<-`(colnames(sp.obj))
        return(pred.cnt)
    }
    if (max.cells == 1) {
        pred.cnt <- rep(1, ncol(sp.obj)) %>% names() <- (colnames(sp.obj))
        return(pred.cnt)
    }
    if (w.L < 0 || w.T < 0 || (w.L + w.T) == 0) {
        stop("Weights 'w.L' and 'w.T' must be non-negative and sum to a positive value.")
    }
    Q <- GetAssayData(sp.obj, slot = "counts", assay = "Spatial")
    expr.genes <- Q >= umi.thresh
    T.j <- Matrix::colSums(expr.genes)
    L.j <- Matrix::colSums(Q)

    winsor <- function(x, p = 0.01) {
        lo <- stats::quantile(x, p, names = FALSE, type = 7)
        hi <- stats::quantile(x, 1 - p, names = FALSE, type = 7)
        pmin(pmax(x, lo), hi)
    }
    if (!is.null(winsor.p) && winsor.p > 0) {
        Lw <- winsor(L.j, winsor.p)
        Tw <- winsor(T.j, winsor.p)
    } else {
        Lw <- L.j
        Tw <- T.j
    }
    T.max <- max(Tw)
    L.max <- max(Lw)

    S_j <- w.L * (Lw / L.max) + w.T * (Tw / T.max)
    j_top <- which(S_j >= stats::quantile(S_j, q.ref, names = FALSE))
    L.j_ref <- round(mean(Lw[j_top]))
    T.j_ref <- round(mean(Tw[j_top]))

    m <- max(L.j_ref / max.cells, .Machine$double.eps)
    t <- max(T.j_ref / max.cells, .Machine$double.eps)

    eps <- 1e-8
    denom <- 1 - t / max(T.max, 1)
    denom <- pmin(pmax(denom, eps), 1 - eps)
    frac <- 1 - Tw / max(T.max, 1)
    frac <- pmin(pmax(frac, eps), 1 - eps)

    C.j <- ifelse(Tw <= 0, 1, log(frac) / log(denom))
    C.j_lib <- as.numeric(Lw / max(m, .Machine$double.eps))

    C1 <- pmin(max.cells, pmax(C.j, 1e-6))
    C2 <- pmin(max.cells, pmax(C.j_lib, 1e-6))

    w_sum <- w.L + w.T
    wL <- w.L / w_sum
    wT <- w.T / w_sum

    pred.cnt <- 1 / (wT / C1 + wL / C2)
    pred.cnt <- pmax(1, pmin(max.cells, round(pred.cnt)))
    names(pred.cnt) <- colnames(Q)
    return(pred.cnt)
}

#' @title adjustScObj
#'
#' @description Adjust single-cell (SC) data by incorporating inferred proportions and cell counts for mapping to spatial transcriptomics (ST) spots
#'
#' @param sc.obj Seurat object of single-cell (SC) data.
#' @param st.prop A data.frame containing cell counts for ST data.
#' @param assay Assay type to use for SC data normalization. Default: SCT.
#' @param nfeatures Number of features used for SC and ST integration analysis. Default: 3000.
#' @return Adjusted Seurat object of SC data.

adjustScObj <- function(sc.obj, st.prop.cnts, assay = "SCT", nfeatures = 3000) {
    ctype.cnt <- colSums(st.prop.cnts)[colSums(st.prop.cnts) > 0]
    sc.idents <- Idents(sc.obj)

    counts.all <- GetAssayData(sc.obj, slot = "counts") # 默认 assay
    sys.list <- future.apply::future_lapply(names(ctype.cnt), function(xx) {
        cells.of.type <- names(sc.idents)[sc.idents == xx]
        need <- ctype.cnt[[xx]]
        have <- length(cells.of.type)
        if (have == 0) {
            return(list(
                EXPR = Matrix::Matrix(0, nrow = nrow(counts.all), ncol = 0, sparse = TRUE),
                CT = character(0)
            ))
        }
        if (have >= need) {
            pick <- sample(cells.of.type, need, replace = FALSE)
            mat <- counts.all[, pick, drop = FALSE]
        } else {
            mat.base <- counts.all[, cells.of.type, drop = FALSE]
            add.n <- need - have
            pick.add <- sample(cells.of.type, add.n, replace = TRUE)
            mat.add <- counts.all[, pick.add, drop = FALSE]
            if (add.n > 0) {
                noise <- Matrix::Matrix(rpois(n = length(mat.add@x), lambda = 0.1), sparse = TRUE)
                mat.add@x <- pmax(0, mat.add@x + noise@x)
                colnames(mat.add) <- paste0("Pseudo_", pick.add, "_XYZ", seq_len(add.n))
            }
            mat <- Matrix::cbind2(mat.base, mat.add)
        }
        list(EXPR = mat, CT = rep(xx, ncol(mat)))
    }, future.seed = TRUE)

    expr.list <- lapply(sys.list, `[[`, "EXPR")
    expr.mat <- if (length(expr.list) == 1) expr.list[[1]] else Reduce(Matrix::cbind2, expr.list)
    cell.types <- unlist(lapply(sys.list, `[[`, "CT"), use.names = FALSE)

    original.ids <- sub("_XYZ\\d+$", "", sub("^Pseudo_", "", colnames(expr.mat)))
    md <- sc.obj@meta.data[original.ids, , drop = FALSE]
    rownames(md) <- colnames(expr.mat)

    sc.new <- CreateSeuratObject(counts = expr.mat, meta.data = md, assay = "RNA")
    if (assay == "SCT") {
        sc.new <- SCTransform(sc.new, verbose = FALSE, method = "glmGamPoi")
    } else {
        sc.new <- sc.new %>%
            NormalizeData(verbose = FALSE) %>%
            FindVariableFeatures(nfeatures = nfeatures, verbose = FALSE) %>%
            ScaleData(verbose = FALSE) %>%
            RunPCA(verbose = FALSE)
    }
    DefaultAssay(sc.new) <- assay
    Idents(sc.new) <- factor(cell.types, levels = names(ctype.cnt))
    return(sc.new)
}

.adjustScObj <- function(...) {
    suppressWarnings({
        sc.new <- adjustScObj(...)
    })
    return(sc.new)
}
