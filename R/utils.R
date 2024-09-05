#' @title println
#'
#' @description Print out message if verbose equals TRUE.
#' @param infos Message that need to be printed.
#' @param status Normal running messages, warn, or error. Default: INFO (normal status).
#' @param verbose Bool value, print out message or not, default: TRUE.
#' @return NULL.

println <- function(X, verbose = TRUE, status = c("INFO", "ERROR", "WARN"), ...) {
    status <- match.arg(status)
    infos <- do.call(sprintf, c(list(X, ...))) %>% paste0("[", status, "] ", .)
    if (verbose || status == "ERROR") {
        cat(paste0(infos, "\n"))
        if (status == "ERROR") stop()
    }
}

#' @title multipleProcess
#'
#' @description Open multiple workers for Seurat processing.
#' @param n.workers Number of cores. Default: 4.
#' @return NULL

multipleProcess <- function(n.workers = 4) {
    options(future.globals.maxSize = 200000 * 1024^2, future.seed = TRUE)
    future::plan("multicore", workers = n.workers)
}

#' @title garbageCollection
#'
#' @return NULL

garbageCollection <- function(...) {
    objects_to_remove <- as.character(match.call(expand.dots = FALSE)$...)
    rm(list = objects_to_remove, envir = parent.frame())
    gc()
}

#' @title downSamplSeurat
#'
#' @description Down-sampling for Seurat obeject.
#' @param obj.seu Seurat object
#' @param cnt Sampling size for each ident. Default: 500.
#' @param percent Sampling percentage of each ident (value range: 0-1). Default: NULL.
#' @param seed Set random seed, default: 123456.
#' @return Subset of seurat object.
#' @export downSamplSeurat

downSamplSeurat <- function(obj.seu, cnt = 500, percent = NULL, seed = 123456) {
    set.seed(seed)
    cells <- Idents(obj.seu) %>% table()
    sub.cells <- sapply(names(cells), function(xx) {
        sub.cells <- Idents(obj.seu)[Idents(obj.seu) == xx] %>% names()
        cnt <- ifelse(is.null(percent), cnt, length(sub.cells) * percent)
        if (length(sub.cells) > cnt) sub.cells <- sample(sub.cells, cnt, replace = FALSE)
        return(sub.cells)
    }) %>% unlist(use.names = FALSE)
    subset(obj.seu, cells = sub.cells)
}

#' @title findScMarkers
#'
#' @description Find markers among cell types based on SC data.
#' @param sc.obj Seurat object of single-cell data.
#' @param group.size Marker size of each subset derived from SC data. Default: 100.
#' @param select.markers Strategy for inferring marker genes. Default: shannon.
#' @param verbose Show running messages or not. Default: TRUE.
#' @return A tibble of identified markers.
#' @export findScMarkers

findScMarkers <- function(sc.obj, group.size, select.markers = c("shannon", "wilcox"), verbose = TRUE) {
    select.markers <- match.arg(select.markers)
    switch(select.markers,
        wilcox = {
            sc.obj.sub <- downSamplSeurat(sc.obj, cnt = 500)
            sc.markers <- sc.obj.sub %>%
                PrepSCTFindMarkers(., verbose = verbose) %>%
                FindAllMarkers(., assay = "SCT", only.pos = TRUE, verbose = verbose) %>%
                .[!grepl("^MT-|^RP[L|S]", .$gene), ] %>%
                mutate(EScore = {
                    .$avg_log2FC * .$pct.1
                }) %>%
                group_by(cluster) %>%
                top_n(wt = EScore, n = group.size) %>%
                {
                    split(.$gene, .$cluster)
                } %>%
                .[lapply(., length) > 0]
        },
        shannon = {
            X <-
                {
                    AverageExpression(sc.obj, assay = "SCT")$SCT + 1e-2
                } %>% as.matrix()
            nrep.sns <- dim(X)[2]
            X.norm <- X / rowSums(X)
            X.mu <- rowMeans(X.norm)
            w.genes <- MatrixGenerics::rowMaxs(X) %>%
                {
                    (. - min(.)) / (max(.) - min(.))
                }
            shannon.dist <- 1 / nrep.sns * rowSums(X.norm / X.mu * log2(X.norm / X.mu))
            group.ctype <- apply(X, 1, which.max) %>%
                colnames(X)[.] %>%
                as.vector()
            score.df <- cbind.data.frame(cluster = group.ctype, Score = shannon.dist * w.genes, gene = rownames(X))
            sc.markers <- score.df %>%
                arrange(desc(Score)) %>%
                arrange(cluster) %>%
                .[!grepl("^MT-|^RP[L|S]", .$gene), ] %>%
                group_by(cluster) %>%
                top_n(wt = Score, n = group.size) %>%
                {
                    split(.$gene, .$cluster)
                } %>%
                .[lapply(., length) > 0]
        }
    )
    return(sc.markers)
}

#' @title coExistIndex
#'
#' @description Calculation of co-exists index.
#' @param sce Seurat or SingleCellExperiment object returned from runMap2Sp function.
#' @param min.cells Cell types with number of cells must be greater than k. Default: 0.
#' @return A data.frame of coexists J-index.
#' @export coExistIndex

coExistIndex <- function(sce, min.cells = 0) {
    if (class(sce) == "Seurat") {
        meta.data <- sce@meta.data
    } else {
        meta.data <- SingleCellExperiment::colData(sce)
    }
    tar.cells <- meta.data$CellType %>%
        table() %>%
        {
            . > min.cells
        } %>%
        names()
    df.mat <- table(meta.data$CellType, meta.data$centerSPOT) %>%
        t() %>%
        as.data.frame.matrix() %>%
        .[, tar.cells]
    df.mat[df.mat > 0] <- 1

    j.index <- lapply(colnames(df.mat), function(xx) {
        ref <- df.mat[, xx]
        idxes <- lapply(colnames(df.mat), function(yy) {
            qry <- df.mat[, yy]
            tmp <- ref + qry
            idx <- sum(tmp == 2) / min(sum(ref > 0) + 1, sum(qry > 0) + 1)
            return(idx)
        }) %>%
            unlist() %>%
            as.vector() %>%
            `names<-`(colnames(df.mat))
    }) %>%
        do.call(rbind, .) %>%
        `rownames<-`(colnames(df.mat))
    return(j.index)
}

#' @title createSpatialObject
#'
#' @description Create Spatial seurat object for ST data from other platforms.
#' @param counts UMI count matrix, rows are genes and columns are barcodes.
#' @param coord.df A data.frame of coordinates for spatial barcodes. The coordinates are placed in the first two columns.
#' @param class Class name of the Spatial image slot. Default: SlideSeq.
#' @return A Seurat object for ST data.

createSpatialObject <- function(counts, coord.df, coord.label = c("x", "y"), meta.data = NULL, class = "SlideSeq") {
    ovp.spots <- intersect(colnames(counts), rownames(coord.df))
    if (!is.null(meta.data)) ovp.spots <- intersect(ovp.spots, rownames(meta.data))
    if (length(ovp.spots) == 0) println("No shared barcodes between count matrix and coordinates", verbose = TRUE, status = "ERROR")
    obj.seu <- CreateSeuratObject(counts = counts[, ovp.spots], assay = "Spatial")
    obj.seu@images$image <- new(
        Class = class,
        assay = "Spatial",
        key = "image_",
        coordinates = coord.df[ovp.spots, coord.label]
    )
    if (!is.null(meta.data)) {
        obj.seu@meta.data <- cbind.data.frame(obj.seu@meta.data, meta.data[ovp.spots, ])
    }
    return(obj.seu)
}
