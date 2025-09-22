#' @title println
#'
#' @description Print out message if verbose equals TRUE.
#'
#' @param infos Message that need to be printed.
#' @param status Normal running messages, warn, or error. Default: INFO (normal status).
#' @param verbose Bool value, print out message or not. Default: TRUE.
#' @return NULL.

println <- function(X, ..., verbose = TRUE, status = c("INFO", "WARN", "ERROR")) {
    status <- match.arg(status)
    msg <- sprintf(X, ...)
    msg <- paste0("[", status, "] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " - ", msg)

    if (verbose || status == "ERROR") {
        cat(msg, "\n")
        if (status == "ERROR") {
            stop(msg, call. = FALSE)
        }
    }
    invisible(NULL)
}

#' @title multipleProcess
#'
#' @description Open multiple workers for Seurat processing.
#'
#' @param n.workers Number of cores. Default: 4.
#' @param max.size Maximum allowed size (in bytes) of global variables.
#' @param use.multicore Logical, use multicore backend if available. Default: TRUE.
#' @return NULL

multipleProcess <- function(n.workers = 4, max.size = 128 * 1024^3, use.multicore = TRUE) {
    if (is.null(n.workers)) {
        n.workers <- parallel::detectCores()
    }
    options(future.globals.maxSize = max.size, future.seed = TRUE)
    if (use.multicore && .Platform$OS.type == "unix") {
        future::plan("multicore", workers = n.workers)
    } else {
        future::plan("multisession", workers = n.workers)
    }
    invisible(NULL)
}

#' @title garbageCollection
#'
#' @return NULL

garbageCollection <- function(...) {
    objects.to.remove <- as.character(match.call(expand.dots = FALSE)$...)
    rm(list = objects.to.remove, envir = parent.frame())
    gc()
}

#' @title downSamplSeurat
#'
#' @description Down-sampling for Seurat obeject.
#'
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

#' @title binarizeMax
#'
#' @description Converts a matrix of scores into a binary matrix where each row contains a 1 at the position of the maximum value in that row.
#'
#' @param sp.score A matrix of scores (e.g., spatial transcriptomics scores). Each row represents a different spot or entity, and each column represents a different gene or feature.
#' @return A binary matrix of the same size as the input matrix, where each row contains a 1 at the position of the maximum value in that row, and 0 elsewhere.

binarizeMax <- function(sp.score) {
    sp.score <- as.matrix(sp.score)
    sp.binary <- matrix(0, nrow = nrow(sp.score), ncol = ncol(sp.score))
    rownames(sp.binary) <- rownames(sp.score)
    colnames(sp.binary) <- colnames(sp.score)

    max_idx <- max.col(sp.score, ties.method = "first")
    sp.binary[cbind(1:nrow(sp.score), max_idx)] <- 1

    return(sp.binary)
}

#' @title prop2counts
#'
#' @description Converts proportions to integer cell counts based on a given number of cells for each spot.
#'
#' @param st.prop A matrix of proportions, where each row represents a different spot (e.g., tissue region or location), and each column represents a different cell type or category.
#' @param num.cells A vector of length equal to the number of rows in `st.prop`, representing the total number of cells available for each corresponding spot.
#' @return A matrix of integer cell counts corresponding to the proportions in `st.prop`, ensuring the total number of cells for each spot matches `num.cells`.

prop2counts <- function(st.prop, num.cells) {
    st.prop <- as.matrix(st.prop)
    stopifnot(nrow(st.prop) == length(num.cells))
    raw <- sweep(st.prop, 1, num.cells, "*")
    counts <- floor(raw)
    remainder <- raw - counts
    need_add <- num.cells - rowSums(counts)

    for (i in seq_len(nrow(st.prop))) {
        if (need_add[i] > 0) {
            ord <- order(remainder[i, ], decreasing = TRUE)[1:need_add[i]]
            counts[i, ord] <- counts[i, ord] + 1
        }
    }
    rownames(counts) <- rownames(st.prop)
    colnames(counts) <- colnames(st.prop)
    return(counts)
}

#' @title roundCountsUseProp
#'
#' @description Round proportional cell type counts so that the total matches the observed number of cells.
#'
#' @param st.prop A matrix or data.frame of cell type proportions per spot. Columns correspond to cell types.
#' @param num.cells A numeric vector of observed cell counts per spot, used to determine the target total number of cells.
#' @return A named integer vector of rounded cell type counts, whose sum equals the sum of `num.cells`.

roundCountsUseProp <- function(st.prop, num.cells) {
    target.total <- sum(num.cells)
    counts.raw <- colSums(st.prop)
    counts.prop <- counts.raw / sum(counts.raw) * target.total
    counts.floor <- floor(counts.prop)

    remainder <- target.total - sum(counts.floor)
    if (remainder > 0) {
        order.frac <- order(counts.prop - counts.floor, decreasing = TRUE)
        counts.floor[order.frac[seq_len(remainder)]] <- counts.floor[order.frac[seq_len(remainder)]] + 1
    }
    return(counts.floor)
}

#' @title createSpatialObject
#'
#' @description Create Spatial seurat object for ST data from other platforms.
#'
#' @param counts UMI count matrix, rows are genes and columns are barcodes.
#' @param coord.df A data.frame of coordinates for spatial barcodes. The coordinates are placed in the first two columns. Default: c("x", "y").
#' @param meta.data A data frame containing metadata to be added to the Seurat object. Default: NULL.
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
        obj.seu@meta.data <- cbind.data.frame(obj.seu@meta.data, meta.data[ovp.spots, , drop = FALSE])
    }
    return(obj.seu)
}
