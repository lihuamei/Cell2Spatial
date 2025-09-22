#' @title weightSimScore
#'
#' @description Applies weighted adjustment to similarity scores between spots and cell types based on provided adjustment weights and hotspot information.
#'
#' @param out.sim Matrix or data frame of similarity scores between spots and cell types.
#' @param adj.w Weight matrix specifying adjustments for each spot and cell type combination.
#' @param spot.name Vector or list of spot names.
#' @param cell.names Vector or list of cell type names.
#' @param hot.pvals A data frame of hotspot p-values indicating hotspot presence in spots and cell types.
#' @param Numeric p-value cutoff for hotspot significance. Default: 0.01.
#' @return Weighted similarity scores matrix where adjustments are applied based on provided weights and hotspot information.

weightSimScore <- function(out.sim, adj.w, spot.name, cell.names, hot.pvals = NULL, p.cut = 0.05) {
    spot.name <- as.vector(spot.name) %>% `names<-`(names(spot.name))
    cell.names <- as.vector(cell.names) %>% `names<-`(names(cell.names))
    out.sim.sub <- out.sim[names(spot.name), names(cell.names)]
    if (ncol(adj.w) != ncol(hot.pvals)) {
        adj.w.sub <- adj.w[names(spot.name), names(cell.names)] %>% `colnames<-`(names(cell.names))
    } else {
        adj.w.sub <- adj.w[spot.name, cell.names] %>%
            `colnames<-`(names(cell.names)) %>%
            `rownames<-`(names(spot.name))
    }
    out.sim.w <-
        {
            out.sim.sub * adj.w.sub
        } %>% t()
    if (!is.null(hot.pvals)) {
        hot.bool <- hot.pvals <= p.cut
        hot.spts.flat <- hot.bool[names(spot.name), cell.names] %>%
            `colnames<-`(names(cell.names)) %>%
            t()
        idxes <- which(hot.spts.flat == FALSE)
        out.sim.w[idxes] <- out.sim.w[idxes] - 1
    }
    return(out.sim.w)
}

#' @title selectByProb
#'
#' @description Computes a similarity score matrix based on spatial and single-cell signature scores, adjusting for probabilities.
#'
#' @param sp.score Matrix or data frame of ST signature scores.
#' @param sc.score Matrix or data frame of SC signature scores.
#' @return Similarity matrix where each row represents a spatial spot and each column represents a single-cell, adjusted based on computed probabilities.

selectByProb <- function(sp.score, sc.score) {
    probs <-
        {
            log2(sc.score + 1e-6) - log2(rowSums(sc.score) + ncol(sc.score) * 1e-6)
        } %>% t()
    testa <- log1p(as.matrix(sp.score + 1e-6)) / log(2)
    out.sim <- {
        as.matrix(testa) %*% as.matrix(probs)
    }
    max.vec <- apply(out.sim, 1, max)
    min.vec <- apply(out.sim, 1, min)
    out.sim <- sweep(out.sim, 1, min.vec, "-") %>% sweep(., 1, max.vec - min.vec, "/")
    return(out.sim)
}

#' @title calcSimilarityDist
#'
#' @description Builds a spot–cell similarity matrix.
#'
#' @param sc.obj Seurat object containing SC data (assay slot must hold a counts/data matrix).
#' @param sp.obj Seurat object containing ST data.
#' @param sc.score Numeric matrix/data frame of SC scores (**cells × features**).
#' @param sp.score Numeric matrix/data frame of ST scores (**spots × features**).
#' @param feature.based Specify whether features for likelihood or correlation calculations between single cells and spots are based on gene expression ('gene.based') or signature scores of cell types ('celltype.based').
#' @param sc.markers List of marker genes per cell type.
#' @return A numeric matrix of shape **spots × cells** with spot–cell similarities.

calcSimlarityDist <- function(sc.obj, sp.obj, sc.score, sp.score, feature.based, sc.markers) {
    if (feature.based == "gene.based") {
        uniq.markers <- intersect(unique(unlist(sc.markers)), rownames(sc.obj))
        sp.score <- GetAssayData(sp.obj) %>%
            .[uniq.markers, ] %>%
            as.matrix() %>%
            t()
        sc.score <- GetAssayData(sc.obj) %>%
            .[uniq.markers, ] %>%
            as.matrix() %>%
            t()
    }
    out.sim <- selectByProb(sp.score, sc.score)
    return(out.sim)
}

#' @title featureSelection
#'
#' @description Performs feature selection and data integration between spatial transcriptomics (ST) and single-cell (SC) data.
#'
#' @param sp.obj Seurat object containing ST data.
#' @param sc.obj Seurat object containing SC data.
#' @param n.features Number of features to select for integration. Default: 3000.
#' @param assay Assay to use for SC data. Default: "SCT".
#' @param reduction Dimensionality reduction technique to align spatial and single-cell data. Default: cca.
#' @param verbose Boolean indicating whether to display verbose messages during processing. Default: TRUE.
#' @return A list containing integrated Seurat objects for SC (`sc.int`) and ST (`st.int`) data.

featureSelection <- function(sp.obj, sc.obj, n.features = 3000, assay = "SCT", reduction = "cca", verbose = TRUE) {
    ifnb.list <- list(SC = sc.obj, ST = sp.obj)
    features <- SelectIntegrationFeatures(object.list = ifnb.list, nfeatures = n.features, verbose = verbose)
    sc.st.anchors <- Seurat::FindTransferAnchors(
        reference = sc.obj,
        query = sp.obj,
        reference.assay = assay,
        query.assay = ifelse(assay == "RNA", "Spatial", "SCT"),
        features = features,
        reduction = reduction,
        verbose = verbose
    )
    st.data.trans <- Seurat::TransferData(
        anchorset = sc.st.anchors,
        refdata = GetAssayData(sc.obj, assay = assay, slot = "data")[features, ],
        weight.reduction = ifelse(reduction == "cca", "cca", "pcaproject"),
        verbose = verbose
    )
    sp.obj@assays$transfer <- st.data.trans

    if (sc.obj@version >= "5.0") {
        st.vv <- data.frame(sp.obj[["transfer"]]$data) %>% `colnames<-`(colnames(sp.obj))
        sc.vv <- data.frame(sc.obj[[assay]]$data[features, ]) %>% `colnames<-`(colnames(sc.obj))
        counts.temp <- cbind(st.vv, sc.vv)
        sc.st.int <- CreateSeuratObject(counts = counts.temp, assay = "traint")
        sc.st.int[["traint"]]$data <- sc.st.int[["traint"]]$counts
    } else {
        st.vv <- data.frame(sp.obj[["transfer"]]@data) %>% `colnames<-`(colnames(sp.obj))
        sc.vv <- data.frame(sc.obj[[assay]]@data[features, ]) %>% `colnames<-`(colnames(sc.obj))
        counts.temp <- cbind(st.vv, sc.vv)
        sc.st.int <- CreateSeuratObject(counts = counts.temp, assay = "traint")
        sc.st.int[["traint"]]@data <- sc.st.int[["traint"]]@counts
    }
    sc.st.int <- ScaleData(sc.st.int, features = features, verbose = verbose) %>% RunPCA(features = features, verbose = verbose)
    sc.int <- subset(sc.st.int, cells = colnames(sc.obj))
    st.int <- subset(sc.st.int, cells = colnames(sp.obj))
    Idents(st.int) <- Idents(sp.obj)
    return(list(sc = sc.int, st = st.int))
}

#' @title getClusterCellCounts
#'
#' @description Sum estimated cell-type counts per ST cluster (from spot-level counts).
#'
#' @param st.prop.cnts Numeric matrix/data frame of estimated counts with \strong{rows = spots}, \strong{cols = cell types}.
#' @param sp.obj Seurat object (ST) whose \code{Idents(sp.obj)} are spot clusters.
#' @return Numeric matrix \strong{clusters × cell types}.

getClusterCellCounts <- function(st.prop.cnts, sp.obj) {
    cnts.cls <- lapply(levels(sp.obj), function(cls) {
        spots <- colnames(sp.obj)[Idents(sp.obj) == cls]
        cnts.sum <-
            {
                st.prop.cnts[spots, , drop = FALSE]
            } %>%
            colSums()
        return(cnts.sum)
    }) %>%
        do.call(rbind, .) %>%
        `rownames<-`(levels(sp.obj)) %>%
        as.matrix()
    return(cnts.cls)
}

#' @title partitionClusters
#'
#' @description Partition single cells into sub-clusters using a deep learning strategy based on integration of spatial transcriptomics (ST) and single-cell (SC) data.
#' @param sp.obj Seurat object containing ST data.
#' @param sc.obj Seurat object containing SC data.
#' @param hot.pvals Data frame of hotspot p-values indicating hotspot presence in spots and cell types.
#' @param st.prop.cnts Estimated spot-level counts matrix (rows=spots, cols=cell types).
#' @param assay Assay to use for SC data. Default: "SCT".
#' @param reduction Dimensionality reduction technique to align spatial and single-cell data. Default: cca.
#' @param npcs Number of PCA components to use. Default: 30.
#' @return A list of single cells partitioned into sub-clusters.
#' @export partitionClusters

partitionClusters <- function(sp.obj, sc.obj, hot.pvals, st.prop.cnts, assay = "SCT", reduction = "cca", npcs = 30) {
    sc.st.int <- featureSelection(sp.obj, sc.obj, assay = assay, reduction = reduction, verbose = FALSE)
    sc.int <- sc.st.int$sc
    sp.int <- sc.st.int$st
    num.cells <- rowSums(st.prop.cnts)
    python.script <- system.file("python/netx.py", package = "Cell2Spatial")
    command <- ifelse(.Platform$OS.type == "windows", "where", "which")
    py.path <- system(sprintf("%s python", command), intern = TRUE)[1]
    reticulate::use_python(py.path)
    source_python(python.script)

    sp.embeddings <- sp.int@reductions$pca@cell.embeddings[, 1:npcs]
    sc.embeddings <- sc.int@reductions$pca@cell.embeddings[, 1:npcs]
    labels <- Idents(sp.int) %>%
        as.vector() %>%
        as.integer()
    netx.pred <- runNetModel(sc.embeddings, sp.embeddings, as.integer(labels), epochs = 1000) %>%
        `rownames<-`(colnames(sc.int)) %>%
        `colnames<-`(levels(sp.obj)) %>%
        as.matrix()
    target.cls.counts <- getClusterCellCounts(st.prop.cnts, sp.obj)
    cell.types <- Idents(sc.obj)[rownames(netx.pred)] %>% as.vector()
    assigned.lst <- assignCellsGreedyClusters(netx.pred, cell.types, target.cls.counts)$assigned_cluster
    return(assigned.lst)
}

#' @title linearSumAssignment
#'
#' @description Assigns single cells to spatial coordinates based on similarity scores between single-cell and spatial transcriptomics (ST) data.
#'
#' @param sp.obj Seurat object containing ST data.
#' @param sc.obj Seurat object containing SC data.
#' @param out.sim Similarity matrix between single cells and spots.
#' @param adj.w Weight matrix for adjusting similarity scores.
#' @param hot.pvals Data frame of hotspot p-values indicating hotspot presence in spots and cell types.
#' @param st.prop.cnts Estimated spot-level counts matrix (rows=spots, cols=cell types).
#' @param num.cells A vector of cell counts for each spot.
#' @param partition Logical indicating whether to split into sub-modules mapped to spatial positions.
#' @param assay Assay to use for SC data. Default: "SCT".
#' @param reduction Dimensionality reduction technique to align spatial and single-cell data. Default: cca.
#' @param thresh.cells If the number of cells exceeds this threshold, the partitioning strategy is executed automatically. Default: 30000.
#' @return A list of assigned cells corresponding to spots.
#' @export linearSumAssignment

linearSumAssignment <- function(sp.obj, sc.obj, out.sim, adj.w, hot.pvals, st.prop.cnts, num.cells, partition, assay, reduction, thresh.cells = 30000) {
    if (partition) {
        index.lst <- partitionClusters(sp.obj, sc.obj, hot.pvals, st.prop.cnts, assay, reduction)
    } else {
        if (nrow(out.sim) > thresh.cells) {
            index.lst <- partitionClusters(sp.obj, sc.obj, hot.pvals, st.prop.cnts, assay, reduction)
        } else {
            index.lst <- list(ENTIRE = NULL)
        }
    }
    tmp.dir <- tempdir()
    assign.res <- future.apply::future_lapply(names(index.lst), function(cls) {
        if (cls != "ENTIRE") {
            spot.sub <- colnames(sp.obj)[Idents(sp.obj) == cls]
            out.sim.sub <- out.sim[spot.sub, index.lst[[cls]]]
        } else {
            out.sim.sub <- out.sim
        }
        spot.name <- Idents(sp.obj)[rownames(out.sim.sub)]
        cell.names <- Idents(sc.obj)[colnames(out.sim.sub)]
        if (!is.null(adj.w)) {
            out.sim.sub <- weightSimScore(out.sim.sub, adj.w, spot.name, cell.names, hot.pvals)
        } else {
            out.sim.sub <- t(out.sim.sub)
        }
        sim.file <- file.path(tmp.dir, sprintf("sim_%s.xls", cls))
        num.file <- file.path(tmp.dir, sprintf("num_%s.xls", cls))
        data.table::fwrite(as.data.frame(out.sim.sub) * (-1), file = sim.file)
        data.table::fwrite(as.data.frame(num.cells[colnames(out.sim.sub)]), file = num.file)
        python.script <- system.file("python/solve.py", package = "Cell2Spatial")
        system(sprintf("python %s %s %s %s %s", python.script, tmp.dir, sim.file, num.file, cls))
        index <- read.table(file.path(tmp.dir, paste0(cls, "_output.txt"))) + 1
        spot.ids <- colnames(out.sim.sub)[index[, 1]]
        cbind.data.frame(CELL = rownames(out.sim.sub), SPOT = spot.ids)
    }, future.seed = TRUE) %>%
        do.call(rbind, .) %>%
        as.data.frame()
    out.sc <- split(assign.res$CELL, assign.res$SPOT)
    out.sc <- out.sc[colnames(sp.obj)]
    return(out.sc)
}

.linearSumAssignment <- function(...) {
    suppressWarnings({
        out.sc <- linearSumAssignment(...)
    })
    return(out.sc)
}
