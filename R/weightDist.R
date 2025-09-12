#' @title pseudoSpotExprUseSC
#'
#' @description Pseudo purity spots expression based on SC data.
#'
#' @param sc.obj Seurat object of SC data.
#' @param sp.obj Seurat object of ST data.
#' @param assay Assay type to use for SC data normalization.
#' @param pseu.cnt Pseudo cell counts for each cell type from SC data. Default: 200.
#' @param lamba Median of cell counts for spots.
#' @return Seurat object of pseudo purity spots data.

pseudoSpotExprUseSC <- function(sc.obj, sp.obj, assay, pseu.cnt = 200, lamba = 5) {
    ctypes <- levels(sc.obj)
    count.mat <- GetAssayData(sc.obj, slot = "counts") %>% as.matrix()
    cell.lst <- split(colnames(sc.obj), Idents(sc.obj) %>% as.vector())
    syn.spot <- future.apply::future_lapply(ctypes, function(ct) {
        syn.spot <- lapply(1:pseu.cnt, function(xx) {
            ct.tmp <- sample(cell.lst[[ct]], lamba, replace = TRUE)
            if (length(ct.tmp) > 1) {
                rowMeans(count.mat[, ct.tmp]) %>% round(.)
            } else {
                count.mat[, ct.tmp]
            }
        }) %>%
            do.call(cbind, .) %>%
            as.data.frame()
        colnames(syn.spot) <- paste0(ct, "_pseu", 1:ncol(syn.spot))
        syn.spot %>% as.data.frame()
    }, future.seed = TRUE) %>% do.call(cbind, .)
    meta.data <- data.frame(CellType = gsub("_pseu.*", "", colnames(syn.spot))) %>% `rownames<-`(colnames(syn.spot))
    syn.spot <- CreateSeuratObject(count = syn.spot, meta.data = meta.data)
    if (assay == "SCT") {
        syn.spot <- syn.spot %>% SCTransform(verbose = FALSE)
    } else {
        syn.spot <- syn.spot %>%
            NormalizeData(verbose = FALSE) %>%
            FindVariableFeatures(., nfeatures = 3000, verbose = FALSE) %>%
            ScaleData(., verbose = FALSE)
    }
    return(syn.spot)
}

#' @title integDataBySeurat
#'
#' @description Integrate SC and ST data normalized by the LogNormalize method using the CCA method.
#'
#' @param sp.obj Seurat object of spatial transcriptome (ST) data.
#' @param sc.obj Seurat object of single-cell data.
#' @param assay Assay type to use for SC data normalization.
#' @param reduction Dimensionality reduction technique to align spatial and single-cell data.
#' @param npcs Number of PCs for running UMAP. Default: 30.
#' @param n.features Number of features set for SelectIntegrationFeatures. Default: 3000.
#' @param dist.based Dimensionality reduction basis used for distance weighting, UMAP or TSNE. Default: UMAP.
#' @param verbose Show running messages or not. Default: TRUE.
#' @return Integrated Seurat object.

integDataBySeurat <- function(sp.obj, sc.obj, assay, reduction, npcs = 30, n.features = 3000, dist.based = c("UMAP", "TSNE"), verbose = TRUE) {
    sc.obj$Batches <- "SC"
    sp.obj$Batches <- "ST"
    if (assay == "RNA") {
        DefaultAssay(sc.obj) <- "RNA"
        DefaultAssay(sp.obj) <- "Spatial"
    }
    normalization.method <- ifelse(assay == "SCT", "SCT", "LogNormalize")
    ifnb.list <- list(SC = sc.obj, ST = sp.obj)
    features <- SelectIntegrationFeatures(object.list = ifnb.list, nfeatures = n.features, verbose = verbose)
    if (reduction == "rpca") {
        ifnb.list <- lapply(ifnb.list, function(obj) {
            obj <- ScaleData(obj, features = features, verbose = verbose)
            obj <- RunPCA(obj, features = features, verbose = verbose)
            return(obj)
        })
    }
    anchors <- FindIntegrationAnchors(
        object.list = ifnb.list,
        normalization.method = normalization.method,
        anchor.features = features,
        reduction = reduction,
        verbose = verbose
    )
    obj <- IntegrateData(anchorset = anchors, normalization.method = normalization.method, k.weight = 30, verbose = verbose)
    obj <- ScaleData(obj, verbose = verbose)
    obj <- RunPCA(obj, verbose = verbose)
    obj <- switch(match.arg(dist.based),
        UMAP = {
            obj <- RunUMAP(obj, reduction = "pca", dims = 1:npcs, reduction.key = "UMAP_", verbose = verbose)
        },
        TSNE = {
            obj <- RunTSNE(obj, reduction = "pca", dims = 1:npcs, reduction.key = "UMAP_", verbose = verbose)
        }
    )
    return(obj)
}

#' @title findClustersForSpData
#'
#' @description Perform clustering on spatial transcriptomics data based on a Seurat object.
#'
#' @param obj.seu A Seurat object containing the spatial transcriptomics (ST) data.
#' @param npcs Number of principal components (PCs) to use for clustering. Default: 30.
#' @param res.start Initial resolution parameter for finding clusters. Higher values result in more clusters. Default: 0.1.
#' @param res.step Step size by which the resolution parameter is incremented in each iteration. Default: 0.05.
#' @param assay The assay to use for the analysis. Default: "SCT" (Seurat’s SCTransform).
#' @param cap The maximum allowed number of cells in a cluster. The clustering resolution will adjust to ensure no cluster exceeds this cap. Default: 20000.
#' @param verbose Logical, whether to print progress messages during the process. Default: TRUE.
#' @return A Seurat object after clustering and UMAP computation.

findClustersForSpData <- function(obj.seu, npcs = 30, res.start = 0.8, res.step = 0.05, assay = "SCT", cap = 20000, verbose = TRUE) {
    if (assay != "SCT") {
        obj.seu <- obj.seu %>%
            FindVariableFeatures(., nfeatures = 3000, verbose = verbose) %>%
            ScaleData(., verbose = verbose) %>%
            RunPCA(., verbose = verbose)
    }
    obj.seu <- obj.seu %>%
        RunPCA(., verbose = verbose) %>%
        FindNeighbors(., dims = 1:npcs, reduction = "pca", verbose = verbose)

    res <- res.start
    repeat {
        obj.seu <- FindClusters(obj.seu, resolution = res, verbose = verbose)
        cluster.sizes <- table(Idents(obj.seu))
        if (all(cluster.sizes <= cap) || res >= res.max) {
            break
        }
        res <- res + res.step
    }
    obj.seu <- obj.seu %>%
        FindClusters(., resolution = res, verbose = verbose) %>%
        RunUMAP(., reduction = "pca", dims = 1:npcs, verbose = verbose)
    return(obj.seu)
}

.findClustersForSpData <- function(...) {
    suppressWarnings({
        sp.obj <- findClustersForSpData(...)
    })
    return(sp.obj)
}

#' @title adjcentScOfSP
#'
#' @description Generating distance matrix between ST spots and SC cells.
#'
#' @param obj Integrated Seurat object.
#' @return A matrix of distance weights between ST and SC.

adjcentScOfSP <- function(obj) {
    umap.coord <- FetchData(obj, vars = c("UMAP_1", "UMAP_2", "Batches"))
    umap.sc <- subset(umap.coord, Batches == "SC")[, 1:2]
    ctype.centroids <- umap.sc %>%
        as.data.frame() %>%
        mutate(celltype = subset(obj, Batches == "SC")$CellType) %>%
        group_by(celltype) %>%
        summarize(UMAP_1 = mean(UMAP_1), UMAP_2 = mean(UMAP_2))

    umap.sp <- subset(umap.coord, Batches == "ST")[, 1:2]
    umap.sp <- umap.sp %>% mutate(spot_name = rownames(umap.sp))
    spot.distances <- umap.sp %>%
        rowwise() %>%
        do({
            spot <- .
            distances <- ctype.centroids %>%
                mutate(distance = sqrt((spot$UMAP_1 - UMAP_1)^2 + (spot$UMAP_2 - UMAP_2)^2)) %>%
                select(celltype, distance)
            spot_name <- spot$spot_name
            distances %>% mutate(spot_name = spot_name)
        }) %>%
        ungroup()

    spot.distances <- spot.distances %>% distinct(spot_name, celltype, .keep_all = TRUE)
    spot.distances.wide <- spot.distances %>%
        pivot_wider(names_from = celltype, values_from = distance) %>%
        as.data.frame() %>%
        `rownames<-`(.$spot_name) %>%
        .[, -1]

    adj.df <- apply(spot.distances.wide, 1, function(res.dist) (max(res.dist) - res.dist) / (1e-10 + max(res.dist) - min(res.dist))) %>%
        t() %>%
        as.matrix()
    return(adj.df)
}

#' @title adjcentScOfSPGlobal
#'
#' @description Generating distance matrix between ST spots and SC cells.
#'
#' @param obj Integrated Seurat object.
#' @param quantile.cut Numeric value specifying the quantile threshold for distance scaling. Default is 0.75, which considers the maximum distance for normalization.
#' @return A matrix of distance weights between ST and SC.

adjcentScOfSPGlobal <- function(obj, quantile.cut = 0.75) {
    umap.coord <- FetchData(obj, vars = c("UMAP_1", "UMAP_2", "Batches"))
    umap.sp <- subset(umap.coord, Batches == "ST")[, 1:2]
    umap.sc <- subset(umap.coord, Batches == "SC")[, 1:2]
    adj.df <- future.apply::future_lapply(1:nrow(umap.sp), function(idx) {
        xx <- umap.sp[idx, ] %>%
            unlist() %>%
            as.vector()
        res.dist <- apply(umap.sc, 1, function(yy) {
            sqrt(sum((xx - yy)^2))
        })
        res.dist <- (quantile(res.dist, quantile.cut) - res.dist) / (1e-10 + quantile(res.dist, quantile.cut) - min(res.dist))
    }, future.seed = TRUE) %>%
        do.call(rbind, .) %>%
        as.data.frame() %>%
        `rownames<-`(rownames(umap.sp))
    adj.df <- as.matrix(adj.df)
    adj.df[adj.df < 0] <- 1e-10
    return(adj.df)
}

#' @title weightDist
#'
#' @description Weight distance between SC and ST data.
#'
#' @param sc.obj Seurat object of SC data.
#' @param sp.obj Seurat object of ST data.
#' @param lamba Median of cell counts for spots.
#' @param assay The assay to use for the normalize analysis.
#' @param reduction Dimensionality reduction technique to align spatial and single-cell data.
#' @param cells.thresh Threshold for cell counts in the SC object. When the cell count exceeds this threshold, the function switches to a faster mode. Default: 50000.
#' @param dist.based Dimensionality reduction basis used for distance weighting. Default: UMAP.
#' @return A matrix of weighted distance matrix. Rows represent spot clusters and columns are cell types.

weightDist <- function(sc.obj, sp.obj, lamba, assay, reduction, cells.thresh = 50000, dist.based = c("UMAP", "TSNE")) {
    use.entire <- ifelse(ncol(sc.obj) < cells.thresh, TRUE, FALSE)
    if (use.entire || length(levels(sc.obj)) == 1) {
        adj.df <- integDataBySeurat(sp.obj, sc.obj, assay, reduction, dist.based = dist.based, verbose = FALSE) %>%
            adjcentScOfSPGlobal(.)
    } else {
        sc.syn <- pseudoSpotExprUseSC(sc.obj, sp.obj, assay, lamba = lamba)
        sp.syn <- downSamplSeurat(sp.obj, cnt = 200)
        adj.df <- integDataBySeurat(sp.syn, sc.syn, assay, reduction, dist.based = dist.based, verbose = FALSE) %>%
            adjcentScOfSP(.) %>%
            as.data.frame(.) %>%
            mutate(CLUSTER = Idents(sp.obj)[rownames(.)]) %>%
            group_by(CLUSTER) %>%
            summarise(across(all_of(levels(sc.obj)), ~ mean(sort(.x, decreasing = TRUE), na.rm = TRUE))) %>%
            data.frame() %>%
            `rownames<-`(.[, 1]) %>%
            {
                .[, -1]
            } %>%
            `colnames<-`(levels(sc.obj))
    }
    return(adj.df)
}

.weightDist <- function(spatial.weight, ...) {
    if (spatial.weight) {
        suppressWarnings({
            adj.w <- weightDist(...)
        })
    } else {
        adj.w <- NULL
    }
    return(adj.w)
}
