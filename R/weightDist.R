#' @title psedoSpotExprUseSC

#' @description Pseudo purity spots expression based on SC data.
#' @param sc.obj Seurat object of SC data.
#' @param sp.obj Seurat object of ST data.
#' @param pseu.cnt Pseudo cell counts for each cell type from SC data. Default: 200.
#' @param lamba Median of cell counts for spots.
#' @param mc.cores Number of cores for parallel running. Default: 4
#' @return Seurat object of pseudo purity spots data.

psedoSpotExprUseSC <- function(sc.obj, sp.obj, pseu.cnt = 200, lamba = 5, mc.cores = 4) {
    ctypes <- levels(sc.obj)
    med.umis <- GetAssayData(sp.obj, slot = "count") %>%
        as.matrix() %>%
        colSums(.) %>%
        median(.)
    count.mat <- GetAssayData(sc.obj, slot = "count") %>% as.matrix()
    cell.lst <- split(colnames(sc.obj), Idents(sc.obj) %>% as.vector())
    syn.spot <- parallel::mclapply(ctypes, function(ct) {
        rand.mat <- lapply(1:pseu.cnt, function(xx) {
            ct.tmp <- sample(cell.lst[[ct]], lamba, replace = TRUE)
            if (length(ct.tmp) > 1) {
                rowSums(count.mat[, ct.tmp])
            } else {
                count.mat[, ct.tmp]
            }
        }) %>%
            do.call(cbind, .) %>%
            as.data.frame()
        syn.spot <- DropletUtils::downsampleMatrix(rand.mat, prop = med.umis / colSums(rand.mat))
        colnames(syn.spot) <- paste0(ct, "_pseu", 1:ncol(syn.spot))
        syn.spot
    }, mc.cores = mc.cores) %>% do.call(cbind, .)
    meta.data <- data.frame(CellType = gsub("_pseu.*", "", colnames(syn.spot))) %>% `rownames<-`(colnames(syn.spot))
    syn.spot <- CreateSeuratObject(count = syn.spot, meta.data = meta.data)
    return(syn.spot)
}

#' @title integDataBySeurat

#' @description Integrate SC and ST data normalized by the LogNormalize method using the CCA method.
#' @param sp.obj Seurat object of spatial transcriptome (ST) data.
#' @param sc.obj Seurat object of single-cell data.
#' @param npcs Number of PCs for running UMAP. Default: 30.
#' @param n.features Number of features set for SelectIntegrationFeatures. Default: 3000.
#' @param verbose Show running messages or not. Default: TRUE.
#' @return Integrated Seurat object.
#' @export integDataBySeurat

integDataBySeurat <- function(sp.obj, sc.obj, npcs = 30, n.features = 3000, verbose = TRUE) {
    sc.obj$Batches <- "SC"
    sp.obj$Batches <- "ST"
    DefaultAssay(sc.obj) <- "RNA"
    DefaultAssay(sp.obj) <- "Spatial"

    options(warn = -1)
    sc.obj <- sc.obj %>%
        NormalizeData(verbose = verbose) %>%
        FindVariableFeatures(verbose = verbose)
    sp.obj <- sp.obj %>%
        NormalizeData(verbose = verbose, assay = "Spatial") %>%
        FindVariableFeatures(verbose = verbose)
    ifnb.list <- list(SC = sc.obj, ST = sp.obj)
    features <- SelectIntegrationFeatures(object.list = ifnb.list, nfeatures = n.features)
    anchors <- FindIntegrationAnchors(
        object.list = ifnb.list,
        normalization.method = "LogNormalize",
        anchor.features = features,
        verbose = verbose
    )
    obj <- IntegrateData(anchorset = anchors, normalization.method = "LogNormalize", verbose = verbose)
    obj <- ScaleData(obj, verbose = verbose)
    obj <- RunPCA(obj, verbose = verbose)
    obj <- RunUMAP(obj, reduction = "pca", dims = 1:npcs, verbose = verbose)
    return(obj)
}

#' @title findClustersForSpData

#' @description Find clusters for single-cell data based on Seurat object.
#' @param obj.seu Seurat object of ST data.
#' @param npcs Number of PCs for clustering. Default: 50.
#' @param res Resolution parameter for finding clusters. Default: 0.8.
#' @return Seurat object after clustering.

findClustersForSpData <- function(obj.seu, npcs = 30, res = 0.8, verbose = TRUE) {
    obj.seu <- obj.seu %>%
        FindNeighbors(., reduction = "pca", verbose = verbose) %>%
        FindClusters(., verbose = verbose) %>%
        RunUMAP(., reduction = "pca", dims = 1:npcs, verbose = verbose)
    return(obj.seu)
}

#' @title adjcentScOfSP

#' @description Generating distance matrix between ST spots and SC cells.
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

    adj.df <- apply(spot.distances.wide, 1, function(res.dist) (max(res.dist) - res.dist) / (max(res.dist) - min(res.dist))) %>%
        t() %>%
        as.matrix()
    return(adj.df)
}

#' @title weightDist

#' @description Weight distance between SC and ST data.

#' @param sc.obj Seurat object of SC data.
#' @param sp.obj Seurat object of ST data.
#' @param lamba Median of cell counts for spots.
#' @param mc.cores Number of cores for parallel running. Default: 4
#' @return A matrix of weighted distance matrix. Rows represent spot clusters and columns are cell types.

weightDist <- function(sc.obj, sp.obj, lamba, mc.cores = 4) {
    sc.syn <- psedoSpotExprUseSC(sc.obj, sp.obj, pseu.cnt = 200, lamba = lamba, mc.cores = mc.cores)
    sp.syn <- downSamplSeurat(sp.obj, cnt = 200)
    adj.df <- integDataBySeurat(sp.syn, sc.syn, verbose = FALSE) %>%
        adjcentScOfSP(.) %>%
        as.data.frame(.) %>%
        mutate(CLUSTER = Idents(sp.obj)[rownames(.)]) %>%
        group_by(CLUSTER) %>%
        summarise(across(all_of(levels(sc.obj)), ~ mean(sort(.x, decreasing = TRUE)[1:10], na.rm = TRUE))) %>%
        data.frame() %>%
        `rownames<-`(.[, 1]) %>%
        {
            .[, -1]
        } %>%
        `colnames<-`(levels(sc.obj))
    message("")
    return(adj.df)
}
