#' @title integDataBySeurat

#' @description Integrate SC and ST data using Harmony method.
#' @param sp.obj Seurat object of spatial transcriptome (ST) data.
#' @param sc.obj Seurat object of single-cell data.
#' @param npcs Number of PCs for running UMAP. Default: 30.
<<<<<<< HEAD
#' @param integ Integrate SC and ST data or not. Default: TRUE.
#' @param verbose Show running messages or not. Default: TRUE.
#' @return A list of information about integration.
#' @export integDataBySeurat

integDataBySeurat <- function(sp.obj, sc.obj, npcs = 30, n.features = 3000, integ = TRUE, verbose = TRUE) {
    sc.obj$Batches <- "SC"
    sp.obj$Batches <- "ST"
    sc.obj <- SCTransform(sc.obj, verbose = verbose) %>% RunPCA(npcs = npcs, verbose = verbose)
    sp.obj <- SCTransform(sp.obj, assay = "Spatial", verbose = verbose) %>% RunPCA(npcs = npcs, verbose = verbose)
    if (integ) {
        ifnb.list <- list(SC = sc.obj, ST = sp.obj)
        features <- SelectIntegrationFeatures(object.list = ifnb.list, nfeatures = n.features)
        ifnb.list <- PrepSCTIntegration(object.list = ifnb.list, anchor.features = features)
        anchors <- FindIntegrationAnchors(
            object.list = ifnb.list,
            normalization.method = "SCT",
            anchor.features = features,
            verbose = verbose
        )
        obj <- IntegrateData(anchorset = anchors, normalization.method = "SCT", verbose = verbose)
        obj <- RunPCA(obj, verbose = verbose)
        obj <- RunUMAP(obj, reduction = "pca", dims = 1:npcs, verbose = verbose)
        adj.w <- adjcentScOfSP(obj)
    } else {
        adj.w <- obj <- NULL
    }
    ovp.genes <- intersect(rownames(sp.obj), rownames(sc.obj))
    return(list(obj = obj, sc.obj = sc.obj[ovp.genes, ], sp.obj = sp.obj[ovp.genes, ], adj.w = adj.w))
=======
#' @param verbose Show running messages or not. Default: TRUE.
#' @return Integrated Seurat object.
#' @export integDataBySeurat

integDataBySeurat <- function(sp.obj, sc.obj, npcs = 30, n.features = 3000, verbose = TRUE) {
  sc.obj$Batches <- "SC"
  sp.obj$Batches <- "ST"
  sc.obj <- SCTransform(sc.obj, verbose = verbose) %>% RunPCA(npcs = npcs, verbose = verbose)
  sp.obj <- SCTransform(sp.obj, assay = "Spatial", verbose = verbose) %>% RunPCA(npcs = npcs, verbose = verbose)
  ifnb.list <- list(SC = sc.obj, ST = sp.obj)
  features <- SelectIntegrationFeatures(object.list = ifnb.list, nfeatures = n.features)
  ifnb.list <- PrepSCTIntegration(object.list = ifnb.list, anchor.features = features)
  anchors <- FindIntegrationAnchors(
    object.list = ifnb.list,
    normalization.method = "SCT",
    anchor.features = features,
    verbose = verbose
  )
  obj <- IntegrateData(anchorset = anchors, normalization.method = "SCT", verbose = verbose)
  obj <- RunPCA(obj, verbose = verbose)
  obj <- RunUMAP(obj, reduction = "pca", dims = 1:npcs, verbose = verbose)
  ovp.genes <- intersect(rownames(sp.obj), rownames(sc.obj))
  return(list(obj = obj, sc.obj = sc.obj[ovp.genes, ], sp.obj = sp.obj[ovp.genes, ]))
>>>>>>> 5859e3527b186cec1e7cc9974128adabbcd95a88
}

#' @title findClustersForSpData

#' @description Find clusters for single-cell data based on Seurat object.
#' @param obj.seu Seurat object of ST data.
#' @param npcs Number of PCs for clustering. Default: 50.
#' @param res Resolution for finding clusters. Default: 0.8.
#' @return Seurat object after clustering.
#' @export findClustersForSpData

findClustersForSpData <- function(obj.seu, npcs = 30, res = 0.8, verbose = TRUE) {
<<<<<<< HEAD
    obj.seu <- obj.seu %>%
        FindNeighbors(., reduction = "pca", verbose = verbose) %>%
        FindClusters(., verbose = verbose) %>%
        RunUMAP(., reduction = "pca", dims = 1:npcs)
    return(obj.seu)
=======
  obj.seu <- obj.seu %>%
    FindNeighbors(., reduction = "pca", verbose = verbose) %>%
    FindClusters(., verbose = verbose) %>%
    RunUMAP(., reduction = "pca", dims = 1:npcs)
  return(obj.seu)
>>>>>>> 5859e3527b186cec1e7cc9974128adabbcd95a88
}

#' @title adjcentScOfSP

#' @description Generating distance matrix between ST spots and SC cells.
#' @param obj Integrated Seurat object.
#' @return A matrix of distance weights between ST and SC.

adjcentScOfSP <- function(obj) {
<<<<<<< HEAD
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
        res.dist <- (max(res.dist) - res.dist) / (max(res.dist) - min(res.dist))
    }, future.seed = TRUE) %>%
        do.call(rbind, .) %>%
        as.data.frame() %>%
        `rownames<-`(rownames(umap.sp))
    adj.df <- as.matrix(adj.df)
    return(adj.df)
=======
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
    res.dist <- (max(res.dist) - res.dist) / (max(res.dist) - min(res.dist))
  }, future.seed = TRUE) %>%
    do.call(rbind, .) %>%
    as.data.frame() %>%
    `rownames<-`(rownames(umap.sp))
  adj.df <- as.matrix(adj.df)
  return(adj.df)
>>>>>>> 5859e3527b186cec1e7cc9974128adabbcd95a88
}

#' @title getGsetScore

#' @description Get singature scores in each cell or spot using a list of markers.
#' @param obj.seu Seurat object.
#' @param sc.markers A list of markers for cell-types.
#' @param assay The assay used to calculate signature scores of cell types. Default: SCT.
#' @return A data.frame of signatrue scores.
#' @export getGsetScore

getGsetScore <- function(obj.seu, sc.markers, assay = "SCT") {
<<<<<<< HEAD
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
=======
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
>>>>>>> 5859e3527b186cec1e7cc9974128adabbcd95a88
}

#' @title estPropInSpots

#' @description Infer whether a cell type belongs to a specific cluster.
#' @param sp.obj Seurat object of spatial transcriptome (ST) data.
#' @param sc.obj Seurat object of single-cell (SC) data.
#' @param st.pvals A data.frame containing p-values to determine whether cell types are significantly present in a spot or not.
#' @param sc.markers A list of markers of cell types.
#' @param pval.cut significant cutoff value. Default: 0.05
#' @param knn K nearest cells used for estimating cellular proportions. Default: 5.
#' @param intercept Using intercept or not for robust linear regression. Default: TRUE.
#' @return A data frame of cellular fractions in each spot.

<<<<<<< HEAD
estPropInSpots <- function(sp.obj, sc.obj, st.pvals, sc.markers, pval.cut = 0.05, knn = 5, intercept = TRUE, chunk.size = 10000) {
    base.ref <-
        {
            AverageExpression(sc.obj, slot = "data", assay = "SCT")$SCT
        } %>%
        {
            .[unlist(sc.markers) %>% unique(), ]
        } %>%
        {
            log10(. + 1)
        }
    image.coord <- GetTissueCoordinates(sp.obj)[colnames(sp.obj), ]
    dist.mat <- dbscan::kNN(na.omit(image.coord[, c(1, 2)]), k = knn)$id

    sp.expr <- GetAssayData(sp.obj, slot = "data", assay = "SCT")[rownames(base.ref), ]
    cell.types <- colnames(base.ref)
    colnames(base.ref) <- gsub(" |/|\\+|-|\\(|\\)", "_", colnames(base.ref))

    if (ncol(sp.obj) > chunk.size) {
        tar.spots <- downSamplSeurat(sp.obj, percent = round(chunk.size / ncol(sp.obj), 1)) %>% colnames()
    } else {
        tar.spots <- colnames(sp.obj)
    }
    w.mat <- ctypesOfClusters(sp.obj, st.pvals)

    spot.prop <- future.apply::future_lapply(tar.spots, function(spot) {
        knn.spots <- dist.mat[spot, ] %>% rownames(dist.mat)[.]
        sp.expr.sub <- sp.expr[, knn.spots] %>%
            {
                10^.
            } %>%
            as.data.frame() %>%
            rowMeans() %>%
            {
                log10(.)
            }
        colnames(base.ref) <- paste0("X", colnames(base.ref))
        dat <- cbind.data.frame(base.ref, Bulk = sp.expr.sub)
        if (intercept) {
            w.coefs <- coef(MASS::rlm(Bulk ~ ., dat))
        } else {
            w.coefs <- coef(MASS::rlm(Bulk ~ . - 1, dat))
        }
        w.coefs[w.coefs < 0] <- 0
        coefs <- w.coefs[(intercept + 1):length(w.coefs)] / sum(w.coefs[(intercept + 1):length(w.coefs)])
        coefs <- coefs * w.mat[as.vector(Idents(sp.obj)[spot]), cell.types]
        coefs <- coefs / sum(coefs)
        names(coefs) <- cell.types
        return(coefs)
    }, future.seed = TRUE) %>%
        do.call(rbind, .) %>%
        as.data.frame() %>%
        `rownames<-`(colnames(tar.spots)) %>%
        `colnames<-`(cell.types)
    spot.prop <- spot.prop[!is.na(rowSums(spot.prop)), ]
    return(spot.prop)
=======
estPropInSpots <- function(sp.obj, sc.obj, st.pvals, sc.markers, pval.cut = 0.05, knn = 5, intercept = TRUE) {
  base.ref <- { AverageExpression(sc.obj, slot = "data", assay = "SCT")$SCT } %>%
    { .[unlist(sc.markers) %>% unique, ] } %>%
    {
      log10(. + 1)
    }
  dist.mat <- calcSpotsDist(sp.obj)
  sp.expr <- GetAssayData(sp.obj, slot = "data", assay = "SCT")[rownames(base.ref), ]
  cell.types <- colnames(base.ref)
  colnames(base.ref) <- gsub(" |/|\\+|-", "_", colnames(base.ref))
  spot.prop <- future.apply::future_lapply(colnames(sp.obj), function(spot) {
    knn.spots <- dist.mat[spot, ] %>%
      order(.) %>%
      {
        rownames(dist.mat)[.[1:(1 + knn)]]
      }
    sp.expr.sub <- sp.expr[, knn.spots] %>%
      {
        10^.
      } %>%
      as.data.frame() %>%
      rowMeans() %>%
      {
        log10(.)
      }
    colnames(base.ref) <- paste0("X", colnames(base.ref))
    dat <- cbind.data.frame(base.ref, Bulk = sp.expr.sub)
    if (intercept) {
      w.coefs <- coef(MASS::rlm(Bulk ~ ., dat))
    } else {
      w.coefs <- coef(MASS::rlm(Bulk ~ . - 1, dat))
    }
    w.coefs[w.coefs < 0] <- 0
    coefs <- w.coefs[(intercept + 1):length(w.coefs)] / sum(w.coefs[(intercept + 1):length(w.coefs)])
    names(coefs) <- cell.types
    adj.ctypes <- colnames(st.pvals)[st.pvals[spot, ] < pval.cut]
    if (length(adj.ctypes) > 0 && sum(coefs[adj.ctypes]) > 0) {
      coefs[setdiff(names(coefs), adj.ctypes)] <- 0
      coefs <- coefs / sum(coefs)
    }
    return(coefs)
  }, future.seed = TRUE) %>%
    do.call(rbind, .) %>%
    as.data.frame() %>%
    `rownames<-`(colnames(sp.obj)) %>%
    `colnames<-`(cell.types)
  return(spot.prop)
>>>>>>> 5859e3527b186cec1e7cc9974128adabbcd95a88
}

#' @title estPvalsOfSpot

#' @description Infer whether signature score in spot significant higher than background.
#' @param sp.obj Seurat object of spatial transcriptome (ST) data.
#' @param sp.score Signature score of cell types.
#' @param knn K nearest cells used. Default: 5.
#' @return A data.frame of p-values.

estPvalsOfSpot <- function(sp.obj, sp.score, knn = 5) {
<<<<<<< HEAD
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
            # t.test(sp.score[knn.spots, ctype], sp.score[-spot.ids[spot, ], ctype], alternative = "greater")$p.value
        }) %>%
            unlist() %>%
            as.vector() %>%
            `names<-`(names(idxes))
        pval.ctype[names(sig.spots)] <- sig.spots
        pval.ctype
    }, future.seed = TRUE) %>%
        do.call(rbind, .) %>%
        `rownames<-`(rownames(sp.score))
    return(spot.pvals)
}

#' @title ctypesOfClusters

ctypesOfClusters <- function(sp.obj, st.pvals) {
    w.mat <- lapply(levels(sp.obj), function(cls) {
        spots <- colnames(sp.obj)[Idents(sp.obj) == cls]
        sub.pvals <- st.pvals[spots, ]
        tmp.cnt <- colSums(sub.pvals < 0.01)
        w.gs <- scales::rescale(tmp.cnt %>% as.numeric(), to = c(0, 1))
    }) %>%
        do.call(rbind, .) %>%
        as.data.frame() %>%
        `rownames<-`(levels(sp.obj)) %>%
        `colnames<-`(colnames(st.pvals))
    return(w.mat)
=======
  dist.mat <- calcSpotsDist(sp.obj)
  dist.lst <- lapply(levels(sp.obj), function(idx) {
    sub.spots <- Idents(sp.obj)[Idents(sp.obj) == idx] %>% names()
    dist.mat[sub.spots, sub.spots]
  })
  rm(dist.mat)
  gc()
  spot.pvals <- future.apply::future_lapply(rownames(sp.score), function(spot) {
    dist.mat.sub <- dist.lst[[Idents(sp.obj)[spot]]]
    knn.spots <- dist.mat.sub[spot, ] %>%
      order(.) %>%
      {
        rownames(dist.mat.sub)[.[1:(1 + knn)]]
      }
    sig.spots <- lapply(colnames(sp.score), function(ctype) {
      t.test(sp.score[knn.spots, ctype], sp.score[, ctype], alternative = "greater")$p.value
    }) %>%
      unlist() %>%
      as.vector() %>%
      `names<-`(colnames(sp.score))
  }) %>%
    do.call(rbind, .) %>%
    `rownames<-`(rownames(sp.score))

  return(spot.pvals)
>>>>>>> 5859e3527b186cec1e7cc9974128adabbcd95a88
}

#' @title selectScByCor

#' @description Select single-cells based on weighted correlation.
#' @param sp.score Score matrix of ST inferred from a list of markers.
#' @param sc.score Score matrix of SC inferred from a list of markers.
#' @param adj.w A weighted matrix infered from distance. Rows represent spot names, and cloumns are single-cells.
#' @return Similarity matrix.

selectScByCor <- function(sp.score, sc.score, adj.w) {
<<<<<<< HEAD
    cor.mat <- cor(t(sp.score), t(sc.score))
    out.sim <-
        {
            if (inherits(adj.w, "matrix")) {
                adj.w <- adj.w[, rownames(sc.score)]
                cor.mat * adj.w
            } else {
                cor.mat
            }
        } %>% t()
    return(out.sim)
=======
  adj.w <- adj.w[, rownames(sc.score)]
  cor.mat <- cor(t(sp.score), t(sc.score))
  out.sim <-
    {
      cor.mat * adj.w
    } %>% t()
  return(out.sim)
>>>>>>> 5859e3527b186cec1e7cc9974128adabbcd95a88
}

#' @title selectByProb

#' @description Measuring the similarity distances using probability model.
#' @param sp.score
#' @param sc.score
#' @param adj.w
#' @return Similarity matrix.

selectByProb <- function(sp.score, sc.score, adj.w) {
<<<<<<< HEAD
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
    out.sim <-
        {
            if (inherits(adj.w, "matrix")) {
                out.sim * adj.w
            } else {
                out.sim
            }
        } %>% t()
    return(out.sim)
}

#' @title featureSelelction

#' @description
#' @param sp.obj Seurat object of ST data.
#' @param sc.obj Seurat object of SC data.
#' @param n.features Common features between SC and ST data.
#' @return
#' @export featureSelelction

featureSelelction <- function(sp.obj, sc.obj, n.features) {
    sc.obj <- SCTransform(sc.obj, verbose = FALSE)
    ifnb.list <- list(SC = sc.obj, ST = sp.obj)
    features <- SelectIntegrationFeatures(object.list = ifnb.list, nfeatures = n.features)
    sc.st.anchors <- Seurat::FindTransferAnchors(
        reference = sc.obj,
        query = sp.obj,
        reference.assay = "SCT",
        query.assay = "SCT",
        features = features,
        reduction = "cca"
    )
    st.data.trans <- Seurat::TransferData(
        anchorset = sc.st.anchors,
        refdata = GetAssayData(sc.obj, assay = "SCT", slot = "data")[features, ],
        weight.reduction = "cca"
    )
    sp.obj@assays$transfer <- st.data.trans

    st.vv <- data.frame(sp.obj[["transfer"]]@data) %>% `colnames<-`(colnames(sp.obj))
    sc.vv <- data.frame(sc.obj[["SCT"]]@data[features, ]) %>% `colnames<-`(colnames(sc.obj))
    counts.temp <- cbind(st.vv, sc.vv)
    sc.st.int <- CreateSeuratObject(counts = counts.temp, assay = "traint")
    sc.st.int[["traint"]]@data <- sc.st.int[["traint"]]@counts
    sc.st.int[["traint"]]@counts <- matrix(NA, nrow = 0, ncol = 0)
    sc.st.int <- ScaleData(sc.st.int, features = features) %>% RunPCA(features = features)
    sc.int <- subset(sc.st.int, cells = colnames(sc.obj))
    st.int <- subset(sc.st.int, cells = colnames(sp.obj))
    Idents(st.int) <- Idents(sp.obj)
    return(list(sc = sc.int, st = st.int))
}

#' @title partionClusters

#' @description Partion single-cells to sub-clusters using deep learning strategy.
#' @param sp.obj Seurat object of ST data.
#' @param sc.obj Seurat object of SC data.
#' @param num.cells A vector of cell counts in spots.
#' @param n.features Number of common features between SC and ST data.
#' @return A list of single-cells in each cluster.
#' @export partionClusters

partionClusters <- function(sp.obj, sc.obj, n.features, num.cells) {
    sc.st.int <- featureSelelction(sp.obj, sc.obj, n.features)
    sc.int <- sc.st.int$sc
    sp.int <- sc.st.int$st
    python.script <- system.file("R/netx.py", package = "Cell2Spatial")
    command <- ifelse(.Platform$OS.type == "windows", "where", "which")
    py.path <- system(sprintf("%s python", command), intern = TRUE)
    reticulate::use_python(py.path)
    source_python(python.script)

    sp.embeddings <- sp.int@reductions$pca@cell.embeddings[, 1:30]
    sc.embeddings <- sc.int@reductions$pca@cell.embeddings[, 1:30]
    labels <- Idents(sp.int) %>%
        as.vector() %>%
        as.integer()
    netx.pred <- runNetModel(sc.embeddings, sp.embeddings, labels, epochs = 1000) %>%
        `colnames<-`(levels(sp.int)) %>%
        `rownames<-`(colnames(sc.int))
    target.counts <- lapply(levels(sp.int), function(xx) {
        spots <- Idents(sp.int)[Idents(sp.int) == xx] %>% names()
        sum(num.cells[spots])
    }) %>%
        unlist() %>%
        `names<-`(levels(sp.int))

    netx.pred.scale <- sweep(netx.pred, MARGIN = 1, apply(netx.pred, 1, max), "/")
    max.idxes <- apply(netx.pred, 1, which.max) %>%
        colnames(netx.pred)[.] %>%
        `names<-`(rownames(netx.pred))
    pred.counts <- rep(0, length(levels(sp.int))) %>% `names<-`(levels(sp.int))
    pred.counts[names(table(max.idxes))] <- table(max.idxes)
    cell.diff <- pred.counts - target.counts

    index.lst <- list()
    flag.id <- flag.id.update <- c()
    while (any(cell.diff) > 0) {
        for (cluster in setdiff(names(target.counts), flag.id.update)) {
            if (cell.diff[cluster] > 0) {
                tmp.index <- netx.pred[, cluster] %>%
                    `names<-`(rownames(netx.pred)) %>%
                    .[order(.) %>% rev()]
                index.lst[[cluster]] <- netx.pred.scale[, cluster] %>%
                    `names<-`(rownames(netx.pred.scale)) %>%
                    .[names(tmp.index)] %>%
                    {
                        .[. == 1]
                    } %>%
                    .[1:target.counts[cluster]] %>%
                    names()
                flag.id <- c(flag.id, cluster)
                pred.counts[cluster] <- target.counts[cluster]
            }
        }
        remain.cells <- setdiff(rownames(netx.pred), unlist(index.lst))
        if (length(setdiff(colnames(netx.pred), flag.id)) == 0) break
        netx.pred <- netx.pred[remain.cells, setdiff(colnames(netx.pred), flag.id), drop = FALSE]
        netx.pred.scale <- sweep(netx.pred, MARGIN = 1, apply(netx.pred, 1, max), "/")

        max.idxes <- apply(netx.pred, 1, which.max) %>%
            colnames(netx.pred)[.] %>%
            `names<-`(rownames(netx.pred))
        pred.counts[names(table(max.idxes))] <- table(max.idxes)
        cell.diff <- pred.counts - target.counts
        flag.id.update <- c(flag.id.update, flag.id)
        flag.id <- c()
    }
    index.lst[[setdiff(names(target.counts), names(index.lst))]] <- setdiff(colnames(sc.int), unlist(index.lst))
    index.lst <- index.lst[names(target.counts)]
    return(index.lst)
}

#' @title linearSumAssignment

#' @description Assign single-cells to spatial coordinates.
#' @param sp.obj Seurat object of ST data.
#' @param sc.obj Seurat object of SC data.
#' @param out.sim Similarity matrix between single-cells and spots.
#' @param num.cells A vector of cell count per spot projected to spots. Default: 10.
#' @param n.features Number of common features between SC and ST data.
#' @param partion Split into sub-modules mapped to spatial positions. Default: TRUE.
#' @return A list of assigned cells corresponding to spots.
#' @export linearSumAssignment

linearSumAssignment <- function(sp.obj, sc.obj, out.sim, num.cells, n.features, partion = TRUE) {
    if (partion) {
        index.lst <- partionClusters(sp.obj, sc.obj, n.features, num.cells)
    } else {
        if (nrow(out.sim) > 30000) {
            index.lst <- partionClusters(sp.obj, sc.obj, n.features, num.cells)
        } else {
            index.lst <- list(ENTIRE = NULL)
        }
    }
    assign.res <- future.apply::future_lapply(names(index.lst), function(cls) {
        if (cls != "ENTIRE") {
            spot.sub <- colnames(sp.obj)[Idents(sp.obj) == cls]
            out.sim.sub <- out.sim[index.lst[[cls]], spot.sub]
        } else {
            out.sim.sub <- out.sim
        }
        tmp.dir <- tempdir()
        sim.file <- file.path(tmp.dir, sprintf("sim_%s.xls", cls))
        num.file <- file.path(tmp.dir, sprintf("num_%s.xls", cls))
        data.table::fwrite(as.data.frame(out.sim.sub) * (-1), file = sim.file)
        data.table::fwrite(as.data.frame(num.cells[colnames(out.sim.sub)]), file = num.file)

        python.script <- system.file("R/solve.py", package = "Cell2Spatial")
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
=======
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
  out.sim <-
    {
      out.sim * adj.w
    } %>% t()
  return(out.sim)
>>>>>>> 5859e3527b186cec1e7cc9974128adabbcd95a88
}

#' @title assignSc2SP

#' @description Assign single-cells to spots using linear assignment method.
#' @param sp.obj Seurat object of ST data.
#' @param sc.obj Seurat object of SC data.
#' @param out.sim Similarity matrix between single-cells and spots.
<<<<<<< HEAD
#' @param num.cells A vector of cell count per spot projected to spots. Default: 10.
#' @param st.pvals A data.frame of p-values of ST spots.
#' @param n.features Number of common features between SC and ST data.
#' @param duplicated Unique assign single-cells to spot or not. Default: FALSE.
#' @param Split into sub-modules mapped to spatial positions. Default: TRUE.
#' @param p.cut Cutoff pvalue to determine negative spots. Default: 0.05
#' @return A list of assigned cells corresponding to spots.

assignSc2SP <- function(sp.obj, sc.obj, out.sim, num.cells, st.pvals, n.features, duplicated = FALSE, partion = TRUE, p.cut = 0.05) {
    idxes <- which(st.pvals > p.cut)
    out.sim[idxes] <- out.sim[idxes] - 1
    if (!duplicated) {
        out.sc <- linearSumAssignment(sp.obj, sc.obj, out.sim, num.cells, n.features, partion)
    } else {
        out.sc <- lapply(1:ncol(out.sim), function(xx) {
            spot <- colnames(out.sim)[xx]
            bb.x <- out.sim[, xx]
            bb.x[order(bb.x) %>% rev()] %>%
                names(.) %>%
                .[1:num.cells[[xx]]] %>%
                .[!is.na(.)]
        }) %>%
            `names<-`(colnames(out.sim)) %>%
            .[(lapply(., length) > 0)]
    }
    return(out.sc)
=======
#' @param num.cells Number of cells projected to spots. Default: 10.
#' @param st.pvals A data.frame of p-values of ST spots.
#' @param duplicated Unique assign single-cells to spot or not. Default: FALSE.
#' @return A list of assigned cells corresponding to spots.

assignSc2SP <- function(sp.obj, sc.obj, out.sim, num.cells, st.pvals, duplicated = FALSE) {
  sc.idents <- Idents(sc.obj)
  for (spot in rownames(st.pvals)) {
    neg.ctypes <- st.pvals[spot, ] %>%
      unlist() %>%
      {
        .[. == 0]
      } %>%
      names()
    neg.cells <- sc.idents[sc.idents %in% neg.ctypes]
    out.sim[names(neg.cells), spot] <- 0
  }
  if (!duplicated) {
    tmp.dir <- tempdir()
    sim.file <- file.path(tmp.dir, "sim.xls")
    num.file <- file.path(tmp.dir, "num.xls")
    data.table::fwrite(as.data.frame(out.sim) * (-1), file = sim.file)
    data.table::fwrite(as.data.frame(num.cells), file = num.file)

    python.script <- system.file("R/solve.py", package = "Cell2Spatial")
    system(sprintf("python %s %s", python.script, tmp.dir))

    index <- read.table(file.path(tmp.dir, "output.txt")) + 1
    spot.ids <- colnames(out.sim)[index[, 1]]
    out.sc <- split(rownames(out.sim), spot.ids)
  } else {
    out.sc <- lapply(1:ncol(out.sim), function(xx) {
      spot <- colnames(out.sim)[xx]
      bb.x <- out.sim[, xx]
      bb.x[order(bb.x) %>% rev()] %>%
        names(.) %>%
        .[1:num.cells[[xx]]] %>%
        .[!is.na(.)]
    }) %>%
      `names<-`(colnames(out.sim)) %>%
      .[(lapply(., length) > 0)]
  }
  return(out.sc)
>>>>>>> 5859e3527b186cec1e7cc9974128adabbcd95a88
}
