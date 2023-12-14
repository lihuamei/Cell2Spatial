#' @title integDataBySeurat

#' @description Integrate SC and ST data using Harmony method.
#' @param sp.obj Seurat object of spatial transcriptome (ST) data.
#' @param sc.obj Seurat object of single-cell data.
#' @param npcs Number of PCs for running UMAP. Default: 30.
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
}

#' @title findClustersForSpData

#' @description Find clusters for single-cell data based on Seurat object.
#' @param obj.seu Seurat object of ST data.
#' @param npcs Number of PCs for clustering. Default: 50.
#' @param res Resolution for finding clusters. Default: 0.8.
#' @return Seurat object after clustering.
#' @export findClustersForSpData

findClustersForSpData <- function(obj.seu, npcs = 30, res = 0.8, verbose = TRUE) {
  obj.seu <- obj.seu %>%
    FindNeighbors(., reduction = "pca", verbose = verbose) %>%
    FindClusters(., verbose = verbose) %>%
    RunUMAP(., reduction = "pca", dims = 1:npcs)
  return(obj.seu)
}

#' @title adjcentScOfSP

#' @description Generating distance matrix between ST spots and SC cells.
#' @param obj Integrated Seurat object.
#' @return A matrix of distance weights between ST and SC.

adjcentScOfSP <- function(obj) {
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
}

#' @title getGsetScore

#' @description Get singature scores in each cell or spot using a list of markers.
#' @param obj.seu Seurat object.
#' @param sc.markers A list of markers for cell-types.
#' @param assay The assay used to calculate signature scores of cell types. Default: SCT.
#' @return A data.frame of signatrue scores.
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
}

#' @title estPvalsOfSpot

#' @description Infer whether signature score in spot significant higher than background.
#' @param sp.obj Seurat object of spatial transcriptome (ST) data.
#' @param sp.score Signature score of cell types.
#' @param knn K nearest cells used. Default: 5.
#' @return A data.frame of p-values.

estPvalsOfSpot <- function(sp.obj, sp.score, knn = 5) {
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
}

#' @title selectScByCor

#' @description Select single-cells based on weighted correlation.
#' @param sp.score Score matrix of ST inferred from a list of markers.
#' @param sc.score Score matrix of SC inferred from a list of markers.
#' @param adj.w A weighted matrix infered from distance. Rows represent spot names, and cloumns are single-cells.
#' @return Similarity matrix.

selectScByCor <- function(sp.score, sc.score, adj.w) {
  adj.w <- adj.w[, rownames(sc.score)]
  cor.mat <- cor(t(sp.score), t(sc.score))
  out.sim <-
    {
      cor.mat * adj.w
    } %>% t()
  return(out.sim)
}

#' @title selectByProb

#' @description Measuring the similarity distances using probability model.
#' @param sp.score
#' @param sc.score
#' @param adj.w
#' @return Similarity matrix.

selectByProb <- function(sp.score, sc.score, adj.w) {
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
}

#' @title assignSc2SP

#' @description Assign single-cells to spots using linear assignment method.
#' @param sp.obj Seurat object of ST data.
#' @param sc.obj Seurat object of SC data.
#' @param out.sim Similarity matrix between single-cells and spots.
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
}
