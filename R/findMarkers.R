#' @title findScMarkersByShannon
#'
#' @description Find markers among cell types based on SC data using Shannon-entropy based strategy.
#'
#' @param sc.obj Seurat object of single-cell data.
#' @param group.size Marker size of each subset derived from SC data.
#' @param assay Assay type to use for SC data normalization. Default: RNA.
#' @param min.cells.of.subset Include cells detected in at least one cell type in the SC data. Default: 5.
#' @param percent.cut Minimum fraction of cells in a cluster expressing a gene. Default: 0.25.
#' @return A list of identified markers.
#' @export findScMarkersByShannon

findScMarkersByShannon <- function(sc.obj, group.size, assay = "RNA", min.cells.of.subset = 5, percent.cut = 0.25) {
    levels(sc.obj) <- intersect(levels(sc.obj), unique(Idents(sc.obj) %>% as.vector()))
    cell.cnts <- table(Idents(sc.obj))
    if (sum(cell.cnts < min.cells.of.subset)) {
        sc.obj <- subset(sc.obj, idents = cell.cnts[cell.cnts >= min.cells.of.subset] %>% names())
    }
    X <-
        {
            AverageExpression(sc.obj, assay = assay)[[assay]] + 1e-2
        } %>% as.matrix()
    nrep.sns <- dim(X)[2]
    X.norm <- X / rowSums(X)
    X.mu <- rowMeans(X.norm)
    w.genes <- MatrixGenerics::rowMaxs(X.norm) %>%
        {
            (. - min(.)) / (max(.) - min(.))
        }
    shannon.dist <- 1 / nrep.sns * rowSums(X.norm / X.mu * log2(X.norm / X.mu))
    group.ctype <- apply(X, 1, which.max) %>%
        colnames(X)[.] %>%
        as.vector()

    cv.filter <- sapply(levels(sc.obj), function(ct) {
        ct.cells <- which(Idents(sc.obj) == ct)
        ct.expr <- GetAssayData(sc.obj, assay = assay, slot = "data")[, ct.cells]
        gene.means <- MatrixGenerics::rowMeans(ct.expr)
        gene.sds <- MatrixGenerics::rowSds(ct.expr)
        cv <- ifelse(gene.means > 0, gene.sds / gene.means, 0)
        all.genes <- rownames(X)
        cv.all <- numeric(length(all.genes))
        names(cv.all) <- all.genes
        cv.all[names(cv)] <- cv
        cv.all
    }, simplify = TRUE)

    score.df <- data.frame(
        cluster = group.ctype,
        Score = shannon.dist * w.genes,
        Shannon = shannon.dist,
        w = w.genes,
        gene = rownames(X),
        stringsAsFactors = FALSE
    )
    cv.bool.df <- future.apply::future_lapply(levels(sc.obj), function(ct) {
        tmp.vec <- rep(FALSE, nrow(score.df))
        threash.cut <- mean(cv.filter[, ct]) + 3 * sd(cv.filter[, ct])
        genes.rm <- which(cv.filter[score.df$gene, ct] < threash.cut)
        tmp.vec[genes.rm] <- TRUE
        return(tmp.vec)
    }, future.seed = TRUE) %>%
        do.call(cbind.data.frame, .) %>%
        `rownames<-`(score.df$gene) %>%
        `colnames<-`(levels(sc.obj))

    keep.ct <- intersect(levels(sc.obj), unique(group.ctype))
    score.df <- score.df[!grepl("^MT-|^RP[L|S]", toupper(score.df$gene)), ]
    q.index <- qualityGenesIndex(sc.obj, assay = assay)
    sc.markers <- future.apply::future_lapply(keep.ct, function(ct) {
        score.df.tmp <- cbind.data.frame(
            score.df,
            logFC = q.index$log2FC[rownames(score.df), ct],
            percent = q.index$Percent.In[rownames(score.df), ct],
            cv = cv.bool.df[rownames(score.df), ct],
            mu.in = q.index$Mu.In[rownames(score.df), ct]
        ) %>% subset(percent >= percent.cut & cv == TRUE & mu.in >= quantile(.$mu.in, 0.25) & cluster == ct)
        logfc.cut <- quantile(score.df.tmp$logFC, 1 - pmin(group.size * 2 / nrow(score.df.tmp), 1))
        score.df.sub <- subset(score.df.tmp, logFC >= logfc.cut)
        if (length(score.df.sub) == 0) {
            return(NULL)
        }
        score.df.sub %>%
            arrange(desc(Score)) %>%
            {
                .$gene[1:min(nrow(.), group.size)]
            }
    }, future.seed = TRUE) %>%
        `names<-`(keep.ct) %>%
        .[lapply(., length) > 0]
    return(sc.markers)
}

#' @title qualityGenesIndex
#'
#' @description Calculation of quality indexes for genes among cell types.
#'
#' @param sc.obj Seurat object of single-cell data.
#' @param assay Assay type to use for SC data normalization. Default: RNA.
#' @return A list of indexes for genes among cell types.

qualityGenesIndex <- function(sc.obj, assay = "RNA") {
    expr <- GetAssayData(sc.obj, assay = assay, slot = "data")
    ids <- Idents(sc.obj)
    cts <- levels(sc.obj)
    genes <- rownames(expr)
    all.cells <- colnames(expr)
    n.all <- length(all.cells)

    ct.cells <- lapply(cts, function(ct) which(ids == ct)) %>% `names<-`(cts)
    n.ct.cells <- vapply(ct.cells, length, 1L)

    mu.in <- AverageExpression(sc.obj, assay = assay)[[assay]] + 1e-2
    sum.in.list <- lapply(cts, function(ct) mu.in[, ct] * n.ct.cells[[ct]])
    total.sum <- do.call(cbind, sum.in.list) %>% rowSums()
    mu.out <- sapply(seq_along(cts), function(j) {
        denom <- (n.all - n.ct.cells[j])
        if (denom <= 0) {
            return(rep(0, length(genes)))
        }
        (total.sum - sum.in.list[[j]]) / denom
    }) %>%
        `colnames<-`(cts) %>%
        `rownames<-`(genes)

    is.expr <- expr > 0
    pct.in <- sapply(cts, function(ct) {
        if (n.ct.cells[[ct]] == 0) {
            return(rep(0, nrow(expr)))
        }
        Matrix::rowMeans(is.expr[, ct.cells[[ct]], drop = FALSE])
    }) %>%
        `colnames<-`(cts) %>%
        `rownames<-`(genes)

    logFC <- log2((mu.in + 1e-8) / (mu.out + 1e-8))
    return(list(log2FC = logFC, Percent.In = pct.in, Mu.In = mu.in))
}

#' @title findMarkersBySeurat
#'
#' @description Find markers among cell types based on SC data using Seurat-inbuilt method.
#'
#' @param sc.obj Seurat object of single-cell data.
#' @param group.size Marker size of each subset derived from SC data.
#' @param assay Assay type to use for SC data normalization. Default: RNA.
#' @return A list of identified markers.
#' @export findMarkersBySeurat

findMarkersBySeurat <- function(sc.obj, group.size, assay = "RNA", verbose = FALSE) {
    sc.obj.sub <- downSamplSeurat(sc.obj, cnt = 500)
    if (assay == "SCT") sc.obj.sub <- PrepSCTFindMarkers(sc.obj.sub, verbose = verbose)
    sc.markers <- sc.obj.sub %>%
        FindAllMarkers(., assay = assay, only.pos = TRUE, verbose = verbose) %>%
        .[!grepl("^MT-|^RP[L|S]", toupper(.$gene)), ] %>%
        mutate(EScore = {
            .$avg_log2FC * .$pct.1
        }) %>%
        group_by(cluster) %>%
        top_n(wt = EScore, n = group.size) %>%
        {
            split(.$gene, .$cluster)
        } %>%
        .[lapply(., length) > 0]

    keep.ct <- intersect(levels(sc.obj), unique(group.ctype))
    sc.markers <- lapply(keep.ct, function(ct) {
        score.df.sub <- subset(score.df, cluster == ct)
        score.df.sub <- score.df.sub[cv.bool.df[score.df.sub$gene, ct], ]
        if (length(score.df.sub) == 0) {
            return(NULL)
        }
        score.df.sub %>%
            arrange(desc(Score)) %>%
            {
                .$gene[1:min(nrow(.), group.size)]
            }
    }) %>%
        `names<-`(keep.ct) %>%
        .[lapply(., length) > 0]
    return(sc.markers)
}


.findScMarkers <- function(select.markers, ...) {
    switch(select.markers,
        wilcox = {
            sc.markers <- findMarkersBySeurat(...)
        },
        shannon = {
            sc.markers <- findScMarkersByShannon(...)
        }
    )
    return(sc.markers)
}
