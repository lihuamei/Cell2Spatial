#' @title findScMarkers
#'
#' @description Find markers among cell types based on SC data.
#'
#' @param sc.obj Seurat object of single-cell data.
#' @param group.size Marker size of each subset derived from SC data.
#' @param assay Assay type to use for SC data normalization. Default: RNA.
#' @param min.cells.of.subset Include cells detected in at least one cell type in the SC data. Default: 5.
#' @return A tibble of identified markers.
#' @export findScMarkers

findScMarkers <- function(sc.obj, group.size, assay = "RNA", min.cells.of.subset = 5) {
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
    w.genes <- MatrixGenerics::rowMaxs(X) %>%
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
    sc.markers <- future.apply::future_lapply(keep.ct, function(ct) {
        score.df.sub <- subset(score.df, cluster == ct)
        score.df.sub <- score.df.sub[cv.bool.df[score.df.sub$gene, ct], ]
        if (length(score.df.sub) == 0) {
            return(NULL)
        }
        score.df.sub %>%
            arrange(desc(Score)) %>%
            .[!grepl("^MT-|^RP[L|S]", .$gene), ] %>%
            {
                .$gene[1:min(nrow(.), group.size)]
            }
    }, future.seed = TRUE) %>%
        `names<-`(keep.ct) %>%
        .[lapply(., length) > 0]
    return(sc.markers)
}
