#' @title estCellPerSpots

#' @description Estimate number of cells per spot.
#' @param sp.obj Seurat object of spatial transcriptome (ST) data.
#' @param max.cells Maximum number of cells in Spots. Default: 10 (10x Genomics).
<<<<<<< HEAD
#' @param fix.cells.in.spot Fixed number of cells assigned to a spot. Default: NULL.
#' @param quantile.cut Spots with a number of UMIs greater than the quantile cutoff are set as a reference. Default: 0.99.
#' @param num.genes Set the number of stable genes. Default: 2000.
#' @return A vector of estimated cells per spot.
#' @export estCellPerSpots

estCellPerSpots <- function(sp.obj, max.cells = 10, fix.cells.in.spot = NULL, quantile.cut = 0.99, num.genes = 2000) {
    if (!is.null(fix.cells.in.spot)) {
        pred.cnt <- rep(fix.cells.in.spot, ncol(sp.obj)) %>% `names<-`(colnames(sp.obj))
    } else {
        count.mat <- GetAssayData(sp.obj, slot = "count", assay = "Spatial") %>% as.data.frame()
        count.mat <- count.mat[!grepl("^MT|^RP[SL]", rownames(count.mat)), ]
        rank.mat <- apply(count.mat, 2, rank)
        rank.sd <- apply(rank.mat, 1, sd)
        genes <- rank.sd[order(rank.sd)][1:num.genes] %>% names()
        count.mat <- count.mat[genes, ]
        cnt.sum <- colSums(count.mat)
        ref.cut <- quantile(cnt.sum, quantile.cut)
        ref.df <- rowMeans(count.mat[, which(cnt.sum >= ref.cut), drop = FALSE])
        pred.cnt <- apply(count.mat, 2, function(gg) ceiling(coef(lm(gg ~ ref.df - 1)) * max.cells))
        pred.cnt[pred.cnt <= 0] <- 1
        pred.cnt[pred.cnt > max.cells] <- max.cells
    }
    return(pred.cnt)
=======
#' @return A vector of estimated cells per spot.
#' @export estCellPerSpots

estCellPerSpots <- function(sp.obj, max.cells = 10, quantile.cut = 0.99, num.genes = 2000) {
  count.mat <- GetAssayData(sp.obj, slot = "count", assay = "Spatial") %>% as.data.frame()
  count.mat <- count.mat[!grepl("^MT|^RP[SL]", rownames(count.mat)), ]
  rank.mat <- apply(count.mat, 2, rank)
  rank.sd <- apply(rank.mat, 1, sd)
  genes <- rank.sd[order(rank.sd)][1:num.genes] %>% names()
  count.mat <- count.mat[genes, ]
  cnt.sum <- colSums(count.mat)
  ref.cut <- quantile(cnt.sum, quantile.cut)
  ref.df <- rowMeans(count.mat[, which(cnt.sum >= ref.cut), drop = FALSE])
  pred.cnt <- apply(count.mat, 2, function(gg) ceiling(coef(lm(gg ~ ref.df - 1)) * max.cells))
  pred.cnt[pred.cnt <= 0] <- 1
  pred.cnt[pred.cnt > max.cells] <- max.cells
  return(pred.cnt)
>>>>>>> 5859e3527b186cec1e7cc9974128adabbcd95a88
}

#' @title adjustScObj

#' @param sc.obj Seurat object of single-cell (SC) data.
#' @param st.prop A data.frame containing proportions for ST data.
#' @param num.cells A vector of number of cells for ST data.
#' @return Adjusted Seurat object of SC data.

adjustScObj <- function(sc.obj, st.prop, num.cells) {
<<<<<<< HEAD
    prop.vec <- colSums(st.prop) / sum(colSums(st.prop))
    ctype.cnt <-
        {
            sum(num.cells) * prop.vec
        } %>% ceiling()
    sc.idents <- Idents(sc.obj)
    ctype.cnt[which.max(ctype.cnt)] <- ctype.cnt[which.max(ctype.cnt)] + sum(num.cells) - sum(ctype.cnt)

    sample.cells <- future.apply::future_lapply(names(ctype.cnt), function(xx) {
        cells.tmp <- sc.idents[sc.idents == xx]
        rep.bool <- ifelse(length(cells.tmp) > ctype.cnt[xx], FALSE, TRUE)
        cells.sample <- sample(cells.tmp, ctype.cnt[xx], replace = rep.bool) %>% names()
    }, future.seed = TRUE) %>%
        unlist() %>%
        as.vector()

    uniq.names <- make.unique(sample.cells)
    cell.types <- sc.idents[sample.cells] %>% `names<-`(uniq.names)
    expr.mat <-
        {
            GetAssayData(sc.obj, assay = "SCT")[, sample.cells]
        } %>% `colnames<-`(uniq.names)
    meta.data <-
        {
            sc.obj@meta.data[sample.cells, ]
        } %>% `rownames<-`(uniq.names)
    sc.obj.new <- CreateSeuratObject(count = expr.mat, meta.data = meta.data, assay = "RNA")
    sc.obj.new@meta.data$RawName <- sample.cells
    Idents(sc.obj.new) <- cell.types
    return(sc.obj.new)
=======
  prop.vec <- colSums(st.prop) / sum(colSums(st.prop))
  ctype.cnt <-
    {
      sum(num.cells) * prop.vec
    } %>% ceiling()
  sc.idents <- Idents(sc.obj)
  ctype.cnt[which.max(ctype.cnt)] <- ctype.cnt[which.max(ctype.cnt)] + sum(num.cells) - sum(ctype.cnt)

  sample.cells <- future.apply::future_lapply(names(ctype.cnt), function(xx) {
    cells.tmp <- sc.idents[sc.idents == xx]
    rep.bool <- ifelse(length(cells.tmp) > ctype.cnt[xx], FALSE, TRUE)
    cells.sample <- sample(cells.tmp, ctype.cnt[xx], replace = rep.bool) %>% names()
  }, future.seed = TRUE) %>%
    unlist() %>%
    as.vector()

  uniq.names <- make.unique(sample.cells)
  cell.types <- sc.idents[sample.cells] %>% `names<-`(uniq.names)
  expr.mat <-
    {
      GetAssayData(sc.obj, assay = "SCT")[, sample.cells]
    } %>% `colnames<-`(uniq.names)
  meta.data <-
    {
      sc.obj@meta.data[sample.cells, ]
    } %>% `rownames<-`(uniq.names)
  sc.obj.new <- CreateSeuratObject(count = expr.mat, meta.data = meta.data)
  sc.obj.new@meta.data$RawName <- sample.cells
  Idents(sc.obj.new) <- cell.types
  return(sc.obj.new)
>>>>>>> 5859e3527b186cec1e7cc9974128adabbcd95a88
}
