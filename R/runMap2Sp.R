#' @title runMap2SP

#' @description Accurate mapping of single cells to spatial coordinates.
#' @param sp.obj Seurat object of spatial transcriptome (ST) data.
#' @param sc.obj Seurat object of single-cell (SC) data.
#' @param sc.markers When the number of cell types identified from SC data is one, markers specific to that cell type must be provided. Default: NULL.
#' @param ctype Specify the column name for the cell type in the meta.data slot of the SC Seurat object. Default: idents.
#' @param group.size Specify the marker size for each subset derived from single-cell data. Default: 30.
#' @param res Resolution for clustering ST spots. Default: 0.8.
#' @param partion Split into sub-modules when mapping SC to ST spots, with 'partion' set to TRUE. Default: TRUE.
#' @param max.cells.in.spot Maximum number of cells in ST spots. Default: 10 (for 10x Visium). For high-resolution ST data (such as Image-based ST technology), set `max.cells.in.spot` to 1.
#' @param fix.cells.in.spot Fixed number of cells assigned to a spot. Default: NULL.
#' @param knn Utilize the k nearest spots to filter out spots that do not contain cell types present in the SC reference. Default: 5.
#' @param sample.size Down-sample a small subset of SC data for mapping to ST coordinates. Default: NULL.
#' @param detect.hotspot.method Detect hot spots based on the signature scores of a specific cell type in ST spots using the Getis-Ord [getis.ord] method or t-test [t.test] framework. Default: getis.ord.
#' @param p.value.threshold P-value threshold for determining hot spots of cell types. Default: 0.3.
#' @param quantile.threshold Quantile threshold for selecting single cells as a reference to estimate cell counts in spots. Default: 0.99.
#' @param select.markers Method for selecting cell-type specific markers when `sc.markers = NULL`: modified Shannon-entropy strategy (shannon) or Wilcoxon test (wilcox) implemented by the Seurat package. Default: shannon.
#' @param dist.method Measure the distance between single cells and spots, using maximum likelihood model (mle) or correlation (cor). Default: mle.
#' @param return.type Assigned results can be of the object type Seurat or SingleCellExperiment. Default: Seurat.
#' @param n.workers Number of cores for parallel processing. Default: 4.
#' @param verbose Show running messages or not. Default: TRUE.
#' @return Assigned results wrapped as a Seurat or SingleCellExperiment object.
#' @export runMap2SP
#'
#' @examples
#' sp.obj <- system.file("data", "Kindney_SP.RDS", package = "Cell2Spatial") %>% readRDS(.)
#' sc.obj <- system.file("data", "Kindney_SC.RDS", package = "Cell2Spatial") %>% readRDS(.)
#' sce <- runMap2SP(sp.obj, sc.obj, ctype = "mainCtype", res = 0.8, group.size = 30)
runMap2SP <- function(sp.obj,
                      sc.obj,
                      ctype = "idents",
                      sc.markers = NULL,
                      group.size = 30,
                      res = 0.8,
                      partion = TRUE,
                      max.cells.in.spot = 10,
                      fix.cells.in.spot = NULL,
                      knn = 5,
                      sample.size = NULL,
                      detect.hotspot.method = c("getis.ord", "t.test"),
                      p.value.threshold = 0.1,
                      quantile.threshold = 0.99,
                      select.markers = c("shannon", "wilcox"),
                      dist.method = c("mle", "cor"),
                      return.type = c("Seurat", "SingleCellExperiment"),
                      n.workers = 4,
                      verbose = TRUE) {
    options(warn = -1)
    options(future.rng.onMisuse = "ignore")
    bl.status <- as.list(environment()) %>% checkParams.runMap2SP(.)
    multipleProcess(n.workers = n.workers)
    println(sprintf("%g of single-cells and %g spots detected", ncol(sc.obj), ncol(sp.obj)), verbose = verbose)

    println(sprintf("Preprocessing and normalizing data using the SCTransform method..."), verbose = verbose)
    obj.lst <- preprocessSeqData(sp.obj, sc.obj, sample.size, ctype, verbose = FALSE)

    if (length(levels(obj.lst$SC)) > 1) {
        println("Finding specific markers across cell types based on SC data", verbose = verbose)
    }
    sc.markers <- selectMakers(obj.lst$SC, sc.markers, match.arg(select.markers), group.size, verbose = verbose)

    println("Detecting hotspot regions for ST data", verbose = verbose)
    sp.score <- getGsetScore(obj.lst$ST, sc.markers, assay = "SCT")
    hot.spts <- switch(match.arg(detect.hotspot.method),
        getis.ord = {
            hot.spts <- dectHotSpotsByGetisOrdGi(obj.lst$ST, sp.score, knn, p.value.threshold)
        },
        t.test = {
            hot.spts <- dectHotSpotsByTtest(obj.lst$ST, sp.score, knn, p.value.threshold)
        }
    )
    keep.spots <- rownames(hot.spts)[rowSums(hot.spts) > 0]
    println("Clustering spots for ST data and estimating the cell counts per spot", verbose = verbose)
    sp.obj <- findClustersForSpData(obj.lst$ST, res = res, verbose = FALSE)
    num.cells <- estCellPerSpots(sp.obj, max.cells.in.spot, fix.cells.in.spot, quantile.cut = quantile.threshold)

    println("Weighting the distance between SC and ST data...", verbose = verbose)
    lamba <- median(num.cells[keep.spots])
    adj.w <- weightDist(obj.lst$SC, sp.obj, lamba, mc.cores = n.workers)

    sc.obj <- obj.lst$SC
    sp.obj <- subset(sp.obj, cells = keep.spots)
    hot.spts <- hot.spts[keep.spots, ]
    num.cells <- num.cells[keep.spots]
    println("Estimating cellular proportions in each spot and adjusting SC data for mapping", verbose = verbose)
    st.prop <- estPropInSpots(sp.obj, sc.obj, hot.spts, sc.markers, intercept = TRUE)
    sc.obj <- adjustScObj(sc.obj, st.prop, num.cells)
    println("Similarity estimation for single-cells and spots", verbose = verbose)
    sp.score <- sp.score[keep.spots, ]
    sc.score <- getGsetScore(sc.obj, sc.markers, assay = "RNA")
    out.sim <- switch(match.arg(dist.method),
        mle = {
            out.sim <- selectByProb(sp.score, sc.score)
        },
        cor = {
            out.sim <- cor(t(sp.score), t(sc.score))
        }
    )
    println(sprintf("Assigning %g single-cells to spots and generating spatial coordinates", sum(num.cells)))
    out.sc <- linearSumAssignment(sp.obj, sc.obj, out.sim, adj.w, num.cells, hot.spts, partion)
    sce <- assignSCcords(sp.obj, sc.obj, out.sc, lapply(out.sc, length), match.arg(return.type), n.workers = n.workers)
    println("Finished!", verbose = verbose)
    return(sce)
}
