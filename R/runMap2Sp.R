#' @title runCell2Spatial
#'
#' @description Accurate mapping of single cells to spatial coordinates.
#' @param sp.obj Seurat object of spatial transcriptome (ST) data.
#' @param sc.obj Seurat object of single-cell (SC) data.
#' @param cell.type.column Specify the column name for the cell type in the meta.data slot of the SC Seurat object. Default: idents.
#' @param sample.size Down-sample a small subset of SC data for mapping to ST coordinates. Default: NULL.
#' @param sc.markers When the number of cell types identified from SC data is one, markers specific to that cell type must be provided. Default: NULL.
#' @param marker.selection.method Method for selecting cell-type specific markers when `sc.markers = NULL`: modified Shannon-entropy strategy (shannon) or Wilcoxon test (wilcox) implemented by the Seurat package. Default: shannon.
#' @param group.size Specify the marker size for each subset derived from single-cell data. Default: 30.
#' @param resolution Resolution for clustering ST spots. Default: 0.8.
#' @param knn Utilize the k nearest spots to filter out spots that do not contain cell types present in the SC reference. Default: 5.
#' @param partition Split into sub-modules when mapping SC to ST spots, with 'partion' set to TRUE. Default: FALSE.
#' @param max.cells.in.spot Maximum number of cells in ST spots. Default: 10 (for 10x Visium). For high-resolution ST data (such as Image-based ST technology), set `max.cells.in.spot` to 1.
#' @param fix.cells.in.spot Fixed number of cells assigned to a spot. Default: NULL.
#' @param signature.scoring.method Method for scoring the signature of cell types in ST data: AddModuleScore, UCell, or AverageExpr. Default: AddModuleScore.
#' @param hotspot.detection.method Detect hot spots based on the signature scores of a specific cell type in ST spots using the Getis-Ord [getis.ord] method or t-test [t.test] framework. Default: getis.ord.
#' @param hotspot.detection.threshold P-value threshold for determining hot spots of cell types, range from 0 to 1. Default: 1.
#' @param integ.entire.dataset Estimate the distance weight between SC and ST using entire or pseudo and downsampled data. Default: FALSE.
#' @param feature.based Specify whether features for likelihood or correlation calculations between single cells and spots are based on gene expression ('gene.based') or signature scores of cell types ('celltype.based'). Default: 'gene.based'.
#' @param dist.method Measure the distance between single cells and spots, using maximum likelihood model (mle) or correlation (cor). Default: mle.
#' @param hclust Logical. If TRUE, hierarchical clustering is applied to group ST clusters. Default: TRUE.
#' @param output.type Assigned results can be of the object type Seurat or SingleCellExperiment. Default: Seurat.
#' @param n.workers Number of cores for parallel processing. Default: 4.
#' @param verbose Show running messages or not. Default: TRUE.
#' @return Assigned results wrapped as a Seurat or SingleCellExperiment object.
#' @export runCell2Spatial
#'
#' @examples
#' sp.obj <- system.file("data", "Kindney_SP.RDS", package = "Cell2Spatial") %>% readRDS(.)
#' sc.obj <- system.file("data", "Kindney_SC.RDS", package = "Cell2Spatial") %>% readRDS(.)
#' sce <- runCell2Spatial(sp.obj, sc.obj, cell.type.column = "mainCtype", res = 0.8, group.size = 100)
runCell2Spatial <- function(sp.obj,
                            sc.obj,
                            cell.type.column = "idents",
                            sample.size = NULL,
                            sc.markers = NULL,
                            marker.selection.method = c("shannon", "wilcox"),
                            group.size = 100,
                            resolution = 0.8,
                            knn = 5,
                            partition = FALSE,
                            max.cells.in.spot = 10,
                            fix.cells.in.spot = NULL,
                            signature.scoring.method = c("AddModuleScore", "UCell", "AverageExpr"),
                            hotspot.detection.method = c("getis.ord", "t.test"),
                            hotspot.detection.threshold = 1,
                            integ.entire.dataset = FALSE,
                            feature.based = c("gene.based", "celltype.based"),
                            dist.method = c("mle", "cor"),
                            hclust = TRUE,
                            output.type = c("Seurat", "SingleCellExperiment"),
                            n.workers = 4,
                            verbose = TRUE) {
    options(warn = -1)
    options(future.rng.onMisuse = "ignore")
    bl.status <- as.list(environment()) %>% checkParams.runMap2SP(.)
    multipleProcess(n.workers = n.workers)
    println(sprintf("%g of single-cells and %g spots detected", ncol(sc.obj), ncol(sp.obj)), verbose = verbose)

    println(sprintf("Preprocessing and normalizing data using the SCTransform method..."), verbose = verbose)
    obj.lst <- preprocessSeqData(sp.obj, sc.obj, sample.size, cell.type.column, verbose = FALSE)
    sp.obj <- obj.lst$ST
    sc.obj <- obj.lst$SC
    garbageCollection(obj.lst)

    if (length(levels(sc.obj)) > 1) {
        println("Finding specific markers across cell types based on SC data", verbose = verbose)
    }
    sc.markers <- selectMakers(sc.obj, sc.markers, match.arg(marker.selection.method), group.size, verbose = verbose)

    println("Detecting hotspot regions for ST data", verbose = verbose)
    sp.score <- getGsetScore(sp.obj, sc.markers, assay = "SCT", match.arg(signature.scoring.method))
    hot.spts <- switch(match.arg(hotspot.detection.method),
        getis.ord = {
            hot.spts.lst <- dectHotSpotsByGetisOrdGi(sp.obj, sp.score, knn, hotspot.detection.threshold)
        },
        t.test = {
            hot.spts.lst <- dectHotSpotsByTtest(sp.obj, sp.score, knn, hotspot.detection.threshold)
        }
    )
    hot.spts <- hot.spts.lst$x
    keep.spots <- rownames(hot.spts)[rowSums(hot.spts) > 0]
    println("Clustering spots for ST data and estimating the cell counts per spot", verbose = verbose)
    suppressWarnings({
        sp.obj <- findClustersForSpData(sp.obj, res = resolution, verbose = FALSE)
    })
    num.cells <- estCellPerSpots(sp.obj, max.cells.in.spot, fix.cells.in.spot)

    println("Weighting the distance between SC and ST data...", verbose = verbose)
    lamba <- median(num.cells[keep.spots])
    suppressWarnings({
        adj.w <- weightDist(sc.obj, sp.obj, lamba, mc.cores = n.workers, use.entire = integ.entire.dataset)
    })
    sp.obj <- subset(sp.obj, cells = keep.spots)
    hot.spts <- hot.spts[keep.spots, ,drop = FALSE]
    num.cells <- num.cells[keep.spots]
    println("Estimating cellular proportions in each spot and adjusting SC data for mapping", verbose = verbose)
    hot.pvals <- hot.spts.lst$y[keep.spots, , drop = FALSE]
    if (length(sc.markers) > 1) {
        st.prop <- estPropInSpots(sp.obj, sc.obj, hot.pvals, sc.markers, intercept = TRUE)
    } else {
        st.prop <- data.frame(rep(1, ncol(sp.obj))) %>%
            `colnames<-`(names(sc.markers)) %>%
            `rownames<-`(colnames(sp.obj))
    }
    sc.obj <- adjustScObj(sc.obj, st.prop, num.cells)
    garbageCollection(st.prop)

    println("Similarity estimation for single-cells and spots", verbose = verbose)
    feature.based <- ifelse(length(sc.obj) == 1, 'gene.based', match.arg(feature.based))
    if (feature.based == "gene.based") {
        sp.score <- GetAssayData(sp.obj) %>%
            .[unique(unlist(sc.markers)), ] %>%
            as.matrix() %>%
            t()
        sc.score <- GetAssayData(sc.obj) %>%
            .[unique(unlist(sc.markers)), ] %>%
            as.matrix() %>%
            t()
    } else {
        sc.score <- getGsetScore(sc.obj, sc.markers, assay = "RNA", match.arg(signature.scoring.method))
    }
    sp.score <- sp.score[keep.spots, ]
    out.sim <- switch(match.arg(dist.method),
        mle = {
            out.sim <- selectByProb(sp.score, sc.score)
        },
        cor = {
            out.sim <- cor(t(sp.score), t(sc.score))
        }
    )
    garbageCollection(sp.score, sc.score)

    println(sprintf("Assigning %g single-cells to spots and generating spatial coordinates", sum(num.cells)))
    suppressWarnings({
        out.sc <- linearSumAssignment(sp.obj, sc.obj, out.sim, adj.w, num.cells, hot.pvals, partition, hclust)
    })
    sce <- assignSCcords(sp.obj, sc.obj, out.sc, lapply(out.sc, length), match.arg(output.type), n.workers = n.workers)
    garbageCollection(sp.obj, sc.obj, out.sim, adj.w, hot.spts, num.cells, out.sc)
    println("Finished!", verbose = verbose)
    return(sce)
}
