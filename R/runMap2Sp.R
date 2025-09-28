#' @title runCell2Spatial
#'
#' @description Accurate mapping of single cells to spatial coordinates.
#'
#' @param sp.obj Seurat object of spatial transcriptome (ST) data.
#' @param sc.obj Seurat object of single-cell (SC) data.
#' @param cell.type.column Specify the column name for the cell type in the meta.data slot of the SC Seurat object. Default: NULL.
#' @param cell.type.markers When the number of cell types identified from SC data is one, markers specific to that cell type must be provided. Default: NULL.
#' @param normalize.method Normalization method for scRNA-seq data. Default: SCTransform.
#' @param marker.selection Method for selecting cell-type specific markers when `sc.markers = NULL`: modified Shannon-entropy strategy (shannon) or Wilcoxon test (wilcox) implemented by the Seurat package, or expression fold changes (logFC). Default: shannon.
#' @param group.size Specify the marker size for each subset derived from single-cell data. Default: 30.
#' @param knn.spots Number of nearest neighbors to consider for spot-based deconvolution analysis. When set to 0, neighbors are not used. Default: 5.
#' @param resolution Specify the resolution for spatial clustering. Default: 0.8.
#' @param max.cells.in.spot Maximum number of cells in ST spots. Default: 10 (for 10x Visium). For high-resolution ST data (such as Image-based ST technology), set `max.cells.in.spot` to 1.
#' @param fix.cells.in.spot Fixed number of cells assigned to a spot or not. Default: FALSE.
#' @param signature.scoring.method Method for scoring the signature of cell types in ST data: AddModuleScore, UCell, or AverageExpr. Default: AddModuleScore.
#' @param hotspot.detection.threshold P-value threshold for determining hot spots of cell types, range from 0 to 1. Default: 1.
#' @param adjust.deconv Specify the strategy for adjusting deconvolution estimation, indicating whether to incorporate hotspot weights into cellular composition estimation; if enabled, a cleaner spatial architecture may be obtained, with regions corresponding to cell types with low confidence potentially removed. Options: "Scales" (rescales hotspot significance within tissue identities), "NONE" (applies no weighting). Default: "Scales".
#' @param feature.based Specify whether features for likelihood or correlation calculations between single cells and spots are based on gene expression ('gene.based') or signature scores of cell types ('celltype.based'). Default: gene.based.
#' @param dist.based Dimensionality reduction basis used for distance weighting, UMAP or TSNE. Default: UMAP.
#' @param spatial.weight Whether to apply spatial weights to the analysis. Default: FALSE.
#' @param reduction Dimensionality reduction technique to align spatial and single-cell data. Default: cca.
#' @param platform.res Resolution of the spatial platform (Low or High). Default: Low.
#' @param partition Split into sub-modules when mapping SC to ST spots, with 'partion' set to TRUE. Default: FALSE.
#' @param n.workers Number of cores for parallel processing. Default: 4.
#' @param verbose Show running messages or not. Default: TRUE.
#' @return Assigned results wrapped as a Seurat object.
#' @export runCell2Spatial
#'
#' @examples
#' sp.obj <- system.file("data", "Kindney_SP.RDS", package = "Cell2Spatial") %>% readRDS(.)
#' sc.obj <- system.file("data", "Kindney_SC.RDS", package = "Cell2Spatial") %>% readRDS(.)
#' sce <- runCell2Spatial(sp.obj, sc.obj, cell.type.column = "mainCtype", group.size = 30, fix.cells.in.spot = TRUE, partition = TRUE)
runCell2Spatial <- function(sp.obj,
                            sc.obj,
                            cell.type.column = NULL,
                            cell.type.markers = NULL,
                            normalize.method = c("SCTransform", "LogNormalize"),
                            marker.selection = c("shannon", "wilcox", "logFC"),
                            group.size = 30,
                            knn.spots = 5,
                            resolution = 0.8,
                            max.cells.in.spot = 10,
                            fix.cells.in.spot = FALSE,
                            signature.scoring.method = c("AddModuleScore", "UCell", "AverageExpr"),
                            hotspot.detection.threshold = 1,
                            adjust.deconv = c("Scales", "NONE"),
                            feature.based = c("gene.based", "celltype.based"),
                            dist.based = c("UMAP", "TSNE"),
                            spatial.weight = FALSE,
                            reduction = c("cca", "rpca"),
                            platform.res = c("Low", "High"),
                            partition = FALSE,
                            n.workers = 4,
                            verbose = TRUE) {
    bl.status <- as.list(environment()) %>% checkParams.runMap2SP(.)
    println(sprintf("%g of single-cells and %g spots detected", ncol(sc.obj), ncol(sp.obj)), verbose = verbose)

    multipleProcess(n.workers = n.workers)
    println(sprintf("Preprocessing and normalizing count data..."), verbose = verbose)
    obj.lst <- preprocessSeqData(sp.obj, sc.obj, cell.type.column, match.arg(normalize.method), verbose = FALSE)
    sp.obj <- obj.lst$ST
    sc.obj <- obj.lst$SC
    garbageCollection(obj.lst)

    if (length(levels(sc.obj)) > 1 & is.null(cell.type.markers)) {
        println("Finding specific markers across cell types based on SC data", verbose = verbose)
    }
    select.markers <- match.arg(marker.selection)
    assay.type <- ifelse(match.arg(normalize.method) == "SCTransform", "SCT", "RNA")
    sc.markers <- selectMakers(sc.obj, cell.type.markers, group.size, assay.type, select.markers)
    sc.obj <- subset(sc.obj, idents = names(sc.markers))

    feature.based <- ifelse(length(levels(sc.obj)) == 1, "gene.based", match.arg(feature.based))
    if (feature.based == "celltype.based" && match.arg(signature.scoring.method) == "AddModuleScore") {
        println("Signature scoring strategy has been reset from AddModuleScore to the UCell method", verbose = TRUE, status = "WARN")
        signature.scoring.method <- "UCell"
    }
    println("Detecting hotspot regions for ST data", verbose = verbose)
    assay.sp <- ifelse(assay.type == "SCT", "SCT", "Spatial")
    sp.score <- getGsetScore(sp.obj, sc.markers, assay = assay.sp, signature.scoring.method)
    hot.spots.lst <- dectHotSpotsByGetisOrdGi(sp.obj, sp.score, p.value.threshold = hotspot.detection.threshold)
    keep.spots <- rownames(hot.spots.lst$x)[rowSums(hot.spots.lst$x) > 0]
    hot.pvals <- hot.spots.lst$y[keep.spots, , drop = FALSE]
    hot.spots <- hot.spots.lst$x[keep.spots, , drop = FALSE]

    if (max.cells.in.spot > 1) {
        println("Estimating the cell counts per spot in low-resolution ST data", verbose = verbose)
    }
    num.cells <- inferCellNumbers(sp.obj[, keep.spots], max.cells = max.cells.in.spot, fix.cells.in.spot = fix.cells.in.spot)

    println("Clustering of spatial spots")
    sp.obj <- .findClustersForSpData(obj.seu = sp.obj, assay = assay.type, res.start = resolution, verbose = FALSE) %>% subset(, cells = keep.spots)

    println("Estimating cellular proportions in each spot and adjusting SC data for mapping", verbose = verbose)
    platform.res <- match.arg(platform.res)
    adjust.deconv <- match.arg(adjust.deconv)
    st.prop.lst <- .estPropInSpots(
        sp.score,
        max.cells.in.spot,
        platform.res,
        partition,
        sp.obj = sp.obj,
        sc.obj = sc.obj,
        sc.markers = sc.markers,
        num.cells = num.cells,
        hot.pvals = hot.pvals,
        weight = adjust.deconv,
        assay = assay.type,
        knn = knn.spots
    )
    sc.obj <- .adjustScObj(sc.obj, st.prop.lst$cnts, assay.type)

    println("Weighting the distance between SC and ST data...", verbose = verbose)
    dist.based <- match.arg(dist.based)
    reduction <- match.arg(reduction)
    adj.w <- .weightDist(
        spatial.weight,
        sc.obj,
        sp.obj,
        lamba = median(num.cells),
        assay.type,
        reduction,
        dist.based = dist.based
    )
    println("Similarity estimation for single cells and spots", verbose = verbose)
    sc.score <- getGsetScore(sc.obj, sc.markers, assay = assay.type, signature.scoring.method)
    sp.score <- sp.score[keep.spots, , drop = FALSE]
    out.sim <- calcSimlarityDist(sc.obj, sp.obj, sc.score, sp.score, feature.based, sc.markers)

    garbageCollection(sp.score, sc.score)
    println(sprintf("Assigning %g single-cells to spots and generating spatial coordinates", sum(num.cells)), verbose = verbose)
    out.sc <- .linearSumAssignment(
        sp.obj,
        sc.obj,
        out.sim,
        adj.w,
        hot.pvals,
        st.prop.lst$cnts,
        num.cells,
        partition,
        assay.type,
        reduction
    )
    sce <- assignSCcords(sp.obj, sc.obj, out.sc, lapply(out.sc, length))
    garbageCollection(sp.obj, sc.obj, out.sim, adj.w, hot.spots, num.cells, out.sc, st.prop.lst, hot.pvals)
    println("Finished!", verbose = verbose)
    return(sce)
}
