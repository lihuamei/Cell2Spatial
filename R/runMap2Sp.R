#' @title runMap2SP

#' @description Accurate mapping of single cells to spatial coordinates.
#' @param sp.obj Seurat object of spatial transcriptome (ST) data.
#' @param sc.obj Seurat object of single-cell (SC) data.
#' @param sc.markers A list of markers for cell types. If not provided, automated inference is performed using the modified Shannon-entropy method. Default: NULL.
#' @param ctype Specify the column name for the cell type in meta.data slot of the SC Seurat object. Default: idents.
#' @param group.size Specify the marker size for each subset derived from single-cell data. This information is crucial for estimating the activity of markers within each cell type. Default: 30.
#' @param res Resolution for clustering ST spot. Default: 0.8.
#' @param integ Integration of SC and ST data or not. For high-resolution ST data, please set FALSE to accelerate the running time. Default: TRUE.
#' @param duplicated Assigning individual cells to ST coordinates with duplicated or not. Default: FALSE.
#' @param partion Split into sub-modules mapped to spatial positions when 'duplicated' set to TRUE. Default: TRUE.
#' @param min.cells.of.subset Include cells detected in at least one cell type in the SC data. Default: 5.
#' @param max.cells.in.spot Maximum number of cells in ST spots. Default: 10 (for 10x Visium).
#' @param fix.cells.in.spot Fixed number of cells assigned to a spot. Default: NULL.
#' @param n.features Top k highly variable features used for integrating SC and ST data. Default: 3000.
#' @param knn Utilizing the k nearest spots to filter out spots that do not contain cell types present in the SC reference. Default: 5.
#' @param sample.size Down-sampling a small subset of SC data for mapping to ST coordinates. Default: NULL.
#' @param select.markers Method for selecting cell-type specific markers when 'sc.markers = NULL': modified Shannon-entropy strategy (shannon) or Wilcoxon test (wilcox) implemented by the Seurat package. Default: shannon.
#' @param dist.method Measure the distance between single-cell and spots, maximum likehood model (mle) or correlation (cor). Default: mle.
#' @param return.type Assigned results can be of the object type Seurat or SingleCellExperiment. Default: Seurat.
#' @param n.wrokers Number of cores for parallel processing. Default: 4.
#' @param verbose Show running messages or not. Default: TRUE.
#' @return Assigning results wrapped as a Seurat or SingleCellExperiment object.
#' @export runMap2SP
#'
#' @examples
#' sp.obj <- system.file("data", "Kindney_SP.RDS", package = "Cell2Spatial") %>% readRDS(.)
#' sc.obj <- system.file("data", "Kindney_SC.RDS", package = "Cell2Spatial") %>% readRDS(.)
#' sce <- runMap2SP(sp.obj, sc.obj, ctype = "mainCtype", res = 0.8, group.size = 30)

runMap2SP <- function(sp.obj,
                      sc.obj,
                      sc.markers = NULL,
                      ctype = "idents",
                      group.size = 30,
                      res = 0.8,
                      integ = TRUE,
                      duplicated = FALSE,
                      partion = TRUE,
                      min.cells.of.subset = 5,
                      max.cells.in.spot = 10,
                      fix.cells.in.spot = NULL,
                      n.features = 3000,
                      knn = 5,
                      sample.size = NULL,
                      select.markers = c("shannon", "wilcox"),
                      dist.method = c("mle", "cor"),
                      return.type = c("Seurat", "SingleCellExperiment"),
                      n.workers = 4,
                      verbose = TRUE) {
    options(warn = -1)
    as.list(environment()) %>% checkParams.runMap2SP(.)
    multipleProcess(n.workers = n.workers)
    println(sprintf("%g of single-cells and %g spots", ncol(sc.obj), ncol(sp.obj)), verbose = verbose)
    sc.obj <- preprocessSCData(sc.obj, sample.size, ctype, min.cells.of.subset)
    sp.obj <- sp.obj[, intersect(colnames(sp.obj), rownames(GetTissueCoordinates(sp.obj)))]
    sp.obj@images[[1]]@coordinates <- sp.obj@images[[1]]@coordinates[colnames(sp.obj), ]

    println("Integrating SC and ST data", verbose = verbose)
    obj.list <- integDataBySeurat(sp.obj, sc.obj, n.features = n.features, integ = integ, verbose = verbose)
    sc.obj <- obj.list$sc.obj
    if (is.null(sc.markers)) {
        println("Finding specific markers across cell types on single-cell data", verbose = verbose)
        sc.markers <- findScMarkers(sc.obj, group.size, select.markers, verbose = verbose)
    }
    println("Clustering spots for ST data and assign the closest clusters to single-cells", verbose = verbose)
    sp.obj <- findClustersForSpData(obj.list$sp.obj, res = res, verbose = verbose)

    println("Calculating signature scores in each cell or spot using a list of markers", verbose = verbose)
    sp.score <- getGsetScore(sp.obj, sc.markers, assay = "SCT")
    st.pvals <- estPvalsOfSpot(sp.obj, sp.score, knn = knn)
    num.cells <- estCellPerSpots(sp.obj, max.cells.in.spot, fix.cells.in.spot, quantile.cut = 0.95)

    if (!duplicated) {
        println("Estimating cellular proportions in each spot", verbose = verbose)
        st.prop <- estPropInSpots(sp.obj, sc.obj, st.pvals, sc.markers, intercept = FALSE)
        sc.obj <- adjustScObj(sc.obj, st.prop, num.cells)
        if (integ) {
            adj.w <- obj.list$adj.w[, sc.obj$RawName] %>% `colnames<-`(colnames(sc.obj))
        } else {
            adj.w <- obj.list$adj.w
        }
    } else {
        adj.w <- obj.list$adj.w
    }
    sc.score <- getGsetScore(sc.obj, sc.markers, assay = "RNA")

    println("Similarity estimation for single-cells and spots", verbose = verbose)
    out.sim <- switch(match.arg(dist.method),
        mle = {
            out.sim <- selectByProb(sp.score, sc.score, adj.w)
        },
        cor = {
            out.sim <- selectScByCor(sp.score, sc.score, adj.w)
        }
    )
    println(sprintf("Assigning %g single-cells to spots...", sum(num.cells)))
    partion <- ifelse(duplicated, FALSE, partion)
    out.sc <- assignSc2SP(sp.obj, sc.obj, out.sim, num.cells, st.pvals, n.features, duplicated, partion, p.cut = 0.05)
    println("Mapping single-cells to spatial coordinates...")
    sce <- assignSCcords(sp.obj, sc.obj, out.sc, lapply(out.sc, length), match.arg(return.type), n.workers = n.workers)
    println("Finished!", verbose = verbose)
    return(sce)
}
