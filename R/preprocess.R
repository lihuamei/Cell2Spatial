#' @title checkParams.runMap2SP
#'
#' @description Check the legality of input arguments.
#'
#' @param argg A list of input arguments.
#' @return A boolean value: TRUE.

checkParams.runMap2SP <- function(argg) {
    stopifnot(
        is(argg[["sp.obj"]], "Seurat"),
        is(argg[["sc.obj"]], "Seurat"),
        is.null(argg[["cell.type.column"]]) || is.character(argg[["cell.type.column"]]),
        if (!is.null(argg[["cell.type.markers"]])) is(argg[["cell.type.markers"]], "list") else TRUE,
        is.numeric(argg[["group.size"]]) && argg[["group.size"]] >= 5,
        is.numeric(argg[["knn.spots"]]) && argg[["knn.spots"]] >= 0,
        if (!is.null(argg[["max.cells.in.spot"]])) is.numeric(argg[["max.cells.in.spot"]]) && argg[["max.cells.in.spot"]] > 0 else TRUE,
        is.logical(argg[["fix.cells.in.spot"]]),
        is.character(argg[["signature.scoring.method"]]),
        is.numeric(argg[["hotspot.detection.threshold"]]) && argg[["hotspot.detection.threshold"]] <= 1 && argg[["hotspot.detection.threshold"]] > 0,
        is.character(argg[["feature.based"]]),
        is.character(argg[["dist.based"]]),
        is.logical(argg[["spatial.weight"]]),
        is.character(argg[["reduction"]]),
        is.character(argg[["platform.res"]]),
        is.logical(argg[["partition"]]),
        is.numeric(argg[["n.workers"]]) && argg[["n.workers"]] >= 1,
        is.logical(argg[["verbose"]])
    )
    return(TRUE)
}

#' @title preprocessSeqData
#'
#' @description Preprocesses single-cell (SC) and spatial transcriptomics (ST) data by down-sampling, filtering low-quality cell types, and normalizing both datasets.
#'
#' @param sp.obj Seurat object containing spatial transcriptomics data.
#' @param sc.obj Seurat object containing single-cell RNA-seq data.
#' @param cell.type.column Column name in the meta.data slot of SC Seurat object specifying cell types.
#' @param normalize.method Normalization method for SC and ST data, "SCTransform" or "LogNormalize".
#' @param min.cells.of.subset Minimum number of cells required per cell type in the SC data. Default: 5.
#' @param verbose Boolean indicating whether to show progress messages. Default: TRUE.
#' @return A list containing preprocessed SC and ST Seurat objects.

preprocessSeqData <- function(sp.obj, sc.obj, cell.type.column, normalize.method, min.cells.of.subset = 5, verbose = TRUE) {
    if (!is.null(cell.type.column)) Idents(sc.obj) <- sc.obj@meta.data[, cell.type.column]
    levels(sc.obj) <- intersect(levels(sc.obj), unique(Idents(sc.obj) %>% as.vector()))
    cell.cnts <- table(Idents(sc.obj))
    if (sum(cell.cnts < min.cells.of.subset)) {
        sc.obj <- subset(sc.obj, idents = cell.cnts[cell.cnts >= min.cells.of.subset] %>% names())
    }
    sp.obj <- sp.obj[, intersect(colnames(sp.obj), rownames(GetTissueCoordinates(sp.obj)))]
    sp.obj@images[[1]]@coordinates <- sp.obj@images[[1]]@coordinates[colnames(sp.obj), ]

    ovp.genes <- intersect(rownames(sc.obj), rownames(sp.obj))
    norm.func <- ifelse(normalize.method == "SCTransform", SCTransform, NormalizeData)
    sc.obj <- sc.obj[ovp.genes, ] %>% norm.func(., verbose = verbose)
    sp.obj <- sp.obj[ovp.genes, ] %>% norm.func(., verbose = verbose, assay = "Spatial")
    ovp.genes <- intersect(rownames(sc.obj), rownames(sp.obj))
    DefaultAssay(sc.obj) <- ifelse(normalize.method == "SCTransform", "SCT", "RNA")
    DefaultAssay(sp.obj) <- ifelse(normalize.method == "SCTransform", "SCT", "Spatial")
    return(list(SC = sc.obj[ovp.genes, ], ST = sp.obj[ovp.genes, ]))
}

#' @title selectMakers
#'
#' @description Select marker genes for each cell type based on single-cell (SC) data.
#'
#' @param sc.obj Seurat object containing SC data.
#' @param cell.type.markers A list of markers for each cell type. Must be provided if SC data contains only one cell type. Default: NULL.
#' @param group.size The size of the subsets used to find markers. Default: 100.
#' @param assay.type Assay type to use for marker finding. Default: RNA.
#' @param verbose Boolean indicating whether to show running messages. Default: TRUE.
#' @return A list of marker genes for each cell type.

selectMakers <- function(sc.obj, cell.type.markers, group.size, assay.type) {
    if (length(levels(sc.obj)) == 1 && is.null(cell.type.markers)) {
        println("When the SC data contains only one cell type, the marker genes for that cell type must be provided", status = "ERROR")
    }
    if (is.null(cell.type.markers)) {
        sc.markers <- findScMarkers(sc.obj, group.size)
    } else {
        bl.status <- identical(names(cell.type.markers) %>% sort(), levels(sc.obj) %>% sort())
        if (!bl.status) {
            println("Specified marker list of cell types does not match the identifiers of SC data", status = "WARN")
            sc.markers <- findScMarkers(sc.obj, group.size, assay = assay.type)
        } else {
            sc.markers <- cell.type.markers
        }
    }
    if (inherits(sc.markers, "character")) sc.markers <- list(XX = sc.markers) %>% `names<-`(levels(sc.obj))
    return(sc.markers)
}
