#' @title checkParams.runMap2SP
#'
#' @description Check the legality of input arguments.
#' @param argg A list of input arguments.
#' @return A boolean value: TRUE.

checkParams.runMap2SP <- function(argg) {
    stopifnot(
        is(argg[["sp.obj"]], "Seurat"),
        is(argg[["sc.obj"]], "Seurat"),
        is.character(argg[["cell.type.column"]]),
        if (!is.null(argg[["sc.markers"]])) is(argg[["sc.markers"]], "list") else TRUE,
        is.character(argg[["marker.selection.method"]]),
        is.numeric(argg[["group.size"]]) && argg[["group.size"]] >= 5,
        is.numeric(argg[["resolution"]]) && argg[["resolution"]] > 0,
        is.logical(argg[["partition"]]),
        if (!is.null(argg[["max.cells.in.spot"]])) is.numeric(argg[["max.cells.in.spot"]]) && argg[["max.cells.in.spot"]] > 0 else TRUE,
        if (!is.null(argg[["fix.cells.in.spot"]])) is.numeric(argg[["fix.cells.in.spot"]]) && argg[["fix.cells.in.spot"]] > 0 else TRUE,
        is.character(argg[["signature.scoring.method"]]),
        is.numeric(argg[["hotspot.detection.threshold"]]) && argg[["hotspot.detection.threshold"]] <= 1 && argg[["hotspot.detection.threshold"]] > 0,
        is.logical(argg[["integ.entire.dataset"]]),
        is.character(argg[["feature.based"]]),
        is.character(argg[["dist.method"]]),
        is.character(argg[["dist.based"]]),
        is.numeric(argg[["dist.quantile.cut"]]) && argg[["dist.quantile.cut"]] <= 1 && argg[["dist.quantile.cut"]] > 0,
        is.character(argg[["output.type"]]),
        is.numeric(argg[["n.workers"]]) && argg[["n.workers"]] >= 1,
        is.logical(argg[["verbose"]])
    )
    return(TRUE)
}

#' @title preprocessSeqData
#'
#' @description The preprocessing of single-cell (SC) data includes several steps such as down-sampling, filtering out low-quality cell types, and normalization.
#' @param sp.obj Seurat object of ST data.
#' @param sc.obj Seurat object of SC data.
#' @param cell.type.column Specify the column name for the cell type in the meta.data slot of the SC Seurat object. Default: idents.
#' @param min.cells.of.subset Include cells detected in at least one cell type in the SC data. Default: 5.
#' @param verbose Boolean indicating whether to show running messages. Default: TRUE.
#' @return A list of Seurat objects of preprocessed data.

preprocessSeqData <- function(sp.obj, sc.obj, cell.type.column, min.cells.of.subset = 5, verbose = TRUE) {
    if (cell.type.column != "idents") Idents(sc.obj) <- sc.obj@meta.data[, cell.type.column]
    levels(sc.obj) <- intersect(levels(sc.obj), unique(Idents(sc.obj) %>% as.vector()))
    cell.cnts <- table(Idents(sc.obj))
    if (sum(cell.cnts < min.cells.of.subset)) {
        sc.obj <- subset(sc.obj, idents = cell.cnts[cell.cnts >= min.cells.of.subset] %>% names())
    }
    sp.obj <- sp.obj[, intersect(colnames(sp.obj), rownames(GetTissueCoordinates(sp.obj)))]
    sp.obj@images[[1]]@coordinates <- sp.obj@images[[1]]@coordinates[colnames(sp.obj), ]

    ovp.genes <- intersect(rownames(sc.obj), rownames(sp.obj))
    sc.obj <- sc.obj[ovp.genes, ] %>% SCTransform(., verbose = verbose)
    sp.obj <- sp.obj[ovp.genes, ] %>% SCTransform(., verbose = verbose, assay = "Spatial")

    ovp.genes <- intersect(rownames(sc.obj), rownames(sp.obj))
    return(list(SC = sc.obj[ovp.genes, ], ST = sp.obj[ovp.genes, ]))
}

#' @title selectMakers
#'
#' @description Select marker genes for each cell type based on single-cell (SC) data.
#' @param sc.obj Seurat object of SC data.
#' @param sc.markers When the number of cell types identified from SC data is one, markers specific to that cell type must be provided. Default: NULL.
#' @param select.markers Method for selecting cell-type specific markers when `sc.markers = NULL`.
#' @param group.size Specify the marker size for each subset derived from single-cell data.
#' @param verbose Boolean indicating whether to show running messages. Default: TRUE.
#' @return A list of markers for cell types.

selectMakers <- function(sc.obj, sc.markers, select.markers, group.size, verbose = TRUE) {
    if (length(levels(sc.obj)) == 1 && is.null(sc.markers)) {
        println("When the SC data contains only one cell type, the marker genes for that cell type must be provided", status = "ERROR")
    }
    if (is.null(sc.markers)) {
        sc.markers <- findScMarkers(sc.obj, group.size, select.markers, verbose = verbose)
    } else {
        bl.status <- identical(names(sc.markers) %>% sort(), levels(sc.obj) %>% sort())
        if (!bl.status) {
            println("Specified marker list of cell types does not match the identifiers of SC data", status = "WARN")
            sc.markers <- findScMarkers(sc.obj, group.size, select.markers, verbose = verbose)
        }
    }
    if (inherits(sc.markers, "character")) sc.markers <- list(XX = sc.markers) %>% `names<-`(levels(sc.obj))
    return(sc.markers)
}
