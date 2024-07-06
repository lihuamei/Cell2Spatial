#' @title checkParams.runMap2SP
#'
#' @description Check the legality of input arguments.
#' @param argg A list of input arguments.
#' @return A boolean value: TRUE.

checkParams.runMap2SP <- function(argg) {
    stopifnot(
        is(argg[["sp.obj"]], "Seurat"),
        is(argg[["sc.obj"]], "Seurat"),
        is.character(argg[["ctype"]]),
        if (!is.null(argg[["sc.markers"]])) is(argg[["sc.markers"]], "list") else TRUE,
        is.numeric(argg[["group.size"]]) && argg[["group.size"]] >= 5,
        is.logical(argg[["partion"]]),
        if (!is.null(argg[["max.cells.in.spot"]])) is.numeric(argg[["max.cells.in.spot"]]) && argg[["max.cells.in.spot"]] > 0 else TRUE,
        if (!is.null(argg[["fix.cells.in.spot"]])) is.numeric(argg[["fix.cells.in.spot"]]) && argg[["fix.cells.in.spot"]] > 0 else TRUE,
        is.numeric(argg[["knn"]]) && argg[["knn"]] >= 1,
        if (!is.null(argg[["sample.size"]])) is.numeric(argg[["sample.size"]]) && argg[["sample.size"]] > 0 else TRUE,
        is.character(argg[["detect.hotspot.method"]]),
        is.numeric(argg[["p.value.threshold"]]) && argg[["p.value.threshold"]] <= 1 && argg[["p.value.threshold"]] > 0,
        is.numeric(argg[["quantile.threshold"]]) && argg[["quantile.threshold"]] <= 1 && argg[["quantile.threshold"]] > 0,
        is.character(argg[["select.markers"]]),
        is.character(argg[["dist.method"]]),
        is.logical(argg[["integ.entire"]]),
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
#' @param sample.size Proportion for downsampling single cells of each cell type.
#' @param ctype Specify the column name for the cell type in the meta.data slot of the SC Seurat object. Default: idents.
#' @param min.cells.of.subset Include cells detected in at least one cell type in the SC data. Default: 5.
#' @param verbose Boolean indicating whether to show running messages. Default: TRUE.
#' @return A list of Seurat objects of preprocessed data.

preprocessSeqData <- function(sp.obj, sc.obj, sample.size, ctype, min.cells.of.subset = 5, verbose = TRUE) {
    if (!is.null(sample.size)) {
        if (sample.size > 1) {
            sc.obj <- downSamplSeurat(sc.obj, cnt = sample.size)
        } else {
            sc.obj <- downSamplSeurat(sc.obj, percent = sample.size)
        }
    }
    if (ctype != "idents") Idents(sc.obj) <- sc.obj@meta.data[, ctype]
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
    return(list(SC = sc.obj, ST = sp.obj))
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
