#' @title checkParams

#' @description Check the legality of input arguments.
#' @param argg A list of input arguments.
#' @return TRUE

checkParams.runMap2SP <- function(argg) {
<<<<<<< HEAD
    stopifnot(
        is(argg[["sp.obj"]], "Seurat"),
        is(argg[["sc.obj"]], "Seurat"),
        if (!is.null(argg[["sc.markers"]])) is(argg[["sc.markers"]], "list") else TRUE,
        is.character(argg[["ctype"]]),
        is.numeric(argg[["group.size"]]) && argg[["group.size"]] >= 5,
        is.logical(argg[["duplicated"]]),
        is.logical(argg[["partion"]]),
        is.logical(argg[["integ"]]),
        is.numeric(argg[["min.cells.of.subset"]]) && argg[["min.cells.of.subset"]] >= 1,
        if (!is.null(argg[["max.cells.in.spot"]])) is.numeric(argg[["max.cells.in.spot"]]) && argg[["max.cells.in.spot"]] > 0 else TRUE,
        if (!is.null(argg[["fix.cells.in.spot"]])) is.numeric(argg[["fix.cells.in.spot"]]) && argg[["fix.cells.in.spot"]] > 0 else TRUE,
        is.numeric(argg[["n.features"]]) && argg[["n.features"]] >= 30,
        is.numeric(argg[["knn"]]) && argg[["knn"]] >= 1,
        if (!is.null(argg[["sample.size"]])) is.numeric(argg[["sample.size"]]) && argg[["sample.size"]] > 0 else TRUE,
        is.character(argg[["select.markers"]]),
        is.character(argg[["dist.method"]]),
        is.numeric(argg[["n.workers"]]) && argg[["n.workers"]] >= 1,
        is.logical(argg[["verbose"]])
    )
    if (argg[["duplicated"]] && argg[["partion"]]) {
        println("Parameters 'duplicated' and 'partion' cannot be set TRUE at the same time, reset partion = FALSE", verbose = TRUE, status = "WARN")
    }
    return(TRUE)
=======
  stopifnot(
    is(argg[["sp.obj"]], "Seurat"),
    is(argg[["sc.obj"]], "Seurat"),
    if (!is.null(argg[["sc.markers"]])) is(argg[["sc.markers"]], "list") else TRUE,
    is.character(argg[["ctype"]]),
    is.numeric(argg[["group.size"]]) && argg[["group.size"]] >= 5,
    is.logical(argg[["duplicated"]]),
    is.numeric(argg[["min.cells.of.subset"]]) && argg[["min.cells.of.subset"]] >= 1,
    if (!is.null(argg[["max.cells.in.spot"]])) is.numeric(argg[["max.cells.in.spot"]]) && argg[["max.cells.in.spot"]] > 0 else TRUE,
    if (!is.null(argg[["fix.cells.in.spot"]])) is.numeric(argg[["fix.cells.in.spot"]]) && argg[["fix.cells.in.spot"]] > 0 else TRUE,
    is.numeric(argg[["n.features"]]) && argg[["n.features"]] >= 30,
    is.numeric(argg[["knn"]]) && argg[["knn"]] >= 1,
    if (!is.null(argg[["sample.size"]])) is.numeric(argg[["sample.size"]]) && argg[["sample.size"]] > 0 else TRUE,
    is.character(argg[["select.markers"]]),
    is.character(argg[["dist.method"]]),
    is.numeric(argg[["n.workers"]]) && argg[["n.workers"]] >= 1,
    is.logical(argg[["verbose"]])
  )
  return(TRUE)
>>>>>>> 5859e3527b186cec1e7cc9974128adabbcd95a88
}

#' @description Preprocessing of single-cell (SC) data involves a series of steps.
#' @return Seurat object of preprocessed single cell data.

preprocessSCData <- function(sc.obj, sample.size, ctype, min.cells.of.subset) {
<<<<<<< HEAD
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
    return(sc.obj)
=======
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
  return(sc.obj)
>>>>>>> 5859e3527b186cec1e7cc9974128adabbcd95a88
}
