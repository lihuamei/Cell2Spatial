#' @title calcSpotsDist

#' @description Identified significant spots located in medulla areas and calculate the distances between spots.
#' @param sp.obj Spatial seurat object data.
#' @return Eulerian distance matrix.

calcSpotsDist <- function(sp.obj) {
  image.coord <- GetTissueCoordinates(sp.obj)
  elu.dist <- dist(image.coord) %>% as.matrix()
  return(elu.dist)
}

#' @title getRandomCords

#' @description Get random coordinates for mapping single-cells to spatial.
#' @param sp.coord Coordinates of ST image data.
#' @param num.cells A list of number of cells for mapping to spot.
#' @param seed Set seed for generating radom coordinates. Default: 123456.
#' @param n.workers Number of cores for parallel processing. Default: 4.
#' @return A data.frame of generated spatial coordinates.
#' @export getRandomCords

getRandomCords <- function(sp.coord, num.cells, seed = 123456, n.workers = 4) {
  set.seed(seed)
  future::plan("multicore", workers = n.workers)
  if (is.numeric(num.cells)) num.cells <- rep(num.cells, nrow(sp.coord)) %>% as.list()
  sp.ED <- fields::rdist(as.matrix(sp.coord))
  min.dist <- apply(sp.ED, 1, function(xx) order(xx)[2] %>% xx[.]) %>%
    {
      median(.) / 2
    }
  sp.coord.new <- future.apply::future_lapply(1:nrow(sp.coord), function(idx) {
    circle <- spatstat.random::runifdisc(
      num.cells[[idx]],
      radius = min.dist,
      centre = c(sp.coord[idx, 1], sp.coord[idx, 2]),
      nsim = 1,
      drop = TRUE
    )
    tmp.cord <- data.frame(x = circle$x, y = circle$y) %>% `colnames<-`(c("x", "y"))
    tmp.cord$centerSPOT <- paste0(sp.coord[idx, 1], "x", sp.coord[idx, 2])
    tmp.cord$centerX <- sp.coord[idx, 1]
    tmp.cord$centerY <- sp.coord[idx, 2]
    tmp.cord$SpotName <- rownames(sp.coord)[idx]
    return(tmp.cord)
  }, future.seed = TRUE) %>%
    do.call(rbind, .) %>%
    as.data.frame()
  sp.coord.new <- sp.coord.new[!duplicated(paste0(sp.coord.new$x, "x", sp.coord.new$y)), ] %>%
    `rownames<-`(paste0(sp.coord.new$x, "x", sp.coord.new$y))
  return(sp.coord.new)
}

#' @title assignSCcords
#'
#' @description Assign single-cells to spatial coordinates.
#' @param sp.obj Seurat object of spatial transcriptome (ST) data.
#' @param sc.obj Seurat object of single-cell data.
#' @param out.sc A list of ordered single-cells correspond to spots.
#' @param num.cells A list of number of cells for mapping to spot.
#' @param return.type Assigned results can be of the object type Seurat or SingleCellExperiment. Default: Seurat.
#' @param n.workers  Number of cores for parallel processing. Default: 4.
#' @return Seurat or SingleCellExperiment object for mapped results.
#' @export assignSCcords

assignSCcords <- function(sp.obj, sc.obj, out.sc, num.cells, n.workers = 4, return.format = "Seurat") {
  sp.cords <- GetTissueCoordinates(sp.obj) %>% .[intersect(rownames(.), names(out.sc)), ]
  sp.obj <- sp.obj[, rownames(sp.cords)]
  out.sc <- out.sc[rownames(sp.cords)]
  cord.new <- getRandomCords(sp.cords, num.cells = num.cells, n.workers = n.workers)

  future::plan("multicore", workers = n.workers)
  map.cords <- future.apply::future_lapply(1:length(out.sc), function(ispot) {
    ispot.cords <- sp.cords[ispot, ]
    sub.cords <- cord.new[cord.new$centerSPOT == paste0(ispot.cords[1], "x", ispot.cords[2]), ]
    sub.cords$EDwithCenter <- fields::rdist(as.matrix(sub.cords[, 1:2]), as.matrix(ispot.cords))
    sub.cords$Cell2Spatial <- Idents(sc.obj)[out.sc[[ispot]]] %>% as.vector()
    sub.cords$Cell <- out.sc[[ispot]]
    sub.cords
  }, future.seed = TRUE) %>%
    do.call(rbind, .) %>%
    as.data.frame() %>%
    `rownames<-`(make.unique(.$Cell)) %>%
    {
      cbind.data.frame(sc.obj@meta.data[.$Cell, ], .)
    }

  count.sc <- GetAssayData(sc.obj, slot = "count")
  count.CT <- as(count.sc[, map.cords$Cell], "sparseMatrix")
  colnames(count.CT) <- rownames(map.cords)

  if (return.format == "Seurat") {
    sce <- CreateSeuratObject(count = count.CT, meta.data = map.cords, project = "Cell2Spatial", assay = "Spatial")
    sce@images <- sp.obj@images
    sce@images[[1]]@assay <- DefaultAssay(sce)
    sce@images[[1]]@coordinates <- data.frame(imagerow = map.cords$x, imagecol = map.cords$y) %>% `rownames<-`(rownames(map.cords))
    sce@images[[1]]@scale.factors <- sp.obj@images[[1]]@scale.factors
    sce@images[[1]]@coordinates <- sce@images[[1]]@coordinates / sce@images[[1]]@scale.factors$lowres
  } else {
    sce <- SingleCellExperiment::SingleCellExperiment(
      list(counts = count.CT),
      colData = as.data.frame(map.cords),
      rowData = as.data.frame(rownames(count.CT))
    )
  }
  return(sce)
}
