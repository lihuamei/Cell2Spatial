#' @title plotSc2Spatial
#'
#' @description Visualizing the distribution of single-cells on spatial coordinates based on SingleCellExperiment object.
#' @param sp.obj Seurat object of spatial transcriptome (ST) data.
#' @param sce SingleCellExperiment object represents the output generated by the Cell2Spatial tool after projecting single-cell data onto spatial coordinates
#' @param tar.ctypes A vector of cell types to show exclusively. Default: NULL.
#' @param pt.size Point size of single-cells on spatial image. Default: 1.0.
#' @param colors A vector colors corresponding to the cell type. Default: NULL.
#' @param shape Point shape. Default: 21.
#' @param percent Down-sampling cell types if percent is not NULL (0 < percent <= 1).
#' @return ggplot object
#' @export plotSc2Spatial

plotSc2Spatial <- function(sp.obj, sce, tar.ctypes = NULL, pt.size = 1.0, colors = NULL, shape = 21, percent = NULL) {
    pseud.coord <- sce@colData %>% as.data.frame()
    img <- sp.obj@images[[names(sp.obj@images)]]@image
    if (!is.null(tar.ctypes)) {
        tar.ctypes.sub <- intersect(tar.ctypes, pseud.coord$Cell2Spatial %>% unique())
        if (length(tar.ctypes.sub) == 0) println("No cell types found in sce data, exiting...", status = "ERROR")
        pseud.coord <- subset(pseud.coord, Cell2Spatial %in% tar.ctypes.sub)
    }
    tar.ctypes <- pseud.coord$Cell2Spatial %>% unique()
    if (!is.null(colors) && is.null(names(colors))) {
        colors <- colors[1:length(tar.ctypes)] %>% `names<-`(tar.ctypes)
    } else {
        colors <- scales::hue_pal()(tar.ctypes %>% length()) %>% `names<-`(tar.ctypes)
    }
    colors <- colors[tar.ctypes]

    if (!is.null(percent) && percent > 0 && percent <= 1) {
        pseud.coord <- lapply(pseud.coord$Cell2Spatial %>% unique(), function(ct) {
            df <- subset(pseud.coord, Cell2Spatial == ct)
            int.idx <- sample(nrow(df), floor(nrow(df) * percent))
            df[int.idx, ]
        }) %>%
            do.call(rbind, .) %>%
            as.data.frame()
    }
    gp <- ggplot(pseud.coord, aes(x = y, y = x, fill = Cell2Spatial)) +
        annotation_raster(img, xmin = 0, xmax = ncol(img), ymin = 0, ymax = -nrow(img)) +
        geom_point(shape = shape, size = pt.size) +
        ggplot2::scale_y_reverse() +
        theme_bw() +
        theme(
            panel.background = element_rect(fill = "black"),
            panel.border = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.ticks = element_blank(),
            axis.text = element_blank(),
            plot.title = element_text(hjust = 0.5)
        ) +
        scale_fill_manual(values = colors) +
        guides(fill = guide_legend(override.aes = list(size = 4))) +
        ggtitle(sp.obj@images %>% names()) +
        xlab("") +
        ylab("")
    return(gp)
}

#' @title plotCoexistCtypes
#'
#' @description Showing coexist indexes of cell-types.
#' @param j.index A data.frame of coexist indexes of cell-types.
#' @param j.cut J-index value below than the specified cutoff will be set to zero. Default: 0.
#' @param cols A vector colors with names, which contains three colors, representing low, medium, and high. Default: c("#91a28c", "white", "#8f2c37").
#' @return ggplot object.
#' @export plotCoexistCtypes

plotCoexistCtypes <- function(j.index, j.cut = 0, cols = c("#91a28c", "white", "#8f2c37")) {
    cols <- cols[1:3]
    j.index[j.index < j.cut] <- 0

    res.p <- ggcorrplot::ggcorrplot(
        j.index %>% as.matrix(),
        hc.order = F,
        outline.color = "white",
        tl.srt = 60,
        tl.cex = 12,
        lab_size = 7,
        colors = cols
    ) +
        theme(
            plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
            panel.background = element_blank(),
            plot.background = element_blank(),
            panel.border = element_rect(colour = "grey89", fill = NA, size = 0.5),
            axis.text = element_text(size = 12),
            axis.ticks = element_blank(),
            axis.title = element_blank(),
            legend.title = element_text(size = 16, face = "bold"),
            legend.text = element_text(size = 16),
            legend.key = element_rect(colour = "transparent", fill = "white"),
            legend.key.size = unit(0.45, "cm")
        ) +
        coord_fixed() +
        ggtitle("J-index") + theme(plot.title = element_text(size = 22, face = "bold"))
    return(res.p)
}
