# Cell2Spatial
Cell2Spatial is a sophisticated tool specifically designed for decoding spatial transcriptomic spots at the individual cell level. Ensuring accurate alignment and maximizing practical applications necessitates a match in tissue origin or major cell type representation between single-cell and spatial transcriptomic data.
<p align="center">
	<img src="vignette_files/overview.png" alt="Resized Image" width="800">
</p>

In this tutorial, we'll showcase the installation and usage of Cell2Spatial, allowing precise interpretation of spatial transcriptomic spots at a single-cell granularity.

## 1. Essential dependencies

The Cell2Spatial's code comprises both R and Python components, necessitating essential dependencies.
* Python3:
	* Numpy
	* Pandas
	* tensorflow
	* keras
	* sklearn
	* lapjv
* R (>= 4.0):
	* reticulate
	* Seurat

## 2. Installation

(1) Python3 must be installed and configured in the environment for use with `reticulate`. Additionally, dependencies for Python libraries can be installed using the following command:

``` bash
pip install -r requirements.txt 

```

(2) Installing Cell2Spatial package

``` r
library(devtools)
install_github("lihuamei/Cell2Spatial")

```
## 3. Loading the packages and datasets (scRNA-seq and ST data)
``` r
library(Cell2Spatial)
library(Seurat)
library(dplyr)
library(randomcoloR)
library(tidydr)
```
``` r
sp.obj <- system.file("data", "Kindney_SP.RDS", package = "Cell2Spatial") %>% readRDS(.)
sc.obj <- system.file("data", "Kindney_SC.RDS", package = "Cell2Spatial") %>% readRDS(.)

```

## 4. Assigning single-cells to spatial coordinates.
```r
sce <- runMap2SP(sp.obj, sc.obj, ctype = "mainCtype", res = 0.8, group.size = 30, fix.cells.in.spot = 10)

```

### 5. Visualization of mapping results.

``` r
sc.obj <- SCTransform(sc.obj, ncells = 3000, verbose = FALSE) %>%
    RunPCA(verbose = FALSE) %>%
    RunUMAP(dims = 1 : 30, verbose = FALSE)
```
``` r
set.seed(2023063)
cell.colors <- randomcoloR::distinctColorPalette(length(unique(sc.obj$mainCtype)))  %>% `names<-`(unique(sc.obj$mainCtype))
gp1 <- SpatialPlot(sce, group.by = 'Cell2Spatial', pt.size.factor=0.6, cols = cell.colors, image.alpha = 0.5, stroke = NA)
gp2 <- DimPlot(sc.obj, label = TRUE, cols = cell.colors) + theme_dr(xlength = 0.2, ylength = 0.2, arrow = grid::arrow(length = unit(0.1, "inches"), ends = 'last', type = "closed")) + theme(panel.grid = element_blank())
gp1 + gp2
```
<p align="center">
	<img src="vignette_files/mapping_results.png" alt="Resized Image" width="800">
</p>

## 6. Setting parameters
|**Parameters**|**Description**                      |
|----------|-----------------------------------------|
|sp.obj	|Seurat object of spatial transcriptome (ST) data.|
|sc.obj	|Seurat object of single-cell (SC) data.|
|sc.markers	|A list of markers for cell types. If not provided, automated inference is performed using the modified Shannon-entropy method. Default: NULL.|
|ctype	|Specify the column name for the cell type in meta.data slot of the SC Seurat object. Default: idents.|
|group.size	|Specify the marker size for each subset derived from single-cell data. This information is crucial for estimating the activity of markers within each cell type. Default: 30.|
|res	|Resolution for clustering ST spot. Default: 0.8.|
|integ	|Integration of SC and ST data or not. For high-resolution ST data, please set FALSE to accelerate the running time. Default: TRUE.|
|duplicated	|Assigning individual cells to ST coordinates with duplicated or not. Default: FALSE.|
|partion	|Split into sub-modules mapped to spatial positions when 'duplicated' set to TRUE. Default: TRUE.|
|min.cells.of.subset	|Include cells detected in at least one cell type in the SC data. Default: 5.|
|max.cells.in.spot	|Maximum number of cells in ST spots. Default: 10 (for 10x Visium).|
|fix.cells.in.spot	|Fixed number of cells assigned to a spot. Default: NULL.|
|n.features	|Top k highly variable features used for integrating SC and ST data. Default: 3000.|
|knn	|Utilizing the k nearest spots to filter out spots that do not contain cell types present in the SC reference. Default: 5.|
|sample.size	|Down-sampling a small subset of SC data for mapping to ST coordinates. Default: NULL.|
|select.markers	|Method for selecting cell-type specific markers when 'sc.markers = NULL': modified Shannon-entropy strategy (shannon) or Wilcoxon test (wilcox) implemented by the Seurat package. Default: shannon.|
|dist.method	|Measure the distance between single-cell and spots, maximum likehood model (mle) or correlation (cor). Default: mle.|
|return.type	|Assigned results can be of the object type Seurat or SingleCellExperiment. Default: Seurat.|
|n.wrokers	|Number of cores for parallel processing. Default: 4.|
|verbose	|Show running messages or not. Default: TRUE.|

## 7. Session infos
```r
> sessionInfo()
R version 4.1.2 (2021-11-01)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 22.04.2 LTS

Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.10.0
LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.10.0

locale:
[1] C

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base

other attached packages:
[1] Cell2Spatial_1.0.0 reticulate_1.30    ggplot2_3.4.2      dplyr_1.1.2
[5] SeuratObject_4.1.3 Seurat_4.3.0.1     devtools_2.4.5     usethis_2.2.2

loaded via a namespace (and not attached):
  [1] spam_2.9-1                  plyr_1.8.8
  [3] igraph_1.5.0.1              lazyeval_0.2.2
  [5] sp_2.0-0                    splines_4.1.2
  [7] listenv_0.9.0               scattermore_1.2
  [9] GenomeInfoDb_1.30.1         digest_0.6.33
 [11] htmltools_0.5.5             viridis_0.6.4
 [13] fansi_1.0.4                 magrittr_2.0.3
 [15] memoise_2.0.1               tensor_1.5
 [17] cluster_2.1.2               ROCR_1.0-11
 [19] remotes_2.4.2.1             globals_0.16.2
 [21] matrixStats_1.0.0           spatstat.sparse_3.0-2
 [23] prettyunits_1.1.1           colorspace_2.1-0
 [25] ggrepel_0.9.3               callr_3.7.3
 [27] crayon_1.5.2                RCurl_1.98-1.12
 [29] jsonlite_1.8.7              progressr_0.13.0
 [31] spatstat.data_3.0-1         survival_3.2-13
 [33] zoo_1.8-12                  glue_1.6.2
 [35] polyclip_1.10-4             gtable_0.3.3
 [37] zlibbioc_1.40.0             XVector_0.34.0
 [39] leiden_0.4.3                DelayedArray_0.20.0
 [41] pkgbuild_1.4.2              future.apply_1.11.0
 [43] SingleCellExperiment_1.16.0 maps_3.4.1
 [45] BiocGenerics_0.40.0         abind_1.4-5
 [47] scales_1.2.1                spatstat.random_3.1-5
 [49] miniUI_0.1.1.1              Rcpp_1.0.11
 [51] viridisLite_0.4.2           xtable_1.8-4
 [53] dotCall64_1.0-2             stats4_4.1.2
 [55] profvis_0.3.8               htmlwidgets_1.6.2
 [57] httr_1.4.6                  RColorBrewer_1.1-3
 [59] ellipsis_0.3.2              ica_1.0-3
 [61] urlchecker_1.0.1            pkgconfig_2.0.3
 [63] uwot_0.1.16                 deldir_1.0-9
 [65] utf8_1.2.3                  ggcorrplot_0.1.4
 [67] tidyselect_1.2.0            rlang_1.1.1
 [69] reshape2_1.4.4              later_1.3.1
 [71] munsell_0.5.0               tools_4.1.2
 [73] cachem_1.0.8                cli_3.6.1
 [75] generics_0.1.3              ggridges_0.5.4
 [77] stringr_1.5.0               fastmap_1.1.1
 [79] goftest_1.2-3               processx_3.8.2
 [81] fs_1.6.3                    fitdistrplus_1.1-11
 [83] purrr_1.0.1                 RANN_2.6.1
 [85] pbapply_1.7-2               future_1.33.0
 [87] nlme_3.1-155                mime_0.12
 [89] compiler_4.1.2              rstudioapi_0.15.0
 [91] plotly_4.10.2               png_0.1-8
 [93] spatstat.utils_3.0-3        tibble_3.2.1
 [95] stringi_1.7.12              ps_1.7.5
 [97] desc_1.4.2                  fields_14.1
 [99] lattice_0.20-45             Matrix_1.6-0
[101] vctrs_0.6.3                 pillar_1.9.0
[103] lifecycle_1.0.3             spatstat.geom_3.2-4
[105] lmtest_0.9-40               RcppAnnoy_0.0.21
[107] data.table_1.14.8           cowplot_1.1.1
[109] bitops_1.0-7                irlba_2.3.5.1
[111] GenomicRanges_1.46.1        httpuv_1.6.11
[113] patchwork_1.1.2             R6_2.5.1
[115] promises_1.2.0.1            KernSmooth_2.23-20
[117] gridExtra_2.3               IRanges_2.28.0
[119] parallelly_1.36.0           sessioninfo_1.2.2
[121] codetools_0.2-18            MASS_7.3-55
[123] pkgload_1.3.2.1             SummarizedExperiment_1.24.0
[125] rprojroot_2.0.3             withr_2.5.0
[127] sctransform_0.3.5           S4Vectors_0.32.4
[129] GenomeInfoDbData_1.2.7      parallel_4.1.2
[131] grid_4.1.2                  tidyr_1.3.0
[133] MatrixGenerics_1.6.0        Rtsne_0.16
[135] spatstat.explore_3.2-1      Biobase_2.54.0
[137] shiny_1.7.4.1


```
