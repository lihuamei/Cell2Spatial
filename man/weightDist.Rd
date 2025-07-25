% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/weightDist.R
\name{weightDist}
\alias{weightDist}
\title{weightDist}
\usage{
weightDist(
  sc.obj,
  sp.obj,
  lamba,
  quantile.cut,
  mc.cores = 4,
  use.entire = TRUE,
  dist.based = c("UMAP", "TSNE")
)
}
\arguments{
\item{sc.obj}{Seurat object of SC data.}

\item{sp.obj}{Seurat object of ST data.}

\item{lamba}{Median of cell counts for spots.}

\item{quantile.cut}{Numeric value specifying the quantile threshold for distance scaling. Default is 1, which considers the maximum distance for normalization.}

\item{mc.cores}{Number of cores for parallel running. Default: 4.}

\item{dist.based}{Dimensionality reduction basis used for distance weighting.}
}
\value{
A matrix of weighted distance matrix. Rows represent spot clusters and columns are cell types.
}
\description{
Weight distance between SC and ST data.
}
