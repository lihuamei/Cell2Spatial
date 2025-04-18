% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{findScMarkers}
\alias{findScMarkers}
\title{findScMarkers}
\usage{
findScMarkers(
  sc.obj,
  group.size,
  select.markers = c("shannon", "wilcox"),
  verbose = TRUE
)
}
\arguments{
\item{sc.obj}{Seurat object of single-cell data.}

\item{group.size}{Marker size of each subset derived from SC data. Default: 100.}

\item{select.markers}{Strategy for inferring marker genes. Default: shannon.}

\item{verbose}{Show running messages or not. Default: TRUE.}
}
\value{
A tibble of identified markers.
}
\description{
Find markers among cell types based on SC data.
}
