% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/detectHotSpots.R
\name{dectHotSpotsByTtest}
\alias{dectHotSpotsByTtest}
\title{dectHotSpotsByTtest}
\usage{
dectHotSpotsByTtest(sp.obj, sp.score, knn = 5, p.value.threshold = 0.05)
}
\arguments{
\item{sp.obj}{Seurat object of spatial transcriptome (ST) data.}

\item{sp.score}{Signature score of cell types.}

\item{knn}{K nearest cells used. Default: 5.}

\item{p.value.threshold}{P-value threadhold for filtering out hot-spots. Default: 0.05.}
}
\value{
A list of data frames containing hotspot information and p-values. A value of TRUE indicates a detected hotspot, while FALSE indicates a non-hotspot.
}
\description{
Infer whether the signature scores in spots significant higher than background.
}
