% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/STICCC_utilities.R
\name{plotGrid}
\alias{plotGrid}
\title{Plot smoothed velocities}
\usage{
plotGrid(
  sce,
  colorVar = NA,
  plotLoadings = F,
  loadingFactor = 3.5,
  colorPalette = NA,
  outputDir = NA,
  plotSuffix = NA,
  scalingFactor = NA,
  minMagnitude = 0.01,
  arrowSize = 3,
  arrowheadSize = 0.3
)
}
\arguments{
\item{sce}{viccc object after computing inferred velocities}

\item{colorVar}{character. Colname of colData(sce) by which to color the plot. May be left NA for grey plot}

\item{plotLoadings}{logical. Whether to plot the loadings for each gene. Only recommended for
very small (<5 gene) circuits. Default FALSE.}

\item{loadingFactor}{numeric. Constant to multiply loadings by to improve visual clarity.}

\item{colorPalette}{vector. Vector containing colors - be sure to include at least as many colors
as there are unique values of colorVar. Default color palette includes 8 values.}

\item{outputDir}{character. File path for the directory to save the plot.}

\item{plotSuffix}{character. Suffix for plot filename.}

\item{scalingFactor}{numeric. Constant to multiply vectors by for visual clarity.}

\item{minMagnitude}{numeric. Smallest magnitude for which a vector will be drawn. Default 0.01.}

\item{arrowSize}{numeric. Size of arrows in plot.}

\item{arrowheadSize}{numeric. Size of arrowheads in plot.}
}
\description{
Saves a pdf plot showing the first 2 PCs of the dataset with grid-smoothed
inferred velocities superimposed.
}
