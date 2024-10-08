% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/STICCC_utilities.R
\name{getClusters}
\alias{getClusters}
\title{Cluster by gene expression.}
\usage{
getClusters(k = 2, matrix, features = NA, clustfun = "ward.D2")
}
\arguments{
\item{k}{numeric. Number of clusters.}

\item{matrix}{data.frame or matrix. Non-normalized expression matrix with rows as features
and columns as cells.}

\item{features}{character. Which features (colnames of matrix) to include in clustering analysis.
Default all features.}

\item{clustfun}{Which clustering algorithm to use. Possible values: "ward.D", "ward.D2",
"single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or
"centroid" (= UPGMC).}
}
\value{
data.frame with cells as rows and metadata characteristics as columns.
}
\description{
Returns a vector of cluster numbers corresponding to each cell in a dataset.
}
