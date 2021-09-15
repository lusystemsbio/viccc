#' @export
#' @title Load topology file
#' @description Returns a data.frame from reading a .tpo file with 3 columns: Source, Target, Type. Type 1 indicates
#' activation and type 2 indicates inhibition.
#' @return data.frame of a network topology.
#' @param topoName character. Filename of a .tpo topology table (do not include file extension).
#' @param topoDir character. Path to the directory containing the .tpo file. Defaults to "input_topos"
#' and will create this directory if needed.
loadTopo <- function(topoName,
                     topoDir=NA) {
  ## Import topology
  if(is.na(topoDir)) {
    topoDir <- file.path(getwd(),"input_topos")
  }

  if(!dir.exists(topoDir)) {
    dir.create(topoDir)
  }

  fName <- file.path(topoDir,paste0(topoName,".tpo"))
  if(!file.exists(file.path(topoDir))) {
    print("Topology file not found")
    return(NULL)
  }
  topo <- read.table(fName, header = T)
  return(topo)
}

#' @export
#' @title Simulate topology
#' @description Simple wrapper function to simulate a topology using sracipeSimulate from the sRACIPE package. Returns
#' a RACIPE object.
#' @return sRACIPE object.
#' @import sRACIPE
#' @param topo data.frame describing a network topology with columns Source, Target, and Type, where type 1
#' indicates activation and type 2 indicates inhibition.
#' @param ... Other parameters for the sracipeSimulate function.
simTopo <- function(topo,
                    ...) {
  ## Set seed for reproducibility
  set.seed(123)

  ## Simulate topology
  racipe <- RacipeSE() # Construct an empty RacipeSE object
  sracipeCircuit(racipe) <- topo
  racipe <- sracipeSimulate(racipe, ...)
  return(racipe)
}

#' @export
#' @title Create VICCC object.
#' @description Returns a SingleCellExperiment object initialized for use with VICCC methods.
#' @return SingleCellExperiment object modified slightly for VICCC methods.
#' @import SingleCellExperiment
#' @param topo data.frame. Topology table.
#' @param exprMat data.frame or matrix. Expression matrix with genes as rows and cells as columns.
#' @param normData data.frame or matrix. Normalized (log2+1) expression matrix.
#' @param topoName character. Name of the topology, used in output file names.
#' @param expName character. Name for the experiment (if you intend to do multiple analyses with the same topology).
#' Default <topoName>_SCE.
#' @param radius numeric. Proportion between 0 and 1 of the maximum pairwise distance to evaluate when computing velocity
#' of a cell. Default 0.15
#' @param minNeighbors numeric. Minimum number of neighbors within specified radius to compute a velocity for a cell. Cells
#' with fewer than this number of neighbors will be ignored, but may be counted as neighbors for other eligible cells.
#' Default 5.
#' @param scalingFactor numeric. Factor to multiply the magnitude of vectors for easier visual interpretation. Default 1.
#' @param gridScalingFactor numeric. Factor to multiply the magnitude of vectors when applying grid-based smoothing. Default 1.
#' @param verbose logical. Whether to print progress statements during analysis with this VICCC object. Default TRUE.
vicSE <- function(topo,
                  exprMat,
                  normData,
                  topoName,
                  expName = NA,
                  radius = 0.15,
                  minNeighbors = 5,
                  scalingFactor = 1,
                  gridScalingFactor = 1,
                  verbose = T
                  ) {
  vic <- SingleCellExperiment(assays = SimpleList(counts=exprMat, normcounts=normData))
  if(is.na(expName)) {
    vic@metadata$experimentName <- paste0(topoName,"_SCE")
  } else {
    vic@metadata$experimentName <- expName
  }
  vic@metadata$topoName <- topoName
  vic@metadata$topo <- topo

  vic@metadata$params <- list(sample_radius=radius, plotScalingFactor=scalingFactor, gridPlotScalingFactor=gridScalingFactor,
                              minNeighbors=minNeighbors, verbose=verbose)
  return(vic)
}


#' @export
#' @title Create metadata table.
#' @description Returns a data.frame containing metadata for each cell in a dataset. This can be
#' substituted with a previously existing data.frame and is only necessary if no metadata table
#' exists. This will also divide the dataset into k clusters using the ward.d2 method.
#' @return data.frame with cells as rows and metadata characteristics as columns.
#' @param exprMat data.frame or matrix. Non-normalized expression matrix with rows as features
#' and columns as cells.
#' @param cluster logical. Whether to perform hierarchical clustering.
#' @param ... Additional arguments to the (internal) clustering function - k (number of clusters),
#' features (which features to cluster by, defaults to all), and clustfun (default "ward.D2")
prepMetadata <- function(exprMat,
                         cluster = T,
                         ...            # Parameters for getClusters function
) {
  metadata <- as.data.frame(t(exprMat))
  clusters <- getClusters(matrix = metadata, ...)
  metadata$SampleID <- rownames(metadata)
  metadata$Cluster <- clusters

  return(metadata)

}

#' @title Cluster by gene expression.
#' @description Returns a vector of cluster numbers corresponding to each cell in a dataset.
#' @return data.frame with cells as rows and metadata characteristics as columns.
#' @importFrom stats hclust
#' @param k numeric. Number of clusters.
#' @param matrix data.frame or matrix. Non-normalized expression matrix with rows as features
#' and columns as cells.
#' @param features character. Which features (colnames of matrix) to include in clustering analysis.
#' Default all features.
#' @param clustfun Which clustering algorithm to use. Possible values: "ward.D", "ward.D2",
#' "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or
#' "centroid" (= UPGMC).
getClusters <- function(k = 2,                  # Number of clusters
                        matrix,                 # Input matrix - rows = samples, cols = features
                        features = NA,          # List of matrix colnames
                        clustfun = "ward.D2"    # Clustering function
) {
  # Check that features are properly entered
  if(is.na(features)) {
    features <- colnames(matrix)
  } else if(any(!features %in% colnames(matrix))) {
    print("Error: specified features not in matrix colnames")
  }

  # Clustering
  hclust_dend <- hclust(distfun(matrix[,c(features)]), method = clustfun)
  clusters <- cutree(hclust_dend, k=k)
  return(clusters)

}

#' @export
#' @title Perform PCA on gene expression data
#' @description Returns a vector of cluster numbers corresponding to each cell in a dataset.
#' @return list with class "prcomp" containing elements sdev, rotation, x, center, and scale.
#' @importFrom stats prcomp
#' @param exprMat data.frame or matrix. Expression matrix with rows as features and columns as
#' cells.
runPCA <- function(exprMat) {
  ## TODO: convert this so the input is a vic object and output is same, but modified
  return(prcomp(t(exprMat)))
}


#' @export
#' @title Compute grid points for later smoothing
#' @description Computes grid points for future steps that will smooth the inferred velocities
#' in PCA space.
#' @return viccc object with defined grid points in PCA space
#' @import SingleCellExperiment
#' @importFrom tidyr expand_grid
#' @param sce viccc object. Must have already executed the function runPCA().
#' @param grid.length numeric. How many grid points should be in each side of the smoothed plot.
#' Default 30.
computeGrid <- function(sce,
                        grid.length = 30
) {
  # Get posMat
  posMat <- reducedDim(sce, "PCA")

  # Define grid points, create metadata structure for later
  ## Create uniform grid
  # Set grid area (2-dimensional?) w uniform points
  xMin <- floor(min(posMat[,1]))
  xMax <- ceiling(max(posMat[,1]))
  yMin <- floor(min(posMat[,2]))
  yMax <- ceiling(max(posMat[,2]))
  x.min <- xMin - 1
  x.max <- xMax + 1
  y.min <- yMin - 1
  y.max <- yMax + 1

  x.points <- seq(x.min, x.max, length.out = grid.length)
  y.points <- seq(y.min, y.max, length.out = grid.length)
  grid.dist <- sqrt((x.points[2] - x.points[1])^2 + (y.points[2] - y.points[1])^2) / 2
  grid.x.dist <- x.points[2] - x.points[1]
  grid.y.dist <- y.points[2] - y.points[1]

  grid.df <- as.data.frame(expand_grid(x.points, y.points))
  grid.df$GridPoint <- 1:nrow(grid.df)
  grid.df$dx <- NA
  grid.df$dy <- NA

  sce@metadata$params$xMin <- x.min
  sce@metadata$params$xMax <- x.max
  sce@metadata$params$yMin <- y.min
  sce@metadata$params$yMax <- y.max

  sce@metadata$grid.df <- grid.df
  sce@metadata$params$grid.dist <- grid.dist
  sce@metadata$params$grid.x.dist <- grid.x.dist
  sce@metadata$params$grid.y.dist <- grid.y.dist

  return(sce)
}


#' @export
#' @title Compute pairwise distance between single cells in PCA space
#' @description Computes grid points for future steps that will smooth the inferred velocities
#' in PCA space.
#' @return viccc object with defined grid points in PCA space
#' @import SingleCellExperiment
#' @param sce viccc object. Must have already executed the function runPCA().
#' @param usePCs logical. Whether to use PCs instead of gene expression values for distance computation.
#' Setting to FALSE may drastically increase runtime in datasets with many genes. Default TRUE.
#' @param numPCs numeric. How many PCs to use for distance computation. Larger values will take
#' longer, especially for large datasets, but also incorporate more information. Default 10.
computeDist <- function(sce,
                        usePCs=TRUE,
                        numPCs=10) {

  if(usePCs) {
    PCAmat <- reducedDim(sce,"PCA")
    sce@metadata$dist_mat <- as.matrix(dist(PCAmat[,1:min(numPCs, ncol(PCAmat))], method = "euclidean"))
  } else {
    exprMat <- assay(sce,"normcounts")
    sce@metadata$dist_mat <- as.matrix(dist(t(exprMat), method = "euclidean"))
  }
  sce@metadata$max_dist <- max(sce@metadata$dist_mat)

  return(sce)
}


#' @export
#' @title Grid-based smoothing of inferred velocities
#' @description Uses inverse distance weighting to compute a smoothed grid of velocity vectors.
#' @return viccc object with smoothed data in sce@metadata$grid.df
#' @param sce viccc object after computing inferred velocities.
#' @param scalingFactor numeric. Constant to multiply vectors by for visual clarity.
computeGridVectors <- function(sce,
                               scalingFactor = NA) {
  if(is.na(scalingFactor)) {
    scalingFactor <- sce@metadata$params$gridPlotScalingFactor
  }

  # Get gridpoints from metadata
  plot_df <- as.data.frame(colData(sce))
  grid.df <- sce@metadata$grid.df
  grid.df$numNeighbors <- NA
  grid.dist <- sce@metadata$params$grid.dist
  x.dist <- sce@metadata$params$grid.x.dist
  y.dist <- sce@metadata$params$grid.y.dist
  arrow.max <- grid.dist * 0.5

  ## Iterate thru each gridpoint & calculate weighted k-means average vector
  for(point in grid.df$GridPoint) {

    ## Find points in this grid space
    query <- grid.df[which(grid.df$GridPoint == point),c(1:2)]
    colnames(query) <- c("X","Y")
    neighbors <- which(plot_df$X <= (query$X + x.dist) &
                         plot_df$X >= (query$X - x.dist) &
                         plot_df$Y <= (query$Y + y.dist) &
                         plot_df$Y >= (query$Y - y.dist))


    subset_data <- plot_df[neighbors,c("X","Y","dX","dY")]
    subset_data <- na.omit(subset_data)
    grid.df[point,"numNeighbors"] <- nrow(subset_data)

    # if fewer than 3 points in selected region, set vector to 0
    if(nrow(subset_data) >= 3) {

    } else if(nrow(subset_data) < 3) {
      ## Add to grid.df as 0
      grid.df[which(grid.df$GridPoint == point),"dx"] <- 0
      grid.df[which(grid.df$GridPoint == point),"dy"] <- 0
      next
    }

    ## Calculate distance & proximity (1/dist) from grid point to each nearby cell
    neighbor.dists <- as.matrix(dist(rbind(query,subset_data[,c(1:2)]),
                                     method = "euclidean"))
    neighbor.dists <- neighbor.dists[which(rownames(neighbor.dists) == point),-which(colnames(neighbor.dists) == point)]
    subset_data$Dist <- neighbor.dists
    subset_data$Proximity <- 1 / subset_data$Dist

    ## Calculate weighted vector (closer = more weight)
    avg_dx <- Avg_IDW(subset_data$dX, subset_data$Proximity)
    avg_dy <- Avg_IDW(subset_data$dY, subset_data$Proximity)

    ## Add to grid.df
    grid.df[which(grid.df$GridPoint == point),"dx"] <- avg_dx
    grid.df[which(grid.df$GridPoint == point),"dy"] <- avg_dy
  }

  ## Scale all vectors to the fixed maximum (fixed magnitude for all?)
  grid.df$Magnitude <- sqrt(grid.df$dx^2 + grid.df$dy^2)
  mean.magnitude <- mean(grid.df$Magnitude[which(grid.df$Magnitude > 0)])
  scaling.factor <- grid.dist / mean.magnitude

  grid.df$dx <- grid.df$dx * scaling.factor * scalingFactor
  grid.df$dy <- grid.df$dy * scaling.factor * scalingFactor
  grid.df$Magnitude <- grid.df$Magnitude * scaling.factor * scalingFactor

  ## Write to object
  sce@metadata$grid.df <- grid.df

  return(sce)
}

#' @title Inverse distance weighting
#' @description Small utility function to
#' @return numeric representing the weighted average of vals by weights
#' @param vals numeric vector. Values to compute average of.
#' @param weights numeric vector. Weights (proximities) for each value.
Avg_IDW <- function(vals, weights) {
  return(sum(vals * weights) / sum(weights))
}

#' @export
#' @title Plot single-cell velocities
#' @description Saves a pdf plot showing the first 2 PCs of the dataset with individual cells'
#' inferred velocities superimposed
#' @return NULL
#' @import ggplot2
#' @importFrom varhandle check.numeric
#' @param sce viccc object after computing inferred velocities
#' @param colorVar character. Colname of colData(sce) by which to color the plot
#' @param plotLoadings logical. Whether to plot the loadings for each gene. Only recommended for
#' very small (<5 gene) circuits. Default FALSE.
#' @param loadingFactor numeric. Constant to multiply loadings by to improve visual clarity.
#' @param colorPalette vector. Vector containing colors - be sure to include at least as many colors
#' as there are unique values of colorVar. Default color palette includes 8 values.
#' @param outputDir character. File path for the directory to save the plot.
#' @param plotSuffix character. Suffix for plot filename.
#' @param scalingFactor numeric. Constant to multiply vectors by for visual clarity.
#' @param plotNoVectors logical. Set to TRUE to plot only the PCA without vectors. Default FALSE.
plotVectors <- function(sce,
                          colorVar,
                          plotLoadings = F,
                          loadingFactor = 3.5,
                          colorPalette = NA,
                          outputDir = NA,
                          plotSuffix = NA,
                          scalingFactor = NA,
                          plotNoVectors = F) {

  ## Create directory for output
  if(is.na(outputDir)) {
    topoName <- sce@metadata$topoName
    expName <- sce@metadata$experimentName
    outputDir <- file.path(getwd(),topoName,expName)
    print(paste0("Using default directory ",outputDir))
    if(!dir.exists(outputDir)) {
      dir.create(outputDir, recursive = T)
    }
  } else if(!dir.exists(outputDir)) {
    print("Specified directory does not exist - creating it..")
    dir.create(outputDir)
  }

  if(is.na(scalingFactor)) {
    scalingFactor <- sce@metadata$params$plotScalingFactor
  }

  if(is.na(plotSuffix)) {
    plotSuffix <- Sys.time()
  }
  fileName <- file.path(outputDir,paste0("PCA_vectors_",expName,"_",plotSuffix,".pdf"))

  # Get utility vars out
  expName <- sce@metadata$experimentName


  # Get 2D plotting coordinates of samples from first 2 cols of position matrix
  plot_df <- as.data.frame(colData(sce))
  if(!colorVar %in% colnames(plot_df)) {
    print("Error: colorVar not found in colnames of colData")
  }

  xMin <- sce@metadata$params$xMin
  xMax <- sce@metadata$params$xMax
  yMin <- sce@metadata$params$yMin
  yMax <- sce@metadata$params$yMax

  # Get proportion of variance for PCs 1 and 2 for axis labels
  pc1_weight <- round(100*sce@metadata$pca_summary$importance[2,1],2)
  pc2_weight <- round(100*sce@metadata$pca_summary$importance[2,2],2)
  plot_xlab <- paste("PC1 (",pc1_weight,"%)",sep="")
  plot_ylab <- paste("PC2 (",pc2_weight,"%)",sep="")

  # get loadings if applicable
  if(plotLoadings) {
    loadingDF <- as.data.frame(sce@metadata$pca_data$rotation)
    loadingDF$gene <- rownames(loadingDF)
    loadingDF$X <- 0
    loadingDF$Y <- 0
  }


  # Add relevant metadata to plotting df
  plot_df$colorVar <- plot_df[,colorVar]

  if(colorVar == "Cluster") {
    plot_df$colorVar <- factor(plot_df$colorVar)
  }



  image <- ggplot(plot_df, aes(x=X,y=Y)) +
    geom_point(mapping=aes(size=1, color=colorVar)) +
    xlab(plot_xlab) +
    ylab(plot_ylab) +
    scale_size(range=c(1.75, 3)) +
    guides(alpha="none", size="none", color=guide_legend(title = colorVar, override.aes = list(size = 5))) +
    xlim(xMin,xMax) +
    ylim(yMin,yMax) +
    ggtitle(paste0("PCA on ",expName," cells"))
  if(all(!check.numeric(plot_df$colorVar)) | colorVar == "Cluster") {
    if(is.na(colorPalette)) {
      colorPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")
    }
    image <- image + scale_color_manual(values=c(colorPalette))
  }
  if(plotLoadings) {
    loadingLabelFactor <- loadingFactor+0.5
    image <- image + geom_segment(data = loadingDF[,],aes(xend=PC1*loadingFactor, yend=PC2*loadingFactor, x = X, y = Y), arrow = arrow(length = unit(0.6,"cm"))) +
      geom_text(data = loadingDF[,], aes(x = PC1*loadingLabelFactor, y = PC2*loadingLabelFactor, label=gene), size=8, fontface="bold")
  }

  # Write to file without vectors
  if(plotNoVectors) {
    noVec_fName <- file.path(outputDir,paste0("PCA_",expName,"_",plotSuffix,".pdf"))
    pdf(file = noVec_fName, width = 10, height = 10)
    print(image)
    dev.off()
  }

  # Write to file with vectors
  pdf(file = fileName, width = 10, height = 10)
  print(image + geom_segment(aes(xend=X+dX*scalingFactor, yend=Y+dY*scalingFactor),
                             arrow = arrow(length = unit(0.2,"cm"))))
  dev.off()

}

#' @export
#' @title Plot smoothed velocities
#' @description Saves a pdf plot showing the first 2 PCs of the dataset with grid-smoothed
#' inferred velocities superimposed.
#' @return NULL
#' @import ggplot2
#' @importFrom varhandle check.numeric
#' @param sce viccc object after computing inferred velocities
#' @param colorVar character. Colname of colData(sce) by which to color the plot
#' @param plotLoadings logical. Whether to plot the loadings for each gene. Only recommended for
#' very small (<5 gene) circuits. Default FALSE.
#' @param loadingFactor numeric. Constant to multiply loadings by to improve visual clarity.
#' @param colorPalette vector. Vector containing colors - be sure to include at least as many colors
#' as there are unique values of colorVar. Default color palette includes 8 values.
#' @param outputDir character. File path for the directory to save the plot.
#' @param plotSuffix character. Suffix for plot filename.
#' @param scalingFactor numeric. Constant to multiply vectors by for visual clarity.
plotGrid <- function(sce,
                       colorVar,
                       plotLoadings = F,
                       loadingFactor = 3.5,
                       colorPalette = NA,
                       outputDir = NA,
                       plotSuffix = NA,
                       scalingFactor = NA) {

  ## Create directory for output
  if(is.na(outputDir)) {
    topoName <- sce@metadata$topoName
    expName <- sce@metadata$experimentName
    outputDir <- file.path(getwd(),topoName,expName)
    print(paste0("Using default directory ",outputDir))
    if(!dir.exists(outputDir)) {
      dir.create(outputDir, recursive = T)
    }
  } else if(!dir.exists(outputDir)) {
    print("Specified directory does not exist - creating it..")
    dir.create(outputDir)
  }

  if(is.na(scalingFactor)) {
    scalingFactor <- sce@metadata$params$gridPlotScalingFactor
  }

  if(is.na(plotSuffix)) {
    plotSuffix <- Sys.time()
  }
  fileName <- file.path(outputDir,paste0("PCA_grid_",expName,"_",plotSuffix,".pdf"))


  # Get 2D plotting coordinates of samples from first 2 cols of position matrix
  plot_df <- as.data.frame(colData(sce))
  if(!colorVar %in% colnames(plot_df)) {
    print("Error: colorVar not found in colnames of colData")
  }

  xMin <- sce@metadata$params$xMin
  xMax <- sce@metadata$params$xMax
  yMin <- sce@metadata$params$yMin
  yMax <- sce@metadata$params$yMax

  # Get proportion of variance for PCs 1 and 2 for axis labels
  pc1_weight <- round(100*sce@metadata$pca_summary$importance[2,1],2)
  pc2_weight <- round(100*sce@metadata$pca_summary$importance[2,2],2)
  plot_xlab <- paste("PC1 (",pc1_weight,"%)",sep="")
  plot_ylab <- paste("PC2 (",pc2_weight,"%)",sep="")

  # get loadings if applicable
  if(plotLoadings) {
    loadingDF <- as.data.frame(sce@metadata$pca_data$rotation)
    loadingDF$gene <- rownames(loadingDF)
    loadingDF$X <- 0
    loadingDF$Y <- 0
  }


  # Add relevant metadata to plotting df
  grid.df <- sce@metadata$grid.df
  x.points <- unique(sce@metadata$grid.df$x.points)
  y.points <- unique(sce@metadata$grid.df$y.points)
  plot_df$colorVar <- plot_df[,colorVar]

  if(colorVar == "Cluster") {
    plot_df$colorVar <- factor(plot_df$colorVar)
  }

  ## Plot grid points
  image <- ggplot(grid.df[which(grid.df$Magnitude > 0.01),], aes(x=x.points,y=y.points)) +
    geom_point(data=plot_df[,],mapping=aes(x=X,y=Y,size=1, color=colorVar)) +
    geom_point(mapping=aes(size=1)) +
    xlab(plot_xlab) +
    ylab(plot_ylab) +
    scale_size(range=c(1.0, 3)) +
    xlim(xMin,xMax) +
    ylim(yMin,yMax) +
    #scale_color_manual(values=c(cbPalette[2:8])) +
    geom_segment(aes(xend=x.points+dx, yend=y.points+dy), arrow = arrow(length = unit(0.2,"cm"))) +
    guides(alpha="none", size="none", color=guide_legend(title = colorVar, override.aes = list(size = 5))) +
    theme(axis.text = element_text(size=20), axis.title = element_text(size=28))
  if(all(!check.numeric(plot_df$colorVar)) | colorVar == "Cluster") {
    if(is.na(colorPalette)) {
      colorPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")
    }
    image <- image + scale_color_manual(values=c(colorPalette))
  }
  if(plotLoadings) {
    loadingLabelFactor <- loadingFactor+0.5
    image <- image + geom_segment(data = loadingDF[,],aes(xend=PC1*loadingFactor, yend=PC2*loadingFactor, x = X, y = Y), arrow = arrow(length = unit(0.6,"cm"))) +
      geom_text(data = loadingDF[,], aes(x = PC1*loadingLabelFactor, y = PC2*loadingLabelFactor, label=gene), size=8, fontface="bold")
  }


  ## Write to file
  pdf(file = fileName, width = 10, height = 10)
  print(image)
  dev.off()
}











