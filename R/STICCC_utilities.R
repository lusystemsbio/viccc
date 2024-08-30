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
sticSE <- function(topo,
                  exprMat = NA,
                  normData,
                  topoName,
                  expName = NA,
                  radius = 0.05,
                  minNeighbors = 5,
                  scalingFactor = 1,
                  gridScalingFactor = 1,
                  verbose = T,
                  useOriginalFeatures = F,
                  nPCs = NA
                  ) {
  if(!all(is.na(exprMat))) {
    stic <- SingleCellExperiment(assays = SimpleList(counts=exprMat, normcounts=normData))  
  } else {
    stic <- SingleCellExperiment(assays = SimpleList(normcounts=normData))
  }
  
  if(is.na(expName)) {
    stic@metadata$experimentName <- paste0(topoName,"_SCE")
  } else {
    stic@metadata$experimentName <- expName
  }
  
  if(is.na(nPCs)) {
    nPCs <- nrow(normData)
  }
  stic@metadata$topoName <- topoName
  stic@metadata$topo <- topo

  stic@metadata$params <- list(sample_radius=radius, plotScalingFactor=scalingFactor, gridPlotScalingFactor=gridScalingFactor,
                              minNeighbors=minNeighbors, verbose=verbose, useOriginalFeatures=useOriginalFeatures, nPCs=nPCs)
  return(stic)
}


#' @export
#' @title Create metadata table.
#' @description Returns a data.frame containing metadata for each cell in a dataset. This can be
#' substituted with a previously existing data.frame and is only necessary if no metadata table
#' exists. This will also divide the dataset into k clusters using the ward.d2 method.
#' @return data.frame with cells as rows and metadata characteristics as columns.
#' @param sce SingleCellExperiment object.
#' @param exprMat data.frame or matrix. Non-normalized expression matrix with rows as features
#' and columns as cells.
#' @param cluster logical. Whether to perform hierarchical clustering.
#' @param ... Additional arguments to the (internal) clustering function - k (number of clusters),
#' features (which features to cluster by, defaults to all), and clustfun (default "ward.D2")
prepMetadata <- function(sce,
                         exprMat,
                         cluster = T,
                         ...            # Parameters for getClusters function
) {
  metadata <- as.data.frame(t(exprMat))
  clusters <- getClusters(matrix = metadata, ...)
  metadata$SampleID <- rownames(metadata)
  metadata$Cluster <- clusters

  colData(sce) <- DataFrame(metadata)
  colnames(sce) <- colData(sce)$SampleID
  return(sce)

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
#' @return modified viccc object containing pca transformed data in reducedDim(obj, "PCA") and pca parameters in obj@metadata$pca_data.
#' @importFrom stats prcomp
#' @param sce viccc object.
#' @param save logical. Whether to save pca results to file
#' @param overwrite logical. If true, pca file will overwrite preivous pca if one exists. Default TRUE. 
#' Set this to FALSE if you intend to use a specific PCA projection for further analysis
#' @param fname character. Where to save pca results. Default will be an Rds file named PCA_res in the data 
#' subdirectory of the topology directory.
runPCA <- function(sce,
                   save=T,
                   overwrite=T,
                   fname=NA) {
  
  if(save && is.na(fname)) {
    fname <- file.path(getwd(),sce@metadata$topoName,"data","PCA_res.Rds")
  }

  exprMat <- assay(sce, "normcounts")
  pca <- prcomp(t(exprMat))
  
  if(save & (!file.exists(fname) | overwrite)) {
    saveRDS(pca, file=fname)
  }

  # add PCA to SCE object
  reducedDim(sce, "PCA") <- pca$x
  sce@metadata$pca_data <- pca[1:4]
  sce@metadata$pca_summary <- summary(pca)

  return(sce)
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
#' @param xMin numeric. Lower bound for plotting on x-axis.
#' @param xMax numeric. Upper bound for plotting on x-axis.
#' @param yMin numeric. Lower bound for plotting on y-axis.
#' @param yMax numeric. Upper bound for plotting on y-axis.
computeGrid <- function(sce,
                        grid.length = 30,
                        xMin = NA,
                        xMax = NA,
                        yMin = NA,
                        yMax = NA
) {
  # Get posMat
  posMat <- reducedDim(sce, "PCA")

  # Define grid points, create metadata structure for later
  ## Create uniform grid
  # Set grid area (2-dimensional?) w uniform points
  if(any(is.na(xMin))) {
    xMin <- floor(min(posMat[,1]))
    xMax <- ceiling(max(posMat[,1]))
    yMin <- floor(min(posMat[,2]))
    yMax <- ceiling(max(posMat[,2]))
  }
  x.min <- xMin
  x.max <- xMax
  y.min <- yMin
  y.max <- yMax

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
#' @param unitVectors logical. Whether to scale grid point vectors to a magnitude of 1. Default FALSE. 
#' @param inVectors logical. Whether to compute grid vectors using v2 predictions instead of v1. Default FALSE
#' @param combine logical. Whether to combine v1 and v2 for either net flow or reversibility in grid output. Default FALSE.
#' @param how character. Describes how to combine v1 and v2 for grid vectors. Must be one of "net", "rev", "v1+v2", or "v1-v2".
computeGridVectors <- function(sce,
                               scalingFactor = NA,
                               unitVectors = TRUE,
                               inVectors = FALSE,
                               combine = FALSE,
                               how = NA) {
  if(is.na(scalingFactor)) {
    scalingFactor <- sce@metadata$params$gridPlotScalingFactor
  }
  
  if(combine & is.na(how)) {
    print("Error: must specify how if combine == T")
    return(NULL)
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

    
    # Decide which vector (or sum) to be plotted
    if(inVectors) {
      subset_data <- plot_df[neighbors,c("X","Y","dX_in","dY_in")]
    } else {
      
      if(combine) {
        
        if(how == "v1+v2") {
          subset_data <- plot_df[neighbors,c("X","Y","dX","dY")]
          subset_data$dX = subset_data$dX + plot_df$dX_in[neighbors]
          subset_data$dY = subset_data$dY + plot_df$dY_in[neighbors]
        } else if(how == "v1-v2") {
          subset_data <- plot_df[neighbors,c("X","Y","dX","dY")]
          subset_data$dX = subset_data$dX - plot_df$dX_in[neighbors]
          subset_data$dY = subset_data$dY - plot_df$dY_in[neighbors]
        } else if(how == "net") {
          subset_data <- plot_df[neighbors,c("X","Y","dX","dY")]
          subset_data$dX = (subset_data$dX + plot_df$dX_in[neighbors]) / 2
          subset_data$dY = (subset_data$dY + plot_df$dY_in[neighbors]) / 2
        } else if(how == "rev") {
          subset_data <- plot_df[neighbors,c("X","Y","dX","dY")]
          subset_data$dX = (subset_data$dX - plot_df$dX_in[neighbors]) / 2
          subset_data$dY = (subset_data$dY - plot_df$dY_in[neighbors]) / 2
        } else {
          print("Error: must specify `how`, if `combine` is set to TRUE")
        }
        
        
      } else {
        subset_data <- plot_df[neighbors,c("X","Y","dX","dY")]
      }
      
    }

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

  if(unitVectors) {
    scaleR <- function(x) {x / sqrt(sum(x^2))}
    scaled_vecs <-  as.data.frame(t(apply(grid.df[,c("dx","dy")], 1, scaleR)))
    grid.df$dx = scaled_vecs[,1] * grid.dist * 0.7
    grid.df$dy = scaled_vecs[,2] * grid.dist * 0.7
    grid.df$Magnitude <- sqrt(grid.df$dx^2 + grid.df$dy^2)
    
  } else {
    ## Scale all vectors to the target mean?
    grid.df$Magnitude <- sqrt(grid.df$dx^2 + grid.df$dy^2)
    
    # mean.magnitude <- mean(grid.df$Magnitude[which(grid.df$Magnitude > 0)])
    # scaling.factor <- grid.dist / mean.magnitude
    # grid.df$dx <- grid.df$dx * scaling.factor * scalingFactor
    # grid.df$dy <- grid.df$dy * scaling.factor * scalingFactor
    # grid.df$Magnitude <- grid.df$Magnitude * scaling.factor * scalingFactor
  }
  

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



#' @title Distance function
#' @description Small utility function to compute dist of a given matrix
#' @return dist
#' @param x numeric input matrix.
distfun=function(x) {
  return(as.dist((1-cor(t(x), method = 's'))/20))
}



#' @title Convert list to numeric
#' @description Small utility function to ensure a list of data is in numeric format
#' @return list converted to numeric data type
#' @param x input list
listToNumeric <- function(x) {
  return(as.numeric(unname(unlist(x))))
}




#' @export
#' @title Plot single-cell velocities
#' @description Saves a pdf plot showing the first 2 PCs of the dataset with individual cells'
#' inferred velocities superimposed
#' @return NULL
#' @import ggplot2
#' @importFrom varhandle check.numeric
#' @param sce viccc object after computing inferred velocities
#' @param colorVar character. Colname of colData(sce) by which to color the plot, or NA to plot all grey
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
                          colorVar = NA,
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
    plotSuffix <- format(Sys.time(), "%b_%d")
  }
  fileName <- file.path(outputDir,paste0("PCA_vectors_",expName,"_",plotSuffix,".pdf"))

  # Get utility vars out
  expName <- sce@metadata$experimentName


  # Get 2D plotting coordinates of samples from first 2 cols of position matrix
  plot_df <- as.data.frame(colData(sce))
  if(!is.na(colorVar)) {
    if(!colorVar %in% colnames(plot_df)) {
      print("Error: colorVar not found in colnames of colData")
    }  
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
  
  ## Plot grid points
  if(!is.na(colorVar)) {
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
      ggtitle(paste0("PCA on ",expName," cells")) +
      theme(axis.text = element_text(size=28), axis.title = element_text(size=36)) 
    if(all(!check.numeric(plot_df$colorVar)) | colorVar == "Cluster") {
      if(is.na(colorPalette)) {
        colorPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")
      }
      image <- image + scale_color_manual(values=c(colorPalette))
    }
    
  } else {
    image <- ggplot(plot_df, aes(x=X,y=Y)) +
      geom_point(mapping=aes(size=1), color="grey", alpha=0.8) +
      xlab(plot_xlab) +
      ylab(plot_ylab) +
      scale_size(range=c(1.75, 3)) +
      xlim(xMin,xMax) +
      ylim(yMin,yMax) +
      ggtitle(paste0("PCA on ",expName," cells")) +
      guides(alpha="none", size="none", color="none") +
      theme(axis.text = element_text(size=28), axis.title = element_text(size=36)) 
    
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
#' @param colorVar character. Colname of colData(sce) by which to color the plot. May be left NA for grey plot
#' @param plotLoadings logical. Whether to plot the loadings for each gene. Only recommended for
#' very small (<5 gene) circuits. Default FALSE.
#' @param loadingFactor numeric. Constant to multiply loadings by to improve visual clarity.
#' @param colorPalette vector. Vector containing colors - be sure to include at least as many colors
#' as there are unique values of colorVar. Default color palette includes 8 values.
#' @param outputDir character. File path for the directory to save the plot.
#' @param plotSuffix character. Suffix for plot filename.
#' @param scalingFactor numeric. Constant to multiply vectors by for visual clarity.
#' @param minMagnitude numeric. Smallest magnitude for which a vector will be drawn. Default 0.01.
#' @param arrowSize numeric. Size of arrows in plot.
#' @param arrowheadSize numeric. Size of arrowheads in plot.
plotGrid <- function(sce,
                     colorVar = NA,
                     plotLoadings = F,
                     loadingFactor = 3.5,
                     colorPalette = NA,
                     outputDir = NA,
                     plotSuffix = NA,
                     scalingFactor = NA,
                     minMagnitude = 0.01,
                     arrowSize = 3,
                     arrowheadSize = 0.3) {

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
    plotSuffix <- format(Sys.time(), "%b_%d")
  }
  fileName <- file.path(outputDir,paste0("PCA_grid_",expName,"_",plotSuffix,".pdf"))


  # Get 2D plotting coordinates of samples from first 2 cols of position matrix
  plot_df <- as.data.frame(colData(sce))
  if(!is.na(colorVar)) {
    if(!colorVar %in% colnames(plot_df)) {
      print("Error: colorVar not found in colnames of colData")
    }  
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
  
  ## Plot grid points
  if(!is.na(colorVar)) {
    plot_df$colorVar <- plot_df[,colorVar]
    
    if(colorVar == "Cluster") {
      plot_df$colorVar <- factor(plot_df$colorVar)
    }
    
    image <- ggplot(grid.df[which(grid.df$Magnitude > minMagnitude),], aes(x=x.points,y=y.points)) +
      geom_point(data=plot_df[,],mapping=aes(x=X,y=Y,size=1, color=colorVar)) +
      geom_point(mapping=aes(size=1)) +
      xlab(plot_xlab) +
      ylab(plot_ylab) +
      scale_size(range=c(1.0, 3)) +
      xlim(xMin,xMax) +
      ylim(yMin,yMax) +
      geom_segment(aes(xend=x.points+dx*scalingFactor, yend=y.points+dy*scalingFactor), 
                   arrow = arrow(length = unit(arrowheadSize,"cm")), size=arrowSize) +
      guides(alpha="none", size="none", color=guide_legend(title = colorVar, override.aes = list(size = 5))) +
      theme(axis.text = element_text(size=28), axis.title = element_text(size=36))
    if(all(!check.numeric(plot_df$colorVar)) | colorVar == "Cluster") {
      if(is.na(colorPalette)) {
        colorPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")
      }
      image <- image + scale_color_manual(values=c(colorPalette))
    }
  } else {
    image <- ggplot(grid.df[which(grid.df$Magnitude > minMagnitude),], aes(x=x.points,y=y.points)) +
      geom_point(data=plot_df[,],mapping=aes(x=X,y=Y,size=1),color="grey",alpha=0.8) + 
      geom_point(mapping=aes(size=1)) + 
      xlab(plot_xlab) +
      ylab(plot_ylab) +
      scale_size(range=c(1.0, 3)) +
      xlim(xMin,xMax) +
      ylim(yMin,yMax) +
      geom_segment(aes(xend=x.points+dx*scalingFactor, yend=y.points+dy*scalingFactor), 
                   arrow = arrow(length = unit(arrowheadSize,"cm")), size=arrowSize) +
      guides(alpha="none", size="none", color="none") +
      theme(axis.text = element_text(size=28), axis.title = element_text(size=36)) 
    
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







#' @export
#' @title Compute single-cell transition vector
#' @description Calculates predicted transition vectors v1 and v2 for a query point.
#' @return named list of results containing X and Y (position on first two dimensions), boolean values for logging purposes, 
#' and predicted vectors.
#' @importFrom MASS ginv
#' @importFrom plyr revalue
#' @importFrom stats cor
#' @param sce sticcc object 
#' @param query_point character. Rowname of posMat indicating which sample to compute a vector for.
#' @param useGinv logical. Whether to use generalized inverse for computing least-squares regression. 
#' Default FALSE, but may be helpful with very homogeneous datasets.
#' @param v2 logical. If TRUE, will calculate v2 (incoming transition) as well as v1 (outgoing). Default TRUE.
#' @param invertV2 logical. If TRUE, will multiply v2 by -1 such that it points "forward in time" as v1. Default FALSE. Note this step
#' may be necessary for downstream analysis. 
#' @param maxNeighbors numeric. Maximum number of neighbors to include in linear regression - if there are more neighbors in range 
#' than this limit, this value will be used to randomly sample a subset.
computeVector <- function(sce, query_point, useGinv=F, v2=T, invertV2=F, maxNeighbors = 100) {
  ## Prepare results
  rs_list <- list("X"=NA,"Y"=NA,"fewNeighbors"=FALSE,"numNeighbors"=NA,"usedGinv"=FALSE,"nonInvertible"=FALSE,"selfCorNA"=FALSE,"selfCor"=NA,"b_vec"=NA,"b_vec_in"=NA,"det"=NA)
  
  ## extract variables from sce
  # Create sign vector
  topo <- sce@metadata$topo
  sign_vector <- revalue(factor(topo$Type), c('1'='1','2'='-1','3'='1','4'='-1','5'='1','6'='-1'), warn_missing = FALSE)
  sign_vector <- as.numeric(levels(sign_vector))[sign_vector]
  
  # Compute sampling radius
  sampling_radius <- sce@metadata$params$sample_radius * sce@metadata$max_dist
  
  # get posMat (position matrix - this should be either identical or closely correlated w/ the data used to compute dist_mat)
  if(is.null(sce@metadata$params$nPCs)) {
    nPCs <- nrow(sce)
    sce@metadata$params$nPCs <- nPCs
    
  } else {
    nPCs <- sce@metadata$params$nPCs
  }

  # use first nPCs PCs as spatial coordinates unless genes are specified
  pcaMat <- reducedDim(sce, "PCA")[colnames(sce),1:nPCs]
  if(sce@metadata$params$useOriginalFeatures) {
    posMat <- as.data.frame(t(assay(sce, "normcounts")))[colnames(sce),]
  } else {
    posMat <- pcaMat
  }
  
  posList <- colnames(posMat)
  rownames(posMat) <- colnames(sce)
  
  # get exprMat (normalized expression matrix - high-dimensional data (rows=features, cols=samples))
  exprMat <- assay(sce,"normcounts")
  
  # Get query point data
  query_data <- posMat[query_point,]
  rs_list[["X"]] <- query_data[1]
  rs_list[["Y"]] <- query_data[2]
  
  ## Find neighbors
  neighbors <- which(sce@metadata$dist_mat[query_point,colnames(sce)] <= sampling_radius)
  #neighbors <- neighbors[which(neighbors %in% colnames(sce))]
  
  ## Skip sample if number of neighbors too small
  if(length(neighbors) <= (sce@metadata$params$minNeighbors+1)) {
    rs_list[["fewNeighbors"]] <- TRUE
    rs_list[["numNeighbors"]] <- length(neighbors)-1
    return(rs_list)
  }
  
  ## Select neighbors
  subset_models <- as.data.frame(posMat[neighbors,])
  subset_models <- subset_models[-which(rownames(subset_models) %in% query_point),]
  rs_list[["numNeighbors"]] <- nrow(subset_models)
  
  ## Subset neighbors if over max
  if(nrow(subset_models) > maxNeighbors) {
    neighborSubset <- sample(rownames(subset_models), maxNeighbors)
    subset_models <- subset_models[neighborSubset,]
  }
  
  ## Multiple regression
  # Create relative positional matrix of neighbor samples
  subset_models[,1:ncol(subset_models)] <- apply(subset_models[,1:ncol(subset_models)], 2, listToNumeric)
  deriv_df <- sweep(subset_models, 2, as.numeric(query_data), "-")
  
  # get relative expression matrix & transpose (neighboring cells only)
  geneex_mat <- as.matrix(deriv_df[,posList])
  geneex_mat_transposed <- t(geneex_mat)
  
  # Multiply relative positional matrix by its transpose
  mat_product <- geneex_mat_transposed %*% geneex_mat
  rs_list[["det"]] <- det(mat_product)
  
  # Take the inverse if possible - if not, skip this sample
  if(useGinv) {
    tol = 10^-6
    if(rcond(mat_product) < tol){
      inverse_mat = ginv(mat_product)
      rs_list[["usedGinv"]] <- TRUE
      
    } else{
      inverse_mat <- tryCatch({
        solve(mat_product)
      }, error = function(e) {
        "non-invertible"
      })
    }
  } else {
    inverse_mat <- tryCatch({
      solve(mat_product)
    }, error = function(e) {
      "non-invertible"
    })
  }
  
  
  if(inverse_mat[1] == "non-invertible") {
    rs_list[["nonInvertible"]] <- TRUE
    return(rs_list)
  }
  
  # Multiply the inverted matrix by the transposed positional matrix
  x_mat <- inverse_mat %*% geneex_mat_transposed
  
  
  ## Calculate delayed correlation
  # Self correlation - skip this sample if NA
  src_activity <- as.numeric(exprMat[topo$Source,query_point] * sign_vector)
  self_target_activity <- as.numeric(exprMat[topo$Target,query_point]) 
  self_cor <- cor(self_target_activity, src_activity)
  if(is.na(self_cor)) {
    rs_list[["selfCorNA"]] <- TRUE
    return(rs_list)
  } 
  rs_list[["selfCor"]] <- self_cor
  
  # For each neighbor:
  deriv_df$DelayedCorr <- NA
  for(model in rownames(deriv_df)) {
    # Correlate neighbor target activity with query source activity * sign_vector
    neighbor_target_activity <- as.numeric(exprMat[topo$Target, model]) 
    if(sd(neighbor_target_activity) > 0 & sd(src_activity)) {
      neighbor_cor <- cor(neighbor_target_activity,src_activity) 
      deriv_df[model,"DelayedCorr"] <- (neighbor_cor - self_cor)
    } else {
      #print(paste("Issue with correlation in models ",query_point, " and ", model))
      deriv_df[model,"DelayedCorr"] <- NA
    }
    
    ## QQ REMOVE THIS
    #plot(src_activity, neighbor_target_activity)
  }
  
  # Remove NA values caused by zero variance vectors
  na_models <- which(is.na(deriv_df$DelayedCorr))
  if(length(na_models) > 0) {
    deriv_df <- deriv_df[-na_models,]
    subset_models <- subset_models[-na_models,]
    x_mat <- x_mat[,-na_models]
  }
  
  
  # Compute product of earlier matrix product and delayed corr vector
  corr_vec <- as.matrix(deriv_df[,"DelayedCorr"])
  b_vec <- x_mat %*% corr_vec
  rs_list[["b_vec"]] <- b_vec
  
  ### Bias didn't seem to make a big difference
  # if(bias) {
  #   ## Split data into train/test 80:20
  #   train_idx <- sample(nrow(deriv_df), round(nrow(deriv_df)*0.8))
  #   df.train <- deriv_df[train_idx,]
  #   df.test <- deriv_df[-train_idx,]
  #   
  #   alpha.guess <- 0
  #   beta.guess <- rep(0,nPCs)
  #   
  #   h <- 100
  #   error <- rep(NA,h)
  #   lambda <- seq(0,1,length.out = h)
  #   
  #   fit <- lm(DelayedCorr~.-1, data=df.train)
  #   yhat.test <- predict(fit, newdata = df.test)
  #   yhat.guess <- alpha.guess + as.matrix(df.test[,-ncol(df.test)]) %*% beta.guess
  #   
  #   error <- sapply(lambda, function(j) sum(((1-j)*yhat.test + j*yhat.guess - df.test$DelayedCorr)^2))
  #   rs_list[["error"]] <- error
  # }
  # 
  
  
  ## Compute v2 as well, if specified
  if(v2) {
    # For each neighbor:
    deriv_df$DelayedCorr_in <- NA
    for(model in rownames(deriv_df)) {
      # Correlate neighbor source activity * sign_vector with query target activity 
      neighbor_src_activity <- as.numeric(exprMat[topo$Source, model]) * sign_vector
      neighbor_cor <- cor(neighbor_src_activity,self_target_activity)
      deriv_df[model,"DelayedCorr_in"] <- (neighbor_cor - self_cor)
    }
    
    # Remove NA values caused by zero variance vectors
    na_models <- which(is.na(deriv_df$DelayedCorr))
    if(length(na_models) > 0) {
      deriv_df <- deriv_df[-na_models,]
      subset_models <- subset_models[-na_models,]
      x_mat <- x_mat[,-na_models]
    }
    
    
    # Compute product of earlier matrix product and delayed corr vector
    corr_vec <- as.matrix(deriv_df[,"DelayedCorr_in"])
    b_vec_in <- x_mat %*% corr_vec
    
    if(invertV2) {
      b_vec_in <- -1*b_vec_in
    }
    
    rs_list[["b_vec_in"]] <- b_vec_in
    
  }
  
  return(rs_list)
}




#' @export
#' @title Plot network to file
#' @description Saves an image of the network topology
#' @return NULL
#' @import visNetwork
#' @param sce sticcc object 
#' @param outputDir character. Directory to save plot. Default is a subdirectory of the working 
#' directory with the same name as sce@metadata$topoName, which will be created if it does not exist
#' @param plotSuffix character. String to append to plot filename. Default is the current time. 
plotNetwork <- function(sce, outputDir=NA, plot_suffix=NA) {
  if(is.na(outputDir)) {
    outputDir <- file.path(getwd(),sce@metadata$topoName)
  }
  if(!dir.exists(outputDir)) {
    dir.create(outputDir)
  }
  
  if(is.na(plot_suffix)) {
    plot_suffix <- Sys.time()
  }
  
  
  net_file <- file.path(outputDir,paste0("networkVis",plot_suffix,".html"))
  
  topo <- sce@metadata$topo
  
  topo[which(as.numeric(topo$Type) %% 2 == 0),"Type"] <- 2
  topo[which(as.numeric(topo$Type) %% 2 == 1),"Type"] <- 1
  
  
  node_list <-
    unique(c(topo[, 1], topo[, 2]))
  
  nodes <-
    data.frame(
      id = node_list,
      label = node_list,
      font.size = 50,
      value = c(rep(1, length(node_list)))
    )
  edge_col <- data.frame(c(1, 2), c("blue", "darkred"))
  arrow_type <- data.frame(c(1, 2), c("arrow", "circle"))
  colnames(arrow_type) <- c("type", "color")
  colnames(edge_col) <- c("type", "color")
  edges <-
    data.frame(
      from = c(topo[, 1]),
      to = c(topo[, 2]),
      arrows.to.type	= arrow_type$color[c(as.numeric(topo[, 3]))],
      width = 3,
      color = edge_col$color[c(as.numeric(topo[, 3]))]
    )
  
  networkPlot <-
    visNetwork::visNetwork(nodes, edges, height = "1000px", width = "100%") %>%
    visEdges(arrows = "to") %>%
    visOptions(manipulation = FALSE) %>%
    visLayout(randomSeed = 123) %>%
    #visNodes(scaling = list(label = list(enabled = T))) %>%
    visPhysics(solver = "forceAtlas2Based", stabilization = FALSE)
  visNetwork::visSave(networkPlot, file = net_file, selfcontained = FALSE)
  
}




#' @export
#' @title Compute a smoothed vector around a specific sample or gene expression state
#' @description Finds neighboring samples around the supplied point, computes vectors for each, and 
#' returns a distance-weighted average representing the smoothed vector for the region as if it 
#' were a grid point from computeGridVectors
#' @return A list containing the following elements: "start": the query point, "v1_x": v1 prediction
#' for the first PC, "v1_y": v1 prediction for the second PC, "v2_x", "v2_y", "net_x", "net_y", "rev_x",
#' "rev_y" denoting the other predictions in PC1 and PC2, "neighborhoodSize": number of nearby samples
#' used to compute the , "smoothingNumber": 
#' @param sce sticcc object 
#' @param queryPoint character or vector. Either a rowname of sce, or a vector to be compared to 
#' the normalized expression data.
#' @param neighborhoodRadius numeric. Proportion of the maximum pairwise distance within sce to find
#' neighbors to compute vectors for.
smoothVector <- function(sce,
                         queryPoint, # 
                         neighborhoodRadius,
                         ...) {
  
  
  
  neighbors <- getNeighbors(sce,
                            queryPoint,
                            neighborhoodRadius)
  
  neighbor.ids <- neighbors[["neighbor.ids"]]
  neighbor.dists <- neighbors[["neighbor.dists"]]
  
  # check that there are sufficient neighbors
  if(length(neighbor.ids) < 5) {
    print(paste0("Error: not enough cells within the given neighborhoodRadius, only found ",length(neighbor.ids)))
    return(NULL)
  }
  
  
  # Compute vectors for each neighbor
  sce@metadata$params$sample_radius <- neighborhoodRadius
  neighborhoodVectors <- data.frame(SampleID = neighbor.ids, Dist=neighbor.dists,
                                    V1_x=NA,  V1_y=NA, 
                                    V2_x=NA, V2_y=NA, 
                                    Net_x=NA, Net_y=NA, 
                                    Rev_x=NA, Rev_y=NA)
  for(neighbor in neighbor.ids) {
    # compute vector
    v_list <- computeVector(sce, 
                              neighbor,
                              v2 = T,
                              ...)
    
    # add to aggregated results
    neighborhoodVectors[which(neighborhoodVectors$SampleID == neighbor), "V1_x"] <- v_list[["b_vec"]][1,]
    neighborhoodVectors[which(neighborhoodVectors$SampleID == neighbor), "V1_y"] <- v_list[["b_vec"]][2,]
    neighborhoodVectors[which(neighborhoodVectors$SampleID == neighbor), "V2_x"] <- v_list[["b_vec_in"]][1,]
    neighborhoodVectors[which(neighborhoodVectors$SampleID == neighbor), "V2_y"] <- v_list[["b_vec_in"]][2,]
    
    neighborhoodVectors[which(neighborhoodVectors$SampleID == neighbor), "Net_x"] <- (v_list[["b_vec"]][1,] + v_list[["b_vec_in"]][1,]) / 2
    neighborhoodVectors[which(neighborhoodVectors$SampleID == neighbor), "Net_y"] <- (v_list[["b_vec"]][2,] + v_list[["b_vec_in"]][2,]) / 2
    neighborhoodVectors[which(neighborhoodVectors$SampleID == neighbor), "Rev_x"] <- (v_list[["b_vec"]][1,] - v_list[["b_vec_in"]][1,]) / 2
    neighborhoodVectors[which(neighborhoodVectors$SampleID == neighbor), "Rev_y"] <- (v_list[["b_vec"]][2,] - v_list[["b_vec_in"]][2,]) / 2
  }
  
  
  # subset neighborhood based on smoothing radius
  # check that there are sufficient neighbors in the smoothing radius
  smoothingDist <- neighborhoodRadius * sce@metadata$max_dist
  neighborhoodVectors <- neighborhoodVectors[which(neighborhoodVectors$Dist < smoothingDist),]
  
  if(nrow(neighborhoodVectors) < 3) {
    print("Warning: not enough vectors to smooth - consider increasing smoothingRadius")
  }
  
  # compute IDW-weighted smoothed vectors
  neighborhoodVectors$Proximity <- 1 / neighborhoodVectors$Dist
  
  avg_dx_v1 <- Avg_IDW(neighborhoodVectors$V1_x, neighborhoodVectors$Proximity)
  avg_dy_v1 <- Avg_IDW(neighborhoodVectors$V1_y, neighborhoodVectors$Proximity)
  
  avg_dx_v2 <- Avg_IDW(neighborhoodVectors$V2_x, neighborhoodVectors$Proximity)
  avg_dy_v2 <- Avg_IDW(neighborhoodVectors$V2_y, neighborhoodVectors$Proximity)
  
  avg_dx_net <- Avg_IDW(neighborhoodVectors$Net_x, neighborhoodVectors$Proximity)
  avg_dy_net <- Avg_IDW(neighborhoodVectors$Net_y, neighborhoodVectors$Proximity)
  
  avg_dx_rev <- Avg_IDW(neighborhoodVectors$Rev_x, neighborhoodVectors$Proximity)
  avg_dy_rev <- Avg_IDW(neighborhoodVectors$Rev_y, neighborhoodVectors$Proximity)
  
  
  # return output vectors
  out_list <- list(start=queryPoint,
                   v1_x=avg_dx_v1,
                   v1_y=avg_dy_v1,
                   v2_x=avg_dx_v2,
                   v2_y=avg_dy_v2,
                   net_x=avg_dx_net,
                   net_y=avg_dy_net,
                   rev_x=avg_dx_rev,
                   rev_y=avg_dy_rev,
                   neighborhoodSize=length(neighbor.ids),
                   smoothingNumber=nrow(neighborhoodVectors))
  
  
  return(out_list)
}


#' @export
#' @title Get nearby samples for a given query sample.
#' @description Finds neighboring samples around the supplied point.
#' @return A list containing the following elements: "neighbor.ids": the names of the nearby samples,
#' "neighbor.dists": distances from each neighbor to the query point.
#' @param sce sticcc object 
#' @param queryPoint character or vector. Either a rowname of sce, or a vector to be compared to 
#' the normalized expression data.
#' @param neighborhoodRadius numeric. Proportion of the maximum pairwise distance within sce to find
#' neighbors.
getNeighbors <- function(sce,
                         queryPoint,
                         neighborhoodRadius) {
  # Find neighborhood members
  if(all(queryPoint %in% colnames(sce))) {
    queryData <- assay(sce)[,queryPoint]
  } else if(length(queryPoint) != nrow(assay(sce)))  {
    print(paste0("Error: query should be a colname of sce or have same rows as sce assay: ",nrow(assay(sce))))
    return(NULL)
  } else {
    queryData <- queryPoint
  }
  
  max_dist <- sce@metadata$max_dist
  lim_dist <- max_dist * neighborhoodRadius
  est_k <- ncol(sce)*neighborhoodRadius*2
  neighbors <- FNN::get.knnx(data=reducedDim(sce,"PCA")[,1:sce@metadata$params$nPCs], query=queryData, k=est_k)
  
  # filter for those below neighborhoodRadius
  neighborIDs <- neighbors$nn.index[which(neighbors$nn.dist <= lim_dist)]
  neighborDists <- neighbors$nn.dist[which(neighbors$nn.dist <= lim_dist)]
  
  
  out_list <- list(neighbor.ids=colnames(sce)[neighborIDs],
                   neighbor.dists=neighborDists)
  return(out_list)
  
}




