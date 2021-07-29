
suppressPackageStartupMessages({
  library(ggplot2)
  library(plyr)
  library(dplyr)
  library(tidyr)
  library(FNN)
  library(ComplexHeatmap) # remove?
  library(RColorBrewer)
  library(ggfortify)
  library(sRACIPE) # remove? 
  library(visNetwork)
  library(varhandle)
  library(SingleCellExperiment)
})


## TODO: fix "-" by simply removing hyphens from topology files
## Simulate using RACIPE - 10k models, default iNoise, 10 noise levels, step size 0.01 and time 100 
simTopo <- function(topo,
                    ...) {
  
  ## Import topology
  topo$Source <- gsub("-","",topo$Source)
  topo$Target <- gsub("-","",topo$Target)
  nGenes <- length(unique(c(topo$Source, topo$Target)))
  #genelist <- unique(c(topo$Source, topo$Target))
  
  ## Set seed for reproducibility
  set.seed(123)
  
  ## Simulate topology (once)
  racipe <- RacipeSE() # Construct an empty RacipeSE object
  sracipeCircuit(racipe) <- topo
  
  racipe <- sracipeSimulate(racipe, numModels = 10000, ...)
  return(racipe)
  #saveRDS(racipe, file = file.path(outputDir,paste0("simData_",topoName,".Rds")))
  ## TODO: Check for extreme low expression before normalizing in some topologies 
  #saveRDS(sracipeNormalize(racipe), file = file.path(outputDir,paste0("simData_",topoName,"_NORM.Rds")))
}



# Load topo in as dataframe
loadTopo <- function(topoName,
                     topoDir=NA) {
  ## Import topology
  if(is.na(topoDir)) {
    topoDir <- file.path(getwd(),"input_topos")
  }
  
  if(!dir.exists(topoDir)) {
    dir.create(topoDir)
  }
  topo <- read.table(file.path(topoDir,paste0(topoName,".tpo")), header = T)
  topo$Source <- gsub("-","",topo$Source)
  topo$Target <- gsub("-","",topo$Target)
  return(topo)
}



runPCA <- function(exprMat) {
  return(prcomp(t(exprMat)))
}






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



computeDist <- function(sce,
                        usePCs=TRUE,
                        numPCs=10) {
  
  if(usePCs) {
    PCAmat <- reducedDim(sce,"PCA")
    sce@metadata$dist_mat <- as.matrix(dist(PCAmat[,1:min(numPCs, ncol(PCAmat))], method = "euclidean"))
    sce@metadata$max_dist <- max(sce@metadata$dist_mat)
  } else {
    print("Error: this function is in development. Try using usePCs=TRUE for now")
  }
  
  return(sce)
}



DCComputeTrajectorySCE <- function(sce
) {
  
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
  posMat <- reducedDim(sce, "PCA")[colnames(sce),1:nPCs]
  posList <- colnames(posMat)
  rownames(posMat) <- colnames(sce)
  
  # get exprMat (normalized expression matrix - high-dimensional data (rows=features, cols=samples))
  exprMat <- assay(sce,"normcounts")
  
  # Begin velocity calculation
  if(sce@metadata$params$verbose) {
    print("Computing inferred velocity")
  }
  
  sample_list <- colData(sce)$SampleID
  colData(sce)$numNeighbors <- NA
  colData(sce)$selfCor <- NA
  colData(sce)$X <- NA
  colData(sce)$Y <- NA
  colData(sce)$dX <- NA
  colData(sce)$dY <- NA
  colData(sce)$fewNeighbors <- FALSE
  colData(sce)$nonInvertible <- FALSE
  colData(sce)$selfCorNA <- FALSE
  
  weighted_vector_list <- vector(mode = "list", length = length(sample_list))
  pb = txtProgressBar(min = 0, max = length(sample_list), 
                      initial = 0, style = 3)
  stepi <- 0
  for(query_point in sample_list) {
    ## Update progress bar
    stepi <- stepi+1
    setTxtProgressBar(pb, stepi)
    
    ## Get query point data
    query_data <- posMat[query_point,]
    colData(sce)[query_point,"X"] <- query_data[1]
    colData(sce)[query_point,"Y"] <- query_data[2]
    
    ## Find neighbors
    neighbors <- which(sce@metadata$dist_mat[query_point,colnames(sce)] <= sampling_radius)
    #neighbors <- neighbors[which(neighbors %in% colnames(sce))]
    
    ## Skip sample if number of neighbors too small
    if(length(neighbors) <= (sce@metadata$params$minNeighbors+1)) {
      colData(sce)[query_point,"fewNeighbors"] <- TRUE
      colData(sce)[query_point,"numNeighbors"] <- length(neighbors)-1
      next
    }
    
    ## Select neighbors
    subset_models <- as.data.frame(posMat[neighbors,])
    subset_models <- subset_models[-which(rownames(subset_models) %in% query_point),]
    colData(sce)[query_point,"numNeighbors"] <- nrow(subset_models)
    
    ## Multiple regression
    # Create relative positional matrix of neighbor samples
    subset_models[,1:ncol(subset_models)] <- apply(subset_models[,1:ncol(subset_models)], 2, listToNumeric)
    deriv_df <- sweep(subset_models, 2, as.numeric(query_data), "-")
    
    # get relative expression matrix & transpose (neighboring cells only)
    geneex_mat <- as.matrix(deriv_df[,posList])
    geneex_mat_transposed <- t(geneex_mat)
    
    # Multiply relative positional matrix by its transpose
    mat_product <- geneex_mat_transposed %*% geneex_mat
    
    # Take the inverse if possible - if not, skip this sample
    inverse_mat <- tryCatch({
      solve(mat_product)
    }, error = function(e) {
      "non-invertible"
    })
    
    if(inverse_mat[1] == "non-invertible") {
      colData(sce)[query_point,"nonInvertible"] <- TRUE
      next
    }
    
    # Multiply the inverted matrix by the transposed positional matrix
    x_mat <- inverse_mat %*% geneex_mat_transposed
    
    
    ## Calculate delayed correlation
    # Self correlation - skip this sample if NA
    src_activity <- as.numeric(exprMat[topo$Source,query_point] * sign_vector)
    self_target_activity <- as.numeric(exprMat[topo$Target,query_point]) 
    self_cor <- cor(self_target_activity, src_activity)
    if(is.na(self_cor)) {
      colData(sce)[query_point,"selfCorNA"] <- TRUE
      next
    } 
    colData(sce)[query_point,"selfCor"] <- self_cor
    
    # For each neighbor:
    deriv_df$DelayedCorr <- NA
    for(model in rownames(deriv_df)) {
      # Correlate neighbor target activity with query source activity * sign_vector
      neighbor_target_activity <- as.numeric(exprMat[topo$Target, model]) 
      neighbor_cor <- cor(neighbor_target_activity,src_activity)
      deriv_df[model,"DelayedCorr"] <- (neighbor_cor - self_cor)
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
    
    ## Save vector to master list
    saved_vector <- as.data.frame(t(b_vec))
    colnames(saved_vector) <- paste0("d",posList)
    weighted_vector_list[[stepi]] <- saved_vector
    names(weighted_vector_list)[[stepi]] <- query_point
    
    ## Update colData
    colData(sce)[query_point,"dX"] <- saved_vector[1]
    colData(sce)[query_point,"dY"] <- saved_vector[2]
    
  }
  
  ## Consolidate vectors and write to file
  weighted_vectors <- do.call(rbind, weighted_vector_list)
  sce@metadata$vectors <- weighted_vectors
  
  return(sce)
  
}





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
    #visEdges(arrows = "to") %>%
    visOptions(manipulation = FALSE) %>%
    visLayout(randomSeed = 123) %>%
    #visNodes(scaling = list(label = list(enabled = T))) %>%
    visPhysics(solver = "forceAtlas2Based", stabilization = FALSE)
  visNetwork::visSave(networkPlot, file = net_file, selfcontained = FALSE)
  
}


DCPlotVectors <- function(sce, 
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
    #labs(colorVar = as.character(pVar)) +
    scale_size(range=c(1.75, 3)) +
    guides(alpha=FALSE, size=FALSE, color=guide_legend(title = colorVar, override.aes = list(size = 5))) +
    #scale_color_manual(values=c(cbPalette[2:8])) +
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


DCPlotGrid <- function(sce,
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
    guides(alpha=FALSE, size=FALSE, color=guide_legend(title = colorVar, override.aes = list(size = 5))) +
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


DCComputeGridVectors <- function(sce,
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
    
    ## Find neighbors
    query <- grid.df[which(grid.df$GridPoint == point),c(1:2)]
    colnames(query) <- c("X","Y")
    neighbors <- which(plot_df$X <= (query$X + x.dist) & 
                         plot_df$X >= (query$X - x.dist) &
                         plot_df$Y <= (query$Y + y.dist) & 
                         plot_df$Y >= (query$Y - y.dist))
    
    
    subset_data <- plot_df[neighbors,c("X","Y","dX","dY")]
    subset_data <- na.omit(subset_data)
    grid.df[point,"numNeighbors"] <- nrow(subset_data)
    
    if(nrow(subset_data) >= 3) {
      
    } else if(nrow(subset_data) < 3) {
      ## Add to grid.df as 0
      grid.df[which(grid.df$GridPoint == point),"dx"] <- 0
      grid.df[which(grid.df$GridPoint == point),"dy"] <- 0
      next
    }
    
    ## Calculate distance & proximity (1/dist) to each neighbor
    neighbor.dists <- as.matrix(dist(rbind(query,subset_data[,c(1:2)]), 
                                     method = "euclidean"))
    neighbor.dists <- neighbor.dists[which(rownames(neighbor.dists) == point),-which(colnames(neighbor.dists) == point)]
    subset_data$Dist <- neighbor.dists
    subset_data$Proximity <- 1 / subset_data$Dist

    
    
    ## Calculate weighted vector (closer = more weight)
    avg_dx <- Avg_IDW(subset_data$dX, subset_data$Proximity)
    avg_dy <- Avg_IDW(subset_data$dY, subset_data$Proximity)
    #avg_dx <- sum(subset_data$dX * subset_data$Proximity) / sum(subset_data$Proximity)
    #avg_dy <- sum(subset_data$dY * subset_data$Proximity) / sum(subset_data$Proximity)
    
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

Avg_IDW <- function(vals, weights) {
  
  return(sum(vals * weights) / sum(weights))
  # avg_dx <- sum(subset_data$dX * subset_data$Proximity) / sum(subset_data$Proximity)
  # avg_dy <- sum(subset_data$dY * subset_data$Proximity) / sum(subset_data$Proximity)
}

DCPLotLocal <- function(sce,
                        id,
                        outputDir = NA,
                        plotSuffix = NA) {
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
  
  

  if(is.na(plotSuffix)) {
    plotSuffix <- Sys.time()
  }
  fileName <- file.path(outputDir,paste0("localPlot_model",id,"_",plotSuffix,".pdf"))
  
  # Get utility vars out
  expName <- sce@metadata$experimentName
  
  # Create sign vector
  topo <- sce@metadata$topo
  sign_vector <- revalue(factor(topo$Type), c('1'='1','2'='-1','3'='1','4'='-1','5'='1','6'='-1'), warn_missing = FALSE)
  sign_vector <- as.numeric(levels(sign_vector))[sign_vector]
  
  # Compute sampling radius
  sampling_radius <- sce@metadata$params$sample_radius * sce@metadata$max_dist
  
  # get posMat (position matrix - this should be either identical or closely correlated w/ the data used to compute dist_mat)
  posMat <- reducedDim(sce, "PCA")[colnames(sce),]
  posList <- colnames(posMat)
  rownames(posMat) <- colnames(sce)
  
  # get exprMat (normalized expression matrix - high-dimensional data (rows=features, cols=samples))
  exprMat <- assay(sce,"normcounts")
  
  
  ## Get query point data
  if(!id %in% rownames(posMat)) {
    print("Error: id not found")
  }
  query_data <- posMat[id,]
  
  ## Find neighbors
  neighbors <- which(sce@metadata$dist_mat[id,colnames(sce)] <= sampling_radius)
  #neighbors <- neighbors[which(neighbors %in% colnames(sce))]
  
  ## Skip sample if number of neighbors too small
  if(length(neighbors) <= (sce@metadata$params$minNeighbors+1)) {
    next
  }
  
  ## Select neighbors
  subset_models <- as.data.frame(posMat[neighbors,])
  subset_models <- subset_models[-which(rownames(subset_models) %in% id),]
  
  ## Multiple regression
  # Create relative positional matrix of neighbor samples
  deriv_df <- sweep(subset_models, 2, query_data, "-")
  
  # get relative expression matrix & transpose (neighboring cells only)
  geneex_mat <- as.matrix(deriv_df[,posList])
  geneex_mat_transposed <- t(geneex_mat)
  
  # Multiply relative positional matrix by its transpose
  mat_product <- geneex_mat_transposed %*% geneex_mat
  
  # Take the inverse if possible - if not, skip this sample
  inverse_mat <- tryCatch({
    solve(mat_product)
  }, error = function(e) {
    "non-invertible"
  })
  
  if(inverse_mat[1] == "non-invertible") {
    next
  }
  
  # Multiply the inverted matrix by the transposed positional matrix
  x_mat <- inverse_mat %*% geneex_mat_transposed
  
  
  ## Calculate delayed correlation
  # Self correlation - skip this sample if NA
  src_activity <- as.numeric(exprMat[topo$Source,id] * sign_vector)
  self_target_activity <- as.numeric(exprMat[topo$Target,id]) 
  self_cor <- cor(self_target_activity, src_activity)
  if(is.na(self_cor)) {
    next
  } 
  
  # For each neighbor:
  deriv_df$DelayedCorr <- NA
  for(model in rownames(deriv_df)) {
    # Correlate neighbor target activity with query source activity * sign_vector
    neighbor_target_activity <- as.numeric(exprMat[topo$Target, model]) 
    neighbor_cor <- cor(neighbor_target_activity,src_activity)
    deriv_df[model,"DelayedCorr"] <- (neighbor_cor - self_cor)
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
  
  ## Save vector to master list
  saved_vector <- as.data.frame(t(b_vec))
  colnames(saved_vector) <- paste0("d",posList)
  
  
  ## Plot local region
  plot_df <- subset_models
  plot_df$DC <- deriv_df$DelayedCorr
  center <- data.frame(X=query_data[1],Y=query_data[2],dX=saved_vector[1],dY=saved_vector[2])
  colnames(center) <- c("X","Y","dX","dY")
  
  
  image <- ggplot(plot_df) +
    geom_point(aes(x=PC1,y=PC2,color=DC), size=15) +
    scale_color_gradient() +
    geom_point(data=center,aes(x=X,y=Y),color="red",size=15) +
    guides(alpha=FALSE, size=FALSE, color=FALSE) +
    theme(panel.background = element_rect(fill = "transparent"),
          panel.border = element_rect(color = "black", fill=NA, size=3),
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank())
  print(image)
  
  
  
  print(image +
          geom_segment(data=center, aes(x=X,y=Y,xend=X+dX, yend=Y+dY), size=4,
                       arrow = arrow(length = unit(1.5,"cm"))))
  
  
  
}





DCPlotVectors_Minimal <- function(sce, 
                          color,
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
  fileName <- file.path(outputDir,paste0("PCA_vectors_",expName,"_",plotSuffix,"_minimal.pdf"))
  
  # Get utility vars out
  expName <- sce@metadata$experimentName
  
  
  # Get 2D plotting coordinates of samples from first 2 cols of position matrix
  plot_df <- as.data.frame(colData(sce))
  # if(length(color) == 1 | length(color) != nrow(plot_df)) {
  #   plot_df$Color <-  rep(color[1], nrow(plot_df))
  # } else {
  #   plot_df$Color <- color
  # }
  
  xMin <- sce@metadata$params$xMin
  xMax <- sce@metadata$params$xMax
  yMin <- sce@metadata$params$yMin
  yMax <- sce@metadata$params$yMax
  
  
  # Add relevant metadata to plotting df
  
  image <- ggplot(plot_df, aes(x=X,y=Y)) +
    geom_point(size=15, color=color) + 
    #xlab(plot_xlab) +
    #ylab(plot_ylab) +
    #labs(colorVar = as.character(pVar)) +
    #scale_size(range=c(1.75, 3)) +
    guides(alpha=FALSE, size=FALSE, color=FALSE) +
    #scale_color_manual(values=c(cbPalette[2:8])) +
    xlim(xMin+2,xMax-2) +
    ylim(yMin+2.2,yMax-2) +
    theme(panel.background = element_rect(fill = "transparent"),
          panel.border = element_rect(color = "black", fill=NA, size=3),
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank())
  
  # Write to file without vectors
  if(plotNoVectors) {
    noVec_fName <- file.path(outputDir,paste0("PCA_",expName,"_",plotSuffix,"_minimal.pdf"))
    pdf(file = noVec_fName, width = 10, height = 10)
    print(image)
    dev.off()
  }
  
  # Write to file with vectors
  pdf(file = fileName, width = 10, height = 10)
  print(image + geom_segment(aes(xend=X+dX*scalingFactor, yend=Y+dY*scalingFactor), size=2,
                             arrow = arrow(length = unit(0.4,"cm"))))
  dev.off()
  
}


DCPlotVectors_Minimal_colVector <- function(sce, 
                                  colors,
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
  fileName <- file.path(outputDir,paste0("PCA_vectors_",expName,"_",plotSuffix,"_minimal.pdf"))
  
  # Get utility vars out
  expName <- sce@metadata$experimentName
  
  
  # Get 2D plotting coordinates of samples from first 2 cols of position matrix
  plot_df <- as.data.frame(colData(sce))
  plot_df$Color <- colors
  # if(length(color) == 1 | length(color) != nrow(plot_df)) {
  #   plot_df$Color <-  rep(color[1], nrow(plot_df))
  # } else {
  #   plot_df$Color <- color
  # }
  
  xMin <- sce@metadata$params$xMin
  xMax <- sce@metadata$params$xMax
  yMin <- sce@metadata$params$yMin
  yMax <- sce@metadata$params$yMax
  
  
  # Add relevant metadata to plotting df
  
  image <- ggplot(plot_df, aes(x=X,y=Y)) +
    geom_point(size=15, aes(color=Color)) + 
    #xlab(plot_xlab) +
    #ylab(plot_ylab) +
    #labs(colorVar = as.character(pVar)) +
    #scale_size(range=c(1.75, 3)) +
    guides(alpha=FALSE, size=FALSE, color=FALSE) +
    scale_color_manual(values=unique(colors)) +
    xlim(xMin+2,xMax-2) +
    ylim(yMin+2.2,yMax-2) +
    theme(panel.background = element_rect(fill = "transparent"),
          panel.border = element_rect(color = "black", fill=NA, size=3),
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank())
  
  # Write to file without vectors
  if(plotNoVectors) {
    noVec_fName <- file.path(outputDir,paste0("PCA_",expName,"_",plotSuffix,"_minimal.pdf"))
    pdf(file = noVec_fName, width = 10, height = 10)
    print(image)
    dev.off()
  }
  
  # Write to file with vectors
  pdf(file = fileName, width = 10, height = 10)
  print(image + geom_segment(aes(xend=X+dX*scalingFactor, yend=Y+dY*scalingFactor), size=2,
                             arrow = arrow(length = unit(0.4,"cm"))))
  dev.off()
  
}



listToNumeric <- function(x) {
  return(as.numeric(unname(unlist(x))))
}

processDataTS <- function(racipe,model) {
  tsData_model <- t(racipe@metadata$timeSeries)
  tsData_model <- data.frame(tsData_model, Time=row.names(tsData_model), Model=model)
  nCols <- ncol(tsData_model)-2
  
  tsData_model[,1:nCols] <- apply(tsData_model[,1:nCols], 2, listToNumeric)
  
  tsData_model[,1:nCols] <- log2(1+tsData_model[,1:nCols]) # Log transform
  tsData_model[,1:nCols] <- sweep(tsData_model[,1:nCols], 2, tmpMeans, FUN = "-") # scale
  tsData_model[,1:nCols] <- sweep(tsData_model[,1:nCols], 2, tmpSds, FUN = "/") # scale
  newPca <- (scale(tsData_model[,1:nCols], pca$center, pca$scale) %*%
               pca$rotation)
  
  tsData_model$PC1 <- newPca[,1]
  tsData_model$PC2 <- newPca[,2]
  
  return(tsData_model)
}




DCPlotGrid_Grey <- function(sce,
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
  image <- ggplot(grid.df[which(grid.df$Magnitude > 0.01),], aes(x=x.points,y=y.points)) +
    geom_point(data=plot_df[,],mapping=aes(x=X,y=Y,size=1),color="grey",alpha=0.8) + 
    geom_point(mapping=aes(size=1)) + 
    xlab(plot_xlab) +
    ylab(plot_ylab) +
    scale_size(range=c(1.0, 3)) +
    xlim(xMin,xMax) +
    ylim(yMin,yMax) +
    #scale_color_manual(values=c(cbPalette[2:8])) +
    geom_segment(aes(xend=x.points+dx, yend=y.points+dy), arrow = arrow(length = unit(0.2,"cm"))) +
    guides(alpha=FALSE, size=FALSE, color=FALSE) +
    theme(axis.text = element_text(size=20), axis.title = element_text(size=28)) 

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



DCPlotGrid_Grey_colVectors <- function(sce,
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
  image <- ggplot(grid.df[which(grid.df$Magnitude > 0.01),], aes(x=x.points,y=y.points)) +
    geom_point(data=plot_df[,],mapping=aes(x=X,y=Y,size=1),color="grey",alpha=0.8) + 
    geom_point(mapping=aes(size=1,color=Diff_Magnitude)) + 
    xlab(plot_xlab) +
    ylab(plot_ylab) +
    scale_size(range=c(1.0, 3)) +
    xlim(xMin,xMax) +
    ylim(yMin,yMax) +
    #scale_color_manual(values=c(cbPalette[2:8])) +
    geom_segment(aes(xend=x.points+dx, yend=y.points+dy, color=Diff_Magnitude), arrow = arrow(length = unit(0.3,"cm")),size=1.25) +
    scale_color_gradient2(low="white",high="red", midpoint=-8) +
    guides(alpha=FALSE, size=FALSE) +
    theme(axis.text = element_text(size=20), axis.title = element_text(size=28)) 
  
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




compareVectorFieldsV2 <- function(grid1, 
                                  grid2,
                                  dimRed,
                                  output_dir = NA,
                                  comparisonName = NA
) {
  
  # create name if none given
  if(is.na(comparisonName)) {
    comparisonName <- paste0("comparison_",Sys.time())
  }
  
  ## Set up directory structure
  ##TODO: remove this, it's hardcoded to make my life easier but breaks things for other users - set expDir to something reasonable by default instead
  if(is.na(output_dir)) {
    output_dir <- file.path(getwd(),"VectorFieldComparisons",comparisonName)
  }

  # Create directory if needed
  if(!dir.exists(output_dir)) {
    dir.create(output_dir)
  } 
  
  
  ## Plot subtracted df by plotVars
  # Find points with non-zero vectors for both experiments
  nonzero.exp1 <- which(exp1.grid.df$Magnitude >= 0.01)
  nonzero.exp2 <- which(exp2.grid.df$Magnitude >= 0.01)
  shared.pts <- nonzero.exp1[which(nonzero.exp1 %in% nonzero.exp2)]
  
  # Subtract relevant values - show exp2 - exp1
  diff.grid <- exp2.grid.df
  diff.grid$dx <- diff.grid$dx - exp1.grid.df$dx
  diff.grid$dy <- diff.grid$dy - exp1.grid.df$dy
  
  diff.grid <- diff.grid[shared.pts,]
  
  # Only include non-zero resulting vectors
  nonzero.after <- which(sqrt(diff.grid$dx^2 + diff.grid$dy^2) > 0.01)
  diff.grid.placeholder <- diff.grid
  diff.grid <- diff.grid[nonzero.after,]
  
  ## Summary statistics on diff.grid
  median.magnitude <- median(diff.grid$Magnitude)
  mean.magnitude <- mean(diff.grid$Magnitude)
  sd.magnitude <- sd(diff.grid$Magnitude)
  summary_list <- list(Median.Magnitude=median.magnitude, Mean.Magnitude=mean.magnitude, SD.Magnitude=sd.magnitude)
  #saveRDS(summary_list, file = file.path(data_dir,paste0(expName1,"_minus_",expName2,"_summary_list.Rds")))
  
  ## Re-scale grid
  grid.dist <- sqrt((unique(diff.grid$x.points)[2] - unique(diff.grid$x.points)[1])^2 +
                      (unique(diff.grid$y.points)[2] - unique(diff.grid$y.points)[1])^2) / 2
  mean.magnitude <- mean(diff.grid$Magnitude[which(diff.grid$Magnitude > 0)])
  scaling.factor <- grid.dist / mean.magnitude
  
  diff.grid$dx <- diff.grid$dx * scaling.factor
  diff.grid$dy <- diff.grid$dy * scaling.factor
  diff.grid$Magnitude <- diff.grid$Magnitude * scaling.factor
  
  
  ## return diff.grid
  return(diff.grid)
  
}






compareVectorFields <- function(expName1, 
                                expName2,
                                metadata,
                                dimRed = NA,
                                plotVars,
                                dimRedType = "PCA",
                                expDir = NA,
                                output_dir = NA,
                                gridPlotScalingFactor = 1,
                                xMin = NA,
                                xMax = NA,
                                yMin = NA,
                                yMax = NA
) {
  
  ## Set up directory structure
  ##TODO: remove this, it's hardcoded to make my life easier but breaks things for other users - set expDir to something reasonable by default instead
  if(is.na(expDir)) {
    expDir <- getwd()
  }
  
  
  if(is.na(output_dir)) {
    output_dir <- file.path(expDir,"VectorFieldComparisons")
  }
  plot_dir <- file.path(output_dir, paste0(expName1,"_minus_",expName2))
  data_dir <- file.path(output_dir, "data")
  # Create directories if needed
  if(!dir.exists(expDir)) {
    dir.create(expDir)
  }
  if(!dir.exists(output_dir)) {
    dir.create(output_dir)
    dir.create(plot_dir)
    dir.create(data_dir)
  } else {
    if(!dir.exists(plot_dir)) {
      dir.create(plot_dir)
    }
    if(!dir.exists(data_dir)) {
      dir.create(data_dir)
    }
  }
  
  
  ## Read in dimRed object if not present already
  ## TODO
  
  ## Read in vector field data for both experiments
  exp1.grid.df <- readRDS(file = file.path(expDir,expName1,"data",paste0(expName1,"_vectorField.Rds")))
  exp2.grid.df <- readRDS(file = file.path(expDir,expName2,"data",paste0(expName2,"_vectorField.Rds")))
  
  ## If dimRed not supplied, pull from exp1
  if(is.na(dimRed)) {
    dimRed <- readRDS(file.path(expDir,expName1,"data",paste0("dimRed_",dimRedType,".Rds")))
  }
  
  ## Summary stats: mean/median ddx, ddy, dMagnitude
  ##TODO
  
  
  ## Plot subtracted df by plotVars
  # Find points with non-zero vectors for both experiments
  nonzero.exp1 <- which(exp1.grid.df$Magnitude >= 0.01)
  nonzero.exp2 <- which(exp2.grid.df$Magnitude >= 0.01)
  shared.pts <- nonzero.exp1[which(nonzero.exp1 %in% nonzero.exp2)]
  
  # Subtract relevant values - show exp2 - exp1
  diff.grid <- exp2.grid.df
  diff.grid$dx <- diff.grid$dx - exp1.grid.df$dx
  diff.grid$dy <- diff.grid$dy - exp1.grid.df$dy
  
  diff.grid <- diff.grid[shared.pts,]
  
  # Only include non-zero resulting vectors
  nonzero.after <- which(sqrt(diff.grid$dx^2 + diff.grid$dy^2) > 0.01)
  diff.grid.placeholder <- diff.grid
  diff.grid <- diff.grid[nonzero.after,]
  
  ## Summary statistics on diff.grid
  median.magnitude <- median(diff.grid$Magnitude)
  mean.magnitude <- mean(diff.grid$Magnitude)
  sd.magnitude <- sd(diff.grid$Magnitude)
  summary_list <- list(Median.Magnitude=median.magnitude, Mean.Magnitude=mean.magnitude, SD.Magnitude=sd.magnitude)
  saveRDS(summary_list, file = file.path(data_dir,paste0(expName1,"_minus_",expName2,"_summary_list.Rds")))
  
  ## Re-scale grid
  grid.dist <- sqrt((unique(diff.grid$x.points)[2] - unique(diff.grid$x.points)[1])^2 +
                      (unique(diff.grid$y.points)[2] - unique(diff.grid$y.points)[1])^2) / 2
  mean.magnitude <- mean(diff.grid$Magnitude[which(diff.grid$Magnitude > 0)])
  scaling.factor <- grid.dist / mean.magnitude
  
  diff.grid$dx <- diff.grid$dx * scaling.factor
  diff.grid$dy <- diff.grid$dy * scaling.factor
  diff.grid$Magnitude <- diff.grid$Magnitude * scaling.factor
  
  
  ## Prep for plotting
  if(dimRedType == "PCA") {
    # Make plotting df
    plot_df <- posMat[,c(1:2)] ## used to be as.data.frame(dimRed$x)[,c(1:2)]
    plot_df$UNIQ_ID <- rownames(plot_df)
    colnames(plot_df) <- c("X","Y","UNIQ_ID")
    
    # Get PCA info
    pc1 <- 1
    pc2 <- 2
    pc1_weight <- round((100*summary(dimRed)$importance[2,pc1]), digits = 2)
    pc2_weight <- round((100*summary(dimRed)$importance[2,pc2]), digits = 2)
    
    # Get proportion of variance for PCs 1 and 2 for axis labels
    plot_xlab <- paste("PC",pc1," (",pc1_weight,"%)",sep="")
    plot_ylab <- paste("PC",pc2," (",pc2_weight,"%)",sep="")
    
  } else if(dimRedType == "UMAP") {
    ## TODO
    print("Error: plotting UMAP needs to be fixed in the compareVectorFields function")
    print("Fix: instantiate plot_df from UMAP object dimRed, fix column names as well and $UNIQ_ID column")
    
    # Set axis labels
    plot_xlab <- "UMAP1"
    plot_ylab <- "UMAP2"
  }
  
  if(is.na(xMin) | is.na(xMax) | is.na(yMin) | is.na(yMax)) {
    print("Using default x and y lims as not all were provided")
    xMin <- min(c(exp1.grid.df$x.points, exp2.grid.df$x.points))
    xMax <- max(c(exp1.grid.df$x.points, exp2.grid.df$x.points))
    yMin <- min(c(exp1.grid.df$y.points, exp2.grid.df$y.points))
    yMax <- max(c(exp1.grid.df$y.points, exp2.grid.df$y.points))
  }
  
  ## Plot grid points w/ real data
  for(pVar in plotVars) {
    # Add relevant metadata to plotting df
    plot_df$colorVar <- metadata[which(metadata$SampleID %in% rownames(plot_df)),pVar]
    
    image <- ggplot(diff.grid[,], aes(x=x.points,y=y.points)) +
      geom_point(data=plot_df[,],mapping=aes(x=X,y=Y,size=1, color=colorVar), alpha=0.5) + 
      geom_point(mapping=aes(size=1)) + 
      xlab(plot_xlab) +
      ylab(plot_ylab) +
      scale_size(range=c(1.0, 3)) +
      xlim(xMin,xMax) +
      ylim(yMin,yMax) +
      geom_segment(aes(xend=x.points+dx, yend=y.points+dy), arrow = arrow(length = unit(0.08,"cm"))) +
      guides(alpha=FALSE, size=FALSE, color=guide_legend(title = pVar, override.aes = list(size = 5))) +
      ggtitle(paste0(expName2, " minus ", expName1))
    
    # Write to file
    fileName <- file.path(plot_dir, paste0(expName1,"_minus_",expName2,"_",dimRedType,"_by",pVar,".pdf"))
    pdf(file = fileName, width = 10, height = 10)
    print(image)
    dev.off()
  }
  
  ## Save diff.grid
  saveRDS(diff.grid, file = file.path(data_dir,paste0(expName1,"_minus_",expName2,"_diff.grid.Rds")))
  
  
  
  
}




distfun=function(x) {
  as.dist((1-cor(t(x), method = 's'))/20)
}




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






prepMetadata <- function(exprMat,       # expression matrix - rows = features, cols = samples
                         cluster = T,   # Logical - whether to perform clustering
                         ...            # Parameters for getClusters function
) {
  metadata <- as.data.frame(t(exprMat))
  clusters <- getClusters(matrix = metadata, ...)
  metadata$SampleID <- rownames(metadata)
  metadata$Cluster <- clusters
  
  return(metadata)
  
}


DCComputeTrajectoryFast <- function(exprMat,
                                    posMat,
                                    topo,
                                    distMat,
                                    expName,
                                    gridTemplate,
                                    outputDir,
                                    plotVars,
                                    metadata,
                                    plotScalingFactor = 1,
                                    gridPlotScalingFactor = 1,
                                    minNeighbors = 5,
                                    radius = 0.1,
                                    verbose = T,
                                    plot = T,
                                    plotAvgTimepoint = T
) {
  
  ## test: transpose exprMat
  exprMat <- as.data.frame(t(exprMat))
  
  
  # Create experiment dir
  expDir <- file.path(outputDir,expName)
  if(!dir.exists(expDir)) {
    dir.create(expDir)
  }
  
  # Save network as visNetwork
  net_file <- file.path(expDir,"networkVis.html")
  
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
    #visEdges(arrows = "to") %>%
    visOptions(manipulation = FALSE) %>%
    visLayout(randomSeed = 123) %>%
    #visNodes(scaling = list(label = list(enabled = T))) %>%
    visPhysics(solver = "forceAtlas2Based", stabilization = FALSE)
  visNetwork::visSave(networkPlot, file = net_file, selfcontained = FALSE)
  
  
  # Create sign vector
  sign_vector <- revalue(factor(topo$Type), c('1'='1','2'='-1','3'='1','4'='-1','5'='1','6'='-1'))
  sign_vector <- as.numeric(levels(sign_vector))[sign_vector]
  
  # Compute sampling radius
  max_dist <- max(distMat)
  sampling_radius <- radius * max_dist
  
  # Begin velocity calculation
  if(verbose) {
    print("Computing inferred velocity")
  }
  posMat$UNIQ_ID <- rownames(posMat)
  sample_list <- posMat$UNIQ_ID
  #query_point <- sample(sample_list, 1)
  output_df <- data.frame(row.names = posMat$UNIQ_ID, UNIQ_ID = posMat$UNIQ_ID,X=NA,Y=NA,dX=NA,dY=NA,selfCor=NA,numNeighbors=NA,
                          minLocalDC=NA,quart1LocalDC=NA,medianLocalDC=NA,meanLocalDC=NA,quart3LocalDC=NA,
                          maxLocalDC=NA)
  weighted_vector_list <- list()
  pb = txtProgressBar(min = 0, max = length(sample_list), 
                      initial = 0, style = 3)
  stepi <- 0
  for(query_point in sample_list) {
    ## Update progress bar
    stepi <- stepi+1
    setTxtProgressBar(pb, stepi)
    
    ## Get query point data
    query_data <- as.data.frame(posMat[query_point,])
    rownames(query_data) <- query_point
    posList <- colnames(posMat)[which(!colnames(posMat) == "UNIQ_ID")]
    
    ## Find neighbors
    query_distances <- distMat[which(rownames(distMat) == query_point),]
    neighbors <- names(query_distances[which(query_distances <= sampling_radius)])
    
    ## Select neighbors
    subset_models <- as.data.frame(posMat[neighbors,])
    subset_models$UNIQ_ID <- rownames(subset_models)
    subset_models <- subset_models[which(subset_models$UNIQ_ID != query_point),]
    
    ## Skip sample if number of neighbors too small
    if(nrow(subset_models) <= minNeighbors) {
      next
    }
    
    ## Multiple regression
    # Create relative positional matrix of neighbor samples
    deriv_df <- subset_models
    for(pos in posList) {
      deriv_df[,pos] <- deriv_df[,pos] - query_data[,pos]
    }
    ## qq what happens if points overlap perfectly?
    ## qq maybe try a slight rotation of PCA, maybe PCA on a subset
    ##  to detect NAs: apply(is.na(expl_matrix1), 2, which)
    #browser()
    geneex_mat <- as.matrix(deriv_df[,posList])
    geneex_mat_transposed <- t(geneex_mat)
    
    # Multiply relative positional matrix by its transpose
    mat_product <- geneex_mat_transposed %*% geneex_mat
    
    # Take the inverse if possible - if not, skip this sample
    inverse_mat <- tryCatch({
      solve(mat_product)
    }, error = function(e) {
      "non-invertible"
    })
    
    if(inverse_mat[1] == "non-invertible") {
      next
    }
    
    # Multiply the inverted matrix by the transposed positional matrix
    x_mat <- inverse_mat %*% geneex_mat_transposed
    
    
    ## Calculate delayed correlation
    # Self correlation - skip this sample if NA
    src_activity <- as.numeric(exprMat[query_point,topo$Source] * sign_vector)
    self_target_activity <- as.numeric(exprMat[query_point,topo$Target]) 
    self_cor <- cor(self_target_activity, src_activity)
    if(is.na(self_cor)) {
      next
    } 
    
    # For each neighbor:
    for(model in deriv_df$UNIQ_ID) {
      # Correlate neighbor target activity with query source activity * sign_vector
      neighbor_target_activity <- as.numeric(exprMat[model, topo$Target]) 
      neighbor_cor <- cor(neighbor_target_activity,src_activity)
      deriv_df[which(deriv_df$UNIQ_ID == model),"DelayedCorr"] <- (neighbor_cor - self_cor)
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
    
    ## Save vector to master list
    saved_vector <- as.data.frame(t(b_vec))
    colnames(saved_vector) <- paste0("d",posList)
    saved_vector$UNIQ_ID <- query_point
    weighted_vector_list[[(length(weighted_vector_list) + 1)]] <- saved_vector
    
    
    ## Update output df
    #output_df <- data.frame(row.names = posMat$UNIQ_ID, UNIQ_ID = posMat$UNIQ_ID,X=NA,Y=NA,dX=NA,dY=NA,selfCor=NA,numNeighbors=NA,
    #                        minLocalDC=NA,quart1LocalDC=NA,medianLocalDC=NA,meanLocalDC=NA,quart3LocalDC=NA,
    #                        maxLocalDC=NA)
    output_df[query_point,"X"] <- query_data[,1]
    output_df[query_point,"Y"] <- query_data[,2]
    output_df[query_point,"dX"] <- saved_vector[,1]
    output_df[query_point,"dY"] <- saved_vector[,2]
    output_df[query_point,"selfCor"] <- self_cor
    output_df[query_point,"numNeighbors"] <- nrow(subset_models)
    localDCStat <- unname(quantile(deriv_df$DelayedCorr))
    output_df[query_point,"minLocalDC"] <- localDCStat[1]
    output_df[query_point,"quart1LocalDC"] <- localDCStat[2]
    output_df[query_point,"medianLocalDC"] <- localDCStat[3]
    output_df[query_point,"quart3LocalDC"] <- localDCStat[4]
    output_df[query_point,"maxLocalDC"] <- localDCStat[5]
    output_df[query_point,"mean"] <- mean(deriv_df$DelayedCorr)
    
  }
  print("Done - plotting and wrapping up")
  
  ## Consolidate vectors and write to file
  weighted_vectors <- do.call(rbind, weighted_vector_list)
  rownames(weighted_vectors) <- weighted_vectors$MODEL_NO
  colnames(weighted_vectors)[1:length(posList)] <- paste0("d",posList)
  saveRDS(weighted_vectors, file = file.path(expDir, paste0("weighted_vectors_",expName,".Rds")))
  saveRDS(output_df, file = file.path(expDir, paste0("output_df_",expName,".Rds")))
  saveRDS(topo, file = file.path(expDir, paste0("topo_",expName,".Rds")))
  
  
  
  ## Plot all vectors together
  if(plot) {
    # Get 2D plotting coordinates of samples from first 2 cols of position matrix
    plot_df <- posMat[weighted_vectors$UNIQ_ID,]
    colnames(plot_df)[1:2] <- c("X","Y")
    xMin <- floor(min(posMat[,1]))
    xMax <- ceiling(max(posMat[,1]))
    yMin <- floor(min(posMat[,2]))
    yMax <- ceiling(max(posMat[,2]))
    
    # Add vector endpoints for each sample to plot df
    plot_df$dX <- weighted_vectors$dPC1
    plot_df$dY <- weighted_vectors$dPC2
    
    # Get proportion of variance for PCs 1 and 2 for axis labels
    #plot_xlab <- paste("PC",pc1," (",pc1_weight,"%)",sep="")
    #plot_ylab <- paste("PC",pc2," (",pc2_weight,"%)",sep="")
    plot_xlab <- "PC1"
    ploy_ylab <- "PC2"
    
    for(pVar in plotVars) {
      # Add relevant metadata to plotting df
      plot_df$colorVar <- metadata[which(metadata$SampleID %in% rownames(plot_df)),pVar]
      
      if(pVar == "Cluster") {
        plot_df$colorVar <- factor(plot_df$colorVar)
      }
      
      image <- ggplot(plot_df, aes(x=X,y=Y)) +
        geom_point(mapping=aes(size=1, color=colorVar)) + 
        #xlab(plot_xlab) +
        #ylab(plot_ylab) +
        #labs(colorVar = as.character(pVar)) +
        scale_size(range=c(1.75, 3)) +
        guides(alpha=FALSE, size=FALSE, color=guide_legend(title = pVar, override.aes = list(size = 5))) +
        #scale_color_manual(values=c(cbPalette[2:8])) +
        xlim(xMin,xMax) +
        ylim(yMin,yMax) +
        ggtitle(paste0("PCA on ",expName," cells"))
      
      # Write to file without vectors
      fileName <- file.path(expDir, paste0("PCAPlot_",expName,"_by",pVar,".pdf"))
      pdf(file = fileName, width = 10, height = 10)
      print(image)
      dev.off()
      
      # Write to file with vectors
      fileName <- file.path(expDir, paste0("PCAPlotVectors_",expName,"_by",pVar,".pdf"))
      pdf(file = fileName, width = 10, height = 10)
      print(image + geom_segment(aes(xend=X+dX*plotScalingFactor, yend=Y+dY*plotScalingFactor), 
                                 arrow = arrow(length = unit(0.08,"cm"))))
      dev.off()
    }
  }
  
  
  
  ## Create uniform grid
  # Set grid area (2-dimensional?) w uniform points
  grid.df <- gridTemplate
  x.points <- unique(gridTemplate$x.points)
  y.points <- unique(gridTemplate$y.points)
  x.dist <- x.points[2] - x.points[1]
  y.dist <- y.points[2] - y.points[1]
  grid.dist <- sqrt((x.points[2] - x.points[1])^2 + (y.points[2] - y.points[1])^2) / 2
  arrow.max <- grid.dist * 0.5
  
  
  ## Iterate thru each gridpoint & calculate weighted k-means average vector
  pca_vectors <- plot_df
  #pb = txtProgressBar(min = 0, max = length(grid.df$GridPoint), 
  #                    initial = 0, style = 3)
  stepi <- 0
  for(point in grid.df$GridPoint) {
    stepi <- stepi+1
    #setTxtProgressBar(pb, stepi)
    ## Find neighbors
    query <- grid.df[which(grid.df$GridPoint == point),c(1:2)]
    colnames(query) <- c("X","Y")
    #knn_list <- get.knnx(data=pca_vectors[,c(1:2)], query=query, k=(nNeighbors))
    
    neighbors <- rownames(pca_vectors[which(pca_vectors$X <= (query$X + x.dist) & 
                                              pca_vectors$X >= (query$X - x.dist) &
                                              pca_vectors$Y <= (query$Y + y.dist) & 
                                              pca_vectors$Y >= (query$Y - y.dist)),])
    
    
    if(length(neighbors) >= 3) {
      #subset_data <- pca_vectors[knn_list$nn.index,c("PC1","PC2","dPC1","dPC2","Time","Model","UniqueID")]
      subset_data <- pca_vectors[neighbors,c("X","Y","dX","dY")]
      
      if(plotAvgTimepoint) {
        avgTimepoints <- mean(as.numeric(metadata[neighbors,"Timepoint"]))
      }
      
    } else if(length(neighbors) < 3) {
      ## Add to grid.df
      grid.df[which(grid.df$GridPoint == point),"dx"] <- 0
      grid.df[which(grid.df$GridPoint == point),"dy"] <- 0
      next
    }
    
    
    ## Calculate distance & proximity (1/dist) to each neighbor
    neighbor.dists <- as.matrix(dist(rbind(query,pca_vectors[neighbors,c(1:2)]), 
                                     method = "euclidean"))
    neighbor.dists <- neighbor.dists[which(rownames(neighbor.dists) == point),-which(colnames(neighbor.dists) == point)]
    subset_data$Dist <- neighbor.dists
    subset_data$Proximity <- 1 / subset_data$Dist
    
    ## Calculate weighted average vector (closer = more weight)
    avg_dx <- sum(subset_data$dX * subset_data$Proximity) / sum(subset_data$Proximity)
    avg_dy <- sum(subset_data$dY * subset_data$Proximity) / sum(subset_data$Proximity)
    
    ## Add to grid.df
    grid.df[which(grid.df$GridPoint == point),"dx"] <- avg_dx
    grid.df[which(grid.df$GridPoint == point),"dy"] <- avg_dy
    
  }
  
  ## Scale all vectors to the fixed maximum (fixed magnitude for all?)
  grid.df$Magnitude <- sqrt(grid.df$dx^2 + grid.df$dy^2)
  mean.magnitude <- mean(grid.df$Magnitude[which(grid.df$Magnitude > 0)])
  scaling.factor <- grid.dist / mean.magnitude
  
  grid.df$dx <- grid.df$dx * scaling.factor * gridPlotScalingFactor
  grid.df$dy <- grid.df$dy * scaling.factor * gridPlotScalingFactor
  grid.df$Magnitude <- grid.df$Magnitude * scaling.factor * gridPlotScalingFactor
  
  
  ## Write to file
  saveRDS(grid.df, file = file.path(expDir, paste0("PCAVectorField_",expName,".Rds")))
  
  for(pVar in plotVars) {
    # Add relevant metadata to plotting df
    pca_vectors$colorVar <- metadata[which(metadata$SampleID %in% rownames(pca_vectors)),pVar]
    
    ## Plot grid points
    image <- ggplot(grid.df[which(grid.df$Magnitude > 0.01),], aes(x=x.points,y=y.points)) +
      geom_point(data=pca_vectors[,],mapping=aes(x=X,y=Y,size=1, color=colorVar)) + 
      geom_point(mapping=aes(size=1)) + 
      #xlab(plot_xlab) +
      #ylab(plot_ylab) +
      scale_size(range=c(1.0, 3)) +
      xlim(min(x.points),max(x.points)) +
      ylim(min(y.points),max(y.points)) +
      #scale_color_manual(values=c(cbPalette[2:8])) +
      geom_segment(aes(xend=x.points+dx, yend=y.points+dy), arrow = arrow(length = unit(0.08,"cm"))) +
      guides(alpha=FALSE, size=FALSE, color=guide_legend(title = pVar, override.aes = list(size = 5))) +
      ggtitle(paste0("PCA on ",expName," cells"))
    
    ## Write to file
    fileName <- file.path(expDir, paste0("PCAVectorField_",expName,"_by",pVar,".pdf"))
    pdf(file = fileName, width = 10, height = 10)
    print(image)
    dev.off()
  }
  
}








