

#' @export
#' @title Infer velocities of single cells
#' @description Computes the inferred velocity of each cell in the dataset according to the topology supplied previously.
#' @return sce object with velocities - see colData(sce) for the first two components (dX and dY),
#' and sce@metadata$vectors for all components.
#' @importFrom plyr revalue
#' @importFrom stats cor
#'
#' @param sce viccc object. Must have already supplied a topology and executed the functions runPCA() and computeDist().
DCComputeTrajectorySCE <- function(sce
) {

  # Create sign vector
  topo <- sce@metadata$topo
  sign_vector <- revalue(factor(topo$Type), c('1'='1','2'='-1','3'='1','4'='-1','5'='1','6'='-1'), warn_missing = FALSE)
  sign_vector <- as.numeric(levels(sign_vector))[sign_vector]

  # Compute sampling radius
  sampling_radius <- sce@metadata$params$sample_radius * sce@metadata$max_dist

  # get posMat (position matrix - this should be the data used to compute dist_mat)
  # if nPCs not specified, use all PCs (this may take a long time!)
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

    ## Find neighbors within sampling radius
    neighbors <- which(sce@metadata$dist_mat[query_point,colnames(sce)] <= sampling_radius)

    ## Skip sample if number of neighbors too small
    if(length(neighbors) <= (sce@metadata$params$minNeighbors+1)) {
      colData(sce)[query_point,"fewNeighbors"] <- TRUE
      colData(sce)[query_point,"numNeighbors"] <- length(neighbors)-1
      next
    }

    ## Select position data for neighbors (remove query point)
    subset_models <- as.data.frame(posMat[neighbors,])
    subset_models <- subset_models[-which(rownames(subset_models) %in% query_point),]
    colData(sce)[query_point,"numNeighbors"] <- nrow(subset_models)

    ## Multiple regression
    # Create relative positional matrix of neighbor samples by subtracting position of query point
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
    # This is equivalent to the least squares estimator
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
