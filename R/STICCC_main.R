

#' @export
#' @title Infer velocities of single cells
#' @description Computes the inferred velocity of each cell in the dataset according to the topology supplied previously.
#' @return sce object with velocities - see colData(sce) for the first two components (dX and dY),
#' and sce@metadata$vectors for all components.
#' @importFrom plyr revalue
#' @param sce sticcc object. Must have already supplied a topology and executed the functions runPCA() and computeDist().
#' @param v2 logical. If TRUE, will calculate v2 (incoming transition) as well as v1 (outgoing). Default TRUE.
#' @param useGinv logical. Whether to use generalized inverse for computing least-squares regression. 
#' Default FALSE, but may be helpful with very homogeneous datasets.
runSTICCC <- function(sce, v2=TRUE, useGinv=FALSE, ...
) {

  # Create sign vector - 1 for activation, -1 for inhibition
  topo <- sce@metadata$topo
  sign_vector <- revalue(factor(topo$Type), c('1'='1','2'='-1','3'='1','4'='-1','5'='1','6'='-1'), warn_missing = FALSE)
  sign_vector <- as.numeric(levels(sign_vector))[sign_vector]

  # Compute sampling radius from parameter & precomputed max dist
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
  colData(sce)$Magnitude <- NA
  colData(sce)$fewNeighbors <- FALSE
  colData(sce)$nonInvertible <- FALSE
  colData(sce)$selfCorNA <- FALSE
  #colData(sce)$rcond <- NA
  #sce@metadata$det_list <- list()

  weighted_vector_list <- vector(mode = "list", length = length(sample_list))
  pb = txtProgressBar(min = 0, max = length(sample_list),
                      initial = 0, style = 3)
  stepi <- 0
  for(query_point in sample_list) {
    ## Update progress bar
    stepi <- stepi+1
    setTxtProgressBar(pb, stepi)

    ## Compute vector
    rs_list <- computeVector(sce, query_point, v2 = v2, useGinv = useGinv, ...) #, bias = bias
    
    
    #"X"=NA,"Y"=NA,"fewNeighbors"=FALSE,"numNeighbors"=NA,"usedGinv"=FALSE,"nonInvertible"=FALSE,"selfCorNA"=FALSE,"selfCor"=NA,"b_vec"=NA
    # Store data
    colData(sce)[query_point,"X"] <- rs_list[["X"]]
    colData(sce)[query_point,"Y"] <- rs_list[["Y"]]
    colData(sce)[query_point,"fewNeighbors"] <- rs_list[["fewNeighbors"]]
    colData(sce)[query_point,"numNeighbors"] <- rs_list[["numNeighbors"]]
    colData(sce)[query_point,"usedGinv"] <- rs_list[["usedGinv"]]
    colData(sce)[query_point,"nonInvertible"] <- rs_list[["nonInvertible"]]
    colData(sce)[query_point,"selfCorNA"] <- rs_list[["selfCorNA"]]
    colData(sce)[query_point,"selfCor"] <- rs_list[["selfCor"]]
    #sce@metadata$det_list[[query_point]] <- rs_list[["det_list"]]
    
    if(length(rs_list[["b_vec"]]) == 1) {
      next
    } 
    
    
    if(bias) {
      error[stepi,] <- rs_list[["error"]]
    }
    
    
    ## Save vector to master list
    saved_vector <- as.data.frame(t(rs_list[["b_vec"]]))
    colnames(saved_vector) <- paste0("d",posList)
    weighted_vector_list[[stepi]] <- saved_vector
    names(weighted_vector_list)[[stepi]] <- query_point
    
    ## Update colData
    colData(sce)[query_point,"dX"] <- saved_vector[1]
    colData(sce)[query_point,"dY"] <- saved_vector[2]
    
    
    ## Add in v2 if necessary
    if(v2) {
      ## Save vector to master list
      saved_vector_in <- as.data.frame(t(rs_list[["b_vec_in"]]))
      colnames(saved_vector_in) <- paste0("d",posList)
      weighted_vector_list_in[[stepi]] <- saved_vector_in
      names(weighted_vector_list_in)[[stepi]] <- query_point
      
      ## Update colData
      colData(sce)[query_point,"dX_in"] <- saved_vector_in[1]
      colData(sce)[query_point,"dY_in"] <- saved_vector_in[2]
    }
    
    
  }

  ## Consolidate vectors and write to file
  weighted_vectors <- do.call(rbind, weighted_vector_list)
  sce@metadata$vectors <- weighted_vectors

  if(v2) {
    weighted_vectors_in <- do.call(rbind, weighted_vector_list_in)
    sce@metadata$vectors_in <- weighted_vectors_in
  }
  
  
  return(sce)

}









