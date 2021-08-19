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
  racipe <- sracipeSimulate(racipe, numModels = 10000, ...)
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



# prepMetadata

# runPCA

# computeGrid

# computeDist

# computeTrajectory - NOT IN THIS FILE

# computeGridVectors

# plotVectors - OTHER FILE

# plotGrid - OTHER FILE











