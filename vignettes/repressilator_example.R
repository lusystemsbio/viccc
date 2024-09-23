rm(list=ls())


# these lines install the prerequisites
#install.packages(c( "ellipsis","tidyr","varhandle", "ggplot2","BiocManager"))
#library(BiocManager)
#BiocManager::install(c("sRACIPE","SingleCellExperiment"))

# these lines install this package if you cloned the git branch - ignore otherwise!
# these will be removed at a later date, when the package is hosted on CRAN or Bioconductor
#install.packages("devtools")
library(roxygen2)
roxygen2::roxygenize()
library(devtools)
devtools::install()


library(sRACIPE)
library(SingleCellExperiment)
library(STICCC)
set.seed(123)


# global params
topoName <- "repressilator"
nSamples <- 2000


# directory setup
dir.create(file.path(getwd(),"input_topos"))
topoDir <- file.path(getwd(),topoName)
outputDir = file.path(topoDir,"data")
if(!dir.exists(topoDir)) {
  dir.create(topoDir)
}
if(!dir.exists(outputDir)) {
  dir.create(outputDir)
}

# Let's make a simple topology to work with: a 3-gene repressilator
# We'll save it to a folder called input_topos, which is where the following function will check by default
topo <- data.frame(Source=c("A","B","C"),Target=c("B","C","A"),Type=c(2,2,2))
write.table(topo, file = file.path(getwd(),"input_topos","repressilator.tpo"), sep = "\t")

# load topology file (if you have topos in a different folder than above, specify using the parameter topoDir)
topo <- loadTopo(topoName)

# simulate topology using sRACIPE - 2000 models by default, can specify with the param numModels
racipe <- simTopo(topo, numModels = nSamples)
saveRDS(racipe, file.path(outputDir, paste0("simData_",topoName,".Rds")))

# normalize data using sRACIPE function
racipe_norm <- sracipeNormalize(racipe)
exprMat <- assay(racipe)
exprMat_norm <- assay(racipe_norm)


# create SCE object
# you can also specify other params radius, minNeighbors, etc - see ?sticSE
stic <- sticSE(topo = topo, exprMat = exprMat, normData = exprMat_norm,
             topoName = topoName, expName = "repressilator_example")


# Cluster expression data
stic <- prepMetadata(sce = stic, exprMat = exprMat, cluster = T, k = 6)

# run PCA
stic <- runPCA(stic)


# compute grid based on PCA
stic <- computeGrid(stic)

# compute pairwise distance between points
stic <- computeDist(stic)

# infer velocities - for a repressilator with 2000 models, expect this to take ~5 min
stic <- runSTICCC(stic)
saveRDS(stic, file.path(topoDir,"stic_repressilator_example.Rds"))
stic <- readRDS(file.path(topoDir,"stic_repressilator_example.Rds"))

# grid-based smoothing of velocities
stic <- computeGridVectors(stic)

# plot individual vectors
plotVectors(stic, scalingFactor = 0.4,
            colorVar = "Cluster")

# plot grid-smoothed vectors
plotGrid(stic, colorVar = "Cluster")
