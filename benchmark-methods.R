###############
## Slingshot ##
###############

runSlingshot <- function(counts){
  library(tictoc)
  library(slingshot)

  tic()
  # processing for Slingshot
  # adapted from slingshot vignette
  # do not do any normalisation

  # dimensionality reduction
  pca <- prcomp(t(log1p(counts)), scale. = FALSE)
  rd1 <- pca$x[,1:2]
  cl1 <- kmeans(rd1, centers = 5)$cluster
  lin1 <- getLineages(rd1, cl1)
  crv1 <- getCurves(lin1)

  time <- toc(quiet=T)
  dt <- time$toc - time$tic
  return(dt)
}


###############
## Monocle 2 ##
###############

runMonocle <- function(counts){
  library(tictoc)
  library(monocle)
  # use sparse matrices
  library(Matrix)
  pd <- new("AnnotatedDataFrame", data = data.frame(cellids=colnames(counts), row.names=colnames(counts)))
  fd <- new("AnnotatedDataFrame", data=data.frame(gene_short_name=rownames(counts), row.names=rownames(counts)))

  SPLT <- newCellDataSet(as(counts, "sparseMatrix"),
                  phenoData = pd,
                  featureData = fd,
                  expressionFamily=negbinomial.size())

  tic()
  # processing for Monocle
  # adapted from Monocle vignette
  SPLT <- estimateSizeFactors(SPLT)
  SPLT <- estimateDispersions(SPLT)

  # select ordering genes based on variance
  # NB there are several alternatives to select genes
  disp_table <- dispersionTable(SPLT)
  ordering_genes <- subset(disp_table,
                    mean_expression >= 0.5 &
                    dispersion_empirical >= 1 * dispersion_fit)$gene_id

  SPLT <- setOrderingFilter(SPLT, ordering_genes)
  # reduce data dimensionality
  SPLT <- reduceDimension(SPLT, max_components = 2,
      method = 'DDRTree')
  # order cells along the trajectory
  SPLT <- orderCells(SPLT)

  time <- toc(quiet=T)
  dt <- time$toc - time$tic
  return(dt)
}


###########
## TSCAN ##
###########

runTSCAN <- function(counts){
  library(tictoc)
  library(TSCAN)

  tic()
  # Preprocess Gene Expression Profiles
  procdata <- preprocess(counts)
  # Constructing Pseudotime
  lpsmclust <- exprmclust(procdata)
  lpsorder <- TSCANorder(lpsmclust)

  time <- toc(quiet=T)
  dt <- time$toc - time$tic
  return(dt)
}


###################
## DTP (destiny) ##
###################
# https://bioconductor.org/packages/release/bioc/html/destiny.html

runDTP <- function(counts){
  library(tictoc)
  library(destiny)
  library(Biobase)

  ct <- ExpressionSet(assayData=counts)

  tic()
  dm <- DiffusionMap(ct)
  dpt <- DPT(dm)

  time <- toc(quiet=T)
  dt <- time$toc - time$tic
  return(dt)
}


###########################
## Import Python Methods ##
###########################

library("reticulate")
use_python("/usr/bin/python3.6")

source_python("code/trajectory-inference-methods.py")


###############
## Benchmark ##
###############

runBenchmark <- function(seed, nGenes, nCells, method){
  library(splatter)
  set.seed(seed)
  params <- readRDS("data/splatter-params.rds")
  params <- setParam(params, "seed", seed)
  counts <- counts(splatSimulate(params, nGenes=nGenes, batchCells=nCells, verbose=F))
  dt <- NA
  if (method == "slingshot"){
    try(dt <- runSlingshot(counts))
  } else if (method == "monocle"){
    try(dt <- runMonocle(counts))
  } else if (method == "TSCAN"){
    try(dt <- runTSCAN(counts))
  } else if (method == "DTP"){
    try(dt <- runDTP(counts))
  } else if (method == "scanpyAGA"){
    try(dt <- runScanpyAGA(counts))
  } else if (method == "scanpyDPT"){
    try(dt <- runScanpyDPT(counts))
  } else if (method == "wishbone"){
    try(dt <- runWishbone(counts))
  } else if (method == "GPfates"){
    try(dt <- runGPfates(counts))
  }
  return(dt)
}


runBatch <- function(method, seeds, nGenesVec, nCellsVec){
  # run batch for single method.
  # add to file after each run
  if (length(nGenesVec) == 1 & length(nCellsVec) > 1){
    fileid <- paste0("data/benchmark/", method, "-cell-runningtime.csv")
  } else if (length(nGenesVec) > 1 & length(nCellsVec) == 1){
    fileid <- paste0("data/benchmark/", method, "-gene-runningtime.csv")
  } else {
    fileid <- paste0("data/benchmark/", method, "-runningtime.csv")
  }

  if (!file.exists(fileid)){
    write(paste(c("genes", "cells", "running_time", "seed"), collapse=","),file=fileid,append=TRUE)
  }
  for (nGenes in nGenesVec){
    for (nCells in nCellsVec){
      for (seed in seeds){
        dtresults <- runBenchmark(seed, nGenes, nCells, method)
        write(paste(c(nGenes, nCells, dtresults, seed) , collapse=","), file=fileid, append=TRUE)
      }
    }
  }
}
