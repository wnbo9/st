library(ggplot2)
library(dplyr)
library(SpatialPCA)
library(bluster)
library(Seurat)
library(MAST)
library(SeuratObject)

# specify data
region = c('forebrain','tract')
method = c('simple','imagej','cellprofiler','scikit','stardist','cellpose')
for (i in 1:2) {
  for (j in 1:6) {
    print(paste0(region[i], " : ", method[j]))
    # load cell info
    info = read.csv(paste0('stereo/label_mask/', region[i], '_info_', method[j], '.csv'))
    rownames(info) = paste0('cell', info$X+1)
    # load countmat
    load(paste0('stereo/CGE/', region[i], '_CGE_', method[j], '.RData'))
    
    # remove missing cells
    info = info[,-1]
    info = info %>% filter(rownames(info) %in% colnames(countmat))
    # get location
    location = info[,c('centroid.1', 'centroid.0')]
    colnames(location) = c('xcoord', 'ycoord')
    location = as.matrix(location)
    
    # run SpatialPCA
    ## create SpatialPCA object
    stereo = CreateSpatialPCAObject(counts=countmat,
                                    location=location,
                                    project = "SpatialPCA",
                                    gene.type="spatial",
                                    sparkversion="sparkx",
                                    numCores_spark=4,
                                    gene.number=3000,
                                    customGenelist=NULL,
                                    min.loctions = 20,
                                    min.features=20)
    dim(stereo@normalized_expr)
    
    # estimate kernel
    stereo = SpatialPCA_buildKernel(stereo,
                                    kerneltype="gaussian",
                                    bandwidthtype="Silverman",
                                    bandwidth.set.by.user=NULL,
                                    sparseKernel=TRUE,sparseKernel_tol=1e-20,
                                    sparseKernel_ncore=4)
    stereo = SpatialPCA_EstimateLoading(stereo,fast=TRUE,SpatialPCnum=20)
    stereo = SpatialPCA_SpatialPCs(stereo, fast=TRUE)
    # create SpatialPCA result
    SpatialPCA_result = list()
    SpatialPCA_result$SpatialPCs  = as.matrix(stereo@SpatialPCs)
    SpatialPCA_result$normalized_expr  = stereo@normalized_expr
    SpatialPCA_result$location = stereo@location
    
    # louvain clustering
    SpatialPCA_louvain = louvain_clustering(6, latent_dat = as.matrix(stereo@SpatialPCs), 100)
    metadata = data.frame(SpatialPCA_louvain)
    
    # Seurat analysis
    Seu <- CreateSeuratObject(counts = countmat, project = "xx", min.cells = 20, min.features = 20)
    Seu = SCTransform(Seu, return.only.var.genes = FALSE, variable.features.n = NULL,  variable.features.rv.th = 1.3)
    Seu <- RunPCA(Seu, features = VariableFeatures(object = Seu))
    Seu <- FindNeighbors(Seu, dims = 1:10)
    Seu <- FindClusters(Seu, resolution = 0.5)
    Seu <- RunUMAP(Seu, reduction = 'pca', dims=1:10)
    Idents(Seu) = paste0("cluster",metadata$SpatialPCA_louvain)
    Seu = AddMetaData(object = Seu, metadata = data.frame(SpatialPCA_result$location))
    Seu@images$image = new(
      Class='SlideSeq',
      assay='Spatial',
      key='image_',
      coordinates=Seu@meta.data[,c("ycoord", "xcoord")]
    )
    
    save(SpatialPCA_result, Seu, file = paste0('stereo/SpatialPCA/SpatialPCA_', region[i], '_', method[j], '.RData'))
    rm(list=ls())
  }
}
