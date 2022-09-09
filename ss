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
i=1
j=1
# load cell info
info = read.csv(paste0('D:/File/Summer_I/Sep/label_mask/', region[i], '/', region[i], '_info_', method[j], '.csv'))
rownames(info) = paste0('cell', info$X+1)
# load countmat
load(paste0('D:/File/Summer_I/Sep/CGE/', region[i], '_CGE_', method[j], '.RData'))

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
                                numCores_spark=1,
                                gene.number=3000,
                                customGenelist=NULL,
                                min.loctions = 20,
                                min.features=20)
dim(stereo@normalized_expr)

# estimate kernel
stereo = SpatialPCA_buildKernel(stereo, kerneltype="gaussian", bandwidthtype="Silverman",bandwidth.set.by.user=NULL,sparseKernel=TRUE,sparseKernel_tol=1e-20,sparseKernel_ncore=1)
stereo = SpatialPCA_EstimateLoading(stereo,fast=TRUE,SpatialPCnum=20)
stereo = SpatialPCA_SpatialPCs(stereo, fast=TRUE)

SpatialPCA_result = list()
SpatialPCA_result$SpatialPCs  = as.matrix(stereo@SpatialPCs)
SpatialPCA_result$normalized_expr  = stereo@normalized_expr
SpatialPCA_result$location = stereo@location
save(SpatialPCA_result, file = paste0('SpatialPCA_', region[i], '_', method[j], '.RData'))



# analysis
## obtain clustering result
SpatialPCA_louvain = louvain_clustering(6,latent_dat=as.matrix(stereo@SpatialPCs), 100) 

cbp_spatialpca = c( "#FD7446" ,"#709AE1", "#31A354","#9EDAE5",
                    "#DE9ED6" ,"#BCBD22", "#CE6DBD" ,"#DADAEB" ,
                    "yellow", "#FF9896","#91D1C2", "#C7E9C0" ,
                    "#6B6ECF", "#7B4173" )
metadata = data.frame(SpatialPCA_louvain)




# region specific gene detection
Seu <- CreateSeuratObject(counts = countmat, project = "xx", min.cells = 20, min.features = 20)
Seu = SCTransform(Seu, return.only.var.genes = FALSE, variable.features.n = NULL,  variable.features.rv.th = 1.3)
Seu <- RunPCA(Seu, features = VariableFeatures(object = Seu))
Seu <- FindNeighbors(Seu, dims = 1:10)
Seu <- FindClusters(Seu, resolution = 0.5)
Seu <- RunUMAP(Seu, reduction = 'pca', dims=1:10)
Idents(Seu) = paste0("cluster",metadata$SpatialPCA_louvain)
Seu <- RenameIdents(object = Seu, 
                    `cluster1` = "Meninges", `cluster2` = "LPall", `cluster3` = "DPall_v",
                    `cluster4` = 'DPall_s', `cluster5` = 'DPall_i', `cluster6`='DPall_p')

# add location metadata
st = AddMetaData(object = Seu, metadata = data.frame(SpatialPCA_result$location))
st@images$image = new(
  Class='SlideSeq',
  assay='Spatial',
  key='image_',
  coordinates=st@meta.data[,c("ycoord", "xcoord")]
)

# umap plot
DimPlot(st, reduction = "umap")

# spatial plot of marker gene
SpatialFeaturePlot(st, features = 'Pbx3', pt.size.factor = 1.5)

plot_cluster(location=SpatialPCA_result$location,
             SpatialPCA_louvain,
             pointsize=1.5,
             title_in=paste0(region[i], '_', method[j]),
             color_in=cbp,legend="right") + coord_fixed(ratio=1)
#
plot_cluster(location=SpatialPCA_result$location,
             Seu@active.ident,
             pointsize=1.5,
             title_in=paste0(region[i], '_', method[j]),
             color_in=cbp,legend="right") + coord_fixed(ratio=1)





DE_gene = list()
DEgene_spatialPCA=c()
each_num = c()
for(cluster in 1:6){
  print(cluster)
  
  DE_gene[[cluster]] = FindMarkers(Seu, ident.1 = paste0("cluster",cluster), ident.2 =NULL, test.use = "MAST")
  each_num[cluster] = dim(DE_gene[[cluster]])[1]
  DEgene_spatialPCA = c(DEgene_spatialPCA, rownames(DE_gene[[cluster]]))
}
DEgene_spatialPCA=unique(DEgene_spatialPCA)
length(DEgene_spatialPCA)
each_num


# maker gene mean expression in each spatial domain (line plot)
make_lineplot = function(genename,clusterlabel,cluster_order){
  counts = countmat[which(rownames(countmat) %in% paste0(genename)),match(rownames(SpatialPCA_result$location), colnames(countmat))]
  cluster = as.character(clusterlabel)
  data=data.frame(counts, cluster)
  dat = data
  dat=dat[order(dat$counts),]
  dat$cellid = factor(paste0(1:dim(dat)[1]),levels=c(paste0(1:dim(dat)[1])),order=T)
  dat$cluster = factor(dat$cluster,levels=c(paste0(clusterorder)),order=T)
  ca <- dat %>%group_by(cluster) %>%summarise(
    mean = mean(counts),
    sd = sd(counts),
    n = n(),
    se = sd / sqrt(n))
  dattt =as.data.frame(ca)
  cbp_spatialpca = c( "#FD7446" ,"#709AE1", "#31A354","#9EDAE5",
                      "#DE9ED6" ,"#BCBD22", "#CE6DBD" ,"#DADAEB" ,
                      "yellow", "#FF9896","#91D1C2", "#C7E9C0" ,
                      "#6B6ECF", "#7B4173" )
  p<- ggplot(dattt, aes(x=cluster, y=mean, color=cluster,group=cluster)) + 
    #geom_bar(stat="identity", color="black", position=position_dodge()) +
    scale_fill_manual(values = cbp_spatialpca)+
    geom_vline(data=dattt, aes(xintercept=cluster, color=cluster),
               linetype="dashed",alpha=0.5)+
    geom_rect(
      aes(xmin = 0.5, xmax = 1.5, fill = cbp_spatialpca[3],alpha = 0.05), colour = "white",
      ymin = -Inf, ymax = Inf, alpha = 0.05) +
    geom_rect(
      aes(xmin = 1.5, xmax = 2.5, fill = cbp_spatialpca[2],alpha = 0.05), colour = "white",
      ymin = -Inf, ymax = Inf, alpha = 0.05) +
    geom_rect(
      aes(xmin = 2.5, xmax = 3.5, fill = cbp_spatialpca[6],alpha = 0.05), colour = "white",
      ymin = -Inf, ymax = Inf, alpha = 0.05) +
    geom_rect(
      aes(xmin = 3.5, xmax = 4.5, fill = cbp_spatialpca[5],alpha = 0.05), colour = "white",
      ymin = -Inf, ymax = Inf, alpha = 0.05) +
    geom_rect(
      aes(xmin = 4.5, xmax = 5.5, fill = cbp_spatialpca[1],alpha = 0.05), colour = "white",
      ymin = -Inf, ymax = Inf, alpha = 0.05) +
    geom_rect(
      aes(xmin = 5.5, xmax = 6.5, fill = cbp_spatialpca[4],alpha = 0.05), colour ="white",
      ymin = -Inf, ymax = Inf, alpha = 0.05) +
    #coord_flip()+
    theme_bw(base_size = 22)+
    theme(plot.title = element_text(size = 22),
          legend.position = "none")+
    geom_line(aes(group=1),color="black", size=1) + ### !!!!!!!!!! aes(group=1) is important
    geom_point( size=3, color="#20854E")+
    geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1,
                  position=position_dodge(.9),color="#20854E") +
    labs(title=paste0(genename),x="Layer", y = "Expression")
  return(p)
}

clusterorder=c(1,2,4,5,6,3)

#pdf("Markergenes_lineplot_SpatialPCA_slideseqV2.pdf",width=3,height=10)
#3
make_lineplot("Fabp7",metadata$SpatialPCA_louvain,clusterorder)
#6
make_lineplot("Sox11",metadata$SpatialPCA_louvain,clusterorder)
# 5
make_lineplot("Tmsb4x",metadata$SpatialPCA_louvain,clusterorder)
# 4
make_lineplot("Map1b",metadata$SpatialPCA_louvain,clusterorder)
# 2
make_lineplot("Pbx3",metadata$SpatialPCA_louvain,clusterorder)
make_lineplot("Nfib",metadata$SpatialPCA_louvain,clusterorder)
# 1
make_lineplot("Col1a2",metadata$SpatialPCA_louvain,clusterorder) # decode the most abundant protein of skin

#dev.off()


