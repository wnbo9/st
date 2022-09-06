library(Seurat)
library(SeuratObject)
library(dplyr)
library(ggplot2)
library(data.table)

# read DGE matrix
dat = fread("ARI/simple/DGE.csv", header = TRUE)
dat = dat[,-1]
dat = as.data.frame(dat)
row.names(dat) = dat[,1]
dat = dat[,-1]

# add cell info
info = read.csv("ARI/simple/info.csv", sep=',', header = T)
info$X = info$X + 1
colnames(info)[1:3] = c('cell', 'Y', 'X')
info = info %>% filter(info$cell %in% colnames(dat))

nCell = dim(dat)[2]
colnames(dat) = c(1:nCell)
rownames(info) = c(1:nCell)

# create Seurat object
st = CreateSeuratObject(dat, assay = "Spatial")
st = AddMetaData(object = st, metadata = info)
st@images$image = new(
  Class='SlideSeq',
  assay='Spatial',
  key='image_',
  coordinates=info[,c(2,3)]
)


# analyze using Seurat
VlnPlot(st, features = 'nCount_Spatial', pt.size = 0) + NoLegend()
SpatialFeaturePlot(st, features = 'nCount_Spatial', pt.size.factor = 1)


#######################
# PCA
st = SCTransform(st, assay = 'Spatial', verbose = FALSE)
st = RunPCA(st, assay='SCT', verbose=FALSE)
st = FindNeighbors(st, reduction = 'pca', dims=1:10)
st = FindClusters(st, verbose = FALSE, resolution = 0.5)
st = RunUMAP(st, reduction = 'pca', dims=1:10)
# decide the PC dimensionality
ElbowPlot(st)

# spatial plot for all
DimPlot(st, reduction = "umap")
SpatialDimPlot(st, pt.size.factor = 1)
SpatialDimPlot(st, 
               cells.highlight = CellsByIdentities(object = st, 
                                                   idents = c(0,1,2,3,4,5)), 
               facet.highlight = TRUE, ncol = 3)
#rename
sn <- RenameIdents(object = st, 
                   `0` = "DPall1", `1` = "DPall2", `2` = "Cavity",
                   `3` = 'DPall3', `4` = 'Vessel', `5`='B2', `6` = 'B1')
# spatial plot for all
DimPlot(sn, reduction = "umap")
SpatialDimPlot(sn, pt.size.factor = 1)
# plot of cluster

SpatialDimPlot(sn, 
               cells.highlight = CellsByIdentities(object = sn, 
                                                   idents = c('DPall1', 'DPall2', 'DPall3')), 
               facet.highlight = TRUE, ncol = 3)
SpatialDimPlot(sn, 
               cells.highlight = CellsByIdentities(object = sn, 
                                                   idents = c('B1', 'B2', 'Vessel')), 
               facet.highlight = TRUE, ncol = 3)
SpatialDimPlot(sn, 
               cells.highlight = CellsByIdentities(object = sn, 
                                                   idents = c('Cavity')), 
               facet.highlight = TRUE, ncol = 3)



# marker gene for cluster
cluster2.markers <- FindMarkers(sn, ident.1 = 0, ident.2 = 1, min.pct = 0.25)
head(cluster2.markers, n = 10)
features = c(rownames(cluster2.markers))

## heatmap
DoHeatmap(st, features = features, size = 3)

# marker gene for all clusters
st.markers <- FindAllMarkers(st, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
st.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(st, features = top10$gene)

## feature plot for cluster
VlnPlot(sn, features = top10$gene,
        pt.size = 0.2, ncol = 4)
## dot plot
DotPlot(sn, features = top10$gene) + RotatedAxis()

`0` = "DPall1", `1` = "DPall2", `2` = "Cavity",
`3` = 'DPall3', `4` = 'Vessel', `5`='B2', `6` = 'B1')
sn@meta.data$iden = 'N'
sn@meta.data$iden[sn@meta.data$seurat_clusters == 0] = 'DPall1'
sn@meta.data$iden[sn@meta.data$seurat_clusters == 1] = 'DPall2'
sn@meta.data$iden[sn@meta.data$seurat_clusters == 2] = 'Cavity'
sn@meta.data$iden[sn@meta.data$seurat_clusters == 3] = 'DPall3'
sn@meta.data$iden[sn@meta.data$seurat_clusters == 4] = 'Vessel'
sn@meta.data$iden[sn@meta.data$seurat_clusters == 5] = 'B2'
sn@meta.data$iden[sn@meta.data$seurat_clusters == 6] = 'B1'
sn@meta.data$iden = factor(sn@meta.data$iden, levels = c("DPall2", "DPall1",  "Vessel",
                                                         'DPall3', 'B2', 'B1', 'Cavity'))
p1 = ggplot(sn@meta.data, aes(x=iden, y=area)) + 
  geom_boxplot(outlier.shape = NA) +
  coord_cartesian(ylim=c(0, 500))
p2 = ggplot(sn@meta.data, aes(x=iden, y=perimeter)) + 
  geom_boxplot(outlier.shape = NA) +
  coord_cartesian(ylim=c(0, 100))
p3 = ggplot(sn@meta.data, aes(x=iden, y=nFeature_Spatial)) + 
  geom_boxplot(outlier.shape = NA) +
  coord_cartesian(ylim=c(0, 300))
p4 = ggplot(sn@meta.data, aes(x=iden, y=nCount_Spatial)) + 
  geom_boxplot(outlier.shape = NA) +
  coord_cartesian(ylim=c(0, 500))
library(patchwork)
p1 + p2 + p3 + p4
