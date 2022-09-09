library(dplyr)
library(ggplot2)
library(Seurat)
library(SeuratObject)

# load SpatialPCA_resul
load('stereo/SpatialPCA/SpatialPCA_forebrain_simple.RData')

# rename the cluster
Seu <- RenameIdents(object = Seu, 
                    `cluster1` = "Meninges", `cluster3` = "LPall", `cluster6` = "DPall_v",
                    `cluster4` = 'DPall_s', `cluster5` = 'DPall_i', `cluster2`='DPall_p')
# plot color
cbp = c("#FD7446" ,"#709AE1", "#9EDAE5",
        "#DE9ED6","#BCBD22" ,"#31A354")

# spatial plot of marker gene
SpatialFeaturePlot(st, features = 'Pbx3', pt.size.factor = 1.5)
# umap plot
DimPlot(Seu, reduction = "umap", cbp)

# spatial plot of domains
plot_cluster(location=SpatialPCA_result$location,
             Seu@active.ident,
             pointsize=1.5,
             title_in=paste0("Forebrian_Simple"),
             color_in=cbp,legend="right") + coord_fixed(ratio=1)


# find domain specific genes
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
  cbp_spatialpca = cbp
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


