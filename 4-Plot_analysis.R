library(dplyr)
library(ggplot2)
library(Seurat)
library(SeuratObject)

# load SpatialPCA_resul
load('stereo/SpatialPCA/SpatialPCA_forebrain_simple.RData')
# load countmat for plot
load('stereo/CGE/forebrain_CGE_simple.csv')

# rename the cluster
Seu <- RenameIdents(object = Seu, 
                    `cluster1` = "Meninges", `cluster3` = "LPall", `cluster6` = "DPall_v",
                    `cluster4` = 'DPall_s', `cluster5` = 'DPall_i', `cluster2`='DPall_p')
# plot color
cbp = c("#FD7446" ,"#709AE1", "#9EDAE5",
        "#DE9ED6","#BCBD22" ,"#31A354")

# spatial plot of marker gene
SpatialFeaturePlot(Seu, features = 'Pbx3', pt.size.factor = 1.5)
# umap plot
DimPlot(Seu, reduction = "umap", cols = cbp)

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

layer = c("Meninges", "LPall", "DPall_v",
        'DPall_s', 'DPall_i', 'DPall_p')
for(i in 1:6){
  cluster = lay[i]
  print(cluster)
  
  DE_gene[[cluster]] = FindMarkers(Seu, ident.1 = cluster, ident.2 =NULL, test.use = "MAST")
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
# DPall_v
make_lineplot("Fabp7",Seu@active.ident,layer)
# DPall_p
make_lineplot("Sox11",Seu@active.ident,layer)
# DPall_i
make_lineplot("Tmsb4x",Seu@active.ident,layer)
# Dpall_s
make_lineplot("Map1b",Seu@active.ident,layer)
# LPall
make_lineplot("Pbx3",Seu@active.ident,layer)
make_lineplot("Nfib",Seu@active.ident,layer)
# Meninges
make_lineplot("Col1a2",Seu@active.ident,layer) # decode the most abundant protein of skin

#dev.off()


#### GI tract
library(dplyr)
library(ggplot2)
library(Seurat)
library(SeuratObject)

# load SpatialPCA_resul
load('stereo/SpatialPCA/SpatialPCA_forebrain_simple.RData')
# load countmat for plot
load('D:/File/Summer_I/Sep/CGE/tract_CGE_simple.RData')

# rename the cluster
Seu <- RenameIdents(object = Seu, 
                    `cluster4` = "Mucosa", `cluster1` = "Submucosa",
                    `cluster2`='Smooth muscle', `cluster3` = 'Connective tissue')
# plot color
cbp = c("#DE9ED6" ,"#709AE1", "#9EDAE5", "#FD7446")

# spatial plot of marker gene
SpatialFeaturePlot(Seu, features = 'Hbb-y', pt.size.factor = 2)
# umap plot
DimPlot(Seu, reduction = "umap", cols = cbp, pt.size = 1)

# spatial plot of domains
plot_cluster(location=SpatialPCA_result$location,
             Seu@active.ident,
             pointsize=1.5,
             title_in=paste0("GI tract - Simple"),
             color_in=cbp,legend="right") + coord_fixed(ratio=1)


# find domain specific genes
DE_gene = list()
DEgene_spatialPCA=c()
each_num = c()

layer = c("Mucosa", "Submucosa", "Smooth muscle", 'Connective tissue')
for(i in 1:4){
  cluster = layer[i]
  print(cluster)
  
  DE_gene[[cluster]] = FindMarkers(Seu, ident.1 = cluster, ident.2 =NULL, test.use = "MAST")
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
  dat$cluster = factor(dat$cluster,levels=c(paste0(cluster_order)),order=T)
  ca <- dat %>%group_by(cluster) %>%summarise(
    mean = mean(counts),
    sd = sd(counts),
    n = n(),
    se = sd / sqrt(n))
  dattt =as.data.frame(ca)
  cbp_spatialpca =  c("#9EDAE5" ,"#DE9ED6", "#709AE1", "#FD7446")
  p<- ggplot(dattt, aes(x=cluster, y=mean, color=cluster,group=cluster)) + 
    #geom_bar(stat="identity", color="black", position=position_dodge()) +
    scale_fill_manual(values = cbp_spatialpca)+
    geom_vline(data=dattt, aes(xintercept=cluster, color=cluster),
               linetype="dashed",alpha=0.5)+
    geom_rect(
      aes(xmin = 0.5, xmax = 1.5, fill = cbp_spatialpca[1], alpha = 0.05), colour = "white",
      ymin = -Inf, ymax = Inf, alpha = 0.1) +
    geom_rect(
      aes(xmin = 1.5, xmax = 2.5, fill = cbp_spatialpca[2],alpha = 0.05), colour = "white",
      ymin = -Inf, ymax = Inf, alpha = 0.1) +
    geom_rect(
      aes(xmin = 2.5, xmax = 3.5, fill = cbp_spatialpca[3],alpha = 0.05), colour = "white",
      ymin = -Inf, ymax = Inf, alpha = 0.1) +
    geom_rect(
      aes(xmin = 3.5, xmax = 4.5, fill = cbp_spatialpca[4],alpha = 0.05), colour = "white",
      ymin = -Inf, ymax = Inf, alpha = 0.1) +
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
# Mucosa
make_lineplot("H19",Seu@active.ident,layer)
# Submucosa
make_lineplot("Hba-a1",Seu@active.ident,layer)
# Smooth muscle
make_lineplot("Acta2",Seu@active.ident,layer)
make_lineplot("Myh11",Seu@active.ident,layer)
# Connective tissue
make_lineplot("Cdk8",Seu@active.ident,layer)
#dev.off()

DE_gene[["Mucosa"]] %>% arrange(desc(avg_log2FC))

DEgene_spatialPCA = c('H19', 'Hsp90ab1', 'Rpl14', 'Oat', 'Rpl41',
                      'Hba-a1', 'Hbb-bs', 'Hba-a2', 'Hbb-bt', 'Hbb-y',
                      'Acta2', 'Myh11', 'Actg2', 'Myl9', 'Mylk',
                      'Cdk8', 'Cmss1', 'Jarid2', 'Col3a1', 'Camk1d')
DoHeatmap(Seu, features = DEgene_spatialPCA)
