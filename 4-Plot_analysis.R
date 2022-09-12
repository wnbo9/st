library(dplyr)
library(ggplot2)
library(Seurat)
library(SeuratObject)
library(ggpubr)
library(SpatialPCA)

# load SpatialPCA_resul
load('stereo/SpatialPCA/SpatialPCA_forebrain_simple.RData')
# load countmat for plot
load('D:/File/Summer_I/Sep/CGE/tract_CGE_cellpose.RData')

SpatialPCA_louvain = louvain_clustering(6, latent_dat = SpatialPCA_result$SpatialPCs, 100)
metadata = data.frame(SpatialPCA_louvain)
Idents(Seu) = paste0("cluster",metadata$SpatialPCA_louvain)
# rename the cluster
Seu <- RenameIdents(object = Seu, 
                    `cluster1` = "Mucosa", `cluster2` = "Submucosa",
                    `cluster5`='Inner muscle', `cluster4` = 'Outer muscle',
                    `cluster6` = 'Gland & Connective tissue', `cluster3` = 'Inter-mucosa')
# plot color
cbp = c( "#DE9ED6", '#7B4173', "#709AE1",
         "#9EDAE5","red", "#FD7446")
# spatial plot of domains
plot_cluster(location=SpatialPCA_result$location,
             Seu@active.ident,
             pointsize=1,
             title_in=paste0("GI tract - Cellpose"),
             color_in=cbp,legend="right") + coord_fixed(ratio=1)

# spatial plot of marker gene
SpatialFeaturePlot(Seu, features = 'Gm19951', pt.size.factor = 2)
# umap plot
DimPlot(Seu, reduction = "umap", cols = cbp, pt.size = 1)

# find domain specific genes
DE_gene = list()
DEgene_spatialPCA=c()
each_num = c()
layer = c("Mucosa", "Submucosa", "Inner muscle", 'Outer muscle', 'Gland', 'Connective tissue')
for(i in 1:6){
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
    geom_rect(
      aes(xmin = 4.5, xmax = 5.5, fill = cbp_spatialpca[4],alpha = 0.05), colour = "white",
      ymin = -Inf, ymax = Inf, alpha = 0.1) +
    geom_rect(
      aes(xmin = 5.5, xmax = 6.5, fill = cbp_spatialpca[4],alpha = 0.05), colour = "white",
      ymin = -Inf, ymax = Inf, alpha = 0.1) +
    #coord_flip()+
    theme_bw(base_size = 10)+
    theme(plot.title = element_text(size = 12),
          legend.position = "none") +  #remove y axis labels)+
    geom_line(aes(group=1),color="black", size=1) + ### !!!!!!!!!! aes(group=1) is important
    geom_point( size=3, color="#20854E")+
    geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1,
                  position=position_dodge(.9),color="#20854E")+
    ylab('') + xlab('Layer')
    #labs(title=paste0(genename),x="", y = "")
  return(p)
}
make_lineplot2 = function(genename,clusterlabel,cluster_order){
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
    theme_bw(base_size = 10)+
    theme(plot.title = element_text(size = 12),
          legend.position = "none") +  #remove y axis labels)+
    geom_line(aes(group=1),color="black", size=1) + ### !!!!!!!!!! aes(group=1) is important
    geom_point( size=3, color="#20854E")+
    geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1,
                  position=position_dodge(.9),color="#20854E")+
    ylab('') + xlab('Layer')
  #labs(title=paste0(genename),x="", y = "")
  return(p)
}
#pdf("Markergenes_lineplot_SpatialPCA_slideseqV2.pdf",width=3,height=10)
# Mucosa
p1 = make_lineplot("H19",Seu@active.ident,layer)
# Submucosa
p2 = make_lineplot("Hba-a1",Seu@active.ident,layer)
# Smooth muscle
p3 = make_lineplot("Gm19951",Seu@active.ident,layer)
p4 = make_lineplot("Myh11",Seu@active.ident,layer)
# Connective tissue
p5 = make_lineplot2("Cdk8",Seu@active.ident,layer)

ggarrange(p1,p2,p3,p4,p5, nrow = 5, labels=c('H19', 'Hba-a1', 'Acta2', 'Myh11', 'Cdk8'))

#dev.off()



# heatmap
DE_gene[["Mucosa"]] %>% arrange(desc(avg_log2FC)) %>% filter(avg_log2FC >0)

DEgene_spatialPCA_simple = c('H19', 'Rpl14', 'Hsp90ab1',
                             'Hba-a1', 'Hbb-bs', 'Hba-a2', 'Hbb-bt', 'Hbb-y', 'Gpc6', 'Adamdec1',
                             'Acta2', 'Actg2', 'Myh11', 'Myl9', 'Mylk', 'Tpm1', 'Flna', 'Cald1', 'Tpm2', 'Myl6', 'Tagln', 'Synpo2', 'Lpp', 'Cnn1', 'Sdk1', 'Prkg1',
                             'Tuba1a', 'Dpysl3',
                             'Cmss1', 'Camk1d', 'Jarid2', 'Cdk8', 'Col1a2', 'Lars2', 'S100a6', 'Zc3h7a', 'Gphn', 'Gm19951', 'D130009I18Rik',
                             'Slit2', 'Col3a1', 'Ptn', 'Pbx1')
DoHeatmap(Seu, features = DEgene_spatialPCA_simple)


# 
