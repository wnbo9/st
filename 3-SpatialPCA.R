library(ggplot2)
library(dplyr)
library(SpatialPCA)
library(bluster)
library(Seurat)
library(MAST)

# load cell info
info = read.csv('D:/File/Summer_I/Sep/label_mask/forebrain/forebrain_info_simple.csv')
rownames(info) = paste0('cell', info$X+1)
# remove missing cells
info = info[,-1]
info = info %>% filter(rownames(info) %in% colnames(countmat))
# get location
location = info[,c('centroid.1', 'centroid.0')]
colnames(location) = c('xcoord', 'ycoord')
location = as.matrix(location)

# run SpatialPCA
simple = CreateSpatialPCAObject(counts=countmat,
                                location=location,
                                project = "SpatialPCA",
                                gene.type="spatial",
                                sparkversion="sparkx",
                                numCores_spark=1,
                                gene.number=3000,
                                customGenelist=NULL,
                                min.loctions = 20,
                                min.features=20)


dim(simple@normalized_expr)

#
simple = SpatialPCA_buildKernel(simple, kerneltype="gaussian", bandwidthtype="Silverman",bandwidth.set.by.user=NULL,sparseKernel=TRUE,sparseKernel_tol=1e-20,sparseKernel_ncore=1)
simple = SpatialPCA_EstimateLoading(simple,fast=TRUE,SpatialPCnum=20)
simple = SpatialPCA_SpatialPCs(simple, fast=TRUE)


SpatialPCA_result = list()
SpatialPCA_result$SpatialPCs  = as.matrix(simple@SpatialPCs)
SpatialPCA_result$normalized_expr  = simple@normalized_expr
SpatialPCA_result$location = simple@location
save(SpatialPCA_result, file = "SlideseqV2_SpatialPCA_result.RData")

# obtain clustering result
SpatialPCA_louvain = louvain_clustering(6,latent_dat=as.matrix(simple@SpatialPCs), 100) 

cbp_spatialpca = c( "#FD7446" ,"#709AE1", "#31A354","#9EDAE5",
                    "#DE9ED6" ,"#BCBD22", "#CE6DBD" ,"#DADAEB" ,
                    "yellow", "#FF9896","#91D1C2", "#C7E9C0" ,
                    "#6B6ECF", "#7B4173" )

loc1 = unlist(simple@location[,1])
loc2 = unlist(simple@location[,2])
cluster = as.character(SpatialPCA_louvain)
datt = data.frame(cluster, loc1, loc2)
p = ggplot(datt, aes(x = loc1, y = loc2, color = cluster)) +
  geom_point( alpha = 1,size=0.5) +
  scale_color_manual(values = cbp_spatialpca)+
  theme_void()+
  theme(plot.title = element_text(size = 20,  face = "bold"),
        text = element_text(size = 20),
        #axis.title = element_text(face="bold"),
        #axis.text.x=element_text(size = 15) ,
        legend.position = "bottom") +
  coord_fixed(ratio=1)
p
ggsave(p, filename ='plot.jpg', width = 8, height = 8, 
       units = c('in'), dpi=300)

metadata = data.frame(SpatialPCA_louvain)




# region specific gene detection
Seu <- CreateSeuratObject(counts = countmat, project = "xx", min.cells = 20, min.features = 20)
Seu = SCTransform(Seu, return.only.var.genes = FALSE, variable.features.n = NULL,  variable.features.rv.th = 1.3)
Seu <- RunPCA(Seu, features = VariableFeatures(object = Seu))
Seu <- FindNeighbors(Seu, dims = 1:10)
Seu <- FindClusters(Seu, resolution = 0.5)

Idents(Seu) = paste0("cluster",metadata$SpatialPCA_louvain)
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



