library(dplyr)
library(data.table)
library(ggplot2)

dat = fread("stereo/L20220116011.export.bin1.txt", header = T)
colnames(dat) = c("geneID", "y", "x", "MIDCounts")
dat$y = max(dat$y) + min(dat$y) - dat$y
dat = cbind(dat, c(1:dim(dat)[1]))
dat = cbind(dat, rep(-1, dim(dat)[1]))
colnames(dat) = c("geneID", "y", "x", "MIDCounts", "ID", "cellID")


dat = dat %>% filter(x <= 4000 & 3000<=y & y<=10000)
dat$y = dat$y - 3000
dat$y[which(dat$y==0)] = 1
dat[1:10000,] %>% ggplot(aes(x=x,y=y)) + geom_point(size=0.1) + coord_fixed(ratio = 1)


info = fread("ARI/Fiji/cell.csv", header = T)
info = info[,-1]
info2 = as.matrix(info)
rm(info)

myFun=function(obj){
  c1 = as.integer(obj[['y']])
  c2 = as.integer(obj[['x']])
  return(info2[c1,c2][[1]])
}

vec = apply(dat, 1, myFun)
dat$cellID = vec


dat_keep = dat %>% filter(cellID != 0)


dat_temp = dat_keep %>% group_by(geneID, cellID) %>% summarise(sum = sum(MIDCounts))



library(reshape2)
DGE = dcast(dat_temp, geneID ~ cellID)
DGE[is.na(DGE)] = 0

