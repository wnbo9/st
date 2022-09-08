# load packages
library(dplyr)
library(data.table)
library(ggplot2)
library(reshape2)

# load barcode file
dat = fread("stereo/L20220116011.export.bin1.txt", header = T)
colnames(dat) = c("geneID", "y", "x", "MIDCounts")
# registration with histology image
dat$x = dat$x - 288
dat$y = dat$y - 861


### forebrain
# select the barcodes within the forebrain region
mat_f = dat %>% filter(x<=4000 & 14000<y & y<=20000)
mat_f$y = mat_f$y - 14000
mat_f = cbind(mat_f, rep(0, dim(mat_f)[1]))
colnames(mat_f) = c("geneID", "y", "x", "MIDCounts", "cellID")
### GI tract
# select the barcodes within the GI tract region
mat_g = dat %>% filter(3500<x & x<=6500 & 4000<y & y<=9000)
mat_g$x = mat_g$x - 3500
mat_g$y = mat_g$y - 4000
mat_g = cbind(mat_g, rep(0, dim(mat_g)[1]))
colnames(mat_g) = c("geneID", "y", "x", "MIDCounts", "cellID")



# load the label csv file
region = c('forebrain','tract')
method = c('simple','imagej','cellprofiler','scikit','stardist','cellpose')

# function of aggregation
myFun = function(obj){
  c1 = as.integer(obj[['y']])
  c2 = as.integer(obj[['x']])
  return(label[c1, c2][[1]])
}

#
for (i in 1:2) {
  for (j in 1:6) {
    print(paste0(region[1], ': ', method[j]))
    # specify file
    file = paste0('stereo/label_mask/', region[i], '_label_', method[j], '.csv')
    # load the file
    label = fread(file, header = FALSE)
    label = label[-1,]
    label = label[,-1]
    label = as.matrix(label)
    
    if(i == 1){
      mat_temp = mat_f
    } else{
      mat_temp = mat_g
    }
    
    # specify the cell labels
    print('specifying')
    vec = apply(mat_temp, 1, myFun)
    mat_temp$cellID = vec
    
    # keep the barcodes within cells
    print('summarising')
    mat_keep = mat_temp %>% filter(cellID != 0)
    # summarise the barcodes within cells
    mat_keep = mat_keep %>% group_by(geneID, cellID) %>% summarise(sum = sum(MIDCounts))
    # reshape
    print('reshaping')
    CGE = dcast(mat_keep, geneID ~ cellID)
    # remove na
    CGE[is.na(CGE)] = 0
    
    # save CGE
    print('writing')
    output = paste0('stereo/CGE/', region[i], '_CGE_', method[j], '.csv')
    write.csv(CGE, output)
    
    # release RAM
    rm(mat_temp)
    rm(mat_keep)
  }
}
