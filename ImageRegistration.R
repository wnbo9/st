library(dplyr)
library(data.table)
library(ggplot2)

# the size of the original image is 15795×22111
# the range of barcode is x: 375~15324; y: 975~22924


# load expression data
dat = fread("stereo/L20220116011.export.bin1.txt", header = T)
colnames(dat) = c("geneID", "y", "x", "MIDCounts")

summary(dat$x)
summary(dat$y)

p1 = dat[1:100000,] %>% ggplot(aes(x=x, y=y)) + 
  geom_point(size=0.02) + 
  coord_fixed(ratio = 1) +
  geom_vline(xintercept = c(0, 17000)) +
  geom_hline(yintercept = c(0, 24000)) + NoLegend()

ggsave(p1, "Barcode.pdf", dpi = 300, width = 16, height = 24, units = "cm")
