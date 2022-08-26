mat_simple = cbind(st_simple@meta.data, method = 'simple')
mat_scikit = cbind(st_scikit@meta.data, method = 'scikit')
mat_imagej = cbind(st_imagej@meta.data, method = 'imagej')
mat_cellprofiler = cbind(st_cellprofiler@meta.data, method = 'cellprofiler')

mat_all = rbind(mat_simple, mat_scikit, mat_imagej, mat_cellprofiler)
mat_all$method = factor(mat_all$method, levels = c("simple", 'imagej', 'scikit', 'cellprofiler'))


p1 = mat_all %>% ggplot(aes(x=method, y=nFeature_Spatial)) +
  geom_boxplot(outlier.shape = NA) + 
  ylim(0, 500)
p2 = mat_all %>% ggplot(aes(x=method, y=nCount_Spatial)) +
  geom_boxplot(outlier.shape = NA) + 
  ylim(0, 500)
p3 = mat_all %>% ggplot(aes(x=method, y=area)) +
  geom_boxplot(outlier.shape = NA) + 
  ylim(0, 1000)
p4 = mat_all %>% ggplot(aes(x=method, y=axis_major_length)) + 
  geom_boxplot(outlier.shape = NA) + 
  ylim(0, 50)
p5 = mat_all %>% ggplot(aes(x=method, y=axis_minor_length)) + 
  geom_boxplot(outlier.shape = NA) + 
  ylim(0, 50)
p6 = mat_all %>% ggplot(aes(x=method, y=solidity)) + 
  geom_boxplot(outlier.shape = NA) + 
  ylim(0.5, 1)
p7 = mat_all %>% ggplot(aes(x=method, y=feret_diameter_max)) + 
  geom_boxplot(outlier.shape = NA) + 
  ylim(0, 50)
p8 = mat_all %>% ggplot(aes(x=method, y=perimeter)) + 
  geom_boxplot(outlier.shape = NA) + 
  ylim(0, 200)
p9 = mat_all %>% ggplot(aes(x=method, y=area_convex)) + 
  geom_boxplot(outlier.shape = NA) + 
  ylim(0, 1000)

p1 + p2 + p3 + p9
p4 + p5 + p6 + p8
s
