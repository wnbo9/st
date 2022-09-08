# load modules
import cv2
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from skimage import measure
from skimage.measure import label, regionprops_table

from cellpose import models, io

# load image
img = cv2.imread('./stereo/forebrain.png', cv2.IMREAD_GRAYSCALE)
plt.imshow(img)

# cell segmentation
# model_type='cyto' or 'nuclei' or 'cyto2'
model = models.Cellpose(model_type='cyto')

# define CHANNELS to run segementation on
# grayscale=0, R=1, G=2, B=3
# channels = [cytoplasm, nucleus]
# if NUCLEUS channel does not exist, set the second channel to 0
channels = [0,0]
# IF ALL YOUR IMAGES ARE THE SAME TYPE, you can give a list with 2 elements
# channels = [0,0] # IF YOU HAVE GRAYSCALE
# channels = [2,3] # IF YOU HAVE G=cytoplasm and B=nucleus
# channels = [2,1] # IF YOU HAVE G=cytoplasm and R=nucleus

# if diameter is set to None, the size of the cells is estimated on a per image basis
# you can set the average cell `diameter` in pixels yourself (recommended)
# diameter can be a list or a single number for all images

masks, flows, styles, diams = model.eval(img, diameter=None, channels=channels)

# mask to label
img_label = measure.label(mm)
plt.imshow(img_label)
pd.DataFrame(img_label).to_csv('./OneDrive/Desktop/label_cellpose.csv')


# extract cell informatipn
props = regionprops_table(img_label, properties=('centroid',
                                          'area',
                                          'perimeter',
                                          'feret_diameter_max',
                                          'solidity',
                                          'axis_major_length',
                                          'axis_minor_length'))
props = pd.DataFrame(props)
props['AR'] = props['axis_major_length'] / props['axis_minor_length']
props.to_csv('./OneDrive/Desktop/info_cellpose.csv')
