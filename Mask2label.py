import pandas as pd
import numpy as np
import skimage


from skimage import data, util, measure
from skimage.filters import threshold_otsu, threshold_local
from skimage.measure import label, regionprops, regionprops_table
from skimage.segmentation import watershed
from skimage.feature import peak_local_max

from scipy import ndimage as ndi
import matplotlib.pyplot as plt
import tifffile as tiff

import cv2
import pandas as pd
from skimage import data, util, measure
from skimage.measure import label, regionprops, regionprops_table

ww = cv2.imread('./DP11.png', cv2.IMREAD_GRAYSCALE)

tt = measure.label(ww)

plt.imshow(ww)

pd.DataFrame(tt).iloc[3000:3100, 1000:1100]

pd.DataFrame(tt).to_csv('./Aug/cellprofiler/cell.csv')

props = regionprops_table(tt, properties=('centroid',
                                          'area', 'area_bbox', 'area_convex',
                                          'orientation',
                                          'axis_major_length',
                                          'axis_minor_length',
                                          'feret_diameter_max',
                                         'solidity',
                                         'perimeter'))

pd.DataFrame(props).to_csv('./Aug/cellprofiler/info.csv')
