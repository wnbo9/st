# load modules
import cv2
import pandas as pd
import matplotlib.pyplot as plt
from skimage import measure
from skimage.measure import label, regionprops, regionprops_table

# load images (labeled mask)
img = cv2.imread('./OneDrive/Desktop/mask.png', cv2.IMREAD_GRAYSCALE)
plt.imshow(img)

# create a img_size dataframe, so that each pixel can be responded to its cell
img_label = measure.label(img)

#save
pd.DataFrame(img_label).to_csv('./OneDrive/Desktop/label.csv')
plt.imshow(img_label)

# extrat the cell information of each cell
props = regionprops_table(img_label, properties=('centroid',
                                          'area',
                                          'perimeter',
                                          'feret_diameter_max',
                                          'solidity',
                                          'axis_major_length',
                                          'axis_minor_length'))
props = pd.DataFrame(props)
props['AR'] = props['axis_major_length'] / props['axis_minor_length']
# save
# Attention: in the mask.csv, '0' means the background, and the numbers represent the cells. Meanwhile, the index of info.csv starts from '0', here it means '1' in the mask.csv
props.to_csv('./OneDrive/Desktop/info.csv')
