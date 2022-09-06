import pandas as pd
import numpy as np
import cv2
import skimage


from skimage import data, util, measure
from skimage.filters import threshold_otsu, threshold_local
from skimage.measure import label, regionprops, regionprops_table
from skimage.segmentation import watershed
from skimage.feature import peak_local_max

from scipy import ndimage as ndi
import matplotlib.pyplot as plt
import tifffile as tiff

img = cv2.imread('./Aug/1.png', cv2.IMREAD_GRAYSCALE)
#img = 255 - img
plt.imshow(img)

# remove the background using a global threshold binarization
global_thresh = threshold_otsu(img)
binary_global = img > global_thresh

# obtain markers
block_size = 41
local_thresh = threshold_local(img, block_size, offset=0.003)
binary_local = img > local_thresh

fig, axes = plt.subplots(ncols=3, figsize=(20, 20))
ax = axes.ravel()
plt.gray()

ax[0].imshow(img)
ax[0].set_title('Original')

ax[1].imshow(binary_global)
ax[1].set_title('Global thresholding')

ax[2].imshow(binary_local)
ax[2].set_title('Local thresholding')

# distance transformation using marker
distance = ndi.distance_transform_edt(binary_local)

fig, axes = plt.subplots(ncols=3, figsize=(20, 20))
ax = axes.ravel()
plt.gray()

ax[0].imshow(img)
ax[0].set_title('Original')

ax[1].imshow(binary_local)
ax[1].set_title('Local thresholding')

ax[2].imshow(distance)
ax[2].set_title('Distance')

local_max = peak_local_max(distance, indices=False, min_distance=11)
markers = ndi.label(local_max, structure=np.ones((3, 3)))[0]
labels = watershed(-distance, markers, mask=binary_global) #here the mask is the foreground object


fig, axes = plt.subplots(ncols=3, figsize=(20, 20), sharex=True, sharey=True)
ax = axes.ravel()

ax[0].imshow(img, cmap=plt.cm.gray)
ax[0].set_title('Overlapping objects')
ax[1].imshow(-distance, cmap=plt.cm.gray)
ax[1].set_title('Distances')
ax[2].imshow(labels, cmap=plt.cm.nipy_spectral)
ax[2].set_title('Separated objects')

tt = measure.label(labels)
#mask = plt.imshow(tt)
#io.imsave('mask.tif', mask)
#mask

