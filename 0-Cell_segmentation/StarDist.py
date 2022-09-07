# install StarDist module
pip install stardist

# load modules
import cv2
import pandas as pd
import matplotlib.pyplot as plt

from stardist.data import test_image_nuclei_2d
from stardist.plot import render_label
from csbdeep.utils import normalize
from skimage import measure
from skimage.measure import label, regionprops_table

from stardist.models import StarDist2D

# prints a list of available models
StarDist2D.from_pretrained()

# creates a pretrained model
model = StarDist2D.from_pretrained('2D_versatile_fluo')

# load H&E image
img = cv2.imread('./OneDrive/Desktop/forebrain_test.png', cv2.IMREAD_GRAYSCALE)
plt.imshow(img)

# cell segmentation
# here nms is overlapping threshold, prob is the confidence threshold
labels, _ = model.predict_instances(normalize(img), nms_thresh=0.3, prob_thresh=0.3)

# show the segmentation result
plt.subplot(1,2,1)
plt.imshow(img, cmap="gray")
plt.axis("off")
plt.title("input image")

plt.subplot(1,2,2)
plt.imshow(render_label(labels, img=img))
plt.axis("off")
plt.title("prediction + input overlay")

# measure the mask to label
img_label = measure.label(labels)
plt.imshow(img_label)
pd.DataFrame(img_label).to_csv('./OneDrive/Desktop/mask.csv')

# extract the cell information
props = regionprops_table(img_label, properties=('centroid',
                                          'area',
                                          'perimeter',
                                          'feret_diameter_max',
                                          'solidity',
                                          'axis_major_length',
                                          'axis_minor_length'))
props = pd.DataFrame(props)
props['AR'] = props['axis_major_length'] / props['axis_minor_length']
props.to_csv('./OneDrive/Desktop/info.csv')
