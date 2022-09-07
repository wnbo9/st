# load modules
import stlearn as st
import pandas as pd
import scipy.io
import numpy as np
from pathlib import Path
st.settings.set_figure_params(dpi=180)

# load cell information
cm = pd.read_csv('D:File/Summer_I/Aug/scikit/info.csv')
cm.iloc[:,1:]

# cell location
cimg = cm.iloc[:,1:3]
cimg.columns=['imagerow', 'imagecol']

# set stlearn image and matrix
ct = pd.DataFrame(np.ones([np.shape(cimg)[0], 300]))
adata = st.create_stlearn(count=ct, spatial=cimg, scale=1, library_id="Sample_test", image_path='D:/File/Summer_I/Aug/HE_image/DP.png')

st.pp.filter_genes(adata,min_cells=1)
st.pp.normalize_total(adata)
st.pp.log1p(adata)

# image tiling
st.pp.tiling(adata, "D:/File/Summer_I/Aug/scikit/stlearn", crop_size=40)
# feature extraction
st.pp.extract_feature(adata)

# save the 2048-dimensional latent vector
pd.DataFrame(adata.obsm["X_tile_feature"]).to_csv('D:/File/Summer_I/Aug/scikit/stlearn.csv')
