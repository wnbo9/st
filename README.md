# **Comparison of segmentation methods for spatial transcriptomic data analysis** #

## **Background**
With the advance of spatail transcriptomics, some recent technologies can achieve single cell even subcellular resolution for the spatail transcriptomic data, including Stereo-seq, HDST, Seq-Scope and etc. However, accurate characterization of cell boundaries and matching spots to cells for these technologies remain as an unsolved problem.

Cell segmentation has been a popular topic these years. Taking advantage of histology images, we can detect cell boundaries and segment these cells with small sizes. There exist some commonly used segmentation methods, and they can be divided into three classes based on their mechanisms: 
- Simple segmentation
  - Morphological segmentation
- Marker-controlled watershed segmentation
  - ImageJ
  - CellProfiler
  - Scikit-image
- Machine learning based
  - Cellpose
  - StarDist

Here, we make use of the mouse embryo tissue (E14.5 E1S3) from Stereo-seq to compare the performance and illustrate the characteristics of these segmentation methods.

## **Workflow**
<p align="left" width="100%">
    <img width="50%" src="https://github.com/wnbo9/st/blob/main/workflow.jpg">
</p>

1. **Image registration.** We first manually register the H&E image and the location image of barcodes. Then we crop the forebrain image from the mouse H&E image.
   - *forebrain.jpg*
   - *GI_tract.jpg*
2. **Cell segmentation.** We then apply the cell segmentation on the cropped forebrian H&E image. We need to get the labeled mask image of these cells.
   - *forebrain_mask_simple.jpg*
3. **Information extraction.** Using Python, we can specify all pixels with their corresponding cell, and the information of cells, including centroid, area and perimeter. These files can be viewed through [Dropbox](https://www.dropbox.com/sh/g8kfbqrtza3u9ex/AADTZNbPNKWDgHV4LsB__RCna?dl=0).
   - *forebrian_label_simple.csv*
   - *forebrain_info_simple.csv*
4. **Data aggregation.** We use *mask.csv* and *info.csv* to aggregate the barcodes within each cell and get a CGE matrix.
   - *forebrain_CGE_simple.csv*
5. **SpatialPCA analysis.** We run SpatialPCA to get the spatial PCs and the normalized expression matrix. Then we use louvain clustering with spatial PCs to get the domain cluster lables.
   - *spatial PCs*
   - *domain labels*
6. **Seurat analysis.** With the domain labels, we can create Seurat object. The following steps include 1) plotting UMAP with the labels, and 2) dectecting the differential genes among these domains. The SpatialPCA result and following Seurat object are stored in the RData file.
   - *SpatialPCA_forebrain_simple.RData*
7. **Ground-truth annotation.** With the plot of some marker genes, we can plot the ground truth of the tissue structure, and use Photoshop to annotate these structures.
   - *forebrain_annotation_simple.jpg*
8. **ARI calculation.** With the *annotation.jpg* and *info.csv*, we can specify the annotation labels for all cells as ground truth. Then with the *domain labels* from SpatialPCA, we canc calculate the ARI.
   - *ARI*
9. ...

