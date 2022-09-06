# **Comparison of segmentation methods for spatial transcriptomic data analysis** #

### **Background**
With the advance of spatail transcriptomics, some recent technologies can achieve single cell even subcellular resolution for the spatail transcriptomic data, including Stereo-seq, HDST, Seq-Scope and etc. However, accurate characterization of cell boundaries and matching spots to cells for these technologies remain as an unsolved problem.

Cell segmentation has been a popular topic these years. Taking advantage of histology images, we can detect cell boundaries and segment these cells with small sizes. There exist some commonly used segmentation methods, and they can be divided into three classes based on their mechanisms: 
- **Simple segmentation**
- **Marker-controlled watershed segmentation**
  - ImageJ
  - Cellprofiler
  - Scikit-image
- **Machine learning based**
  - Cellpose
  - StarDist

Here, we make use of the mouse embryo tissue (E14.5 E1S3) from Stereo-seq to compare the performance and illustrate the characteristics of these segmentation methods.

### **Method**

We first 
