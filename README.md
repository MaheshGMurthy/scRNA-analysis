GSE126321 dataset was downloaded from NCBI-GEO database. The dataset had 3 samples GM18502, GM12878 and GSM3596322. 
scRNA-seq data was analsyed using the Seurat package. Raw data was imported using Read10X function. 
Seurat objective was created having a threshold of atleast 200 genes per cell and genes detected in at least 3 cells.
As for the QC mitochondrial gene content was calculated to identify the low-quality cells and cells with high mitochondrial RNA content (dead cells or cell stress). 
Cells with lesser 5% mitochondrial RNA was retained and between 200 to 2500 detected genes resulting in a clean dataset for further analysis. 
Data was then normalised using the lognormalise function from Seurat’s, followed by the identification of highly variable features. 
Data was scaled up for linear dimensionality reduction (PCA) and clustering followed by nonlinear dimensionality reduction (UMAP)
