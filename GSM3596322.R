install.packages("Seurat")
install.packages("ggplot2")
install.packages("dplyr")
install.packages("patchwork")

#load the librabries
library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)

# Set the correct directory path where the files are located
data_GSM3596322 <- Read10X(data.dir = "C:/Users/manum/Downloads/GSM3596322")

# Check the dimensions of the matrix (rows = genes, columns = cells)
dim(data_GSM3596322)

# Check the first few rows and columns
data_GSM3596322[1:10, 1:10]

#creating the seurat_object
seurat_object <- CreateSeuratObject(counts = data_GSM3596322, project = "CancerData", min.cells = 3, min.features = 200)
head(seurat_object)

#read the celllineage data
lineage_data <- read.table("C:/Users/manum/Downloads/GSM3596322/cellLineage.tsv.gz", header = TRUE, sep = "\t")

# Set the cell barcodes in lineage_data as rownames to match Seurat object barcodes
rownames(lineage_data) <- lineage_data[, 1]  

# Add lineage data as metadata
seurat_object$lineage <- lineage_data[, 2]  
# Check the first few rows of metadata in the Seurat object
head(seurat_object@meta.data)


#calculate the percentage of mitochondrial genes
seurat_object[["percent.mt"]]<-PercentageFeatureSet(seurat_object, pattern = "^MT-")
head(seurat_object)
dim(seurat_object)

#visulaise the data before QC
VlnPlot(seurat_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),ncol = 3)

#QC for the data
seurat_object<- subset(seurat_object, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
dim(seurat_object)

#visualise after QC
VlnPlot(seurat_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")+
  geom_smooth(method = 'lm')
FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "percent.mt")

#normalization
seurat_object <- NormalizeData(seurat_object, normalization.method = "LogNormalize")

#find the highly variable features
seurat_object <- FindVariableFeatures(seurat_object, selection.method = "vst", nfeatures = 2000)

#identify the top 10 highly variable features
top10<- head(VariableFeatures(seurat_object), 10)

#plot the highly variable features
plot_high<- VariableFeaturePlot(seurat_object) 
LabelPoints(plot = plot_high, points = top10, repel = TRUE)

#scaling the data
allgenes <- rownames(seurat_object)
seurat_object <- ScaleData(seurat_object, features = allgenes)

#Perform the linear dimensionality reduction
seurat_object <- RunPCA(seurat_object, features = VariableFeatures(object = seurat_object))

#visualise the PCA results
print(seurat_object[["pca"]], dims = 1:5, nfeatures = 5)
DimHeatmap(seurat_object, dims = 1, cells = 500, balanced = TRUE)

#determine the dimensionality of the data
ElbowPlot(seurat_object)

#clustering the cells
seurat_object <- FindNeighbors(seurat_object, dims = 1:10)

#understanding resolution
seurat_object<- FindClusters(seurat_object, resolution = c(0.1, 0.3, 0.5, 0.7, 1))
View(seurat_object@meta.data)

DimPlot(seurat_object, group.by = "RNA_snn_res.0.3", label = TRUE)

#setting identity of clusters
Idents(seurat_object)<- "RNA_snn_res.0.3"


#non-linear dimensionality reduction
seurat_object <- RunUMAP(seurat_object, dims = 1:10)
#individual clusters
DimPlot(seurat_object, reduction = "umap")

DimPlot(seurat_object, group.by = "lineage")
