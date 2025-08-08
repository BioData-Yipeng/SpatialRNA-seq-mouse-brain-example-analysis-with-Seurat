getwd()
setwd('F:/data science/R/spatialRNA-seq')
library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)
install.packages("Rfast2")

#creat seurat object
seurat <- Load10X_Spatial(data.dir = "F:/data science/R/spatialRNA-seq/data",
  filename = "CytAssist_FFPE_Mouse_Brain_Rep1_filtered_feature_bc_matrix.h5",
  assay = "spatial",
)

View(seurat@meta.data)

#check expression numbers from each position/dots
plot1 <- VlnPlot(seurat, features = "nCount_spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(seurat, features = "nCount_spatial", pt.size.factor = 3) + theme(legend.position = "right")
wrap_plots(plot1, plot2)

#data normalization and scaling
seurat <- SCTransform(seurat, assay = "spatial", verbose = TRUE)

#data visualization
#expression of specific genes on slides
SpatialFeaturePlot(seurat, pt.size.factor = 2, features = c("Prkcd", "Mbp"))
SpatialFeaturePlot(seurat, pt.size.factor = 3, alpha = c(0.1,0.5), features = c("Enpp2", "Plp1"))
#pt.size.factor- This will scale the size of the spots. Default is 1.6
#alpha - minimum and maximum transparency. Default is c(1, 1).
plot <- SpatialFeaturePlot(seurat, pt.size.factor = 3, features = c("Enpp2")) + theme(legend.text = element_text(size = 15),
                                                                  legend.title = element_text(size = 15),
                                                                  legend.key.size = unit(0.5, "cm"))
#save as jpeg files
jpeg(filename = "Enpp2.jpg", height = 700, width = 1200, quality = 50)
print(plot)
dev.off()

#data clustering
seurat <- RunPCA(seurat, assay = "SCT", verbose = TRUE)
seurat <- FindNeighbors(seurat, reduction = "pca", dims = 1:30)
seurat <- FindClusters(seurat, verbose = TRUE)
seurat <- RunUMAP(seurat, reduction = "pca", dims = 1:30)

#data visualization show both umap graph and sptial group label with clusters
p1 <- DimPlot(seurat, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(seurat, label = TRUE, label.size = 3)
p1 + p2

#highlight clusters
SpatialDimPlot(seurat, pt.size.factor=2, cells.highlight = CellsByIdentities(object = seurat, idents = c(2, 1, 4)),
               facet.highlight = TRUE, ncol = 3)

#identification spatial variant features
de_markers <- FindMarkers(seurat, ident.1 = 5, ident.2 = 6)
SpatialFeaturePlot(object = seurat, pt.size.factor=3, features = rownames(de_markers)[1:3], alpha = c(1, 1), ncol = 3)

#search for features exhibiting spatial patterning in the absence of pre-annotation
seurat <- FindSpatiallyVariableFeatures(seurat, assay = "SCT", features = VariableFeatures(seurat)[1:1000],
                                       selection.method = "moransi")
top.features <- head(SpatiallyVariableFeatures(seurat, selection.method = "moransi"), 6)
SpatialFeaturePlot(seurat, features = top.features, ncol = 3, alpha = c(1, 1))


#subset to cluster of interest
seurat_subset <- subset(seurat, idents = c(1, 3))
SpatialDimPlot(seurat_subset, pt.size.factor = 4, label = TRUE, label.size = 3)

