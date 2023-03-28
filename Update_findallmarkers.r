#Obtaining genes that make up the leiden clusters identified by cell2location
#Section: 2.5.4 Seurat FindAllmarkers
#Extract metadata from the cell2location object using python
 
import scanpy as sc

adata = sc.read("/scratch/mshruv003/Updated_integrated_MTG_annotation/cell2location_map/sp_with_clusters.h5ad")

adata.obs.to_csv("/scratch/mshruv003/Updated_integrated_MTG_annotation/metadata_leiden_n21_res08.csv")

#R

library(ggplot2)
library(Matrix)
library(rjson)
library(cowplot)
library(RColorBrewer)
library(grid)
library(readbitmap)
library(Seurat)
library(dplyr)
library(hdf5r)
library(data.table)
library(tidyverse)

#Creating seurat objects
#Code adapted from https://satijalab.org/seurat/articles/spatial_vignette.html

DH2_1 <- Load10X_Spatial(
  data.dir = "/scratch/mshruv003/visium_18_01_22/X201SC21111924-Z01-F001_2/raw_data/run_spaceranger_count_DH2/outs",
  filename = "/filtered_feature_bc_matrix.h5",
  assay = "Spatial",
  slice = "DH2_2",
  filter.matrix = TRUE,
  to.upper = FALSE,
  image = NULL,
)

mt.genes <- rownames(DH2_1)[grep("^MT-",rownames(DH2_1))]
C<-GetAssayData(object = DH2_1, slot = "counts")
percent.mito <- colSums(C[mt.genes,])/Matrix::colSums(C)*100
DH2_1 <- AddMetaData(DH2_1, percent.mito, col.name = "percent.mito")
DH2_1 <- DH2_1[-(which(rownames(DH2_1) %in% as.factor(mt.genes))),]

DH2_1 <- SCTransform(DH2_1, assay = "Spatial", verbose = FALSE)

DH2_1 <- RunPCA(DH2_1, assay = "SCT", verbose = FALSE)
DH2_1 <- FindNeighbors(DH2_1, reduction = "pca", dims = 1:30)
DH2_1 <- FindClusters(DH2_1, verbose = FALSE, resolution = 0.8)
DH2_1 <- RunUMAP(DH2_1, reduction = "pca", dims = 1:30)

DimPlot(DH2_1, reduction = "umap", label = TRUE)
SpatialDimPlot(DH2_1, image.alpha = 0, label = FALSE, pt.size.factor = 1.4)

allen_reference <- readRDS("/scratch/mshruv003/visium_08_10_21/Allen_Brain_Data/allen.rds")
library(dplyr)
allen_reference <- SCTransform(allen_reference, ncells = 3000, verbose = FALSE) %>%
    RunPCA(verbose = FALSE) %>%
    RunUMAP(dims = 1:30)
	
anchors <- FindTransferAnchors(reference = allen_reference, query = DH2_1, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = allen_reference$subclass_label, prediction.assay = TRUE,
    weight.reduction = DH2_1[["pca"]], dims = 1:30)
DH2_1[["predictions"]] <- predictions.assay

DefaultAssay(DH2_1) <- "predictions"
SpatialFeaturePlot(DH2_1, features = c("L2/3 IT", "L4"), pt.size.factor = 1.6, ncol = 2, crop = TRUE)

saveRDS <- (DH2_1, file = "juvenile_4yr_DH2_1_seurat.rds")

DH3_1 <- Load10X_Spatial(
  data.dir = "/scratch/mshruv003/visium_18_01_22/X201SC21111924-Z01-F001_3/raw_data/run_spaceranger_count_DH3/outs",
  filename = "/filtered_feature_bc_matrix.h5",
  assay = "Spatial",
  slice = "DH3_1",
  filter.matrix = TRUE,
  to.upper = FALSE,
  image = NULL,
)

mt.genes <- rownames(DH3_1)[grep("^MT-",rownames(DH3_1))]
C<-GetAssayData(object = DH3_1, slot = "counts")
percent.mito <- colSums(C[mt.genes,])/Matrix::colSums(C)*100
DH3_1 <- AddMetaData(DH3_1, percent.mito, col.name = "percent.mito")
DH3_1 <- DH3_1[-(which(rownames(DH3_1) %in% as.factor(mt.genes))),]

DH3_1 <- SCTransform(DH3_1, assay = "Spatial", verbose = FALSE)

DH3_1 <- RunPCA(DH3_1, assay = "SCT", verbose = FALSE)
DH3_1 <- FindNeighbors(DH3_1, reduction = "pca", dims = 1:30)
DH3_1 <- FindClusters(DH3_1, verbose = FALSE, resolution = 0.8)
DH3_1 <- RunUMAP(DH3_1, reduction = "pca", dims = 1:30)

DimPlot(DH3_1, reduction = "umap", label = TRUE)
SpatialDimPlot(DH3_1, image.alpha = 0, label = FALSE, pt.size.factor = 1.4)

allen_reference <- readRDS("/scratch/mshruv003/visium_08_10_21/Allen_Brain_Data/allen.rds")
library(dplyr)
allen_reference <- SCTransform(allen_reference, ncells = 3000, verbose = FALSE) %>%
    RunPCA(verbose = FALSE) %>%
    RunUMAP(dims = 1:30)
	
anchors <- FindTransferAnchors(reference = allen_reference, query = DH3_1, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = allen_reference$subclass_label, prediction.assay = TRUE,
    weight.reduction = DH3_1[["pca"]], dims = 1:30)
DH3_1[["predictions"]] <- predictions.assay

DefaultAssay(DH3_1) <- "predictions"
SpatialFeaturePlot(DH3_1, features = c("L2/3 IT", "L4"), pt.size.factor = 1.6, ncol = 2, crop = TRUE)

saveRDS(DH3_1, file = "juvenile_2nd_15yr_DH3_1_seurat.rds")


DH4_1 <- Load10X_Spatial(
  data.dir = "/scratch/mshruv003/visium_18_01_22/X201SC21111924-Z01-F001_2/raw_data/run_spaceranger_count_DH4/outs",
  filename = "/filtered_feature_bc_matrix.h5",
  assay = "Spatial",
  slice = "DH4_2",
  filter.matrix = TRUE,
  to.upper = FALSE,
  image = NULL,
)

mt.genes <- rownames(DH4_1)[grep("^MT-",rownames(DH4_1))]
C<-GetAssayData(object = DH4_1, slot = "counts")
percent.mito <- colSums(C[mt.genes,])/Matrix::colSums(C)*100
DH4_1 <- AddMetaData(DH4_1, percent.mito, col.name = "percent.mito")
DH4_1 <- DH4_1[-(which(rownames(DH4_1) %in% as.factor(mt.genes))),]

DH4_1 <- SCTransform(DH4_1, assay = "Spatial", verbose = FALSE)

DH4_1 <- RunPCA(DH4_1, assay = "SCT", verbose = FALSE)
DH4_1 <- FindNeighbors(DH4_1, reduction = "pca", dims = 1:30)
DH4_1 <- FindClusters(DH4_1, verbose = FALSE, resolution = 0.8)
DH4_1 <- RunUMAP(DH4_1, reduction = "pca", dims = 1:30)

DimPlot(DH4_1, reduction = "umap", label = TRUE)
SpatialDimPlot(DH4_1, image.alpha = 0, label = FALSE, pt.size.factor = 1.4)

allen_reference <- readRDS("/scratch/mshruv003/visium_08_10_21/Allen_Brain_Data/allen.rds")
library(dplyr)
allen_reference <- SCTransform(allen_reference, ncells = 3000, verbose = FALSE) %>%
    RunPCA(verbose = FALSE) %>%
    RunUMAP(dims = 1:30)
	
anchors <- FindTransferAnchors(reference = allen_reference, query = DH4_1, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = allen_reference$subclass_label, prediction.assay = TRUE,
    weight.reduction = DH4_1[["pca"]], dims = 1:30)
DH4_1[["predictions"]] <- predictions.assay

DefaultAssay(DH4_1) <- "predictions"
SpatialFeaturePlot(DH4_1, features = c("L2/3 IT", "L4"), pt.size.factor = 1.6, ncol = 2, crop = TRUE)

saveRDS(DH4_1, file = "juvenile_2nd_15yr_DH4_1_seurat.rds")

DH1 <- Load10X_Spatial(
  data.dir = "/scratch/mshruv003/visium_08_10_21/Hockman_Visium/X201SC21082979-Z01-F001_01/raw_data/run_spaceranger_count_DH1a/outs",
  filename = "/filtered_feature_bc_matrix.h5",
  assay = "Spatial",
  slice = "DH1_1",
  filter.matrix = TRUE,
  to.upper = FALSE,
  image = NULL,
)

mt.genes <- rownames(DH1)[grep("^MT-",rownames(DH1))]
C<-GetAssayData(object = DH1, slot = "counts")
percent.mito <- colSums(C[mt.genes,])/Matrix::colSums(C)*100
DH1 <- AddMetaData(DH1, percent.mito, col.name = "percent.mito")
DH1 <- DH1[-(which(rownames(DH1) %in% as.factor(mt.genes))),]

DH1 <- SCTransform(DH1, assay = "Spatial", verbose = FALSE)

DH1 <- RunPCA(DH1, assay = "SCT", verbose = FALSE)
DH1 <- FindNeighbors(DH1, reduction = "pca", dims = 1:30)
DH1 <- FindClusters(DH1, verbose = FALSE, resolution = 0.8)
DH1 <- RunUMAP(DH1, reduction = "pca", dims = 1:30)

DimPlot(DH1, reduction = "umap", label = TRUE)
SpatialDimPlot(DH1, image.alpha = 0, label = FALSE, pt.size.factor = 1.4)

allen_reference <- readRDS("/scratch/mshruv003/visium_08_10_21/Allen_Brain_Data/allen.rds")
library(dplyr)
allen_reference <- SCTransform(allen_reference, ncells = 3000, verbose = FALSE) %>%
    RunPCA(verbose = FALSE) %>%
    RunUMAP(dims = 1:30)
	
anchors <- FindTransferAnchors(reference = allen_reference, query = DH1, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = allen_reference$subclass_label, prediction.assay = TRUE,
    weight.reduction = DH1[["pca"]], dims = 1:30)
DH1[["predictions"]] <- predictions.assay

DefaultAssay(DH1) <- "predictions"
SpatialFeaturePlot(DH1, features = c("L2/3 IT", "L4"), pt.size.factor = 1.6, ncol = 2, crop = TRUE)

saveRDS(DH1, file = "juvenile_first15yr_DH1_seurat.rds")


DH2 <- Load10X_Spatial(
  data.dir = "/scratch/mshruv003/visium_08_10_21/Hockman_Visium/X201SC21082979-Z01-F001_01/raw_data/run_spaceranger_count_DH2a/outs",
  filename = "/filtered_feature_bc_matrix.h5",
  assay = "Spatial",
  slice = "DH2_1",
  filter.matrix = TRUE,
  to.upper = FALSE,
  image = NULL,
)

mt.genes <- rownames(DH2)[grep("^MT-",rownames(DH2))]
C<-GetAssayData(object = DH2, slot = "counts")
percent.mito <- colSums(C[mt.genes,])/Matrix::colSums(C)*100
DH2 <- AddMetaData(DH2, percent.mito, col.name = "percent.mito")
DH2 <- DH2[-(which(rownames(DH2) %in% as.factor(mt.genes))),]

DH2 <- SCTransform(DH2, assay = "Spatial", verbose = FALSE)

DH2 <- RunPCA(DH2, assay = "SCT", verbose = FALSE)
DH2 <- FindNeighbors(DH2, reduction = "pca", dims = 1:30)
DH2 <- FindClusters(DH2, verbose = FALSE, resolution = 0.8)
DH2 <- RunUMAP(DH2, reduction = "pca", dims = 1:30)

DimPlot(DH2, reduction = "umap", label = TRUE)
SpatialDimPlot(DH2, image.alpha = 0, label = FALSE, pt.size.factor = 1.4)

allen_reference <- readRDS("/scratch/mshruv003/visium_08_10_21/Allen_Brain_Data/allen.rds")
library(dplyr)
allen_reference <- SCTransform(allen_reference, ncells = 3000, verbose = FALSE) %>%
    RunPCA(verbose = FALSE) %>%
    RunUMAP(dims = 1:30)
	
anchors <- FindTransferAnchors(reference = allen_reference, query = DH2, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = allen_reference$subclass_label, prediction.assay = TRUE,
    weight.reduction = DH2[["pca"]], dims = 1:30)
DH2[["predictions"]] <- predictions.assay

DefaultAssay(DH2) <- "predictions"
SpatialFeaturePlot(DH2, features = c("L2/3 IT", "L4"), pt.size.factor = 1.6, ncol = 2, crop = TRUE)

saveRDS(DH2, file = "juvenile_first15yr_DH2_seurat.rds")


DH3 <- Load10X_Spatial(
  data.dir = "/scratch/mshruv003/visium_08_10_21/Hockman_Visium/X201SC21082979-Z01-F001_02/raw_data/run_spaceranger_count_DH3a/outs",
  filename = "/filtered_feature_bc_matrix.h5",
  assay = "Spatial",
  slice = "DH3_1",
  filter.matrix = TRUE,
  to.upper = FALSE,
  image = NULL,
)

mt.genes <- rownames(DH3)[grep("^MT-",rownames(DH3))]
C<-GetAssayData(object = DH3, slot = "counts")
percent.mito <- colSums(C[mt.genes,])/Matrix::colSums(C)*100
DH3 <- AddMetaData(DH3, percent.mito, col.name = "percent.mito")
DH3 <- DH3[-(which(rownames(DH3) %in% as.factor(mt.genes))),]

DH3 <- SCTransform(DH3, assay = "Spatial", verbose = FALSE)

DH3 <- RunPCA(DH3, assay = "SCT", verbose = FALSE)
DH3 <- FindNeighbors(DH3, reduction = "pca", dims = 1:30)
DH3 <- FindClusters(DH3, verbose = FALSE, resolution = 0.8)
DH3 <- RunUMAP(DH3, reduction = "pca", dims = 1:30)

DimPlot(DH3, reduction = "umap", label = TRUE)
SpatialDimPlot(DH3, image.alpha = 0, label = FALSE, pt.size.factor = 1.4)

allen_reference <- readRDS("/scratch/mshruv003/visium_08_10_21/Allen_Brain_Data/allen.rds")
library(dplyr)
allen_reference <- SCTransform(allen_reference, ncells = 3000, verbose = FALSE) %>%
    RunPCA(verbose = FALSE) %>%
    RunUMAP(dims = 1:30)
	
anchors <- FindTransferAnchors(reference = allen_reference, query = DH3, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = allen_reference$subclass_label, prediction.assay = TRUE,
    weight.reduction = DH3[["pca"]], dims = 1:30)
DH3[["predictions"]] <- predictions.assay

DefaultAssay(DH3) <- "predictions"
SpatialFeaturePlot(DH3, features = c("L2/3 IT", "L4"), pt.size.factor = 1.6, ncol = 2, crop = TRUE)

saveRDS(DH3, file = "adult_31yr_DH3_seurat.rds")


DH4 <- Load10X_Spatial(
  data.dir = "/scratch/mshruv003/visium_08_10_21/Hockman_Visium/X201SC21082979-Z01-F001_02/raw_data/run_spaceranger_count_DH4a/outs",
  filename = "/filtered_feature_bc_matrix.h5",
  assay = "Spatial",
  slice = "DH4_1",
  filter.matrix = TRUE,
  to.upper = FALSE,
  image = NULL,
)

mt.genes <- rownames(DH4)[grep("^MT-",rownames(DH4))]
C<-GetAssayData(object = DH4, slot = "counts")
percent.mito <- colSums(C[mt.genes,])/Matrix::colSums(C)*100
DH4 <- AddMetaData(DH4, percent.mito, col.name = "percent.mito")
DH4 <- DH4[-(which(rownames(DH4) %in% as.factor(mt.genes))),]

DH4 <- SCTransform(DH4, assay = "Spatial", verbose = FALSE)

DH4 <- RunPCA(DH4, assay = "SCT", verbose = FALSE)
DH4 <- FindNeighbors(DH4, reduction = "pca", dims = 1:30)
DH4 <- FindClusters(DH4, verbose = FALSE, resolution = 0.8)
DH4 <- RunUMAP(DH4, reduction = "pca", dims = 1:30)

DimPlot(DH4, reduction = "umap", label = TRUE)
SpatialDimPlot(DH4, image.alpha = 0, label = FALSE, pt.size.factor = 1.4)

allen_reference <- readRDS("/scratch/mshruv003/visium_08_10_21/Allen_Brain_Data/allen.rds")
library(dplyr)
allen_reference <- SCTransform(allen_reference, ncells = 3000, verbose = FALSE) %>%
    RunPCA(verbose = FALSE) %>%
    RunUMAP(dims = 1:30)
	
anchors <- FindTransferAnchors(reference = allen_reference, query = DH4, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = allen_reference$subclass_label, prediction.assay = TRUE,
    weight.reduction = DH4[["pca"]], dims = 1:30)
DH4[["predictions"]] <- predictions.assay

DefaultAssay(DH4) <- "predictions"
SpatialFeaturePlot(DH4, features = c("L2/3 IT", "L4"), pt.size.factor = 1.6, ncol = 2, crop = TRUE)

saveRDS(DH4, file = "adult_31yr_DH4_seurat.rds")


#Adding leiden clusters into seurat objects
juvenile_4yr_DH1 <- readRDS("/scratch/mshruv003/Seurat_objects_converted_from_cell2loc/updated_visium_4yr_seurat.rds")
scanpy_4yr_DH1 <- read.csv("/scratch/mshruv003/Updated_integrated_MTG_annotation/metadata_leiden_n21_res08.csv")
leiden <- as.data.frame(scanpy_4yr_DH1$leiden)
rownames(leiden) <- scanpy_4yr_DH1$spot_id
leiden$barcode_id <- rownames(leiden)
leiden_new <- leiden[-grep('DH1a', rownames(leiden)), ]
subset_DH1 <- leiden[grep('DH1', rownames(leiden_new)), ]
clean_subset_DH1 <- word(rownames(subset_DH1), 2, sep="_")
rownames(subset_DH1) <- clean_subset_DH1
juvenile_4yr_DH1 <- AddMetaData(juvenile_4yr_DH1, subset_DH1, col.name = NULL)
saveRDS(juvenile_4yr_DH1, file = "juvenile_4yr_DH1_seurat_leiden.rds")

juvenile_4yr_DH2 <- readRDS("/scratch/mshruv003/Seurat_objects_converted_from_cell2loc/juvenile_4yr_DH2_1_seurat.rds")
scanpy_4yr_DH2 <- read.csv("/scratch/mshruv003/Updated_integrated_MTG_annotation/metadata_leiden_n21_res08.csv")
leiden <- as.data.frame(scanpy_4yr_DH2$leiden)
rownames(leiden) <- scanpy_4yr_DH2$spot_id
leiden$barcode_id <- rownames(leiden)
leiden_new <- leiden[-grep('DH2a', rownames(leiden)), ]
subset_DH2 <- leiden[grep('DH2', rownames(leiden_new)), ]
clean_subset_DH2 <- word(rownames(subset_DH2), 2, sep="_")
rownames(subset_DH2) <- clean_subset_DH2
juvenile_4yr_DH2 <- AddMetaData(juvenile_4yr_DH2, subset_DH2, col.name = NULL)
saveRDS(juvenile_4yr_DH2, file = "juvenile_4yr_DH2_seurat_leiden.rds")

juvenile_2nd_15yr_DH3 <- readRDS("/scratch/mshruv003/Seurat_objects_converted_from_cell2loc/juvenile_2nd_15yr_DH3_1_seurat.rds")
scanpy_15yr_DH3 <- read.csv("/scratch/mshruv003/Updated_integrated_MTG_annotation/metadata_leiden_n21_res08.csv")
leiden <- as.data.frame(scanpy_15yr_DH3$leiden)
rownames(leiden) <- scanpy_15yr_DH3$spot_id
leiden$barcode_id <- rownames(leiden)
leiden_new <- leiden[-grep('DH3a', rownames(leiden)), ]
subset_DH3 <- leiden[grep('DH3', rownames(leiden_new)), ]
clean_subset_DH3 <- word(rownames(subset_DH3), 2, sep="_")
rownames(subset_DH3) <- clean_subset_DH3
juvenile_2nd_15yr_DH3 <- AddMetaData(juvenile_2nd_15yr_DH3, subset_DH3, col.name = NULL)
saveRDS(juvenile_2nd_15yr_DH3, file = "juvenile_2nd_15yr_DH3_seurat_leiden.rds")

juvenile_2nd_15yr_DH4 <- readRDS("/scratch/mshruv003/Seurat_objects_converted_from_cell2loc/juvenile_2nd_15yr_DH4_1_seurat.rds")
scanpy_15yr_DH4 <- read.csv("/scratch/mshruv003/Updated_integrated_MTG_annotation/metadata_leiden_n21_res08.csv")
leiden <- as.data.frame(scanpy_15yr_DH4$leiden)
rownames(leiden) <- scanpy_15yr_DH4$spot_id
leiden$barcode_id <- rownames(leiden)
leiden_new <- leiden[-grep('DH4a', rownames(leiden)), ]
subset_DH4 <- leiden[grep('DH4', rownames(leiden_new)), ]
clean_subset_DH4 <- word(rownames(subset_DH4), 2, sep="_")
rownames(subset_DH4) <- clean_subset_DH4
juvenile_2nd_15yr_DH4 <- AddMetaData(juvenile_2nd_15yr_DH4, subset_DH4, col.name = NULL)
saveRDS(juvenile_2nd_15yr_DH4, file = "juvenile_2nd_15yr_DH4_seurat_leiden.rds")


juvenile_first15yr_DH1 <- readRDS("/scratch/mshruv003/Seurat_objects_converted_from_cell2loc/juvenile_first15yr_DH1_seurat.rds")
scanpy_15yr_DH1 <- read.csv("/scratch/mshruv003/Updated_integrated_MTG_annotation/metadata_leiden_n21_res08.csv")
leiden <- as.data.frame(scanpy_15yr_DH1$leiden)
rownames(leiden) <- scanpy_15yr_DH1$spot_id
leiden$barcode_id <- rownames(leiden)
subset_DH1a <- leiden[grep('DH1a', rownames(leiden)), ]
clean_subset_DH1a <- word(rownames(subset_DH1a), 2, sep="_")
rownames(subset_DH1a) <- clean_subset_DH1a
juvenile_first15yr_DH1 <- AddMetaData(juvenile_first15yr_DH1, subset_DH1a, col.name = NULL)
saveRDS(juvenile_first15yr_DH1, file = "juvenile_first15yr_DH1_seurat_leiden.rds")


juvenile_first15yr_DH2 <- readRDS("/scratch/mshruv003/Seurat_objects_converted_from_cell2loc/juvenile_first15yr_DH2_seurat.rds")
scanpy_15yr_DH2 <- read.csv("/scratch/mshruv003/Updated_integrated_MTG_annotation/metadata_leiden_n21_res08.csv")
leiden <- as.data.frame(scanpy_15yr_DH2$leiden)
rownames(leiden) <- scanpy_15yr_DH2$spot_id
leiden$barcode_id <- rownames(leiden)
subset_DH2a <- leiden[grep('DH2a', rownames(leiden)), ]
clean_subset_DH2a <- word(rownames(subset_DH2a), 2, sep="_")
rownames(subset_DH2a) <- clean_subset_DH2a
juvenile_first15yr_DH2 <- AddMetaData(juvenile_first15yr_DH2, subset_DH2a, col.name = NULL)
saveRDS(juvenile_first15yr_DH2, file = "juvenile_first15yr_DH2_seurat_leiden.rds")


adult_31yr_DH3 <- readRDS("/scratch/mshruv003/Seurat_objects_converted_from_cell2loc/adult_31yr_DH3_seurat.rds")
scanpy_31yr_DH3 <- read.csv("/scratch/mshruv003/Updated_integrated_MTG_annotation/metadata_leiden_n21_res08.csv")
leiden <- as.data.frame(scanpy_31yr_DH3$leiden)
rownames(leiden) <- scanpy_31yr_DH3$spot_id
leiden$barcode_id <- rownames(leiden)
subset_DH3a <- leiden[grep('DH3a', rownames(leiden)), ]
clean_subset_DH3a <- word(rownames(subset_DH3a), 2, sep="_")
rownames(subset_DH3a) <- clean_subset_DH3a
adult_31yr_DH3 <- AddMetaData(adult_31yr_DH3, subset_DH3a, col.name = NULL)
saveRDS(adult_31yr_DH3, file = "adult_31yr_DH3_seurat_leiden.rds")

adult_31yr_DH4 <- readRDS("/scratch/mshruv003/Seurat_objects_converted_from_cell2loc/adult_31yr_DH4_seurat.rds")
scanpy_31yr_DH4 <- read.csv("/scratch/mshruv003/Updated_integrated_MTG_annotation/metadata_leiden_n21_res08.csv")
leiden <- as.data.frame(scanpy_31yr_DH4$leiden)
rownames(leiden) <- scanpy_31yr_DH4$spot_id
leiden$barcode_id <- rownames(leiden)
subset_DH4a <- leiden[grep('DH4a', rownames(leiden)), ]
clean_subset_DH4a <- word(rownames(subset_DH4a), 2, sep="_")
rownames(subset_DH4a) <- clean_subset_DH4a
adult_31yr_DH4 <- AddMetaData(adult_31yr_DH4, subset_DH4a, col.name = NULL)
saveRDS(adult_31yr_DH4, file = "adult_31yr_DH4_seurat_leiden.rds")


juvenile_4yr_DH1 <- readRDS("/scratch/mshruv003/Seurat_objects_converted_from_cell2loc/updated_visium_4yr_seurat.rds")
juvenile_4yr_DH2 <- readRDS("/scratch/mshruv003/Seurat_objects_converted_from_cell2loc/juvenile_4yr_DH2_1_seurat.rds")
juvenile_2nd_15yr_DH3 <- readRDS("/scratch/mshruv003/Seurat_objects_converted_from_cell2loc/juvenile_2nd_15yr_DH3_1_seurat.rds")
juvenile_2nd_15yr_DH4 <- readRDS("/scratch/mshruv003/Seurat_objects_converted_from_cell2loc/juvenile_2nd_15yr_DH4_1_seurat.rds")
juvenile_first15yr_DH1 <- readRDS("/scratch/mshruv003/Seurat_objects_converted_from_cell2loc/juvenile_first15yr_DH1_seurat.rds")
juvenile_first15yr_DH2 <- readRDS("/scratch/mshruv003/Seurat_objects_converted_from_cell2loc/juvenile_first15yr_DH2_seurat.rds")
adult_31yr_DH3 <- readRDS("/scratch/mshruv003/Seurat_objects_converted_from_cell2loc/adult_31yr_DH3_seurat.rds")
adult_31yr_DH4 <- readRDS("/scratch/mshruv003/Seurat_objects_converted_from_cell2loc/adult_31yr_DH4_seurat.rds")

#Merging the seurat objects
brain.merge3 <- merge(juvenile_4yr_DH1, y = c(juvenile_4yr_DH2, juvenile_2nd_15yr_DH3, juvenile_2nd_15yr_DH4, juvenile_first15yr_DH1, juvenile_first15yr_DH2, adult_31yr_DH3, adult_31yr_DH4), add.cell.ids = c("DH1", "DH2", "DH3", "DH4", "DH1a", "DH2a", "DH3a", "DH4a"), project = "Visium") 
leiden_clusters<- read.csv("/scratch/mshruv003/Updated_integrated_MTG_annotation/metadata_leiden_n21_res08.csv")
leiden_cell2loc <- as.data.frame(leiden_clusters$leiden)
rownames(leiden_cell2loc) <- leiden_clusters$spot_id
leiden_cell2loc$barcode_id <- rownames(leiden_cell2loc)
brain.merge3 <- AddMetaData(brain.merge3, leiden_cell2loc, col.name = NULL)

brain.merge3@meta.data$sample <- NA
brain.merge3@meta.data$sample[which(str_detect(brain.merge3@meta.data$barcode_id,"^DH1"))] <- "4_year_old_T2_B2"
brain.merge3@meta.data$sample[which(str_detect(brain.merge3@meta.data$barcode_id,"^DH2"))] <- "4_year_old_T1_B2"
brain.merge3@meta.data$sample[which(str_detect(brain.merge3@meta.data$barcode_id,"^DH3"))] <- "15_year_old_T2_B2"
brain.merge3@meta.data$sample[which(str_detect(brain.merge3@meta.data$barcode_id,"^DH4"))] <- "15_year_old_T1_B2"
brain.merge3@meta.data$sample[which(str_detect(brain.merge3@meta.data$barcode_id,"^DH1a"))] <- "15_year_old_T1_B1"
brain.merge3@meta.data$sample[which(str_detect(brain.merge3@meta.data$barcode_id,"^DH2a"))] <- "15_year_old_T2_B1"
brain.merge3@meta.data$sample[which(str_detect(brain.merge3@meta.data$barcode_id,"^DH3a"))] <- "31_year_old_T1_B1"
brain.merge3@meta.data$sample[which(str_detect(brain.merge3@meta.data$barcode_id,"^DH4a"))] <- "31_year_old_T2_B1"


brain.merge3 <- SCTransform(brain.merge3, assay = "Spatial", verbose = FALSE)
Idents(brain.merge3) <- "leiden_clusters.leiden"

brain.merge3$Region_cluster <- brain.merge3$leiden_clusters.leiden

WM <- subset(brain.merge3, idents = c(6, 9, 10, 12, 15, 16))
WM <- SCTransform(WM, assay = "Spatial", verbose = FALSE)

markers <- FindAllMarkers(object = WM, 
                          slot = "data",
                          only.pos = TRUE,
                          logfc.threshold = 0.25) 	
						  
markers %>%
    group_by(cluster) %>%
    slice_max(n = 20, order_by = avg_log2FC)
top20 <- markers %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC)						 

Ave_exp <- AverageExpression(
                  object = WM,
                  return.seurat = TRUE,
                  group.by = c("leiden_clusters.leiden", "sample"),
                  slot = "data")		

						  
 
png("WM_dittoheatmap.png")
pdf("WM_dittoheatmap2.pdf")

Ave_exp@meta.data$sample <- NA
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^6_15_year_old_T1_B1"))] <- "15_year_old_T1_B1"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^6_15_year_old_T2_B1"))] <- "15_year_old_T2_B1"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^6_15_year_old_T1_B2"))] <- "15_year_old_T1_B2"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^6_15_year_old_T2_B2"))] <- "15_year_old_T2_B2"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^6_31_year_old_T1_B1"))] <- "31_year_old_T1_B1"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^6_31_year_old_T2_B1"))] <- "31_year_old_T2_B1"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^6_4_year_old_T1_B2"))] <- "4_year_old_T1_B2"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^6_4_year_old_T2_B2"))] <- "4_year_old_T2_B2"


Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^9_15_year_old_T1_B1"))] <- "15_year_old_T1_B1"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^9_15_year_old_T2_B1"))] <- "15_year_old_T2_B1"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^9_15_year_old_T1_B2"))] <- "15_year_old_T1_B2"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^9_15_year_old_T2_B2"))] <- "15_year_old_T2_B2"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^9_31_year_old_T1_B1"))] <- "31_year_old_T1_B1"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^9_31_year_old_T2_B1"))] <- "31_year_old_T2_B1"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^9_4_year_old_T1_B2"))] <- "4_year_old_T1_B2"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^9_4_year_old_T2_B2"))] <- "4_year_old_T2_B2"

Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^10_15_year_old_T1_B1"))] <- "15_year_old_T1_B1"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^10_15_year_old_T2_B1"))] <- "15_year_old_T2_B1"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^10_15_year_old_T1_B2"))] <- "15_year_old_T1_B2"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^10_15_year_old_T2_B2"))] <- "15_year_old_T2_B2"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^10_31_year_old_T1_B1"))] <- "31_year_old_T1_B1"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^10_31_year_old_T2_B1"))] <- "31_year_old_T2_B1"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^10_4_year_old_T1_B2"))] <- "4_year_old_T1_B2"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^10_4_year_old_T2_B2"))] <- "4_year_old_T2_B2"

Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^12_15_year_old_T1_B1"))] <- "15_year_old_T1_B1"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^12_15_year_old_T2_B1"))] <- "15_year_old_T2_B1"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^12_15_year_old_T1_B2"))] <- "15_year_old_T1_B2"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^12_15_year_old_T2_B2"))] <- "15_year_old_T2_B2"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^12_31_year_old_T1_B1"))] <- "31_year_old_T1_B1"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^12_31_year_old_T2_B1"))] <- "31_year_old_T2_B1"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^12_4_year_old_T1_B2"))] <- "4_year_old_T1_B2"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^12_4_year_old_T2_B2"))] <- "4_year_old_T2_B2"

Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^15_15_year_old_T1_B1"))] <- "15_year_old_T1_B1"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^15_15_year_old_T2_B1"))] <- "15_year_old_T2_B1"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^15_15_year_old_T1_B2"))] <- "15_year_old_T1_B2"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^15_15_year_old_T2_B2"))] <- "15_year_old_T2_B2"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^15_31_year_old_T1_B1"))] <- "31_year_old_T1_B1"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^15_31_year_old_T2_B1"))] <- "31_year_old_T2_B1"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^15_4_year_old_T1_B2"))] <- "4_year_old_T1_B2"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^15_4_year_old_T2_B2"))] <- "4_year_old_T2_B2"

Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^16_15_year_old_T1_B1"))] <- "15_year_old_T1_B1"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^16_15_year_old_T2_B1"))] <- "15_year_old_T2_B1"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^16_15_year_old_T1_B2"))] <- "15_year_old_T1_B2"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^16_15_year_old_T2_B2"))] <- "15_year_old_T2_B2"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^16_31_year_old_T1_B1"))] <- "31_year_old_T1_B1"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^16_31_year_old_T2_B1"))] <- "31_year_old_T2_B1"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^16_4_year_old_T1_B2"))] <- "4_year_old_T1_B2"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^16_4_year_old_T2_B2"))] <- "4_year_old_T2_B2"

Ave_exp$Region_cluster <- Ave_exp$orig.ident

dittoHeatmap(Ave_exp, top20$gene,
    annot.by = c("sample", "Region_cluster"),
	order.by = "Region_cluster", slot = "data", scale = "row", 
	fontsize_row = 4.0)
	
	
cortex <- subset(brain.merge3, idents = c(7, 14))
cortex <- SCTransform(cortex, assay = "Spatial", verbose = FALSE)

markers <- FindAllMarkers(object = cortex, 
                          slot = "data",
                          only.pos = TRUE,
                          logfc.threshold = 0.25) 	

markers %>%
    slice_max(n = 20, order_by = avg_log2FC)
top10 <- markers %>%
  top_n(n = 20, wt = avg_log2FC)
 						  
						 

Ave_exp <- AverageExpression(
                  object = cortex,
                  return.seurat = TRUE,
                  group.by = c("leiden_clusters.leiden", "sample"),
                  slot = "data")		
						  

pdf("cortex714.dittoHeatmap2.pdf")

Ave_exp@meta.data$sample <- NA
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^7_15_year_old_T1_B1"))] <- "15_year_old_T1_B1"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^7_15_year_old_T2_B1"))] <- "15_year_old_T2_B1"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^7_15_year_old_T1_B2"))] <- "15_year_old_T1_B2"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^7_15_year_old_T2_B2"))] <- "15_year_old_T2_B2"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^7_31_year_old_T1_B1"))] <- "31_year_old_T1_B1"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^7_31_year_old_T2_B1"))] <- "31_year_old_T2_B1"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^7_4_year_old_T1_B2"))] <- "4_year_old_T1_B2"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^7_4_year_old_T2_B2"))] <- "4_year_old_T2_B2"


Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^14_15_year_old_T1_B1"))] <- "15_year_old_T1_B1"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^14_15_year_old_T2_B1"))] <- "15_year_old_T2_B1"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^14_15_year_old_T1_B2"))] <- "15_year_old_T1_B2"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^14_15_year_old_T2_B2"))] <- "15_year_old_T2_B2"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^14_31_year_old_T1_B1"))] <- "31_year_old_T1_B1"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^14_31_year_old_T2_B1"))] <- "31_year_old_T2_B1"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^14_4_year_old_T1_B2"))] <- "4_year_old_T1_B2"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^14_4_year_old_T2_B2"))] <- "4_year_old_T2_B2"

Ave_exp$Region_cluster <- Ave_exp$orig.ident


dittoHeatmap(Ave_exp, top10$gene,
    annot.by = c("sample", "Region_cluster"),
	order.by = "Region_cluster", slot = "data", scale = "row", 
	fontsize_row = 5.0)
	
cortex <- subset(brain.merge3, idents = c(17, 18, 20))
cortex <- SCTransform(cortex, assay = "Spatial", verbose = FALSE)

markers <- FindAllMarkers(object = cortex, 
                          slot = "data",
                          only.pos = TRUE,
                          logfc.threshold = 0.25) 	

markers %>%
    slice_max(n = 20, order_by = avg_log2FC)
top10 <- markers %>%
  top_n(n = 20, wt = avg_log2FC)
 						  
						 

Ave_exp <- AverageExpression(
                  object = cortex,
                  return.seurat = TRUE,
                  group.by = c("leiden_clusters.leiden", "sample"),
                  slot = "data")		
						  

pdf("17_18_20cortex.dittoHeatmap2.pdf")

Ave_exp@meta.data$sample <- NA
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^7_15_year_old_T1_B1"))] <- "15_year_old_T1_B1"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^7_15_year_old_T2_B1"))] <- "15_year_old_T2_B1"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^7_15_year_old_T1_B2"))] <- "15_year_old_T1_B2"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^7_15_year_old_T2_B2"))] <- "15_year_old_T2_B2"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^7_31_year_old_T1_B1"))] <- "31_year_old_T1_B1"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^7_31_year_old_T2_B1"))] <- "31_year_old_T2_B1"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^7_4_year_old_T1_B2"))] <- "4_year_old_T1_B2"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^7_4_year_old_T2_B2"))] <- "4_year_old_T2_B2"

Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^1_15_year_old_T1_B1"))] <- "15_year_old_T1_B1"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^1_15_year_old_T2_B1"))] <- "15_year_old_T2_B1"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^1_15_year_old_T1_B2"))] <- "15_year_old_T1_B2"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^1_15_year_old_T2_B2"))] <- "15_year_old_T2_B2"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^1_31_year_old_T1_B1"))] <- "31_year_old_T1_B1"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^1_31_year_old_T2_B1"))] <- "31_year_old_T2_B1"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^1_4_year_old_T1_B2"))] <- "4_year_old_T1_B2"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^1_4_year_old_T2_B2"))] <- "4_year_old_T2_B2"

Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^2_15_year_old_T1_B1"))] <- "15_year_old_T1_B1"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^2_15_year_old_T2_B1"))] <- "15_year_old_T2_B1"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^2_15_year_old_T1_B2"))] <- "15_year_old_T1_B2"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^2_15_year_old_T2_B2"))] <- "15_year_old_T2_B2"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^2_31_year_old_T1_B1"))] <- "31_year_old_T1_B1"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^2_31_year_old_T2_B1"))] <- "31_year_old_T2_B1"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^2_4_year_old_T1_B2"))] <- "4_year_old_T1_B2"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^2_4_year_old_T2_B2"))] <- "4_year_old_T2_B2"

Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^3_15_year_old_T1_B1"))] <- "15_year_old_T1_B1"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^3_15_year_old_T2_B1"))] <- "15_year_old_T2_B1"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^3_15_year_old_T1_B2"))] <- "15_year_old_T1_B2"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^3_15_year_old_T2_B2"))] <- "15_year_old_T2_B2"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^3_31_year_old_T1_B1"))] <- "31_year_old_T1_B1"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^3_31_year_old_T2_B1"))] <- "31_year_old_T2_B1"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^3_4_year_old_T1_B2"))] <- "4_year_old_T1_B2"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^3_4_year_old_T2_B2"))] <- "4_year_old_T2_B2"

Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^0_15_year_old_T1_B1"))] <- "15_year_old_T1_B1"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^0_15_year_old_T2_B1"))] <- "15_year_old_T2_B1"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^0_15_year_old_T1_B2"))] <- "15_year_old_T1_B2"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^0_15_year_old_T2_B2"))] <- "15_year_old_T2_B2"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^0_31_year_old_T1_B1"))] <- "31_year_old_T1_B1"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^0_31_year_old_T2_B1"))] <- "31_year_old_T2_B1"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^0_4_year_old_T1_B2"))] <- "4_year_old_T1_B2"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^0_4_year_old_T2_B2"))] <- "4_year_old_T2_B2"

Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^4_15_year_old_T1_B1"))] <- "15_year_old_T1_B1"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^4_15_year_old_T2_B1"))] <- "15_year_old_T2_B1"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^4_15_year_old_T1_B2"))] <- "15_year_old_T1_B2"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^4_15_year_old_T2_B2"))] <- "15_year_old_T2_B2"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^4_31_year_old_T1_B1"))] <- "31_year_old_T1_B1"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^4_31_year_old_T2_B1"))] <- "31_year_old_T2_B1"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^4_4_year_old_T1_B2"))] <- "4_year_old_T1_B2"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^4_4_year_old_T2_B2"))] <- "4_year_old_T2_B2"

Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^5_15_year_old_T1_B1"))] <- "15_year_old_T1_B1"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^5_15_year_old_T2_B1"))] <- "15_year_old_T2_B1"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^5_15_year_old_T1_B2"))] <- "15_year_old_T1_B2"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^5_15_year_old_T2_B2"))] <- "15_year_old_T2_B2"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^5_31_year_old_T1_B1"))] <- "31_year_old_T1_B1"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^5_31_year_old_T2_B1"))] <- "31_year_old_T2_B1"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^5_4_year_old_T1_B2"))] <- "4_year_old_T1_B2"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^5_4_year_old_T2_B2"))] <- "4_year_old_T2_B2"

Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^8_15_year_old_T1_B1"))] <- "15_year_old_T1_B1"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^8_15_year_old_T2_B1"))] <- "15_year_old_T2_B1"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^8_15_year_old_T1_B2"))] <- "15_year_old_T1_B2"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^8_15_year_old_T2_B2"))] <- "15_year_old_T2_B2"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^8_31_year_old_T1_B1"))] <- "31_year_old_T1_B1"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^8_31_year_old_T2_B1"))] <- "31_year_old_T2_B1"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^8_4_year_old_T1_B2"))] <- "4_year_old_T1_B2"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^8_4_year_old_T2_B2"))] <- "4_year_old_T2_B2"

Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^11_15_year_old_T1_B1"))] <- "15_year_old_T1_B1"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^11_15_year_old_T2_B1"))] <- "15_year_old_T2_B1"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^11_15_year_old_T1_B2"))] <- "15_year_old_T1_B2"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^11_15_year_old_T2_B2"))] <- "15_year_old_T2_B2"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^11_31_year_old_T1_B1"))] <- "31_year_old_T1_B1"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^11_31_year_old_T2_B1"))] <- "31_year_old_T2_B1"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^11_4_year_old_T1_B2"))] <- "4_year_old_T1_B2"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^11_4_year_old_T2_B2"))] <- "4_year_old_T2_B2"

Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^13_15_year_old_T1_B1"))] <- "15_year_old_T1_B1"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^13_15_year_old_T2_B1"))] <- "15_year_old_T2_B1"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^13_15_year_old_T1_B2"))] <- "15_year_old_T1_B2"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^13_15_year_old_T2_B2"))] <- "15_year_old_T2_B2"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^13_31_year_old_T1_B1"))] <- "31_year_old_T1_B1"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^13_31_year_old_T2_B1"))] <- "31_year_old_T2_B1"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^13_4_year_old_T1_B2"))] <- "4_year_old_T1_B2"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^13_4_year_old_T2_B2"))] <- "4_year_old_T2_B2"


Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^14_15_year_old_T1_B1"))] <- "15_year_old_T1_B1"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^14_15_year_old_T2_B1"))] <- "15_year_old_T2_B1"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^14_15_year_old_T1_B2"))] <- "15_year_old_T1_B2"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^14_15_year_old_T2_B2"))] <- "15_year_old_T2_B2"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^14_31_year_old_T1_B1"))] <- "31_year_old_T1_B1"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^14_31_year_old_T2_B1"))] <- "31_year_old_T2_B1"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^14_4_year_old_T1_B2"))] <- "4_year_old_T1_B2"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^14_4_year_old_T2_B2"))] <- "4_year_old_T2_B2"

Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^17_15_year_old_T1_B1"))] <- "15_year_old_T1_B1"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^17_15_year_old_T2_B1"))] <- "15_year_old_T2_B1"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^17_15_year_old_T1_B2"))] <- "15_year_old_T1_B2"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^17_15_year_old_T2_B2"))] <- "15_year_old_T2_B2"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^17_31_year_old_T1_B1"))] <- "31_year_old_T1_B1"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^17_31_year_old_T2_B1"))] <- "31_year_old_T2_B1"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^17_4_year_old_T1_B2"))] <- "4_year_old_T1_B2"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^17_4_year_old_T2_B2"))] <- "4_year_old_T2_B2"

Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^18_15_year_old_T1_B1"))] <- "15_year_old_T1_B1"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^18_15_year_old_T2_B1"))] <- "15_year_old_T2_B1"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^18_15_year_old_T1_B2"))] <- "15_year_old_T1_B2"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^18_15_year_old_T2_B2"))] <- "15_year_old_T2_B2"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^18_31_year_old_T1_B1"))] <- "31_year_old_T1_B1"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^18_31_year_old_T2_B1"))] <- "31_year_old_T2_B1"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^18_4_year_old_T1_B2"))] <- "4_year_old_T1_B2"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^18_4_year_old_T2_B2"))] <- "4_year_old_T2_B2"

Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^19_15_year_old_T1_B1"))] <- "15_year_old_T1_B1"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^19_15_year_old_T2_B1"))] <- "15_year_old_T2_B1"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^19_15_year_old_T1_B2"))] <- "15_year_old_T1_B2"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^19_15_year_old_T2_B2"))] <- "15_year_old_T2_B2"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^19_31_year_old_T1_B1"))] <- "31_year_old_T1_B1"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^19_31_year_old_T2_B1"))] <- "31_year_old_T2_B1"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^19_4_year_old_T1_B2"))] <- "4_year_old_T1_B2"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^19_4_year_old_T2_B2"))] <- "4_year_old_T2_B2"


Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^20_15_year_old_T1_B1"))] <- "15_year_old_T1_B1"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^20_15_year_old_T2_B1"))] <- "15_year_old_T2_B1"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^20_15_year_old_T1_B2"))] <- "15_year_old_T1_B2"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^20_15_year_old_T2_B2"))] <- "15_year_old_T2_B2"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^20_31_year_old_T1_B1"))] <- "31_year_old_T1_B1"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^20_31_year_old_T2_B1"))] <- "31_year_old_T2_B1"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^20_4_year_old_T1_B2"))] <- "4_year_old_T1_B2"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^20_4_year_old_T2_B2"))] <- "4_year_old_T2_B2"


Ave_exp$Region_cluster <- Ave_exp$orig.ident


dittoHeatmap(Ave_exp, top10$gene,
    annot.by = c("sample", "Region_cluster"),
	order.by = "Region_cluster", slot = "data", scale = "row", 
	fontsize_row = 2.5)
	
