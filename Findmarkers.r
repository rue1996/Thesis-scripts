import scanpy as sc

adata = sc.read("/scratch/mshruv003/Updated_integrated_MTG_annotation/cell2location_map/sp_with_clusters.h5ad")

adata.obs.to_csv("/scratch/mshruv003/Updated_integrated_MTG_annotation/metadata_leiden_n21_res08.csv")

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




adult_31yr_DH4_seurat <- readRDS("/scratch/mshruv003/Updated_integrated_MTG_annotation/adult_31yr_DH4_seurat_leiden.rds")
adult_31yr_DH4_seurat <- SCTransform(adult_31yr_DH4_seurat, assay = "Spatial", verbose = FALSE)
Idents(adult_31yr_DH4_seurat) <- "scanpy_31yr_DH4.leiden"
cortex <- subset(adult_31yr_DH4_seurat, idents = c(1, 2, 3, 0, 4, 5, 7, 8, 11, 14, 12, 17, 18, 19, 20))
cortex <- SCTransform(cortex, assay = "Spatial", verbose = FALSE)
markers <- FindAllMarkers(object = cortex, 
                          slot = "data",
                          only.pos = TRUE,
                          logfc.threshold = 0.01) 	

write.csv(markers, 
          file = "cortex_markers_31yr_dh4.csv", 
          quote = FALSE, 
          row.names = FALSE)
						  
markers %>%
    group_by(cluster) %>%
    slice_max(n = 2, order_by = avg_log2FC)
top20 <- markers %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC)
 
pdf("doheatmap_31yr_DH4a_cortex.pdf")

DoHeatmap(object = cortex, features = top20$gene, slot = "data", size = 2.5) + theme(axis.text.y = element_text(size = 2))

WM <- subset(adult_31yr_DH4_seurat, idents = c(6, 9, 10, 15, 16, 13))
WM <- SCTransform(WM, assay = "Spatial", verbose = FALSE)


markers1 <- FindAllMarkers(object = WM, 
                          slot = "data",
                          only.pos = TRUE,
                          logfc.threshold = 0.01) 
	
write.csv(markers1, 
          file = "WM_markers_31yr_dh4.csv", 
          quote = FALSE, 
          row.names = FALSE)
	  
markers1 %>%
    group_by(cluster) %>%
    slice_max(n = 2, order_by = avg_log2FC)
top20 <- markers1 %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC)
 
pdf("doheatmap_31yr_DH4a_WM.pdf")

DoHeatmap(object = WM, features = top20$gene, slot = "data", size = 2.5) + theme(axis.text.y = element_text(size = 2))




adult_31yr_DH3_seurat <- readRDS("/scratch/mshruv003/Updated_integrated_MTG_annotation/adult_31yr_DH3_seurat_leiden.rds")
adult_31yr_DH3_seurat <- SCTransform(adult_31yr_DH3_seurat, assay = "Spatial", verbose = FALSE)
Idents(adult_31yr_DH3_seurat) <- "scanpy_31yr_DH3.leiden"
cortex <- subset(adult_31yr_DH3_seurat, idents = c(1, 2, 3, 0, 4, 5, 7, 8, 11, 14, 17, 18))
cortex <- SCTransform(cortex, assay = "Spatial", verbose = FALSE)
markers <- FindAllMarkers(object = cortex, 
                          slot = "data",
                          only.pos = TRUE,
                          logfc.threshold = 0.01) 	

write.csv(markers, 
          file = "cortex_markers_31yr_dh3.csv", 
          quote = FALSE, 
          row.names = FALSE)
		  
markers %>%
    group_by(cluster) %>%
    slice_max(n = 2, order_by = avg_log2FC)
top20 <- markers %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC)
 
pdf("doheatmap_31yr_DH3a_cortex.pdf")

DoHeatmap(object = cortex, features = top20$gene, slot = "data", size = 2.5) + theme(axis.text.y = element_text(size = 2))

WM <- subset(adult_31yr_DH3_seurat, idents = c(6, 9, 10, 12, 15, 16, 13))
WM <- SCTransform(WM, assay = "Spatial", verbose = FALSE)
markers1 <- FindAllMarkers(object = WM, 
                          slot = "data",
                          only.pos = TRUE,
                          logfc.threshold = 0.01) 
		  
		  
write.csv(markers1, 
          file = "WM_markers_31yr_dh3.csv", 
          quote = FALSE, 
          row.names = FALSE)
		  
markers1 %>%
    group_by(cluster) %>%
    slice_max(n = 2, order_by = avg_log2FC)
top20 <- markers1 %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC)
 
pdf("doheatmap_31yr_DH3a_WM.pdf")

DoHeatmap(object = WM, features = top20$gene, slot = "data", size = 2.5) + theme(axis.text.y = element_text(size = 2))


juvenile_first15yr_DH1_seurat <- readRDS("/scratch/mshruv003/Updated_integrated_MTG_annotation/juvenile_first15yr_DH1_seurat_leiden.rds")
juvenile_first15yr_DH1_seurat <- SCTransform(juvenile_first15yr_DH1_seurat, assay = "Spatial", verbose = FALSE)
Idents(juvenile_first15yr_DH1_seurat) <- "scanpy_15yr_DH1.leiden"
cortex <- subset(juvenile_first15yr_DH1_seurat, idents = c(1, 2, 3, 0, 4, 5, 7, 8, 11, 14, 17, 18, 19, 20))
cortex <- SCTransform(cortex, assay = "Spatial", verbose = FALSE)
markers <- FindAllMarkers(object = cortex, 
                          slot = "data",
                          only.pos = TRUE,
                          logfc.threshold = 0.01) 	

write.csv(markers, 
          file = "cortex_markers_first15yr_dh1.csv", 
          quote = FALSE, 
          row.names = FALSE)
						  
markers %>%
    group_by(cluster) %>%
    slice_max(n = 2, order_by = avg_log2FC)
top20 <- markers %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC)
 
pdf("doheatmap_first15yr_DH1a_cortex.pdf")

DoHeatmap(object = cortex, features = top20$gene, slot = "data", size = 2.5) + theme(axis.text.y = element_text(size = 2))

WM <- subset(juvenile_first15yr_DH1_seurat, idents = c(6, 9, 10, 12, 15, 16, 13))
WM <- SCTransform(WM, assay = "Spatial", verbose = FALSE)


markers1 <- FindAllMarkers(object = WM, 
                          slot = "data",
                          only.pos = TRUE,
                          logfc.threshold = 0.01) 
	

write.csv(markers1, 
          file = "WM_markers_first15yr_dh1.csv", 
          quote = FALSE, 
          row.names = FALSE)
		  
markers1 %>%
    group_by(cluster) %>%
    slice_max(n = 2, order_by = avg_log2FC)
top20 <- markers1 %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC)
 
pdf("doheatmap_first15yr_DH1a_WM.pdf")

DoHeatmap(object = WM, features = top20$gene, slot = "data", size = 2.5) + theme(axis.text.y = element_text(size = 2))


juvenile_first15yr_DH2_seurat <- readRDS("/scratch/mshruv003/Updated_integrated_MTG_annotation/juvenile_first15yr_DH2_seurat_leiden.rds")
juvenile_first15yr_DH2_seurat <- SCTransform(juvenile_first15yr_DH2_seurat, assay = "Spatial", verbose = FALSE)
Idents(juvenile_first15yr_DH2_seurat) <- "scanpy_15yr_DH2.leiden"
cortex <- subset(juvenile_first15yr_DH2_seurat, idents = c(1, 2, 3, 0, 4, 5, 7, 8, 11, 14, 17, 18, 19, 20))
cortex <- SCTransform(cortex, assay = "Spatial", verbose = FALSE)
markers <- FindAllMarkers(object = cortex, 
                          slot = "data",
                          only.pos = TRUE,
                          logfc.threshold = 0.01) 	  
						  
write.csv(markers, 
          file = "cortex_markers_first15yr_dh2.csv", 
          quote = FALSE, 
          row.names = FALSE)
		  
markers %>%
    group_by(cluster) %>%
    slice_max(n = 2, order_by = avg_log2FC)
top20 <- markers %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC)
 
pdf("doheatmap_first15yr_DH2a_cortex.pdf")

DoHeatmap(object = cortex, features = top20$gene, slot = "data", size = 2.5) + theme(axis.text.y = element_text(size = 2))

WM <- subset(juvenile_first15yr_DH2_seurat, idents = c(6, 9, 10, 12, 13, 15, 16))
WM <- SCTransform(WM, assay = "Spatial", verbose = FALSE)


markers1 <- FindAllMarkers(object = WM, 
                          slot = "data",
                          only.pos = TRUE,
                          logfc.threshold = 0.01) 
	

write.csv(markers1, 
          file = "WM_markers_first15yr_dh2.csv", 
          quote = FALSE, 
          row.names = FALSE)
		  
markers1 %>%
    group_by(cluster) %>%
    slice_max(n = 2, order_by = avg_log2FC)
top20 <- markers1 %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC)
 
pdf("doheatmap_first15yr_DH2a_WM.pdf")

DoHeatmap(object = WM, features = top20$gene, slot = "data", size = 2.5) + theme(axis.text.y = element_text(size = 2))


juvenile_2nd_15yr_DH3_seurat <- readRDS("/scratch/mshruv003/Updated_integrated_MTG_annotation/juvenile_2nd_15yr_DH3_seurat_leiden.rds")
juvenile_2nd_15yr_DH3_seurat <- SCTransform(juvenile_2nd_15yr_DH3_seurat, assay = "Spatial", verbose = FALSE)
Idents(juvenile_2nd_15yr_DH3_seurat) <- "scanpy_15yr_DH3.leiden"
cortex <- subset(juvenile_2nd_15yr_DH3_seurat, idents = c(1, 2, 3, 0, 4, 5, 7, 8, 11, 14, 13, 17, 18, 19))
cortex <- SCTransform(cortex, assay = "Spatial", verbose = FALSE)
markers <- FindAllMarkers(object = cortex, 
                          slot = "data",
                          only.pos = TRUE,
                          logfc.threshold = 0.01) 

write.csv(markers, 
          file = "cortex_markers_2nd_15yr_dh3.csv", 
          quote = FALSE, 
          row.names = FALSE)
		  
markers %>%
    group_by(cluster) %>%
    slice_max(n = 2, order_by = avg_log2FC)
top20 <- markers %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC)
 
pdf("doheatmap_2nd_15yr_DH3_cortex.pdf")

DoHeatmap(object = cortex, features = top20$gene, slot = "data", size = 2.5) + theme(axis.text.y = element_text(size = 2))

WM <- subset(juvenile_2nd_15yr_DH3_seurat, idents = c(6, 9, 10, 12, 15, 16))
WM <- SCTransform(WM, assay = "Spatial", verbose = FALSE)


markers1 <- FindAllMarkers(object = WM, 
                          slot = "data",
                          only.pos = TRUE,
                          logfc.threshold = 0.01) 
	
write.csv(markers1, 
          file = "WM_markers_2nd_15yr_dh3.csv", 
          quote = FALSE, 
          row.names = FALSE)
		  
markers1 %>%
    group_by(cluster) %>%
    slice_max(n = 2, order_by = avg_log2FC)
top20 <- markers1 %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC)
 
pdf("doheatmap_2nd_15yr_DH3_WM.pdf")

DoHeatmap(object = WM, features = top20$gene, slot = "data", size = 2.5) + theme(axis.text.y = element_text(size = 2))


juvenile_2nd_15yr_DH4_seurat <- readRDS("/scratch/mshruv003/Updated_integrated_MTG_annotation/juvenile_2nd_15yr_DH4_seurat_leiden.rds")
juvenile_2nd_15yr_DH4_seurat <- SCTransform(juvenile_2nd_15yr_DH4_seurat, assay = "Spatial", verbose = FALSE)
Idents(juvenile_2nd_15yr_DH4_seurat) <- "scanpy_15yr_DH4.leiden"
cortex <- subset(juvenile_2nd_15yr_DH4_seurat, idents = c(1, 2, 3, 0, 4, 5, 7, 8, 11, 14, 13, 17, 18))
cortex <- SCTransform(cortex, assay = "Spatial", verbose = FALSE)
markers <- FindAllMarkers(object = cortex, 
                          slot = "data",
                          only.pos = TRUE,
                          logfc.threshold = 0.01) 

write.csv(markers, 
          file = "cortex_markers_2nd_15yr_dh4.csv", 
          quote = FALSE, 
          row.names = FALSE)
		  
markers %>%
    group_by(cluster) %>%
    slice_max(n = 2, order_by = avg_log2FC)
top20 <- markers %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC)
 
pdf("doheatmap_2nd_15yr_DH4_cortex.pdf")

DoHeatmap(object = cortex, features = top20$gene, slot = "data", size = 2.5) + theme(axis.text.y = element_text(size = 2))

WM <- subset(juvenile_2nd_15yr_DH4_seurat, idents = c(6, 9, 10, 12, 15, 13, 16))
WM <- SCTransform(WM, assay = "Spatial", verbose = FALSE)


markers1 <- FindAllMarkers(object = WM, 
                          slot = "data",
                          only.pos = TRUE,
                          logfc.threshold = 0.01) 
	

write.csv(markers1, 
          file = "WM_markers_2nd_15yr_dh4.csv", 
          quote = FALSE, 
          row.names = FALSE)
		  
markers1 %>%
    group_by(cluster) %>%
    slice_max(n = 2, order_by = avg_log2FC)
top20 <- markers1 %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC)
 
pdf("doheatmap_2nd_15yr_DH4_WM.pdf")

DoHeatmap(object = WM, features = top20$gene, slot = "data", size = 2.5) + theme(axis.text.y = element_text(size = 2))

juvenile_4yr_DH1_seurat <- readRDS("/scratch/mshruv003/Updated_integrated_MTG_annotation/juvenile_4yr_DH1_seurat_leiden.rds")
juvenile_4yr_DH1_seurat <- SCTransform(juvenile_4yr_DH1_seurat, assay = "Spatial", verbose = FALSE)
Idents(juvenile_4yr_DH1_seurat) <- "scanpy_4yr_DH1.leiden"
cortex <- subset(juvenile_4yr_DH1_seurat, idents = c(1, 2, 3, 0, 5, 7, 8, 14, 13, 18, 4, 11))
cortex <- SCTransform(cortex, assay = "Spatial", verbose = FALSE)
markers <- FindAllMarkers(object = cortex, 
                          slot = "data",
                          only.pos = TRUE,
                          logfc.threshold = 0.01) 

write.csv(markers, 
          file = "cortex_markers_4yr_dh1.csv", 
          quote = FALSE, 
          row.names = FALSE)
		  
markers %>%
    group_by(cluster) %>%
    slice_max(n = 2, order_by = avg_log2FC)
top20 <- markers %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC)
 
pdf("doheatmap_4yr_DH1_cortex.pdf")

DoHeatmap(object = cortex, features = top20$gene, slot = "data", size = 2.5) + theme(axis.text.y = element_text(size = 2))

WM <- subset(juvenile_4yr_DH1_seurat, idents = c(6, 9, 10, 12, 12, 15, 16, 13))
WM <- SCTransform(WM, assay = "Spatial", verbose = FALSE)


markers1 <- FindAllMarkers(object = WM, 
                          slot = "data",
                          only.pos = TRUE,
                          logfc.threshold = 0.01) 
	
write.csv(markers1, 
          file = "WM_markers_4yr_dh1.csv", 
          quote = FALSE, 
          row.names = FALSE)
		  
markers1 %>%
    group_by(cluster) %>%
    slice_max(n = 2, order_by = avg_log2FC)
top20 <- markers1 %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC)
 
pdf("doheatmap_4yr_DH1_WM.pdf")

DoHeatmap(object = WM, features = top20$gene, slot = "data", size = 2.5) + theme(axis.text.y = element_text(size = 2))

juvenile_4yr_DH2_seurat <- readRDS("/scratch/mshruv003/Updated_integrated_MTG_annotation/juvenile_4yr_DH2_seurat_leiden.rds")
juvenile_4yr_DH2_seurat <- SCTransform(juvenile_4yr_DH2_seurat, assay = "Spatial", verbose = FALSE)
Idents(juvenile_4yr_DH2_seurat) <- "scanpy_4yr_DH2.leiden"
cortex <- subset(juvenile_4yr_DH2_seurat, idents = c(1, 2, 3, 0, 5, 4, 7, 8, 11, 14, 13, 19, 18))
cortex <- SCTransform(cortex, assay = "Spatial", verbose = FALSE)
markers <- FindAllMarkers(object = cortex, 
                          slot = "data",
                          only.pos = TRUE,
                          logfc.threshold = 0.01) 	

write.csv(markers, 
          file = "cortex_markers_4yr_dh2.csv", 
          quote = FALSE, 
          row.names = FALSE)
		  
markers %>%
    group_by(cluster) %>%
    slice_max(n = 2, order_by = avg_log2FC)
top20 <- markers %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC)
 
pdf("doheatmap_4yr_DH2_cortex.pdf")

DoHeatmap(object = cortex, features = top20$gene, slot = "data", size = 2.5) + theme(axis.text.y = element_text(size = 2))

WM <- subset(juvenile_4yr_DH2_seurat, idents = c(6, 9, 10, 12, 13, 15, 16))
WM <- SCTransform(WM, assay = "Spatial", verbose = FALSE)


markers1 <- FindAllMarkers(object = WM, 
                          slot = "data",
                          only.pos = TRUE,
                          logfc.threshold = 0.01) 
	
write.csv(markers1, 
          file = "WM_markers_4yr_dh2.csv", 
          quote = FALSE, 
          row.names = FALSE)
		  
markers1 %>%
    group_by(cluster) %>%
    slice_max(n = 2, order_by = avg_log2FC)
top20 <- markers1 %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC)
 
pdf("doheatmap_4yr_DH2_WM.pdf")

DoHeatmap(object = WM, features = top20$gene, slot = "data", size = 2.5) + theme(axis.text.y = element_text(size = 2))

