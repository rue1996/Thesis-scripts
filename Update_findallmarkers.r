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
	