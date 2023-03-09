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
library(dittoSeq)

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

WM <- subset(brain.merge3, idents = c(6, 9, 10, 12, 15, 16))
WM <- SCTransform(WM, assay = "Spatial", verbose = FALSE)

markers <- FindAllMarkers(object = WM, 
                          slot = "data",
                          only.pos = TRUE,
                          logfc.threshold = 0.25) 	

write.csv(markers, 
          file = "WM.csv", 
          quote = FALSE, 
          row.names = FALSE)	
		  
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
pdf("WM_dittoheatmap4.pdf")


Ave_exp@meta.data$sample <- NA
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^6_4_year_old_T1_B2"))] <- "4_year_old_T1_B2"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^6_4_year_old_T2_B2"))] <- "4_year_old_T2_B2"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^6_15_year_old_T1_B1"))] <- "15_year_old_T1_B1"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^6_15_year_old_T2_B1"))] <- "15_year_old_T2_B1"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^6_15_year_old_T1_B2"))] <- "15_year_old_T1_B2"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^6_15_year_old_T2_B2"))] <- "15_year_old_T2_B2"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^6_31_year_old_T1_B1"))] <- "31_year_old_T1_B1"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^6_31_year_old_T2_B1"))] <- "31_year_old_T2_B1"

Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^9_4_year_old_T1_B2"))] <- "4_year_old_T1_B2"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^9_4_year_old_T2_B2"))] <- "4_year_old_T2_B2"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^9_15_year_old_T1_B1"))] <- "15_year_old_T1_B1"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^9_15_year_old_T2_B1"))] <- "15_year_old_T2_B1"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^9_15_year_old_T1_B2"))] <- "15_year_old_T1_B2"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^9_15_year_old_T2_B2"))] <- "15_year_old_T2_B2"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^9_31_year_old_T1_B1"))] <- "31_year_old_T1_B1"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^9_31_year_old_T2_B1"))] <- "31_year_old_T2_B1"


Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^10_4_year_old_T1_B2"))] <- "4_year_old_T1_B2"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^10_4_year_old_T2_B2"))] <- "4_year_old_T2_B2"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^10_15_year_old_T1_B1"))] <- "15_year_old_T1_B1"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^10_15_year_old_T2_B1"))] <- "15_year_old_T2_B1"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^10_15_year_old_T1_B2"))] <- "15_year_old_T1_B2"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^10_15_year_old_T2_B2"))] <- "15_year_old_T2_B2"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^10_31_year_old_T1_B1"))] <- "31_year_old_T1_B1"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^10_31_year_old_T2_B1"))] <- "31_year_old_T2_B1"

Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^12_4_year_old_T1_B2"))] <- "4_year_old_T1_B2"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^12_4_year_old_T2_B2"))] <- "4_year_old_T2_B2"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^12_15_year_old_T1_B1"))] <- "15_year_old_T1_B1"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^12_15_year_old_T2_B1"))] <- "15_year_old_T2_B1"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^12_15_year_old_T1_B2"))] <- "15_year_old_T1_B2"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^12_15_year_old_T2_B2"))] <- "15_year_old_T2_B2"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^12_31_year_old_T1_B1"))] <- "31_year_old_T1_B1"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^12_31_year_old_T2_B1"))] <- "31_year_old_T2_B1"

Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^15_4_year_old_T1_B2"))] <- "4_year_old_T1_B2"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^15_4_year_old_T2_B2"))] <- "4_year_old_T2_B2"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^15_15_year_old_T1_B1"))] <- "15_year_old_T1_B1"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^15_15_year_old_T2_B1"))] <- "15_year_old_T2_B1"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^15_15_year_old_T1_B2"))] <- "15_year_old_T1_B2"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^15_15_year_old_T2_B2"))] <- "15_year_old_T2_B2"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^15_31_year_old_T1_B1"))] <- "31_year_old_T1_B1"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^15_31_year_old_T2_B1"))] <- "31_year_old_T2_B1"

Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^16_4_year_old_T1_B2"))] <- "4_year_old_T1_B2"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^16_4_year_old_T2_B2"))] <- "4_year_old_T2_B2"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^16_15_year_old_T1_B1"))] <- "15_year_old_T1_B1"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^16_15_year_old_T2_B1"))] <- "15_year_old_T2_B1"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^16_15_year_old_T1_B2"))] <- "15_year_old_T1_B2"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^16_15_year_old_T2_B2"))] <- "15_year_old_T2_B2"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^16_31_year_old_T1_B1"))] <- "31_year_old_T1_B1"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^16_31_year_old_T2_B1"))] <- "31_year_old_T2_B1"

Ave_exp$sample <- factor(Ave_exp$sample, levels = c("4_year_old_T1_B2", "4_year_old_T2_B2", "15_year_old_T1_B2", "15_year_old_T2_B2", "15_year_old_T1_B1", "15_year_old_T2_B1", "31_year_old_T1_B1", "31_year_old_T2_B1"))

Ave_exp$Region_cluster <- Ave_exp$orig.ident

dittoHeatmap(Ave_exp, top20$gene,
    annot.by = c("sample", "Region_cluster"),
	order.by = c("Region_cluster", "sample"), slot = "data", scale = "row", 
	fontsize_row = 4.0)
	
dev.off()

cortex <- subset(brain.merge3, idents = c(7, 14))
cortex <- SCTransform(cortex, assay = "Spatial", verbose = FALSE)

markers <- FindAllMarkers(object = cortex, 
                          slot = "data",
                          only.pos = TRUE,
                          logfc.threshold = 0.25) 

write.csv(markers, 
          file = "cortex_7_14.csv", 
          quote = FALSE, 
          row.names = FALSE)						  

markers %>%
    slice_max(n = 20, order_by = avg_log2FC)
top10 <- markers %>%
  top_n(n = 20, wt = avg_log2FC)
 						  
						 

Ave_exp <- AverageExpression(
                  object = cortex,
                  return.seurat = TRUE,
                  group.by = c("leiden_clusters.leiden", "sample"),
                  slot = "data")		
						  

pdf("cortex714.dittoHeatmap3.pdf")

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

Ave_exp$sample <- factor(Ave_exp$sample, levels = c("4_year_old_T1_B2", "4_year_old_T2_B2", "15_year_old_T1_B2", "15_year_old_T2_B2", "15_year_old_T1_B1", "15_year_old_T2_B1", "31_year_old_T1_B1", "31_year_old_T2_B1"))

Ave_exp$Region_cluster <- Ave_exp$orig.ident


dittoHeatmap(Ave_exp, top10$gene,
    annot.by = c("sample", "Region_cluster"),
	order.by = c("Region_cluster", "sample"), slot = "data", scale = "row", 
	fontsize_row = 5.0)
	
cortex <- subset(brain.merge3, idents = c(17, 18, 20))
cortex <- SCTransform(cortex, assay = "Spatial", verbose = FALSE)

markers <- FindAllMarkers(object = cortex, 
                          slot = "data",
                          only.pos = TRUE,
                          logfc.threshold = 0.25) 	
						  
write.csv(markers, 
          file = "cortex_17_18_20.csv", 
          quote = FALSE, 
          row.names = FALSE)	

markers %>%
    slice_max(n = 20, order_by = avg_log2FC)
top10 <- markers %>%
  top_n(n = 20, wt = avg_log2FC)
 						  
						 

Ave_exp <- AverageExpression(
                  object = cortex,
                  return.seurat = TRUE,
                  group.by = c("leiden_clusters.leiden", "sample"),
                  slot = "data")		
						  

pdf("17_18_20cortex.dittoHeatmap3.pdf")

Ave_exp@meta.data$sample <- NA
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^17_15_year_old_T1_B1"))] <- "15_year_old_T1_B1"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^17_15_year_old_T2_B1"))] <- "15_year_old_T2_B1"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^17_15_year_old_T1_B2"))] <- "15_year_old_T1_B2"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^17_15_year_old_T2_B2"))] <- "15_year_old_T2_B2"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^17_31_year_old_T1_B1"))] <- "31_year_old_T1_B1"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^17_31_year_old_T2_B1"))] <- "31_year_old_T2_B1"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^17_4_year_old_T1_B2"))] <- "4_year_old_T1_B2"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^17_4_year_old_T2_B2"))] <- "4_year_old_T2_B2"


Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^20_15_year_old_T1_B1"))] <- "15_year_old_T1_B1"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^20_15_year_old_T2_B1"))] <- "15_year_old_T2_B1"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^20_15_year_old_T1_B2"))] <- "15_year_old_T1_B2"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^20_15_year_old_T2_B2"))] <- "15_year_old_T2_B2"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^20_31_year_old_T1_B1"))] <- "31_year_old_T1_B1"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^20_31_year_old_T2_B1"))] <- "31_year_old_T2_B1"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^20_4_year_old_T1_B2"))] <- "4_year_old_T1_B2"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^20_4_year_old_T2_B2"))] <- "4_year_old_T2_B2"

Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^18_15_year_old_T1_B1"))] <- "15_year_old_T1_B1"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^18_15_year_old_T2_B1"))] <- "15_year_old_T2_B1"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^18_15_year_old_T1_B2"))] <- "15_year_old_T1_B2"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^18_15_year_old_T2_B2"))] <- "15_year_old_T2_B2"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^18_31_year_old_T1_B1"))] <- "31_year_old_T1_B1"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^18_31_year_old_T2_B1"))] <- "31_year_old_T2_B1"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^18_4_year_old_T1_B2"))] <- "4_year_old_T1_B2"
Ave_exp@meta.data$sample[which(str_detect(rownames(x = Ave_exp@meta.data),"^18_4_year_old_T2_B2"))] <- "4_year_old_T2_B2"

Ave_exp$sample <- factor(Ave_exp$sample, levels = c("4_year_old_T1_B2", "4_year_old_T2_B2", "15_year_old_T1_B2", "15_year_old_T2_B2", "15_year_old_T1_B1", "15_year_old_T2_B1", "31_year_old_T1_B1", "31_year_old_T2_B1"))

Ave_exp$Region_cluster <- Ave_exp$orig.ident


dittoHeatmap(Ave_exp, top10$gene,
    annot.by = c("sample", "Region_cluster"),
	order.by = c("Region_cluster", "sample"), slot = "data", scale = "row", 
	fontsize_row = 4.0)
	