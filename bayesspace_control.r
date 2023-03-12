#Bayesspace figures for the genes used as controls
#Section: 2.5.5 Enhancement of spatial transcriptomic data visualisation (BayesSpace)
#Code adapted from: https://edward130603.github.io/BayesSpace/articles/BayesSpace.html

library(BayesSpace)
library(ggplot2)
library(SingleCellExperiment)


Juvenile15yr_DH3_2 <- readVisium("/scratch/mshruv003/visium_18_01_22/X201SC21111924-Z01-F001_3/raw_data/run_spaceranger_count_DH3/outs/")
set.seed(102)
Juvenile15yr_DH3_2  <- spatialPreprocess(Juvenile15yr_DH3_2, platform="Visium", 
                              n.PCs=7, n.HVGs=2000, skip.PCA = FALSE, log.normalize=TRUE)
							  
Juvenile15yr_DH3_2  <- qTune(Juvenile15yr_DH3_2, qs=seq(2, 10), platform="Visium", d=7)
qPlot(Juvenile15yr_DH3_2 )
set.seed(149)

pdf("Juvenile15yr_DH3_2_nrep10000_d7_q_8_control.pdf")
Juvenile15yr_DH3_2  <- spatialCluster(Juvenile15yr_DH3_2, q=8, platform="Visium", d=7,
                           init.method="mclust", model="t", gamma=3,
                           nrep=10000, burn.in=1000,
                           save.chain=TRUE)

clusterPlot(Juvenile15yr_DH3_2 )

featurePlot(Juvenile15yr_DH3_2, "SLC17A7") 
featurePlot(Juvenile15yr_DH3_2, "AQP4")
featurePlot(Juvenile15yr_DH3_2, "GAD1")

Juvenile15yr_DH3_2.enhanced <- spatialEnhance(Juvenile15yr_DH3_2, q=8, platform="Visium", d=7,
                                    model="t", gamma=3,
                                    jitter_prior=0.3, jitter_scale=3.5,
                                    nrep=10000, burn.in=1000,
                                    save.chain=TRUE)  

clusterPlot(Juvenile15yr_DH3_2.enhanced)									
markers <- c('SLC17A7', 'AQP4', 'GAD1')
Juvenile15yr_DH3_2.enhanced <- enhanceFeatures(Juvenile15yr_DH3_2.enhanced, Juvenile15yr_DH3_2,
                                     feature_names=markers,
                                     nrounds=0)
									 
logcounts(Juvenile15yr_DH3_2.enhanced)[markers, 1:5]	
rowData(Juvenile15yr_DH3_2.enhanced)[markers, ]

featurePlot(Juvenile15yr_DH3_2.enhanced, "SLC17A7") + scale_fill_gradientn(colors = c(low = 'white', high = 'brown'), limits = c(0.0,4.0), breaks = c(0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0))
featurePlot(Juvenile15yr_DH3_2.enhanced, "AQP4") + scale_fill_gradientn(colors = c(low = 'white', high = 'brown'), limits = c(0.0,4.5), breaks = c(0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 4.5))
featurePlot(Juvenile15yr_DH3_2.enhanced, "GAD1") + scale_fill_gradientn(colors = c(low = 'white', high = 'brown'), limits = c(0.0,3.0), breaks = c(0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0))

								 
enhanced.plots <- purrr::map(markers, function(x) featurePlot(Juvenile15yr_DH3_2.enhanced, x))
patchwork::wrap_plots(enhanced.plots, ncol=2)	

spot.plots <- purrr::map(markers, function(x) featurePlot(Juvenile15yr_DH3_2, x))
patchwork::wrap_plots(c(enhanced.plots, spot.plots), ncol=2)	

dev.off()

Juvenile4yr_DH1 <- readVisium("/scratch/mshruv003/visium_18_01_22/X201SC21111924-Z01-F001_1/raw_data/run_spaceranger_count_DH1/outs/")
set.seed(102)
Juvenile4yr_DH1 <- spatialPreprocess(Juvenile4yr_DH1, platform="Visium", 
                              n.PCs=7, n.HVGs=2000, skip.PCA = FALSE, log.normalize=TRUE)
							  
Juvenile4yr_DH1 <- qTune(Juvenile4yr_DH1, qs=seq(2, 10), platform="Visium", d=7)
qPlot(Juvenile4yr_DH1)
set.seed(149)

pdf("Juvenile4yr_DH1_nrep10000_d7_q8_control.pdf")
Juvenile4yr_DH1 <- spatialCluster(Juvenile4yr_DH1, q=8, platform="Visium", d=7,
                           init.method="mclust", model="t", gamma=3,
                           nrep=10000, burn.in=1000,
                           save.chain=TRUE)

clusterPlot(Juvenile4yr_DH1)

featurePlot(Juvenile4yr_DH1, "SLC17A7")
featurePlot(Juvenile4yr_DH1, "AQP4")
featurePlot(Juvenile4yr_DH1, "GAD1")

Juvenile4yr_DH1.enhanced <- spatialEnhance(Juvenile4yr_DH1, q=8, platform="Visium", d=7,
                                    model="t", gamma=3,
                                    jitter_prior=0.3, jitter_scale=3.5,
                                    nrep=10000, burn.in=1000,
                                    save.chain=TRUE)  

clusterPlot(Juvenile4yr_DH1.enhanced)									
markers <- c('SLC17A7', 'AQP4', 'GAD1')
Juvenile4yr_DH1.enhanced <- enhanceFeatures(Juvenile4yr_DH1.enhanced, Juvenile4yr_DH1,
                                     feature_names=markers,
                                     nrounds=0)
									 
logcounts(Juvenile4yr_DH1.enhanced)[markers, 1:5]	
rowData(Juvenile4yr_DH1.enhanced)[markers, ]


featurePlot(Juvenile4yr_DH1.enhanced, "SLC17A7") + scale_fill_gradientn(colors = c(low = 'white', high = 'brown'), limits = c(0.0,4.0), breaks = c(0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0))
featurePlot(Juvenile4yr_DH1.enhanced, "AQP4") + scale_fill_gradientn(colors = c(low = 'white', high = 'brown'), limits = c(0.0,4.5), breaks = c(0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 4.5))
featurePlot(Juvenile4yr_DH1.enhanced, "GAD1") + scale_fill_gradientn(colors = c(low = 'white', high = 'brown'), limits = c(0.0,3.0), breaks = c(0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0))


enhanced.plots <- purrr::map(markers, function(x) featurePlot(Juvenile4yr_DH1.enhanced, x))
patchwork::wrap_plots(enhanced.plots, ncol=2)	

spot.plots <- purrr::map(markers, function(x) featurePlot(Juvenile4yr_DH1, x))
patchwork::wrap_plots(c(enhanced.plots, spot.plots), ncol=2)	

dev.off()

Juvenile4yr_DH2 <- readVisium("/scratch/mshruv003/visium_18_01_22/X201SC21111924-Z01-F001_2/raw_data/run_spaceranger_count_DH2/outs/")
set.seed(102)
Juvenile4yr_DH2 <- spatialPreprocess(Juvenile4yr_DH2, platform="Visium", 
                              n.PCs=7, n.HVGs=2000, skip.PCA = FALSE, log.normalize=TRUE)
							  
Juvenile4yr_DH2 <- qTune(Juvenile4yr_DH2, qs=seq(2, 10), platform="Visium", d=7)
qPlot(Juvenile4yr_DH2)
set.seed(149)

pdf("Juvenile4yr_DH2_nrep10000_d7_q8_control.pdf")
Juvenile4yr_DH2 <- spatialCluster(Juvenile4yr_DH2, q=8, platform="Visium", d=7,
                           init.method="mclust", model="t", gamma=3,
                           nrep=10000, burn.in=1000,
                           save.chain=TRUE)

clusterPlot(Juvenile4yr_DH2)

featurePlot(Juvenile4yr_DH2, "SLC17A7")
featurePlot(Juvenile4yr_DH2, "AQP4")
featurePlot(Juvenile4yr_DH2, "GAD1")

Juvenile4yr_DH2.enhanced <- spatialEnhance(Juvenile4yr_DH2, q=8, platform="Visium", d=7,
                                    model="t", gamma=3,
                                    jitter_prior=0.3, jitter_scale=3.5,
                                    nrep=10000, burn.in=1000,
                                    save.chain=TRUE)  

clusterPlot(Juvenile4yr_DH2.enhanced)									
markers <- c('SLC17A7', 'AQP4', 'GAD1')
Juvenile4yr_DH2.enhanced <- enhanceFeatures(Juvenile4yr_DH2.enhanced, Juvenile4yr_DH2,
                                     feature_names=markers,
                                     nrounds=0)
									 
logcounts(Juvenile4yr_DH2.enhanced)[markers, 1:5]	
rowData(Juvenile4yr_DH2.enhanced)[markers, ]

featurePlot(Juvenile4yr_DH2.enhanced, "SLC17A7") + scale_fill_gradientn(colors = c(low = 'white', high = 'brown'), limits = c(0.0,4.0), breaks = c(0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0))
featurePlot(Juvenile4yr_DH2.enhanced, "AQP4") + scale_fill_gradientn(colors = c(low = 'white', high = 'brown'), limits = c(0.0,4.5), breaks = c(0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 4.5))
featurePlot(Juvenile4yr_DH2.enhanced, "GAD1") + scale_fill_gradientn(colors = c(low = 'white', high = 'brown'), limits = c(0.0,3.0), breaks = c(0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0))
								 
enhanced.plots <- purrr::map(markers, function(x) featurePlot(Juvenile4yr_DH2.enhanced, x))
patchwork::wrap_plots(enhanced.plots, ncol=2)	

spot.plots <- purrr::map(markers, function(x) featurePlot(Juvenile4yr_DH2, x))
patchwork::wrap_plots(c(enhanced.plots, spot.plots), ncol=2)	

dev.off()

Juvenile15yr_DH4_2 <- readVisium("/scratch/mshruv003/visium_18_01_22/X201SC21111924-Z01-F001_2/raw_data/run_spaceranger_count_DH4/outs/")
set.seed(102)
Juvenile15yr_DH4_2 <- spatialPreprocess(Juvenile15yr_DH4_2, platform="Visium", 
                              n.PCs=7, n.HVGs=2000, skip.PCA = FALSE, log.normalize=TRUE)
							  
Juvenile15yr_DH4_2 <- qTune(Juvenile15yr_DH4_2, qs=seq(2, 10), platform="Visium", d=7)
qPlot(Juvenile15yr_DH4_2)
set.seed(149)

pdf("Juvenile15yr_DH4_2_nrep10000_d7_q8_control.pdf")
Juvenile15yr_DH4_2 <- spatialCluster(Juvenile15yr_DH4_2, q=8, platform="Visium", d=7,
                           init.method="mclust", model="t", gamma=3,
                           nrep=10000, burn.in=1000,
                           save.chain=TRUE)

clusterPlot(Juvenile15yr_DH4_2)

featurePlot(Juvenile15yr_DH4_2, "SLC17A7")
featurePlot(Juvenile15yr_DH4_2, "AQP4")
featurePlot(Juvenile15yr_DH4_2, "GAD1")

Juvenile15yr_DH4_2.enhanced <- spatialEnhance(Juvenile15yr_DH4_2, q=8, platform="Visium", d=7,
                                    model="t", gamma=3,
                                    jitter_prior=0.3, jitter_scale=3.5,
                                    nrep=10000, burn.in=1000,
                                    save.chain=TRUE)  

clusterPlot(Juvenile15yr_DH4_2.enhanced)
markers <- c('SLC17A7', 'AQP4', 'GAD1')
Juvenile15yr_DH4_2.enhanced <- enhanceFeatures(Juvenile15yr_DH4_2.enhanced, Juvenile15yr_DH4_2,
                                     feature_names=markers,
                                     nrounds=0)
									 
logcounts(Juvenile15yr_DH4_2.enhanced)[markers, 1:5]	
rowData(Juvenile15yr_DH4_2.enhanced)[markers, ]

featurePlot(Juvenile15yr_DH4_2.enhanced, "SLC17A7") + scale_fill_gradientn(colors = c(low = 'white', high = 'brown'), limits = c(0.0,4.0), breaks = c(0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0))
featurePlot(Juvenile15yr_DH4_2.enhanced, "AQP4") + scale_fill_gradientn(colors = c(low = 'white', high = 'brown'), limits = c(0.0,4.5), breaks = c(0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 4.5))
featurePlot(Juvenile15yr_DH4_2.enhanced, "GAD1") + scale_fill_gradientn(colors = c(low = 'white', high = 'brown'), limits = c(0.0,3.0), breaks = c(0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0))

								 
enhanced.plots <- purrr::map(markers, function(x) featurePlot(Juvenile15yr_DH4_2.enhanced, x))
patchwork::wrap_plots(enhanced.plots, ncol=2)	

spot.plots <- purrr::map(markers, function(x) featurePlot(Juvenile15yr_DH4_2, x))
patchwork::wrap_plots(c(enhanced.plots, spot.plots), ncol=2)	

dev.off()

Juvenile15yr_DH1_1 <- readVisium("/scratch/mshruv003/visium_08_10_21/Hockman_Visium/X201SC21082979-Z01-F001_01/raw_data/run_spaceranger_count_DH1a/outs")
set.seed(102)
Juvenile15yr_DH1_1  <- spatialPreprocess(Juvenile15yr_DH1_1, platform="Visium", 
                              n.PCs=7, n.HVGs=2000, skip.PCA = FALSE, log.normalize=TRUE)
							  
Juvenile15yr_DH1_1  <- qTune(Juvenile15yr_DH1_1, qs=seq(2, 10), platform="Visium", d=7)
qPlot(Juvenile15yr_DH1_1)
set.seed(149)

pdf("Juvenile15yr_DH1_1_nrep10000_d7_q_8_control.pdf")
Juvenile15yr_DH1_1  <- spatialCluster(Juvenile15yr_DH1_1, q=8, platform="Visium", d=7,
                           init.method="mclust", model="t", gamma=3,
                           nrep=10000, burn.in=1000,
                           save.chain=TRUE)

clusterPlot(Juvenile15yr_DH1_1)


featurePlot(Juvenile15yr_DH1_1, "SLC17A7") 
featurePlot(Juvenile15yr_DH1_1, "AQP4")
featurePlot(Juvenile15yr_DH1_1, "GAD1")


Juvenile15yr_DH1_1.enhanced <- spatialEnhance(Juvenile15yr_DH1_1, q=8, platform="Visium", d=7,
                                    model="t", gamma=3,
                                    jitter_prior=0.3, jitter_scale=3.5,
                                    nrep=10000, burn.in=1000,
                                    save.chain=TRUE)  

clusterPlot(Juvenile15yr_DH1_1.enhanced)									
markers <- c('SLC17A7', 'AQP4', 'GAD1')
Juvenile15yr_DH1_1.enhanced <- enhanceFeatures(Juvenile15yr_DH1_1.enhanced, Juvenile15yr_DH1_1,
                                     feature_names=markers,
                                     nrounds=0)
									 
logcounts(Juvenile15yr_DH1_1.enhanced)[markers, 1:5]	
rowData(Juvenile15yr_DH1_1.enhanced)[markers, ]

featurePlot(Juvenile15yr_DH1_1.enhanced, "SLC17A7") + scale_fill_gradientn(colors = c(low = 'white', high = 'brown'), limits = c(0.0,4.0), breaks = c(0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0))
featurePlot(Juvenile15yr_DH1_1.enhanced, "AQP4") + scale_fill_gradientn(colors = c(low = 'white', high = 'brown'), limits = c(0.0,4.5), breaks = c(0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 4.5))
featurePlot(Juvenile15yr_DH1_1.enhanced, "GAD1") + scale_fill_gradientn(colors = c(low = 'white', high = 'brown'), limits = c(0.0,3.0), breaks = c(0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0))

								 
enhanced.plots <- purrr::map(markers, function(x) featurePlot(Juvenile15yr_DH1_1.enhanced, x))
patchwork::wrap_plots(enhanced.plots, ncol=2)	

spot.plots <- purrr::map(markers, function(x) featurePlot(Juvenile15yr_DH1_1, x))
patchwork::wrap_plots(c(enhanced.plots, spot.plots), ncol=2)	

dev.off()
						
Juvenile15yr_DH2_1 <- readVisium("/scratch/mshruv003/visium_08_10_21/Hockman_Visium/X201SC21082979-Z01-F001_01/raw_data/run_spaceranger_count_DH2a/outs")
set.seed(102)
Juvenile15yr_DH2_1  <- spatialPreprocess(Juvenile15yr_DH2_1, platform="Visium", 
                              n.PCs=7, n.HVGs=2000, skip.PCA = FALSE, log.normalize=TRUE)
							  
Juvenile15yr_DH2_1  <- qTune(Juvenile15yr_DH2_1, qs=seq(2, 10), platform="Visium", d=7)
qPlot(Juvenile15yr_DH2_1)
set.seed(149)

pdf("Juvenile15yr_DH2_1_nrep10000_d7_q_8_control.pdf")
Juvenile15yr_DH2_1  <- spatialCluster(Juvenile15yr_DH2_1, q=8, platform="Visium", d=7,
                           init.method="mclust", model="t", gamma=3,
                           nrep=10000, burn.in=1000,
                           save.chain=TRUE)

clusterPlot(Juvenile15yr_DH2_1)

featurePlot(Juvenile15yr_DH2_1, "SLC17A7")
featurePlot(Juvenile15yr_DH2_1, "AQP4")
featurePlot(Juvenile15yr_DH2_1, "GAD1")

Juvenile15yr_DH2_1.enhanced <- spatialEnhance(Juvenile15yr_DH2_1, q=8, platform="Visium", d=7,
                                    model="t", gamma=3,
                                    jitter_prior=0.3, jitter_scale=3.5,
                                    nrep=10000, burn.in=1000,
                                    save.chain=TRUE)  

clusterPlot(Juvenile15yr_DH2_1.enhanced)
markers <- c('SLC17A7', 'AQP4', 'GAD1')
Juvenile15yr_DH2_1.enhanced <- enhanceFeatures(Juvenile15yr_DH2_1.enhanced, Juvenile15yr_DH2_1,
                                     feature_names=markers,
                                     nrounds=0)
									 
logcounts(Juvenile15yr_DH2_1.enhanced)[markers, 1:5]	
rowData(Juvenile15yr_DH2_1.enhanced)[markers, ]

featurePlot(Juvenile15yr_DH2_1.enhanced, "SLC17A7") + scale_fill_gradientn(colors = c(low = 'white', high = 'brown'), limits = c(0.0,4.0), breaks = c(0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0))
featurePlot(Juvenile15yr_DH2_1.enhanced, "AQP4") + scale_fill_gradientn(colors = c(low = 'white', high = 'brown'), limits = c(0.0,4.5), breaks = c(0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 4.5))
featurePlot(Juvenile15yr_DH2_1.enhanced, "GAD1") + scale_fill_gradientn(colors = c(low = 'white', high = 'brown'), limits = c(0.0,3.0), breaks = c(0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0))
								 
enhanced.plots <- purrr::map(markers, function(x) featurePlot(Juvenile15yr_DH2_1.enhanced, x))
patchwork::wrap_plots(enhanced.plots, ncol=2)	

spot.plots <- purrr::map(markers, function(x) featurePlot(Juvenile15yr_DH2_1, x))
patchwork::wrap_plots(c(enhanced.plots, spot.plots), ncol=2)	

dev.off()

Adult31yr_DH3 <- readVisium("/scratch/mshruv003/visium_08_10_21/Hockman_Visium/X201SC21082979-Z01-F001_02/raw_data/run_spaceranger_count_DH3a/outs")
set.seed(102)
Adult31yr_DH3  <- spatialPreprocess(Adult31yr_DH3, platform="Visium", 
                              n.PCs=7, n.HVGs=2000, skip.PCA = FALSE, log.normalize=TRUE)
							  
Adult31yr_DH3  <- qTune(Adult31yr_DH3, qs=seq(2, 10), platform="Visium", d=7)
qPlot(Adult31yr_DH3)
set.seed(149)

pdf("Adult31yr_DH3_nrep10000_d7_q_8_control.pdf")
Adult31yr_DH3  <- spatialCluster(Adult31yr_DH3, q=8, platform="Visium", d=7,
                           init.method="mclust", model="t", gamma=3,
                           nrep=10000, burn.in=1000,
                           save.chain=TRUE)

clusterPlot(Adult31yr_DH3)

featurePlot(Adult31yr_DH3, "SLC17A7")
featurePlot(Adult31yr_DH3, "AQP4")
featurePlot(Adult31yr_DH3, "GAD1")

Adult31yr_DH3.enhanced <- spatialEnhance(Adult31yr_DH3, q=8, platform="Visium", d=7,
                                    model="t", gamma=3,
                                    jitter_prior=0.3, jitter_scale=3.5,
                                    nrep=10000, burn.in=1000,
                                    save.chain=TRUE)  

clusterPlot(Adult31yr_DH3.enhanced)									
markers <- c('SLC17A7', 'AQP4', 'GAD1')
Adult31yr_DH3.enhanced <- enhanceFeatures(Adult31yr_DH3.enhanced, Adult31yr_DH3,
                                     feature_names=markers,
                                     nrounds=0)
									 
logcounts(Adult31yr_DH3.enhanced)[markers, 1:5]	
rowData(Adult31yr_DH3.enhanced)[markers, ]

featurePlot(Adult31yr_DH3.enhanced, "SLC17A7") + scale_fill_gradientn(colors = c(low = 'white', high = 'brown'), limits = c(0.0,4.0), breaks = c(0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0))
featurePlot(Adult31yr_DH3.enhanced, "AQP4") + scale_fill_gradientn(colors = c(low = 'white', high = 'brown'), limits = c(0.0,4.5), breaks = c(0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 4.5))
featurePlot(Adult31yr_DH3.enhanced, "GAD1") + scale_fill_gradientn(colors = c(low = 'white', high = 'brown'), limits = c(0.0,3.0), breaks = c(0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0))

								 
enhanced.plots <- purrr::map(markers, function(x) featurePlot(Adult31yr_DH3.enhanced, x))
patchwork::wrap_plots(enhanced.plots, ncol=2)	

spot.plots <- purrr::map(markers, function(x) featurePlot(Adult31yr_DH3, x))
patchwork::wrap_plots(c(enhanced.plots, spot.plots), ncol=2)	

dev.off()

Adult31yr_DH4 <- readVisium("/scratch/mshruv003/visium_08_10_21/Hockman_Visium/X201SC21082979-Z01-F001_02/raw_data/run_spaceranger_count_DH4a/outs")
set.seed(102)
Adult31yr_DH4  <- spatialPreprocess(Adult31yr_DH4, platform="Visium", 
                              n.PCs=7, n.HVGs=2000, skip.PCA = FALSE, log.normalize=TRUE)
							  
Adult31yr_DH4  <- qTune(Adult31yr_DH4, qs=seq(2, 10), platform="Visium", d=7)
qPlot(Adult31yr_DH4)
set.seed(149)

pdf("Adult31yr_DH4_nrep10000_d7_q_8_control.pdf")
Adult31yr_DH4  <- spatialCluster(Adult31yr_DH4, q=8, platform="Visium", d=7,
                           init.method="mclust", model="t", gamma=3,
                           nrep=10000, burn.in=1000,
                           save.chain=TRUE)

clusterPlot(Adult31yr_DH4)

featurePlot(Adult31yr_DH4, "SLC17A7")
featurePlot(Adult31yr_DH4, "AQP4")
featurePlot(Adult31yr_DH4, "GAD1")


Adult31yr_DH4.enhanced <- spatialEnhance(Adult31yr_DH4, q=8, platform="Visium", d=7,
                                    model="t", gamma=3,
                                    jitter_prior=0.3, jitter_scale=3.5,
                                    nrep=10000, burn.in=1000,
                                    save.chain=TRUE)  

clusterPlot(Adult31yr_DH4.enhanced)									
markers <- c('SLC17A7', 'AQP4', 'GAD1')
Adult31yr_DH4.enhanced <- enhanceFeatures(Adult31yr_DH4.enhanced, Adult31yr_DH4,
                                     feature_names=markers,
                                     nrounds=0)
									 
logcounts(Adult31yr_DH4.enhanced)[markers, 1:5]	
rowData(Adult31yr_DH4.enhanced)[markers, ]

featurePlot(Adult31yr_DH4.enhanced, "SLC17A7") + scale_fill_gradientn(colors = c(low = 'white', high = 'brown'), limits = c(0.0,4.0), breaks = c(0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0))
featurePlot(Adult31yr_DH4.enhanced, "AQP4") + scale_fill_gradientn(colors = c(low = 'white', high = 'brown'), limits = c(0.0,4.5), breaks = c(0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 4.5))
featurePlot(Adult31yr_DH4.enhanced, "GAD1") + scale_fill_gradientn(colors = c(low = 'white', high = 'brown'), limits = c(0.0,3.0), breaks = c(0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0))

enhanced.plots <- purrr::map(markers, function(x) featurePlot(Adult31yr_DH4.enhanced, x))
patchwork::wrap_plots(enhanced.plots, ncol=2)	

spot.plots <- purrr::map(markers, function(x) featurePlot(Adult31yr_DH4, x))
patchwork::wrap_plots(c(enhanced.plots, spot.plots), ncol=2)	

dev.off()
