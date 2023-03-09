import scanpy as sc
import anndata
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import pickle as pickle

import cell2location
import scvi

from matplotlib import rcParams
rcParams['pdf.fonttype'] = 42 # enables correct plotting of text
import seaborn as sns

run_name = '/scratch/mshruv003/Updated_integrated_MTG_annotation/cell2location_map'

sp_data_file = '/scratch/mshruv003/Updated_integrated_MTG_annotation/cell2location_map/sp_with_clusters.h5ad'
adata_vis = anndata.read(sp_data_file)

mod = cell2location.models.Cell2location.load(f"{run_name}", adata_vis)

clust_labels = ['means_per_cluster_mu_fg_Inh L1-2 SST BAGE2',
				  'means_per_cluster_mu_fg_Inh L1-4 LAMP5 LCP2', 'means_per_cluster_mu_fg_Inh L3-5 SST ADGRG6', 'means_per_cluster_mu_fg_Inh L4-6 SST B3GAT2', 'means_per_cluster_mu_fg_Inh L5-6 PVALB LGR5'
                        ]
clust_col = ['' + str(i) for i in clust_labels] # in case column names differ from labels

from cell2location.utils import select_slide
slide = select_slide(adata_vis, 'DH4a')

with mpl.rc_context({'figure.figsize': (20, 20)}):
    fig = plot_spatial(
        adata=slide,
        # labels to show on a plot
        color=clust_col, labels=clust_labels,
		coords=slide.obsm['spatial'] \
                          * list(slide.uns['spatial'].values())[0]['scalefactors']['tissue_hires_scalef'],
		alpha_scaling=1.0,
        show_img=True,
		img_alpha=0.0,
		img=list(slide.uns['spatial'].values())[0]['images']['hires'],
		adjust_text=True,
        # 'fast' (white background) or 'dark_background'
        style='dark_background',
        # limit color scale at 99.2% quantile of cell abundance
        max_color_quantile=0.992,
        # size of locations (adjust depending on figure size)
        circle_diameter=9,
        colorbar_position='right',
		reorder_cmap=[6, 1, 2, 3, 4, 0],
    ), plt.savefig("DH4a_31yr_Inh_clusters_thesis2.pdf")
	
clust_labels = ['means_per_cluster_mu_fg_Inh L1-4 LAMP5 LCP2', 'means_per_cluster_mu_fg_Inh L4-6 SST B3GAT2', 'means_per_cluster_mu_fg_Inh L5-6 PVALB LGR5'
                        ]
clust_col = ['' + str(i) for i in clust_labels] # in case column names differ from labels

from cell2location.utils import select_slide
slide = select_slide(adata_vis, 'DH4a')

with mpl.rc_context({'figure.figsize': (20, 20)}):
    fig = plot_spatial(
        adata=slide,
        # labels to show on a plot
        color=clust_col, labels=clust_labels,
		coords=slide.obsm['spatial'] \
                          * list(slide.uns['spatial'].values())[0]['scalefactors']['tissue_hires_scalef'],
		alpha_scaling=1.0,
        show_img=True,
		img_alpha=0.0,
		img=list(slide.uns['spatial'].values())[0]['images']['hires'],
		adjust_text=True,
        # 'fast' (white background) or 'dark_background'
        style='dark_background',
        # limit color scale at 99.2% quantile of cell abundance
        max_color_quantile=0.992,
        # size of locations (adjust depending on figure size)
        circle_diameter=9,
        colorbar_position='right',
		reorder_cmap=[2, 1, 0],
    ), plt.savefig("DH4a_31yr_Inh_clusters_thesis3.pdf")
	
	
clust_labels = ['means_per_cluster_mu_fg_Astro L1-2 FGFR3 GFAP', 'means_per_cluster_mu_fg_Exc L5-6 THEMIS C1QL3', 'means_per_cluster_mu_fg_Exc L3-5 RORB ESR1',
				  'means_per_cluster_mu_fg_Exc L2 LAMP5 LTK', 'means_per_cluster_mu_fg_Exc L2-3 LINC00507 FREM3'
                        ]
clust_col = ['' + str(i) for i in clust_labels] # in case column names differ from labels

from cell2location.utils import select_slide
slide = select_slide(adata_vis, 'DH4a')

with mpl.rc_context({'figure.figsize': (20, 20)}):
    fig = plot_spatial(
        adata=slide,
        # labels to show on a plot
        color=clust_col, labels=clust_labels,
		coords=slide.obsm['spatial'] \
                          * list(slide.uns['spatial'].values())[0]['scalefactors']['tissue_hires_scalef'],
		alpha_scaling=1.0,
        show_img=True,
		img_alpha=0.0,
		img=list(slide.uns['spatial'].values())[0]['images']['hires'],
		adjust_text=True,
        # 'fast' (white background) or 'dark_background'
        style='dark_background',
        # limit color scale at 99.2% quantile of cell abundance
        max_color_quantile=0.992,
        # size of locations (adjust depending on figure size)
        circle_diameter=9,
        colorbar_position='right',
		reorder_cmap=[6, 1, 2, 0, 4],
    ), plt.savefig("DH4a_31yr_Exc_clusters_thesis3a.pdf")