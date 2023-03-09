import sys
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

# read a previously trained LocationModel
results_folder = '/scratch/mshruv003/Updated_integrated_MTG_annotation/'
sp_results_folder = f'{results_folder}'
sc_results_folder = f'{results_folder}reference_signature_cluster_label/'

run_name = '/scratch/mshruv003/Updated_integrated_MTG_annotation/cell2location_map'
sp_data_file = '/scratch/mshruv003/Updated_integrated_MTG_annotation/cell2location_map/sp_with_clusters.h5ad'
adata_vis = anndata.read(sp_data_file)


def unpickle_model(path, mod_name):
    r""" Unpickle model
    """
    file = path + 'model_' + mod_name + ".p"
    
    mod1_ann = pickle.load(file = open(file, "rb"))
    return mod1_ann['mod']

n_fact = 15
mod_path = '/scratch/mshruv003/Updated_integrated_MTG_annotation/cell2location_map/CoLocatedComb/CoLocatedGroupsSklearnNMF_27703locations_57factors/models/'
adata_file = '/scratch/mshruv003/Updated_integrated_MTG_annotation/cell2location_map/CoLocatedComb/CoLocatedGroupsSklearnNMF_27703locations_57factors/anndata/sp.h5ad'


#mod_sk = unpickle_model(mod_path, f'n_fact{n_fact}')

adata_vis_sk = anndata.read(adata_file)

adata_snrna_raw = sc.read_h5ad('/scratch/mshruv003/Updated_integrated_MTG_annotation/reference_signature_cluster_label_batch/sc.h5ad')

adata_vis_sk.uns['mod']['fact_names']

location_factors_df = adata_vis_sk.uns['mod_coloc_n_fact15']['post_sample_means']['location_factors']
location_factors_df = pd.DataFrame(location_factors_df, 
                                  index= adata_vis_sk.uns['mod_coloc_n_fact15']['obs_names'],
                                  columns=adata_vis_sk.uns['mod_coloc_n_fact15']['fact_names'])
adata_vis_sk.obsm['mod_coloc_n_fact15'] = location_factors_df
adata_vis_sk.obs[location_factors_df.columns] = location_factors_df

from cell2location.utils import select_slide
slide = select_slide(adata_vis_sk, 'DH4a')

mpl.rcParams['pdf.fonttype'] = 42 # enables correct plotting of text
with mpl.rc_context({'figure.figsize': (10, 10), 'axes.facecolor': "black"}):
    sc.pl.spatial(slide, cmap='magma',
                  color=location_factors_df.columns, ncols=4, 
                  size=1.3, img_key='hires', alpha_img=0,
                  frameon=True, legend_fontsize=20,
                  vmin=0, vmax='p99.5', #save='Fig4D_suppl_spatial_clusters.pdf',
                  palette=sc.pl.palettes.default_102
                 ), plt.savefig("factor15_dh4a.jpeg")
				 
from cell2location.plt.mapping_video import plot_spatial

# select up to 6 clusters 
sel_clust = ['fact_9', 'fact_6', 'fact_10', 'fact_14', 'fact_12', 'fact_2']
sel_clust_col = ['' + str(i) for i in sel_clust]

slide = select_slide(adata_vis_sk, 'DH4a')

# identify spot locations to crop near tissue
crop_max = (slide.obsm['spatial'] \
            * list(slide.uns['spatial'].values())[0]['scalefactors']['tissue_hires_scalef']).max(axis=0)
crop_min = (slide.obsm['spatial'] \
            * list(slide.uns['spatial'].values())[0]['scalefactors']['tissue_hires_scalef']).min(axis=0)

crop_x = [crop_min[0]-50, crop_max[0]+50]
crop_y = [crop_max[1]+20, crop_min[1]-20]
					 
					 
mpl.rcParams['pdf.fonttype'] = 42 # enables correct plotting of text
with mpl.rc_context({'figure.figsize': (10, 10)}):
    fig = plot_spatial(slide.obs[sel_clust_col], 
                  coords=slide.obsm['spatial'] \
                          * list(slide.uns['spatial'].values())[0]['scalefactors']['tissue_hires_scalef'], 
                  show_img=True, img_alpha=0,
                  img=list(slide.uns['spatial'].values())[0]['images']['hires'],
                  max_color_quantile=0.99, 
                  crop_x=crop_x, crop_y=[crop_y[1]] + [crop_y[0]],
                  circle_diameter=5, labels=sel_clust, colorbar_position='right',
                 style='dark_background', 
                     colorbar_shape={'vertical_gaps': 2, 'horizontal_gaps': 0.18}), plt.savefig("DH4a_fact15_2.pdf")
                     
slide = select_slide(adata_vis_sk, 'DH1a')
 
 mpl.rcParams['pdf.fonttype'] = 42 # enables correct plotting of text
with mpl.rc_context({'figure.figsize': (10, 10)}):
    fig = plot_spatial(slide.obs[sel_clust_col], 
                  coords=slide.obsm['spatial'] \
                          * list(slide.uns['spatial'].values())[0]['scalefactors']['tissue_hires_scalef'], 
                  show_img=True, img_alpha=0,
                  img=list(slide.uns['spatial'].values())[0]['images']['hires'],
                  max_color_quantile=0.99, 
                  crop_x=crop_x, crop_y=[crop_y[1]] + [crop_y[0]],
                  circle_diameter=5, labels=sel_clust, colorbar_position='right',
                 style='dark_background', 
                     colorbar_shape={'vertical_gaps': 2, 'horizontal_gaps': 0.18}), plt.savefig("DH1a_15yr_fact15.pdf")                   
                     
mpl.rcParams['pdf.fonttype'] = 42 # enables correct plotting of text
with mpl.rc_context({'figure.figsize': (10, 10)}):
    fig = plot_spatial(slide.obs[sel_clust_col], 
                  coords=slide.obsm['spatial'] \
                          * list(slide.uns['spatial'].values())[0]['scalefactors']['tissue_hires_scalef'], 
                  show_img=True, img_alpha=0,
                  img=list(slide.uns['spatial'].values())[0]['images']['hires'],
                  max_color_quantile=0.99, 
                  crop_x=crop_x, crop_y=[crop_y[1]] + [crop_y[0]],
                  circle_diameter=5, labels=sel_clust, colorbar_position='right',
                 style='dark_background', 
                     colorbar_shape={'vertical_gaps': 2, 'horizontal_gaps': 0.18}), plt.savefig("DH1a_15yr_fact15.jpeg") 

slide = select_slide(adata_vis_sk, 'DH2a')
 
 mpl.rcParams['pdf.fonttype'] = 42 # enables correct plotting of text
with mpl.rc_context({'figure.figsize': (10, 10)}):
    fig = plot_spatial(slide.obs[sel_clust_col], 
                  coords=slide.obsm['spatial'] \
                          * list(slide.uns['spatial'].values())[0]['scalefactors']['tissue_hires_scalef'], 
                  show_img=True, img_alpha=0,
                  img=list(slide.uns['spatial'].values())[0]['images']['hires'],
                  max_color_quantile=0.99, 
                  crop_x=crop_x, crop_y=[crop_y[1]] + [crop_y[0]],
                  circle_diameter=5, labels=sel_clust, colorbar_position='right',
                 style='dark_background', 
                     colorbar_shape={'vertical_gaps': 2, 'horizontal_gaps': 0.18}), plt.savefig("DH2a_15yr_fact15.pdf")                   
                     
mpl.rcParams['pdf.fonttype'] = 42 # enables correct plotting of text
with mpl.rc_context({'figure.figsize': (10, 10)}):
    fig = plot_spatial(slide.obs[sel_clust_col], 
                  coords=slide.obsm['spatial'] \
                          * list(slide.uns['spatial'].values())[0]['scalefactors']['tissue_hires_scalef'], 
                  show_img=True, img_alpha=0,
                  img=list(slide.uns['spatial'].values())[0]['images']['hires'],
                  max_color_quantile=0.99, 
                  crop_x=crop_x, crop_y=[crop_y[1]] + [crop_y[0]],
                  circle_diameter=5, labels=sel_clust, colorbar_position='right',
                 style='dark_background', 
                     colorbar_shape={'vertical_gaps': 2, 'horizontal_gaps': 0.18}), plt.savefig("DH2a_15yr_fact15.jpeg")  

slide = select_slide(adata_vis_sk, 'DH3a')
 
 mpl.rcParams['pdf.fonttype'] = 42 # enables correct plotting of text
with mpl.rc_context({'figure.figsize': (10, 10)}):
    fig = plot_spatial(slide.obs[sel_clust_col], 
                  coords=slide.obsm['spatial'] \
                          * list(slide.uns['spatial'].values())[0]['scalefactors']['tissue_hires_scalef'], 
                  show_img=True, img_alpha=0,
                  img=list(slide.uns['spatial'].values())[0]['images']['hires'],
                  max_color_quantile=0.99, 
                  crop_x=crop_x, crop_y=[crop_y[1]] + [crop_y[0]],
                  circle_diameter=5, labels=sel_clust, colorbar_position='right',
                 style='dark_background', 
                     colorbar_shape={'vertical_gaps': 2, 'horizontal_gaps': 0.18}), plt.savefig("DH3a_31yr_fact15.pdf")                   
                     
mpl.rcParams['pdf.fonttype'] = 42 # enables correct plotting of text
with mpl.rc_context({'figure.figsize': (10, 10)}):
    fig = plot_spatial(slide.obs[sel_clust_col], 
                  coords=slide.obsm['spatial'] \
                          * list(slide.uns['spatial'].values())[0]['scalefactors']['tissue_hires_scalef'], 
                  show_img=True, img_alpha=0,
                  img=list(slide.uns['spatial'].values())[0]['images']['hires'],
                  max_color_quantile=0.99, 
                  crop_x=crop_x, crop_y=[crop_y[1]] + [crop_y[0]],
                  circle_diameter=5, labels=sel_clust, colorbar_position='right',
                 style='dark_background', 
                     colorbar_shape={'vertical_gaps': 2, 'horizontal_gaps': 0.18}), plt.savefig("DH3a_31yr_fact15.jpeg")    

slide = select_slide(adata_vis_sk, 'DH1')
 
 mpl.rcParams['pdf.fonttype'] = 42 # enables correct plotting of text
with mpl.rc_context({'figure.figsize': (10, 10)}):
    fig = plot_spatial(slide.obs[sel_clust_col], 
                  coords=slide.obsm['spatial'] \
                          * list(slide.uns['spatial'].values())[0]['scalefactors']['tissue_hires_scalef'], 
                  show_img=True, img_alpha=0,
                  img=list(slide.uns['spatial'].values())[0]['images']['hires'],
                  max_color_quantile=0.99, 
                  crop_x=crop_x, crop_y=[crop_y[1]] + [crop_y[0]],
                  circle_diameter=5, labels=sel_clust, colorbar_position='right',
                 style='dark_background', 
                     colorbar_shape={'vertical_gaps': 2, 'horizontal_gaps': 0.18}), plt.savefig("DH1_4yr_fact15.pdf")                   
                     
mpl.rcParams['pdf.fonttype'] = 42 # enables correct plotting of text
with mpl.rc_context({'figure.figsize': (10, 10)}):
    fig = plot_spatial(slide.obs[sel_clust_col], 
                  coords=slide.obsm['spatial'] \
                          * list(slide.uns['spatial'].values())[0]['scalefactors']['tissue_hires_scalef'], 
                  show_img=True, img_alpha=0,
                  img=list(slide.uns['spatial'].values())[0]['images']['hires'],
                  max_color_quantile=0.99, 
                  crop_x=crop_x, crop_y=[crop_y[1]] + [crop_y[0]],
                  circle_diameter=5, labels=sel_clust, colorbar_position='right',
                 style='dark_background', 
                     colorbar_shape={'vertical_gaps': 2, 'horizontal_gaps': 0.18}), plt.savefig("DH1_4yr_fact15.jpeg")     

slide = select_slide(adata_vis_sk, 'DH2')
 
mpl.rcParams['pdf.fonttype'] = 42 # enables correct plotting of text
with mpl.rc_context({'figure.figsize': (10, 10)}):
    fig = plot_spatial(slide.obs[sel_clust_col], 
                  coords=slide.obsm['spatial'] \
                          * list(slide.uns['spatial'].values())[0]['scalefactors']['tissue_hires_scalef'], 
                  show_img=True, img_alpha=0,
                  img=list(slide.uns['spatial'].values())[0]['images']['hires'],
                  max_color_quantile=0.99, 
                  crop_x=crop_x, crop_y=[crop_y[1]] + [crop_y[0]],
                  circle_diameter=5, labels=sel_clust, colorbar_position='right',
                 style='dark_background', 
                     colorbar_shape={'vertical_gaps': 2, 'horizontal_gaps': 0.18}), plt.savefig("DH2_4yr_fact15.pdf")                   
                     
mpl.rcParams['pdf.fonttype'] = 42 # enables correct plotting of text
with mpl.rc_context({'figure.figsize': (10, 10)}):
    fig = plot_spatial(slide.obs[sel_clust_col], 
                  coords=slide.obsm['spatial'] \
                          * list(slide.uns['spatial'].values())[0]['scalefactors']['tissue_hires_scalef'], 
                  show_img=True, img_alpha=0,
                  img=list(slide.uns['spatial'].values())[0]['images']['hires'],
                  max_color_quantile=0.99, 
                  crop_x=crop_x, crop_y=[crop_y[1]] + [crop_y[0]],
                  circle_diameter=5, labels=sel_clust, colorbar_position='right',
                 style='dark_background', 
                     colorbar_shape={'vertical_gaps': 2, 'horizontal_gaps': 0.18}), plt.savefig("DH2_4yr_fact15.jpeg")    

slide = select_slide(adata_vis_sk, 'DH3')
 
mpl.rcParams['pdf.fonttype'] = 42 # enables correct plotting of text
with mpl.rc_context({'figure.figsize': (10, 10)}):
    fig = plot_spatial(slide.obs[sel_clust_col], 
                  coords=slide.obsm['spatial'] \
                          * list(slide.uns['spatial'].values())[0]['scalefactors']['tissue_hires_scalef'], 
                  show_img=True, img_alpha=0,
                  img=list(slide.uns['spatial'].values())[0]['images']['hires'],
                  max_color_quantile=0.99, 
                  crop_x=crop_x, crop_y=[crop_y[1]] + [crop_y[0]],
                  circle_diameter=5, labels=sel_clust, colorbar_position='right',
                 style='dark_background', 
                     colorbar_shape={'vertical_gaps': 2, 'horizontal_gaps': 0.18}), plt.savefig("DH3_15yr_fact15.pdf")                   
                     
mpl.rcParams['pdf.fonttype'] = 42 # enables correct plotting of text
with mpl.rc_context({'figure.figsize': (10, 10)}):
    fig = plot_spatial(slide.obs[sel_clust_col], 
                  coords=slide.obsm['spatial'] \
                          * list(slide.uns['spatial'].values())[0]['scalefactors']['tissue_hires_scalef'], 
                  show_img=True, img_alpha=0,
                  img=list(slide.uns['spatial'].values())[0]['images']['hires'],
                  max_color_quantile=0.99, 
                  crop_x=crop_x, crop_y=[crop_y[1]] + [crop_y[0]],
                  circle_diameter=5, labels=sel_clust, colorbar_position='right',
                 style='dark_background', 
                     colorbar_shape={'vertical_gaps': 2, 'horizontal_gaps': 0.18}), plt.savefig("DH3_15yr_fact15.jpeg") 

slide = select_slide(adata_vis_sk, 'DH4')
 
 mpl.rcParams['pdf.fonttype'] = 42 # enables correct plotting of text
with mpl.rc_context({'figure.figsize': (10, 10)}):
    fig = plot_spatial(slide.obs[sel_clust_col], 
                  coords=slide.obsm['spatial'] \
                          * list(slide.uns['spatial'].values())[0]['scalefactors']['tissue_hires_scalef'], 
                  show_img=True, img_alpha=0,
                  img=list(slide.uns['spatial'].values())[0]['images']['hires'],
                  max_color_quantile=0.99, 
                  crop_x=crop_x, crop_y=[crop_y[1]] + [crop_y[0]],
                  circle_diameter=5, labels=sel_clust, colorbar_position='right',
                 style='dark_background', 
                     colorbar_shape={'vertical_gaps': 2, 'horizontal_gaps': 0.18}), plt.savefig("DH4_15yr_fact15.pdf")                   
                     
mpl.rcParams['pdf.fonttype'] = 42 # enables correct plotting of text
with mpl.rc_context({'figure.figsize': (10, 10)}):
    fig = plot_spatial(slide.obs[sel_clust_col], 
                  coords=slide.obsm['spatial'] \
                          * list(slide.uns['spatial'].values())[0]['scalefactors']['tissue_hires_scalef'], 
                  show_img=True, img_alpha=0,
                  img=list(slide.uns['spatial'].values())[0]['images']['hires'],
                  max_color_quantile=0.99, 
                  crop_x=crop_x, crop_y=[crop_y[1]] + [crop_y[0]],
                  circle_diameter=5, labels=sel_clust, colorbar_position='right',
                 style='dark_background', 
                     colorbar_shape={'vertical_gaps': 2, 'horizontal_gaps': 0.18}), plt.savefig("DH4_15yr_fact15.jpeg")                                        