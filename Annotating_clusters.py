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

run_name = '/scratch/mshruv003/cell2loc_integrated_with_Allen_MTG/cell2location_map'
sp_data_file = '/scratch/mshruv003/cell2loc_integrated_with_Allen_MTG/cell2location_map/sp_with_clusters_n21_res08.h5ad'
adata= anndata.read(sp_data_file)

mod = cell2location.models.Cell2location.load(f"{run_name}", adata)

from cell2location.utils import select_slide

adata.obs['new_clusters_DH4a'] = (
    adata.obs["region_cluster"]
    .map(lambda x: {"0":"A",
					"1":"B",
					"2":"B",
					"3":"B",
					"4":"B",
					"5":"B",
					"6":"B",
					"7":"B",
					"8":"B",
					"9":"B",
					"10":"B",
					"11":"B",
					"12":"B",
					"13":"B",
					"14":"B",
					"15":"B",
					"16":"B",
					"17":"B",
					"18":"B",
					"19":"B",
					"20":"B",
					"21":"B",
					"22":"B",
					"23":"B"}.get(x, x))
    .astype("category")
)

slide = select_slide(adata, 'DH4')
with mpl.rc_context({'axes.facecolor':  'black',
                            'figure.figsize': [10, 10]}):
    sc.pl.spatial(slide, #cmap='magma',
                  color=['new_clusters_DH4a'], ncols=4, 
                  #library_id=s,
                  size=1.3, img_key='hires', alpha_img=0,
                  frameon=True, legend_fontsize=20,
                  #crop_coord=crop_x + [crop_y[0]] + [crop_y[1]],
                  vmin=0, vmax='p99.5', save='DH4_cluster_0.png',
                  palette=sc.pl.palettes.default_102
                 )
				 
adata.obs['new_clusters_DH4a'] = (
    adata.obs["region_cluster"]
    .map(lambda x: {"0":"B",
					"1":"A",
					"2":"B",
					"3":"B",
					"4":"B",
					"5":"B",
					"6":"B",
					"7":"B",
					"8":"B",
					"9":"B",
					"10":"B",
					"11":"B",
					"12":"B",
					"13":"B",
					"14":"B",
					"15":"B",
					"16":"B",
					"17":"B",
					"18":"B",
					"19":"B",
					"20":"B",
					"21":"B",
					"22":"B",
					"23":"B"}.get(x, x))
    .astype("category")
)

slide = select_slide(adata, 'DH4')
with mpl.rc_context({'axes.facecolor':  'black',
                            'figure.figsize': [10, 10]}):
    sc.pl.spatial(slide, #cmap='magma',
                  color=['new_clusters_DH4a'], ncols=4, 
                  #library_id=s,
                  size=1.3, img_key='hires', alpha_img=0,
                  frameon=True, legend_fontsize=20,
                  #crop_coord=crop_x + [crop_y[0]] + [crop_y[1]],
                  vmin=0, vmax='p99.5', save='DH4_cluster_1.png',
                  palette=sc.pl.palettes.default_102
                 )
				 
adata.obs['new_clusters_DH4a'] = (
    adata.obs["region_cluster"]
    .map(lambda x: {"0":"B",
					"1":"B",
					"2":"A",
					"3":"B",
					"4":"B",
					"5":"B",
					"6":"B",
					"7":"B",
					"8":"B",
					"9":"B",
					"10":"B",
					"11":"B",
					"12":"B",
					"13":"B",
					"14":"B",
					"15":"B",
					"16":"B",
					"17":"B",
					"18":"B",
					"19":"B",
					"20":"B",
					"21":"B",
					"22":"B",
					"23":"B"}.get(x, x))
    .astype("category")
)

slide = select_slide(adata, 'DH4')
with mpl.rc_context({'axes.facecolor':  'black',
                            'figure.figsize': [10, 10]}):
    sc.pl.spatial(slide, #cmap='magma',
                  color=['new_clusters_DH4a'], ncols=4, 
                  #library_id=s,
                  size=1.3, img_key='hires', alpha_img=0,
                  frameon=True, legend_fontsize=20,
                  #crop_coord=crop_x + [crop_y[0]] + [crop_y[1]],
                  vmin=0, vmax='p99.5', save='DH4_cluster_2.png',
                  palette=sc.pl.palettes.default_102
                 )
				 
adata.obs['new_clusters_DH4a'] = (
    adata.obs["region_cluster"]
    .map(lambda x: {"0":"B",
					"1":"B",
					"2":"B",
					"3":"A",
					"4":"B",
					"5":"B",
					"6":"B",
					"7":"B",
					"8":"B",
					"9":"B",
					"10":"B",
					"11":"B",
					"12":"B",
					"13":"B",
					"14":"B",
					"15":"B",
					"16":"B",
					"17":"B",
					"18":"B",
					"19":"B",
					"20":"B",
					"21":"B",
					"22":"B",
					"23":"B"}.get(x, x))
    .astype("category")
)

slide = select_slide(adata, 'DH4')
with mpl.rc_context({'axes.facecolor':  'black',
                            'figure.figsize': [10, 10]}):
    sc.pl.spatial(slide, #cmap='magma',
                  color=['new_clusters_DH4a'], ncols=4, 
                  #library_id=s,
                  size=1.3, img_key='hires', alpha_img=0,
                  frameon=True, legend_fontsize=20,
                  #crop_coord=crop_x + [crop_y[0]] + [crop_y[1]],
                  vmin=0, vmax='p99.5', save='DH4_cluster_3.png',
                  palette=sc.pl.palettes.default_102
                 )
				 
				 
adata.obs['new_clusters_DH4a'] = (
    adata.obs["region_cluster"]
    .map(lambda x: {"0":"B",
					"1":"B",
					"2":"B",
					"3":"B",
					"4":"A",
					"5":"B",
					"6":"B",
					"7":"B",
					"8":"B",
					"9":"B",
					"10":"B",
					"11":"B",
					"12":"B",
					"13":"B",
					"14":"B",
					"15":"B",
					"16":"B",
					"17":"B",
					"18":"B",
					"19":"B",
					"20":"B",
					"21":"B",
					"22":"B",
					"23":"B"}.get(x, x))
    .astype("category")
)

slide = select_slide(adata, 'DH4')
with mpl.rc_context({'axes.facecolor':  'black',
                            'figure.figsize': [10, 10]}):
    sc.pl.spatial(slide, #cmap='magma',
                  color=['new_clusters_DH4a'], ncols=4, 
                  #library_id=s,
                  size=1.3, img_key='hires', alpha_img=0,
                  frameon=True, legend_fontsize=20,
                  #crop_coord=crop_x + [crop_y[0]] + [crop_y[1]],
                  vmin=0, vmax='p99.5', save='DH4_cluster_4.png',
                  palette=sc.pl.palettes.default_102
                 )
				 
adata.obs['new_clusters_DH4a'] = (
    adata.obs["region_cluster"]
    .map(lambda x: {"0":"B",
					"1":"B",
					"2":"B",
					"3":"B",
					"4":"B",
					"5":"A",
					"6":"B",
					"7":"B",
					"8":"B",
					"9":"B",
					"10":"B",
					"11":"B",
					"12":"B",
					"13":"B",
					"14":"B",
					"15":"B",
					"16":"B",
					"17":"B",
					"18":"B",
					"19":"B",
					"20":"B",
					"21":"B",
					"22":"B",
					"23":"B"}.get(x, x))
    .astype("category")
)

slide = select_slide(adata, 'DH4')
with mpl.rc_context({'axes.facecolor':  'black',
                            'figure.figsize': [10, 10]}):
    sc.pl.spatial(slide, #cmap='magma',
                  color=['new_clusters_DH4a'], ncols=4, 
                  #library_id=s,
                  size=1.3, img_key='hires', alpha_img=0,
                  frameon=True, legend_fontsize=20,
                  #crop_coord=crop_x + [crop_y[0]] + [crop_y[1]],
                  vmin=0, vmax='p99.5', save='DH4_cluster_5.png',
                  palette=sc.pl.palettes.default_102
                 )
				 
adata.obs['new_clusters_DH4a'] = (
    adata.obs["region_cluster"]
    .map(lambda x: {"0":"B",
					"1":"B",
					"2":"B",
					"3":"B",
					"4":"B",
					"5":"B",
					"6":"A",
					"7":"B",
					"8":"B",
					"9":"B",
					"10":"B",
					"11":"B",
					"12":"B",
					"13":"B",
					"14":"B",
					"15":"B",
					"16":"B",
					"17":"B",
					"18":"B",
					"19":"B",
					"20":"B",
					"21":"B",
					"22":"B",
					"23":"B"}.get(x, x))
    .astype("category")
)

slide = select_slide(adata, 'DH4')
with mpl.rc_context({'axes.facecolor':  'black',
                            'figure.figsize': [10, 10]}):
    sc.pl.spatial(slide, #cmap='magma',
                  color=['new_clusters_DH4a'], ncols=4, 
                  #library_id=s,
                  size=1.3, img_key='hires', alpha_img=0,
                  frameon=True, legend_fontsize=20,
                  #crop_coord=crop_x + [crop_y[0]] + [crop_y[1]],
                  vmin=0, vmax='p99.5', save='DH4_cluster_6.png',
                  palette=sc.pl.palettes.default_102
                 )
				 
adata.obs['new_clusters_DH4a'] = (
    adata.obs["region_cluster"]
    .map(lambda x: {"0":"B",
					"1":"B",
					"2":"B",
					"3":"B",
					"4":"B",
					"5":"B",
					"6":"B",
					"7":"A",
					"8":"B",
					"9":"B",
					"10":"B",
					"11":"B",
					"12":"B",
					"13":"B",
					"14":"B",
					"15":"B",
					"16":"B",
					"17":"B",
					"18":"B",
					"19":"B",
					"20":"B",
					"21":"B",
					"22":"B",
					"23":"B"}.get(x, x))
    .astype("category")
)

slide = select_slide(adata, 'DH4')
with mpl.rc_context({'axes.facecolor':  'black',
                            'figure.figsize': [10, 10]}):
    sc.pl.spatial(slide, #cmap='magma',
                  color=['new_clusters_DH4a'], ncols=4, 
                  #library_id=s,
                  size=1.3, img_key='hires', alpha_img=0,
                  frameon=True, legend_fontsize=20,
                  #crop_coord=crop_x + [crop_y[0]] + [crop_y[1]],
                  vmin=0, vmax='p99.5', save='DH4_cluster_7.png',
                  palette=sc.pl.palettes.default_102
                 )
				 
adata.obs['new_clusters_DH4a'] = (
    adata.obs["region_cluster"]
    .map(lambda x: {"0":"B",
					"1":"B",
					"2":"B",
					"3":"B",
					"4":"B",
					"5":"B",
					"6":"B",
					"7":"B",
					"8":"A",
					"9":"B",
					"10":"B",
					"11":"B",
					"12":"B",
					"13":"B",
					"14":"B",
					"15":"B",
					"16":"B",
					"17":"B",
					"18":"B",
					"19":"B",
					"20":"B",
					"21":"B",
					"22":"B",
					"23":"B"}.get(x, x))
    .astype("category")
)

slide = select_slide(adata, 'DH4')
with mpl.rc_context({'axes.facecolor':  'black',
                            'figure.figsize': [10, 10]}):
    sc.pl.spatial(slide, #cmap='magma',
                  color=['new_clusters_DH4a'], ncols=4, 
                  #library_id=s,
                  size=1.3, img_key='hires', alpha_img=0,
                  frameon=True, legend_fontsize=20,
                  #crop_coord=crop_x + [crop_y[0]] + [crop_y[1]],
                  vmin=0, vmax='p99.5', save='DH4_cluster_8.png',
                  palette=sc.pl.palettes.default_102
                 )
				 
adata.obs['new_clusters_DH4a'] = (
    adata.obs["region_cluster"]
    .map(lambda x: {"0":"B",
					"1":"B",
					"2":"B",
					"3":"B",
					"4":"B",
					"5":"B",
					"6":"B",
					"7":"B",
					"8":"B",
					"9":"A",
					"10":"B",
					"11":"B",
					"12":"B",
					"13":"B",
					"14":"B",
					"15":"B",
					"16":"B",
					"17":"B",
					"18":"B",
					"19":"B",
					"20":"B",
					"21":"B",
					"22":"B",
					"23":"B"}.get(x, x))
    .astype("category")
)

slide = select_slide(adata, 'DH4')
with mpl.rc_context({'axes.facecolor':  'black',
                            'figure.figsize': [10, 10]}):
    sc.pl.spatial(slide, #cmap='magma',
                  color=['new_clusters_DH4a'], ncols=4, 
                  #library_id=s,
                  size=1.3, img_key='hires', alpha_img=0,
                  frameon=True, legend_fontsize=20,
                  #crop_coord=crop_x + [crop_y[0]] + [crop_y[1]],
                  vmin=0, vmax='p99.5', save='DH4_cluster_9.png',
                  palette=sc.pl.palettes.default_102
                 )
				 
adata.obs['new_clusters_DH4a'] = (
    adata.obs["region_cluster"]
    .map(lambda x: {"0":"B",
					"1":"B",
					"2":"B",
					"3":"B",
					"4":"B",
					"5":"B",
					"6":"B",
					"7":"B",
					"8":"B",
					"9":"B",
					"10":"A",
					"11":"B",
					"12":"B",
					"13":"B",
					"14":"B",
					"15":"B",
					"16":"B",
					"17":"B",
					"18":"B",
					"19":"B",
					"20":"B",
					"21":"B",
					"22":"B",
					"23":"B"}.get(x, x))
    .astype("category")
)

slide = select_slide(adata, 'DH4')
with mpl.rc_context({'axes.facecolor':  'black',
                            'figure.figsize': [10, 10]}):
    sc.pl.spatial(slide, #cmap='magma',
                  color=['new_clusters_DH4a'], ncols=4, 
                  #library_id=s,
                  size=1.3, img_key='hires', alpha_img=0,
                  frameon=True, legend_fontsize=20,
                  #crop_coord=crop_x + [crop_y[0]] + [crop_y[1]],
                  vmin=0, vmax='p99.5', save='DH4_cluster_10.png',
                  palette=sc.pl.palettes.default_102
                 )
				 
adata.obs['new_clusters_DH4a'] = (
    adata.obs["region_cluster"]
    .map(lambda x: {"0":"B",
					"1":"B",
					"2":"B",
					"3":"B",
					"4":"B",
					"5":"B",
					"6":"B",
					"7":"B",
					"8":"B",
					"9":"B",
					"10":"B",
					"11":"A",
					"12":"B",
					"13":"B",
					"14":"B",
					"15":"B",
					"16":"B",
					"17":"B",
					"18":"B",
					"19":"B",
					"20":"B",
					"21":"B",
					"22":"B",
					"23":"B"}.get(x, x))
    .astype("category")
)

slide = select_slide(adata, 'DH4')
with mpl.rc_context({'axes.facecolor':  'black',
                            'figure.figsize': [10, 10]}):
    sc.pl.spatial(slide, #cmap='magma',
                  color=['new_clusters_DH4a'], ncols=4, 
                  #library_id=s,
                  size=1.3, img_key='hires', alpha_img=0,
                  frameon=True, legend_fontsize=20,
                  #crop_coord=crop_x + [crop_y[0]] + [crop_y[1]],
                  vmin=0, vmax='p99.5', save='DH4_cluster_11.png',
                  palette=sc.pl.palettes.default_102
                 )
				 
adata.obs['new_clusters_DH4a'] = (
    adata.obs["region_cluster"]
    .map(lambda x: {"0":"B",
					"1":"B",
					"2":"B",
					"3":"B",
					"4":"B",
					"5":"B",
					"6":"B",
					"7":"B",
					"8":"B",
					"9":"B",
					"10":"B",
					"11":"B",
					"12":"A",
					"13":"B",
					"14":"B",
					"15":"B",
					"16":"B",
					"17":"B",
					"18":"B",
					"19":"B",
					"20":"B",
					"21":"B",
					"22":"B",
					"23":"B"}.get(x, x))
    .astype("category")
)

slide = select_slide(adata, 'DH4')
with mpl.rc_context({'axes.facecolor':  'black',
                            'figure.figsize': [10, 10]}):
    sc.pl.spatial(slide, #cmap='magma',
                  color=['new_clusters_DH4a'], ncols=4, 
                  #library_id=s,
                  size=1.3, img_key='hires', alpha_img=0,
                  frameon=True, legend_fontsize=20,
                  #crop_coord=crop_x + [crop_y[0]] + [crop_y[1]],
                  vmin=0, vmax='p99.5', save='DH4_cluster_12.png',
                  palette=sc.pl.palettes.default_102
                 )
				 
adata.obs['new_clusters_DH4a'] = (
    adata.obs["region_cluster"]
    .map(lambda x: {"0":"B",
					"1":"B",
					"2":"B",
					"3":"B",
					"4":"B",
					"5":"B",
					"6":"B",
					"7":"B",
					"8":"B",
					"9":"B",
					"10":"B",
					"11":"B",
					"12":"B",
					"13":"A",
					"14":"B",
					"15":"B",
					"16":"B",
					"17":"B",
					"18":"B",
					"19":"B",
					"20":"B",
					"21":"B",
					"22":"B",
					"23":"B"}.get(x, x))
    .astype("category")
)

slide = select_slide(adata, 'DH4')
with mpl.rc_context({'axes.facecolor':  'black',
                            'figure.figsize': [10, 10]}):
    sc.pl.spatial(slide, #cmap='magma',
                  color=['new_clusters_DH4a'], ncols=4, 
                  #library_id=s,
                  size=1.3, img_key='hires', alpha_img=0,
                  frameon=True, legend_fontsize=20,
                  #crop_coord=crop_x + [crop_y[0]] + [crop_y[1]],
                  vmin=0, vmax='p99.5', save='DH4_cluster_13.png',
                  palette=sc.pl.palettes.default_102
                 )
				 
adata.obs['new_clusters_DH4a'] = (
    adata.obs["region_cluster"]
    .map(lambda x: {"0":"B",
					"1":"B",
					"2":"B",
					"3":"B",
					"4":"B",
					"5":"B",
					"6":"B",
					"7":"B",
					"8":"B",
					"9":"B",
					"10":"B",
					"11":"B",
					"12":"B",
					"13":"B",
					"14":"A",
					"15":"B",
					"16":"B",
					"17":"B",
					"18":"B",
					"19":"B",
					"20":"B",
					"21":"B",
					"22":"B",
					"23":"B"}.get(x, x))
    .astype("category")
)

slide = select_slide(adata, 'DH4')
with mpl.rc_context({'axes.facecolor':  'black',
                            'figure.figsize': [10, 10]}):
    sc.pl.spatial(slide, #cmap='magma',
                  color=['new_clusters_DH4a'], ncols=4, 
                  #library_id=s,
                  size=1.3, img_key='hires', alpha_img=0,
                  frameon=True, legend_fontsize=20,
                  #crop_coord=crop_x + [crop_y[0]] + [crop_y[1]],
                  vmin=0, vmax='p99.5', save='DH4_cluster_14.png',
                  palette=sc.pl.palettes.default_102
                 )
				 
adata.obs['new_clusters_DH4a'] = (
    adata.obs["region_cluster"]
    .map(lambda x: {"0":"B",
					"1":"B",
					"2":"B",
					"3":"B",
					"4":"B",
					"5":"B",
					"6":"B",
					"7":"B",
					"8":"B",
					"9":"B",
					"10":"B",
					"11":"B",
					"12":"B",
					"13":"B",
					"14":"B",
					"15":"A",
					"16":"B",
					"17":"B",
					"18":"B",
					"19":"B",
					"20":"B",
					"21":"B",
					"22":"B",
					"23":"B"}.get(x, x))
    .astype("category")
)

slide = select_slide(adata, 'DH4')
with mpl.rc_context({'axes.facecolor':  'black',
                            'figure.figsize': [10, 10]}):
    sc.pl.spatial(slide, #cmap='magma',
                  color=['new_clusters_DH4a'], ncols=4, 
                  #library_id=s,
                  size=1.3, img_key='hires', alpha_img=0,
                  frameon=True, legend_fontsize=20,
                  #crop_coord=crop_x + [crop_y[0]] + [crop_y[1]],
                  vmin=0, vmax='p99.5', save='DH4_cluster_15.png',
                  palette=sc.pl.palettes.default_102
                 )
				 
adata.obs['new_clusters_DH4a'] = (
    adata.obs["region_cluster"]
    .map(lambda x: {"0":"B",
					"1":"B",
					"2":"B",
					"3":"B",
					"4":"B",
					"5":"B",
					"6":"B",
					"7":"B",
					"8":"B",
					"9":"B",
					"10":"B",
					"11":"B",
					"12":"B",
					"13":"B",
					"14":"B",
					"15":"B",
					"16":"A",
					"17":"B",
					"18":"B",
					"19":"B",
					"20":"B",
					"21":"B",
					"22":"B",
					"23":"B"}.get(x, x))
    .astype("category")
)

slide = select_slide(adata, 'DH4')
with mpl.rc_context({'axes.facecolor':  'black',
                            'figure.figsize': [10, 10]}):
    sc.pl.spatial(slide, #cmap='magma',
                  color=['new_clusters_DH4a'], ncols=4, 
                  #library_id=s,
                  size=1.3, img_key='hires', alpha_img=0,
                  frameon=True, legend_fontsize=20,
                  #crop_coord=crop_x + [crop_y[0]] + [crop_y[1]],
                  vmin=0, vmax='p99.5', save='DH4_cluster_16.png',
                  palette=sc.pl.palettes.default_102
                 )
				 
adata.obs['new_clusters_DH4a'] = (
    adata.obs["region_cluster"]
    .map(lambda x: {"0":"B",
					"1":"B",
					"2":"B",
					"3":"B",
					"4":"B",
					"5":"B",
					"6":"B",
					"7":"B",
					"8":"B",
					"9":"B",
					"10":"B",
					"11":"B",
					"12":"B",
					"13":"B",
					"14":"B",
					"15":"B",
					"16":"B",
					"17":"A",
					"18":"B",
					"19":"B",
					"20":"B",
					"21":"B",
					"22":"B",
					"23":"B"}.get(x, x))
    .astype("category")
)

slide = select_slide(adata, 'DH4')
with mpl.rc_context({'axes.facecolor':  'black',
                            'figure.figsize': [10, 10]}):
    sc.pl.spatial(slide, #cmap='magma',
                  color=['new_clusters_DH4a'], ncols=4, 
                  #library_id=s,
                  size=1.3, img_key='hires', alpha_img=0,
                  frameon=True, legend_fontsize=20,
                  #crop_coord=crop_x + [crop_y[0]] + [crop_y[1]],
                  vmin=0, vmax='p99.5', save='DH4_cluster_17.png',
                  palette=sc.pl.palettes.default_102
                 )
				 
adata.obs['new_clusters_DH4a'] = (
    adata.obs["region_cluster"]
    .map(lambda x: {"0":"B",
					"1":"B",
					"2":"B",
					"3":"B",
					"4":"B",
					"5":"B",
					"6":"B",
					"7":"B",
					"8":"B",
					"9":"B",
					"10":"B",
					"11":"B",
					"12":"B",
					"13":"B",
					"14":"B",
					"15":"B",
					"16":"B",
					"17":"B",
					"18":"A",
					"19":"B",
					"20":"B",
					"21":"B",
					"22":"B",
					"23":"B"}.get(x, x))
    .astype("category")
)

slide = select_slide(adata, 'DH4')
with mpl.rc_context({'axes.facecolor':  'black',
                            'figure.figsize': [10, 10]}):
    sc.pl.spatial(slide, #cmap='magma',
                  color=['new_clusters_DH4a'], ncols=4, 
                  #library_id=s,
                  size=1.3, img_key='hires', alpha_img=0,
                  frameon=True, legend_fontsize=20,
                  #crop_coord=crop_x + [crop_y[0]] + [crop_y[1]],
                  vmin=0, vmax='p99.5', save='DH4_cluster_18.png',
                  palette=sc.pl.palettes.default_102
                 )
				 
adata.obs['new_clusters_DH4a'] = (
    adata.obs["region_cluster"]
    .map(lambda x: {"0":"B",
					"1":"B",
					"2":"B",
					"3":"B",
					"4":"B",
					"5":"B",
					"6":"B",
					"7":"B",
					"8":"B",
					"9":"B",
					"10":"B",
					"11":"B",
					"12":"B",
					"13":"B",
					"14":"B",
					"15":"B",
					"16":"B",
					"17":"B",
					"18":"B",
					"19":"A",
					"20":"B",
					"21":"B",
					"22":"B",
					"23":"B"}.get(x, x))
    .astype("category")
)

slide = select_slide(adata, 'DH4')
with mpl.rc_context({'axes.facecolor':  'black',
                            'figure.figsize': [10, 10]}):
    sc.pl.spatial(slide, #cmap='magma',
                  color=['new_clusters_DH4a'], ncols=4, 
                  #library_id=s,
                  size=1.3, img_key='hires', alpha_img=0,
                  frameon=True, legend_fontsize=20,
                  #crop_coord=crop_x + [crop_y[0]] + [crop_y[1]],
                  vmin=0, vmax='p99.5', save='DH4_cluster_19.png',
                  palette=sc.pl.palettes.default_102
                 )
				 
adata.obs['new_clusters_DH4a'] = (
    adata.obs["region_cluster"]
    .map(lambda x: {"0":"B",
					"1":"B",
					"2":"B",
					"3":"B",
					"4":"B",
					"5":"B",
					"6":"B",
					"7":"B",
					"8":"B",
					"9":"B",
					"10":"B",
					"11":"B",
					"12":"B",
					"13":"B",
					"14":"B",
					"15":"B",
					"16":"B",
					"17":"B",
					"18":"B",
					"19":"B",
					"20":"A",
					"21":"B",
					"22":"B",
					"23":"B"}.get(x, x))
    .astype("category")
)

slide = select_slide(adata, 'DH4')
with mpl.rc_context({'axes.facecolor':  'black',
                            'figure.figsize': [10, 10]}):
    sc.pl.spatial(slide, #cmap='magma',
                  color=['new_clusters_DH4a'], ncols=4, 
                  #library_id=s,
                  size=1.3, img_key='hires', alpha_img=0,
                  frameon=True, legend_fontsize=20,
                  #crop_coord=crop_x + [crop_y[0]] + [crop_y[1]],
                  vmin=0, vmax='p99.5', save='DH4_cluster_20.png',
                  palette=sc.pl.palettes.default_102
                 )
				 
adata.obs['new_clusters_DH4a'] = (
    adata.obs["region_cluster"]
    .map(lambda x: {"0":"B",
					"1":"B",
					"2":"B",
					"3":"B",
					"4":"B",
					"5":"B",
					"6":"B",
					"7":"B",
					"8":"B",
					"9":"B",
					"10":"B",
					"11":"B",
					"12":"B",
					"13":"B",
					"14":"B",
					"15":"B",
					"16":"B",
					"17":"B",
					"18":"B",
					"19":"B",
					"20":"B",
					"21":"A",
					"22":"B",
					"23":"B"}.get(x, x))
    .astype("category")
)

slide = select_slide(adata, 'DH4')
with mpl.rc_context({'axes.facecolor':  'black',
                            'figure.figsize': [10, 10]}):
    sc.pl.spatial(slide, #cmap='magma',
                  color=['new_clusters_DH4a'], ncols=4, 
                  #library_id=s,
                  size=1.3, img_key='hires', alpha_img=0,
                  frameon=True, legend_fontsize=20,
                  #crop_coord=crop_x + [crop_y[0]] + [crop_y[1]],
                  vmin=0, vmax='p99.5', save='DH4_cluster_21.png',
                  palette=sc.pl.palettes.default_102
                 )
				 
adata.obs['new_clusters_DH4a'] = (
    adata.obs["region_cluster"]
    .map(lambda x: {"0":"B",
					"1":"B",
					"2":"B",
					"3":"B",
					"4":"B",
					"5":"B",
					"6":"B",
					"7":"B",
					"8":"B",
					"9":"B",
					"10":"B",
					"11":"B",
					"12":"B",
					"13":"B",
					"14":"B",
					"15":"B",
					"16":"B",
					"17":"B",
					"18":"B",
					"19":"B",
					"20":"B",
					"21":"B",
					"22":"A",
					"23":"B"}.get(x, x))
    .astype("category")
)

slide = select_slide(adata, 'DH4')
with mpl.rc_context({'axes.facecolor':  'black',
                            'figure.figsize': [10, 10]}):
    sc.pl.spatial(slide, #cmap='magma',
                  color=['new_clusters_DH4a'], ncols=4, 
                  #library_id=s,
                  size=1.3, img_key='hires', alpha_img=0,
                  frameon=True, legend_fontsize=20,
                  #crop_coord=crop_x + [crop_y[0]] + [crop_y[1]],
                  vmin=0, vmax='p99.5', save='DH4_cluster_22.png',
                  palette=sc.pl.palettes.default_102
                 )
				 
adata.obs['new_clusters_DH4a'] = (
    adata.obs["region_cluster"]
    .map(lambda x: {"0":"B",
					"1":"B",
					"2":"B",
					"3":"B",
					"4":"B",
					"5":"B",
					"6":"B",
					"7":"B",
					"8":"B",
					"9":"B",
					"10":"B",
					"11":"B",
					"12":"B",
					"13":"B",
					"14":"B",
					"15":"B",
					"16":"B",
					"17":"B",
					"18":"B",
					"19":"B",
					"20":"B",
					"21":"B",
					"22":"B",
					"23":"A"}.get(x, x))
    .astype("category")
)

slide = select_slide(adata, 'DH4')
with mpl.rc_context({'axes.facecolor':  'black',
                            'figure.figsize': [10, 10]}):
    sc.pl.spatial(slide, #cmap='magma',
                  color=['new_clusters_DH4a'], ncols=4, 
                  #library_id=s,
                  size=1.3, img_key='hires', alpha_img=0,
                  frameon=True, legend_fontsize=20,
                  #crop_coord=crop_x + [crop_y[0]] + [crop_y[1]],
                  vmin=0, vmax='p99.5', save='DH4_cluster_23.png',
                  palette=sc.pl.palettes.default_102
                 )