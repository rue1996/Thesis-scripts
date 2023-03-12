#6th-step
#Creating figures
#Section: 2.5.3 Annotation and spatial mapping of cell types (Cell2location)
#Code adapted from: https://cell2location.readthedocs.io/en/latest/notebooks/cell2location_tutorial.html

import sys
import scanpy as sc
import anndata
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

import cell2location
import scvi

from matplotlib import rcParams
rcParams['pdf.fonttype'] = 42 # enables correct plotting of text
import seaborn as sns

fpath = "/scratch/mshruv003/Updated_integrated_MTG_annotation"
ref_run_name = f'{fpath}/reference_signature_cluster_label_batch'
run_name = f'{fpath}/cell2location_map'

adata_ref = sc.read("/scratch/mshruv003/Updated_integrated_MTG_annotation/Allen_MTG_seurat.h5ad")

from cell2location.utils.filtering import filter_genes
selected = filter_genes(adata_ref, cell_count_cutoff=5, cell_percentage_cutoff2=0.03, nonz_mean_cutoff=1.12)

# filter the object
adata_ref = adata_ref[~adata_ref.obs['Allen_high_res'].isin(['NA']),:]
adata_ref = adata_ref[:, selected].copy()

mod = cell2location.models.RegressionModel.load(f"{ref_run_name}", adata_ref)
adata_file = f"{ref_run_name}/sc.h5ad"
adata_ref = sc.read_h5ad(adata_file)

# export estimated expression in each cluster
if 'means_per_cluster_mu_fg' in adata_ref.varm.keys():
    inf_aver = adata_ref.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}'
                                    for i in adata_ref.uns['mod']['factor_names']]].copy()
else:
    inf_aver = adata_ref.var[[f'means_per_cluster_mu_fg_{i}'
                                    for i in adata_ref.uns['mod']['factor_names']]].copy()
    inf_aver.columns = adata_ref.uns['mod']['factor_names']
    inf_aver.iloc[0:5, 0:5]
	
adata_vis = sc.read("/scratch/mshruv003/Visium_datasets_combined/sp_all_datasets.h5ad")

# find mitochondria-encoded (MT) genes
adata_vis.var['MT_gene'] = [gene.startswith('MT-') for gene in adata_vis.var['SYMBOL']]

# remove MT genes for spatial mapping (keeping their counts in the object)
adata_vis.obsm['MT'] = adata_vis[:, adata_vis.var['MT_gene'].values].X.toarray()
adata_vis = adata_vis[:, ~adata_vis.var['MT_gene'].values]

# find shared genes and subset both anndata and reference signatures
intersect = np.intersect1d(adata_vis.var_names, inf_aver.index)
adata_vis = adata_vis[:, intersect].copy()
inf_aver = inf_aver.loc[intersect, :].copy()

mod = cell2location.models.Cell2location.load(f"{run_name}", adata_vis)
adata_file = f"{run_name}/sp.h5ad"
adata_vis = sc.read_h5ad(adata_file)

mod.export_posterior(adata_vis)
mod.plot_QC()

adata_vis.obs[adata_vis.uns['mod']['factor_names']] = adata_vis.obsm['q05_cell_abundance_w_sf']

from cell2location.utils import select_slide
slide = select_slide(adata_vis, 'DH1')

with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):
    sc.pl.spatial(slide, cmap='magma',
                  color=['means_per_cluster_mu_fg_Astro L1-2 FGFR3 GFAP', 'means_per_cluster_mu_fg_Endo L2-6 NOSTRIN', 'means_per_cluster_mu_fg_Micro L1-3 TYROBP', 'means_per_cluster_mu_fg_OPC L1-6 PDGFRA', 'means_per_cluster_mu_fg_Oligo L1-6 OPALIN','means_per_cluster_mu_fg_Exc L2 LAMP5 LTK', 'means_per_cluster_mu_fg_Exc L2-3 LINC00507 FREM3', 'means_per_cluster_mu_fg_Exc L4-6 RORB SEMA3E',
                        'means_per_cluster_mu_fg_Exc L5-6 FEZF2 ABO', 'means_per_cluster_mu_fg_Inh L4-5 SST STK32A', 'means_per_cluster_mu_fg_Inh L2-5 PVALB SCUBE3', 'means_per_cluster_mu_fg_Inh L2-5 VIP SERPINF1', 'means_per_cluster_mu_fg_Inh L2-6 LAMP5 CA1'],
                  ncols=6, size=1.3,
                  img_key='hires',
                  vmin=0, vmax='p99.2', save='DH1_4yr_celltypes_thesis.pdf'
                 )
				 
with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):
    sc.pl.spatial(slide, cmap='magma',
                  color=['means_per_cluster_mu_fg_Astro L1-2 FGFR3 GFAP', 'means_per_cluster_mu_fg_Endo L2-6 NOSTRIN', 'means_per_cluster_mu_fg_Micro L1-3 TYROBP', 'means_per_cluster_mu_fg_OPC L1-6 PDGFRA', 'means_per_cluster_mu_fg_Oligo L1-6 OPALIN','means_per_cluster_mu_fg_Exc L2 LAMP5 LTK',
                        'means_per_cluster_mu_fg_Inh L4-5 SST STK32A'],
                  ncols=6, size=1.3,
                  img_key='hires',
                  vmin=0, vmax='p99.2', save='DH1_4yr_celltypes_thesis2.pdf'
                 )				 

				 
with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):
    sc.pl.spatial(slide, cmap='magma',
                  color=['means_per_cluster_mu_fg_Astro L1-6 FGFR3 SLC14A1', 'means_per_cluster_mu_fg_Exc L2 LAMP5 LTK', 'means_per_cluster_mu_fg_Exc L2-3 LINC00507 FREM3', 'means_per_cluster_mu_fg_Exc L3-5 RORB COL22A1',
                        'means_per_cluster_mu_fg_Exc L3-5 RORB ESR1', 'means_per_cluster_mu_fg_Exc L3-5 RORB FILIP1L', 'means_per_cluster_mu_fg_Exc L3-5 RORB TWIST2'],
                  ncols=4, size=1.3,
                  img_key='hires',
                  vmin=0, vmax='p99.2', save='DH1_4yr_Exc_L3_5.pdf'
                 )
				 
with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):
    sc.pl.spatial(slide, cmap='magma',
                  color=['means_per_cluster_mu_fg_Exc L4-5 RORB DAPK2', 'means_per_cluster_mu_fg_Exc L4-5 RORB FOLH1B', 'means_per_cluster_mu_fg_Exc L4-6 FEZF2 IL26', 'means_per_cluster_mu_fg_Exc L4-6 RORB SEMA3E'], 
				  ncols=4, size=1.3,
                  img_key='hires',
                  vmin=0, vmax='p99.2', save='DH1_4yr_Exc_L4_6.pdf'
				 )
				 
with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):
    sc.pl.spatial(slide, cmap='magma',
                  color=['means_per_cluster_mu_fg_Exc L5-6 FEZF2 ABO', 'means_per_cluster_mu_fg_Exc L5-6 THEMIS C1QL3', 'means_per_cluster_mu_fg_Exc L5-6 THEMIS FGF10'
				        ],
                  ncols=4, size=1.3,
                  img_key='hires',
                  vmin=0, vmax='p99.2', save='DH1_4yr_Exc_L5_6.pdf'		
                  )


with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):
    sc.pl.spatial(slide, cmap='magma',
                  color=['means_per_cluster_mu_fg_Inh L1 SST NMBR', 'means_per_cluster_mu_fg_Inh L1-2 GAD1 MC4R', 'means_per_cluster_mu_fg_Inh L1-2 PAX6 CDH12', 
				  'means_per_cluster_mu_fg_Inh L1-2 SST BAGE2', 'means_per_cluster_mu_fg_Inh L1-2 VIP LBH', 'means_per_cluster_mu_fg_Inh L1-2 VIP PCDH20', 'means_per_cluster_mu_fg_Inh L1-2 VIP TSPAN12'],
                  ncols=4, size=1.3,
                  img_key='hires',
                  vmin=0, vmax='p99.2', save='DH1_4yr_Inh_L1_2.pdf'	
				  )
				  

with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):
    sc.pl.spatial(slide, cmap='magma',
                  color=['means_per_cluster_mu_fg_Inh L1-3 SST CALB1', 'means_per_cluster_mu_fg_Inh L1-3 VIP ADAMTSL1', 'means_per_cluster_mu_fg_Inh L1-3 VIP CCDC184', 'means_per_cluster_mu_fg_Inh L1-3 VIP CHRM2', 
				  'means_per_cluster_mu_fg_Inh L1-3 VIP GGH'],
                  ncols=4, size=1.3,
                  img_key='hires',
                  vmin=0, vmax='p99.2', save='DH1_4yr_Inh_L1_3.pdf'	
				  )
				  
				  
with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):
    sc.pl.spatial(slide, cmap='magma',
                  color=['means_per_cluster_mu_fg_Inh L1-4 LAMP5 LCP2', 'means_per_cluster_mu_fg_Inh L1-4 VIP CHRNA6', 'means_per_cluster_mu_fg_Inh L1-4 VIP OPRM1' 
				  ],
                  ncols=4, size=1.3,
                  img_key='hires',
                  vmin=0, vmax='p99.2', save='DH1_4yr_Inh_L1_4.pdf'  
                  )				  

with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):
    sc.pl.spatial(slide, cmap='magma',
                  color=['means_per_cluster_mu_fg_Inh L2-4 PVALB WFDC2', 'means_per_cluster_mu_fg_Inh L2-4 SST FRZB', 'means_per_cluster_mu_fg_Inh L2-4 VIP CBLN1', 'means_per_cluster_mu_fg_Inh L2-4 VIP SPAG17'],
                  ncols=4, size=1.3,
                  img_key='hires',
                  vmin=0, vmax='p99.2', save='DH1_4yr_Inh_L2_4.pdf'	
				  )
				  				 

with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):
    sc.pl.spatial(slide, cmap='magma',
                  color=['means_per_cluster_mu_fg_Inh L2-5 PVALB SCUBE3', 'means_per_cluster_mu_fg_Inh L2-5 VIP SERPINF1', 'means_per_cluster_mu_fg_Inh L2-5 VIP TYR', 'means_per_cluster_mu_fg_Inh L2-6 LAMP5 CA1', 
				  'means_per_cluster_mu_fg_Inh L2-6 VIP QPCT', 'means_per_cluster_mu_fg_Inh L3-5 SST ADGRG6', 'means_per_cluster_mu_fg_Inh L3-6 VIP HS3ST3A1'],
                  ncols=4, size=1.3,
                  img_key='hires',
                  vmin=0, vmax='p99.2', save='DH1_4yr_Inh_L2_6.pdf'	
				  )


with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):
    sc.pl.spatial(slide, cmap='magma',
                  color=['means_per_cluster_mu_fg_Inh L4-5 SST STK32A', 'means_per_cluster_mu_fg_Inh L4-6 PVALB SULF1', 'means_per_cluster_mu_fg_Inh L4-6 SST B3GAT2', 'means_per_cluster_mu_fg_Inh L4-6 SST GXYLT2', 
				  'means_per_cluster_mu_fg_Inh L5-6 PVALB LGR5', 'means_per_cluster_mu_fg_Inh L5-6 SST MIR548F2', 'means_per_cluster_mu_fg_Inh L5-6 SST NPM1P10', 'means_per_cluster_mu_fg_Micro L1-3 TYROBP'],
                  ncols=4, size=1.3,
                  img_key='hires',
                  vmin=0, vmax='p99.2', save='DH1_4yr_Inh_L5_6.pdf'  
				  )
			  

with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):
    sc.pl.spatial(slide, cmap='magma',
                  color=['means_per_cluster_mu_fg_Endo L2-6 NOSTRIN'],
                  ncols=4, size=1.3,
                  img_key='hires',
                  vmin=0, vmax='p99.2', save='DH1_4yr_Endo.pdf' 
				  )
				  
with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):
    sc.pl.spatial(slide, cmap='magma',
                  color=['means_per_cluster_mu_fg_OPC L1-6 PDGFRA', 'means_per_cluster_mu_fg_Oligo L1-6 OPALIN', 'means_per_cluster_mu_fg_Inh L3-6 SST NPY'],
                  ncols=4, size=1.3,
                  img_key='hires',
                  vmin=0, vmax='p99.2', save='DH1_4yr_Oligo.pdf' 
				  )
				  				  
				  
from cell2location.utils import select_slide
slide = select_slide(adata_vis, 'DH2')

with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):
    sc.pl.spatial(slide, cmap='magma',
                  color=['means_per_cluster_mu_fg_Astro L1-2 FGFR3 GFAP', 'means_per_cluster_mu_fg_Endo L2-6 NOSTRIN', 'means_per_cluster_mu_fg_Micro L1-3 TYROBP', 'means_per_cluster_mu_fg_OPC L1-6 PDGFRA', 'means_per_cluster_mu_fg_Oligo L1-6 OPALIN','means_per_cluster_mu_fg_Exc L2 LAMP5 LTK', 'means_per_cluster_mu_fg_Exc L2-3 LINC00507 FREM3', 'means_per_cluster_mu_fg_Exc L4-6 RORB SEMA3E',
                        'means_per_cluster_mu_fg_Exc L5-6 FEZF2 ABO', 'means_per_cluster_mu_fg_Inh L4-5 SST STK32A', 'means_per_cluster_mu_fg_Inh L2-5 PVALB SCUBE3', 'means_per_cluster_mu_fg_Inh L2-5 VIP SERPINF1', 'means_per_cluster_mu_fg_Inh L2-6 LAMP5 CA1'],
                  ncols=6, size=1.3,
                  img_key='hires',
                  vmin=0, vmax='p99.2', save='DH2_4yr_celltypes_thesis.pdf'
                 )
				 
with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):
    sc.pl.spatial(slide, cmap='magma',
                  color=['means_per_cluster_mu_fg_Astro L1-2 FGFR3 GFAP', 'means_per_cluster_mu_fg_Endo L2-6 NOSTRIN', 'means_per_cluster_mu_fg_Micro L1-3 TYROBP', 'means_per_cluster_mu_fg_OPC L1-6 PDGFRA', 'means_per_cluster_mu_fg_Oligo L1-6 OPALIN','means_per_cluster_mu_fg_Exc L2 LAMP5 LTK',
                        'means_per_cluster_mu_fg_Inh L4-5 SST STK32A'],
                  ncols=6, size=1.3,
                  img_key='hires',
                  vmin=0, vmax='p99.2', save='DH2_4yr_celltypes_thesis2.pdf'
                 )		

				 
with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):
    sc.pl.spatial(slide, cmap='magma',
                  color=['means_per_cluster_mu_fg_Astro L1-6 FGFR3 SLC14A1', 'means_per_cluster_mu_fg_Exc L2 LAMP5 LTK', 'means_per_cluster_mu_fg_Exc L2-3 LINC00507 FREM3', 'means_per_cluster_mu_fg_Exc L3-5 RORB COL22A1',
                        'means_per_cluster_mu_fg_Exc L3-5 RORB ESR1', 'means_per_cluster_mu_fg_Exc L3-5 RORB FILIP1L', 'means_per_cluster_mu_fg_Exc L3-5 RORB TWIST2'],
                  ncols=4, size=1.3,
                  img_key='hires',
                  vmin=0, vmax='p99.2', save='DH2_4yr_Exc_L3_5.pdf'
                 )
				 
with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):
    sc.pl.spatial(slide, cmap='magma',
                  color=['means_per_cluster_mu_fg_Exc L4-5 RORB DAPK2', 'means_per_cluster_mu_fg_Exc L4-5 RORB FOLH1B', 'means_per_cluster_mu_fg_Exc L4-6 FEZF2 IL26', 'means_per_cluster_mu_fg_Exc L4-6 RORB SEMA3E'], 
				  ncols=4, size=1.3,
                  img_key='hires',
                  vmin=0, vmax='p99.2', save='DH2_4yr_Exc_L4_6.pdf'
				 )
				 
with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):
    sc.pl.spatial(slide, cmap='magma',
                  color=['means_per_cluster_mu_fg_Exc L5-6 FEZF2 ABO', 'means_per_cluster_mu_fg_Exc L5-6 THEMIS C1QL3', 'means_per_cluster_mu_fg_Exc L5-6 THEMIS FGF10'
				        ],
                  ncols=4, size=1.3,
                  img_key='hires',
                  vmin=0, vmax='p99.2', save='DH2_4yr_Exc_L5_6.pdf'		
                  )


with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):
    sc.pl.spatial(slide, cmap='magma',
                  color=['means_per_cluster_mu_fg_Inh L1 SST NMBR', 'means_per_cluster_mu_fg_Inh L1-2 GAD1 MC4R', 'means_per_cluster_mu_fg_Inh L1-2 PAX6 CDH12',
				  'means_per_cluster_mu_fg_Inh L1-2 SST BAGE2', 'means_per_cluster_mu_fg_Inh L1-2 VIP LBH', 'means_per_cluster_mu_fg_Inh L1-2 VIP PCDH20', 'means_per_cluster_mu_fg_Inh L1-2 VIP TSPAN12'],
                  ncols=4, size=1.3,
                  img_key='hires',
                  vmin=0, vmax='p99.2', save='DH2_4yr_Inh_L1_2.pdf'	
				  )
				  

with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):
    sc.pl.spatial(slide, cmap='magma',
                  color=['means_per_cluster_mu_fg_Inh L1-3 SST CALB1', 'means_per_cluster_mu_fg_Inh L1-3 VIP ADAMTSL1', 'means_per_cluster_mu_fg_Inh L1-3 VIP CCDC184', 'means_per_cluster_mu_fg_Inh L1-3 VIP CHRM2', 
				  'means_per_cluster_mu_fg_Inh L1-3 VIP GGH'],
                  ncols=4, size=1.3,
                  img_key='hires',
                  vmin=0, vmax='p99.2', save='DH2_4yr_Inh_L1_3.pdf'	
				  )
				  
				  
with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):
    sc.pl.spatial(slide, cmap='magma',
                  color=['means_per_cluster_mu_fg_Inh L1-4 LAMP5 LCP2', 'means_per_cluster_mu_fg_Inh L1-4 VIP CHRNA6', 'means_per_cluster_mu_fg_Inh L1-4 VIP OPRM1' 
				  ],
                  ncols=4, size=1.3,
                  img_key='hires',
                  vmin=0, vmax='p99.2', save='DH2_4yr_Inh_L1_4.pdf'  
                  )				  

with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):
    sc.pl.spatial(slide, cmap='magma',
                  color=['means_per_cluster_mu_fg_Inh L2-4 PVALB WFDC2', 'means_per_cluster_mu_fg_Inh L2-4 SST FRZB', 'means_per_cluster_mu_fg_Inh L2-4 VIP CBLN1', 'means_per_cluster_mu_fg_Inh L2-4 VIP SPAG17'],
                  ncols=4, size=1.3,
                  img_key='hires',
                  vmin=0, vmax='p99.2', save='DH2_4yr_Inh_L2_4.pdf'	
				  )
				  
				  
				 

with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):
    sc.pl.spatial(slide, cmap='magma',
                  color=['means_per_cluster_mu_fg_Inh L2-5 PVALB SCUBE3', 'means_per_cluster_mu_fg_Inh L2-5 VIP SERPINF1', 'means_per_cluster_mu_fg_Inh L2-5 VIP TYR', 'means_per_cluster_mu_fg_Inh L2-6 LAMP5 CA1', 
				  'means_per_cluster_mu_fg_Inh L2-6 VIP QPCT', 'means_per_cluster_mu_fg_Inh L3-5 SST ADGRG6', 'means_per_cluster_mu_fg_Inh L3-6 VIP HS3ST3A1'],
                  ncols=4, size=1.3,
                  img_key='hires',
                  vmin=0, vmax='p99.2', save='DH2_4yr_Inh_L2_6.pdf'	
				  )


with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):
    sc.pl.spatial(slide, cmap='magma',
                  color=['means_per_cluster_mu_fg_Inh L4-5 SST STK32A', 'means_per_cluster_mu_fg_Inh L4-6 PVALB SULF1', 'means_per_cluster_mu_fg_Inh L4-6 SST B3GAT2', 'means_per_cluster_mu_fg_Inh L4-6 SST GXYLT2', 
				  'means_per_cluster_mu_fg_Inh L5-6 PVALB LGR5', 'means_per_cluster_mu_fg_Inh L5-6 SST MIR548F2', 'means_per_cluster_mu_fg_Inh L5-6 SST NPM1P10', 'means_per_cluster_mu_fg_Micro L1-3 TYROBP'],
                  ncols=4, size=1.3,
                  img_key='hires',
                  vmin=0, vmax='p99.2', save='DH2_4yr_Inh_L5_6.pdf'  
				  )
			  

with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):
    sc.pl.spatial(slide, cmap='magma',
                  color=['means_per_cluster_mu_fg_Endo L2-6 NOSTRIN'],
                  ncols=4, size=1.3,
                  img_key='hires',
                  vmin=0, vmax='p99.2', save='DH2_4yr_Endo.pdf' 
				  )

with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):
    sc.pl.spatial(slide, cmap='magma',
                  color=['means_per_cluster_mu_fg_OPC L1-6 PDGFRA', 'means_per_cluster_mu_fg_Oligo L1-6 OPALIN', 'means_per_cluster_mu_fg_Inh L3-6 SST NPY'],
                  ncols=4, size=1.3,
                  img_key='hires',
                  vmin=0, vmax='p99.2', save='DH2_4yr_Oligo.pdf' 
				  )
						  
from cell2location.utils import select_slide
slide = select_slide(adata_vis, 'DH3')

with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):
    sc.pl.spatial(slide, cmap='magma',
                  color=['means_per_cluster_mu_fg_Astro L1-2 FGFR3 GFAP', 'means_per_cluster_mu_fg_Endo L2-6 NOSTRIN', 'means_per_cluster_mu_fg_Micro L1-3 TYROBP', 'means_per_cluster_mu_fg_OPC L1-6 PDGFRA', 'means_per_cluster_mu_fg_Oligo L1-6 OPALIN','means_per_cluster_mu_fg_Exc L2 LAMP5 LTK', 'means_per_cluster_mu_fg_Exc L2-3 LINC00507 FREM3', 'means_per_cluster_mu_fg_Exc L4-6 RORB SEMA3E',
                        'means_per_cluster_mu_fg_Exc L5-6 FEZF2 ABO', 'means_per_cluster_mu_fg_Inh L4-5 SST STK32A', 'means_per_cluster_mu_fg_Inh L2-5 PVALB SCUBE3', 'means_per_cluster_mu_fg_Inh L2-5 VIP SERPINF1', 'means_per_cluster_mu_fg_Inh L2-6 LAMP5 CA1'],
                  ncols=6, size=1.3,
                  img_key='hires',
                  vmin=0, vmax='p99.2', save='DH3_2nd_15yr_celltypes_thesis.pdf'
                 )
				 
with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):
    sc.pl.spatial(slide, cmap='magma',
                  color=['means_per_cluster_mu_fg_Astro L1-2 FGFR3 GFAP', 'means_per_cluster_mu_fg_Endo L2-6 NOSTRIN', 'means_per_cluster_mu_fg_Micro L1-3 TYROBP', 'means_per_cluster_mu_fg_OPC L1-6 PDGFRA', 'means_per_cluster_mu_fg_Oligo L1-6 OPALIN','means_per_cluster_mu_fg_Exc L2 LAMP5 LTK',
                        'means_per_cluster_mu_fg_Inh L4-5 SST STK32A'],
                  ncols=6, size=1.3,
                  img_key='hires',
                  vmin=0, vmax='p99.2', save='DH3_2nd_15yr_celltypes_thesis2.pdf'
                 )		

				 
with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):
    sc.pl.spatial(slide, cmap='magma',
                  color=['means_per_cluster_mu_fg_Astro L1-6 FGFR3 SLC14A1', 'means_per_cluster_mu_fg_Exc L2 LAMP5 LTK', 'means_per_cluster_mu_fg_Exc L2-3 LINC00507 FREM3', 'means_per_cluster_mu_fg_Exc L3-5 RORB COL22A1',
                        'means_per_cluster_mu_fg_Exc L3-5 RORB ESR1', 'means_per_cluster_mu_fg_Exc L3-5 RORB FILIP1L', 'means_per_cluster_mu_fg_Exc L3-5 RORB TWIST2'],
                  ncols=4, size=1.3,
                  img_key='hires',
                  vmin=0, vmax='p99.2', save='DH3_2nd_15yr_Exc_L3_5.pdf'
                 )
				 
with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):
    sc.pl.spatial(slide, cmap='magma',
                  color=['means_per_cluster_mu_fg_Exc L4-5 RORB DAPK2', 'means_per_cluster_mu_fg_Exc L4-5 RORB FOLH1B', 'means_per_cluster_mu_fg_Exc L4-6 FEZF2 IL26', 'means_per_cluster_mu_fg_Exc L4-6 RORB SEMA3E'], 
				  ncols=4, size=1.3,
                  img_key='hires',
                  vmin=0, vmax='p99.2', save='DH3_2nd_15yr_Exc_L4_6.pdf'
				 )
				 
with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):
    sc.pl.spatial(slide, cmap='magma',
                  color=['means_per_cluster_mu_fg_Exc L5-6 FEZF2 ABO', 'means_per_cluster_mu_fg_Exc L5-6 THEMIS C1QL3', 'means_per_cluster_mu_fg_Exc L5-6 THEMIS FGF10'
				        ],
                  ncols=4, size=1.3,
                  img_key='hires',
                  vmin=0, vmax='p99.2', save='DH3_2nd_15yr_Exc_L5_6.pdf'		
                  )


with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):
    sc.pl.spatial(slide, cmap='magma',
                  color=['means_per_cluster_mu_fg_Inh L1 SST NMBR', 'means_per_cluster_mu_fg_Inh L1-2 GAD1 MC4R', 'means_per_cluster_mu_fg_Inh L1-2 PAX6 CDH12',
                         'means_per_cluster_mu_fg_Inh L1-2 SST BAGE2', 'means_per_cluster_mu_fg_Inh L1-2 VIP LBH', 'means_per_cluster_mu_fg_Inh L1-2 VIP PCDH20', 'means_per_cluster_mu_fg_Inh L1-2 VIP TSPAN12'],
                  ncols=4, size=1.3,
                  img_key='hires',
                  vmin=0, vmax='p99.2', save='DH3_2nd_15yr_Inh_L1_2.pdf'	
				  )
				  

with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):
    sc.pl.spatial(slide, cmap='magma',
                  color=['means_per_cluster_mu_fg_Inh L1-3 SST CALB1', 'means_per_cluster_mu_fg_Inh L1-3 VIP ADAMTSL1', 'means_per_cluster_mu_fg_Inh L1-3 VIP CCDC184', 'means_per_cluster_mu_fg_Inh L1-3 VIP CHRM2', 
				  'means_per_cluster_mu_fg_Inh L1-3 VIP GGH'],
                  ncols=4, size=1.3,
                  img_key='hires',
                  vmin=0, vmax='p99.2', save='DH3_2nd_15yr_Inh_L1_3.pdf'	
				  )
				  
				  
with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):
    sc.pl.spatial(slide, cmap='magma',
                  color=['means_per_cluster_mu_fg_Inh L1-4 LAMP5 LCP2', 'means_per_cluster_mu_fg_Inh L1-4 VIP CHRNA6', 'means_per_cluster_mu_fg_Inh L1-4 VIP OPRM1' 
				  ],
                  ncols=4, size=1.3,
                  img_key='hires',
                  vmin=0, vmax='p99.2', save='DH3_2nd_15yr_Inh_L1_4.pdf'  
                  )				  

with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):
    sc.pl.spatial(slide, cmap='magma',
                  color=['means_per_cluster_mu_fg_Inh L2-4 PVALB WFDC2', 'means_per_cluster_mu_fg_Inh L2-4 SST FRZB', 'means_per_cluster_mu_fg_Inh L2-4 VIP CBLN1', 'means_per_cluster_mu_fg_Inh L2-4 VIP SPAG17'],
                  ncols=4, size=1.3,
                  img_key='hires',
                  vmin=0, vmax='p99.2', save='DH3_2nd_15yr_Inh_L2_4.pdf'	
				  )
				  
				  
				 

with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):
    sc.pl.spatial(slide, cmap='magma',
                  color=['means_per_cluster_mu_fg_Inh L2-5 PVALB SCUBE3', 'means_per_cluster_mu_fg_Inh L2-5 VIP SERPINF1', 'means_per_cluster_mu_fg_Inh L2-5 VIP TYR', 'means_per_cluster_mu_fg_Inh L2-6 LAMP5 CA1', 
				  'means_per_cluster_mu_fg_Inh L2-6 VIP QPCT', 'means_per_cluster_mu_fg_Inh L3-5 SST ADGRG6', 'means_per_cluster_mu_fg_Inh L3-6 VIP HS3ST3A1'],
                  ncols=4, size=1.3,
                  img_key='hires',
                  vmin=0, vmax='p99.2', save='DH3_2nd_15yr_Inh_L2_6.pdf'	
				  )			  


with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):
    sc.pl.spatial(slide, cmap='magma',
                  color=['means_per_cluster_mu_fg_Inh L4-5 SST STK32A', 'means_per_cluster_mu_fg_Inh L4-6 PVALB SULF1', 'means_per_cluster_mu_fg_Inh L4-6 SST B3GAT2', 'means_per_cluster_mu_fg_Inh L4-6 SST GXYLT2', 
				  'means_per_cluster_mu_fg_Inh L5-6 PVALB LGR5', 'means_per_cluster_mu_fg_Inh L5-6 SST MIR548F2', 'means_per_cluster_mu_fg_Inh L5-6 SST NPM1P10', 'means_per_cluster_mu_fg_Micro L1-3 TYROBP'],
                  ncols=4, size=1.3,
                  img_key='hires',
                  vmin=0, vmax='p99.2', save='DH3_2nd_15yr_Inh_L5_6.pdf'  
				  )
			  

with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):
    sc.pl.spatial(slide, cmap='magma',
                  color=['means_per_cluster_mu_fg_Endo L2-6 NOSTRIN'],
                  ncols=4, size=1.3,
                  img_key='hires',
                  vmin=0, vmax='p99.2', save='DH3_2nd_15yr_Endo.pdf' 
				  )

with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):
    sc.pl.spatial(slide, cmap='magma',
                  color=['means_per_cluster_mu_fg_OPC L1-6 PDGFRA', 'means_per_cluster_mu_fg_Oligo L1-6 OPALIN', 'means_per_cluster_mu_fg_Inh L3-6 SST NPY'],
                  ncols=4, size=1.3,
                  img_key='hires',
                  vmin=0, vmax='p99.2', save='DH3_2nd_15yr_Oligo.pdf' 
				  )
				  
from cell2location.utils import select_slide
slide = select_slide(adata_vis, 'DH4')

with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):
    sc.pl.spatial(slide, cmap='magma',
                  color=['means_per_cluster_mu_fg_Astro L1-2 FGFR3 GFAP', 'means_per_cluster_mu_fg_Endo L2-6 NOSTRIN', 'means_per_cluster_mu_fg_Micro L1-3 TYROBP', 'means_per_cluster_mu_fg_OPC L1-6 PDGFRA', 'means_per_cluster_mu_fg_Oligo L1-6 OPALIN','means_per_cluster_mu_fg_Exc L2 LAMP5 LTK', 'means_per_cluster_mu_fg_Exc L2-3 LINC00507 FREM3', 'means_per_cluster_mu_fg_Exc L4-6 RORB SEMA3E',
                        'means_per_cluster_mu_fg_Exc L5-6 FEZF2 ABO', 'means_per_cluster_mu_fg_Inh L4-5 SST STK32A', 'means_per_cluster_mu_fg_Inh L2-5 PVALB SCUBE3', 'means_per_cluster_mu_fg_Inh L2-5 VIP SERPINF1', 'means_per_cluster_mu_fg_Inh L2-6 LAMP5 CA1'],
                  ncols=6, size=1.3,
                  img_key='hires',
                  vmin=0, vmax='p99.2', save='DH4_2nd_15yr_celltypes_thesis.pdf'
                 )
				 
with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):
    sc.pl.spatial(slide, cmap='magma',
                  color=['means_per_cluster_mu_fg_Astro L1-2 FGFR3 GFAP', 'means_per_cluster_mu_fg_Endo L2-6 NOSTRIN', 'means_per_cluster_mu_fg_Micro L1-3 TYROBP', 'means_per_cluster_mu_fg_OPC L1-6 PDGFRA', 'means_per_cluster_mu_fg_Oligo L1-6 OPALIN','means_per_cluster_mu_fg_Exc L2 LAMP5 LTK',
                        'means_per_cluster_mu_fg_Inh L4-5 SST STK32A'],
                  ncols=6, size=1.3,
                  img_key='hires',
                  vmin=0, vmax='p99.2', save='DH4_2nd_15yr_celltypes_thesis2.pdf'
                 )		

				 
with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):
    sc.pl.spatial(slide, cmap='magma',
                  color=['means_per_cluster_mu_fg_Astro L1-6 FGFR3 SLC14A1', 'means_per_cluster_mu_fg_Exc L2 LAMP5 LTK', 'means_per_cluster_mu_fg_Exc L2-3 LINC00507 FREM3', 'means_per_cluster_mu_fg_Exc L3-5 RORB COL22A1',
                        'means_per_cluster_mu_fg_Exc L3-5 RORB ESR1', 'means_per_cluster_mu_fg_Exc L3-5 RORB FILIP1L', 'means_per_cluster_mu_fg_Exc L3-5 RORB TWIST2'],
                  ncols=4, size=1.3,
                  img_key='hires',
                  vmin=0, vmax='p99.2', save='DH4_2nd_15yr_Exc_L3_5.pdf'
                 )
				 
with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):
    sc.pl.spatial(slide, cmap='magma',
                  color=['means_per_cluster_mu_fg_Exc L4-5 RORB DAPK2', 'means_per_cluster_mu_fg_Exc L4-5 RORB FOLH1B', 'means_per_cluster_mu_fg_Exc L4-6 FEZF2 IL26', 'means_per_cluster_mu_fg_Exc L4-6 RORB SEMA3E'], 
				  ncols=4, size=1.3,
                  img_key='hires',
                  vmin=0, vmax='p99.2', save='DH4_2nd_15yr_Exc_L4_6.pdf'
				 )
				 
with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):
    sc.pl.spatial(slide, cmap='magma',
                  color=['means_per_cluster_mu_fg_Exc L5-6 FEZF2 ABO', 'means_per_cluster_mu_fg_Exc L5-6 THEMIS C1QL3', 'means_per_cluster_mu_fg_Exc L5-6 THEMIS FGF10'
				        ],
                  ncols=4, size=1.3,
                  img_key='hires',
                  vmin=0, vmax='p99.2', save='DH4_2nd_15yr_Exc_L5_6.pdf'		
                  )


with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):
    sc.pl.spatial(slide, cmap='magma',
                  color=['means_per_cluster_mu_fg_Inh L1 SST NMBR', 'means_per_cluster_mu_fg_Inh L1-2 GAD1 MC4R', 'means_per_cluster_mu_fg_Inh L1-2 PAX6 CDH12',
				  'means_per_cluster_mu_fg_Inh L1-2 SST BAGE2', 'means_per_cluster_mu_fg_Inh L1-2 VIP LBH', 'means_per_cluster_mu_fg_Inh L1-2 VIP PCDH20', 'means_per_cluster_mu_fg_Inh L1-2 VIP TSPAN12'],
                  ncols=4, size=1.3,
                  img_key='hires',
                  vmin=0, vmax='p99.2', save='DH4_2nd_15yr_Inh_L1_2.pdf'	
				  )
				  

with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):
    sc.pl.spatial(slide, cmap='magma',
                  color=['means_per_cluster_mu_fg_Inh L1-3 SST CALB1', 'means_per_cluster_mu_fg_Inh L1-3 VIP ADAMTSL1', 'means_per_cluster_mu_fg_Inh L1-3 VIP CCDC184', 'means_per_cluster_mu_fg_Inh L1-3 VIP CHRM2', 
				  'means_per_cluster_mu_fg_Inh L1-3 VIP GGH'],
                  ncols=4, size=1.3,
                  img_key='hires',
                  vmin=0, vmax='p99.2', save='DH4_2nd_15yr_Inh_L1_3.pdf'	
				  )
				  
				  
with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):
    sc.pl.spatial(slide, cmap='magma',
                  color=['means_per_cluster_mu_fg_Inh L1-4 LAMP5 LCP2', 'means_per_cluster_mu_fg_Inh L1-4 VIP CHRNA6', 'means_per_cluster_mu_fg_Inh L1-4 VIP OPRM1' 
				  ],
                  ncols=4, size=1.3,
                  img_key='hires',
                  vmin=0, vmax='p99.2', save='DH4_2nd_15yr_Inh_L1_4.pdf'  
                  )				  

with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):
    sc.pl.spatial(slide, cmap='magma',
                  color=['means_per_cluster_mu_fg_Inh L2-4 PVALB WFDC2', 'means_per_cluster_mu_fg_Inh L2-4 SST FRZB', 'means_per_cluster_mu_fg_Inh L2-4 VIP CBLN1', 'means_per_cluster_mu_fg_Inh L2-4 VIP SPAG17'],
                  ncols=4, size=1.3,
                  img_key='hires',
                  vmin=0, vmax='p99.2', save='DH4_2nd_15yr_Inh_L2_4.pdf'	
				  )
				  
				  
				 

with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):
    sc.pl.spatial(slide, cmap='magma',
                  color=['means_per_cluster_mu_fg_Inh L2-5 PVALB SCUBE3', 'means_per_cluster_mu_fg_Inh L2-5 VIP SERPINF1', 'means_per_cluster_mu_fg_Inh L2-5 VIP TYR', 'means_per_cluster_mu_fg_Inh L2-6 LAMP5 CA1', 
				  'means_per_cluster_mu_fg_Inh L2-6 VIP QPCT', 'means_per_cluster_mu_fg_Inh L3-5 SST ADGRG6', 'means_per_cluster_mu_fg_Inh L3-6 VIP HS3ST3A1'],
                  ncols=4, size=1.3,
                  img_key='hires',
                  vmin=0, vmax='p99.2', save='DH4_2nd_15yr_Inh_L2_6.pdf'	
				  )


with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):
    sc.pl.spatial(slide, cmap='magma',
                  color=['means_per_cluster_mu_fg_Inh L4-5 SST STK32A', 'means_per_cluster_mu_fg_Inh L4-6 PVALB SULF1', 'means_per_cluster_mu_fg_Inh L4-6 SST B3GAT2', 'means_per_cluster_mu_fg_Inh L4-6 SST GXYLT2', 
				  'means_per_cluster_mu_fg_Inh L5-6 PVALB LGR5', 'means_per_cluster_mu_fg_Inh L5-6 SST MIR548F2', 'means_per_cluster_mu_fg_Inh L5-6 SST NPM1P10', 'means_per_cluster_mu_fg_Micro L1-3 TYROBP'],
                  ncols=4, size=1.3,
                  img_key='hires',
                  vmin=0, vmax='p99.2', save='DH4_2nd_15yr_Inh_L5_6.pdf'  
				  )
			  

with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):
    sc.pl.spatial(slide, cmap='magma',
                  color=['means_per_cluster_mu_fg_Endo L2-6 NOSTRIN'],
                  ncols=4, size=1.3,
                  img_key='hires',
                  vmin=0, vmax='p99.2', save='DH4_2nd_15yr_Endo.pdf' 
				  )

with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):
    sc.pl.spatial(slide, cmap='magma',
                  color=['means_per_cluster_mu_fg_OPC L1-6 PDGFRA', 'means_per_cluster_mu_fg_Oligo L1-6 OPALIN', 'means_per_cluster_mu_fg_Inh L3-6 SST NPY'],
                  ncols=4, size=1.3,
                  img_key='hires',
                  vmin=0, vmax='p99.2', save='DH4_2nd_15yr_Oligo.pdf' 
				  )
				  
from cell2location.utils import select_slide
slide = select_slide(adata_vis, 'DH2a')

with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):
    sc.pl.spatial(slide, cmap='magma',
                  color=['means_per_cluster_mu_fg_Astro L1-2 FGFR3 GFAP', 'means_per_cluster_mu_fg_Endo L2-6 NOSTRIN', 'means_per_cluster_mu_fg_Micro L1-3 TYROBP', 'means_per_cluster_mu_fg_OPC L1-6 PDGFRA', 'means_per_cluster_mu_fg_Oligo L1-6 OPALIN','means_per_cluster_mu_fg_Exc L2 LAMP5 LTK', 'means_per_cluster_mu_fg_Exc L2-3 LINC00507 FREM3', 'means_per_cluster_mu_fg_Exc L4-6 RORB SEMA3E',
                        'means_per_cluster_mu_fg_Exc L5-6 FEZF2 ABO', 'means_per_cluster_mu_fg_Inh L4-5 SST STK32A', 'means_per_cluster_mu_fg_Inh L2-5 PVALB SCUBE3', 'means_per_cluster_mu_fg_Inh L2-5 VIP SERPINF1', 'means_per_cluster_mu_fg_Inh L2-6 LAMP5 CA1'],
                  ncols=6, size=1.3,
                  img_key='hires',
                  vmin=0, vmax='p99.2', save='DH2a_first15yr_celltypes_thesis.pdf'
                 )
				 
with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):
    sc.pl.spatial(slide, cmap='magma',
                  color=['means_per_cluster_mu_fg_Astro L1-2 FGFR3 GFAP', 'means_per_cluster_mu_fg_Endo L2-6 NOSTRIN', 'means_per_cluster_mu_fg_Micro L1-3 TYROBP', 'means_per_cluster_mu_fg_OPC L1-6 PDGFRA', 'means_per_cluster_mu_fg_Oligo L1-6 OPALIN','means_per_cluster_mu_fg_Exc L2 LAMP5 LTK',
                        'means_per_cluster_mu_fg_Inh L4-5 SST STK32A'],
                  ncols=6, size=1.3,
                  img_key='hires',
                  vmin=0, vmax='p99.2', save='DH2a_first15yr_celltypes_thesis2.pdf'
                 )		

				 
with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):
    sc.pl.spatial(slide, cmap='magma',
                  color=['means_per_cluster_mu_fg_Astro L1-6 FGFR3 SLC14A1', 'means_per_cluster_mu_fg_Exc L2 LAMP5 LTK', 'means_per_cluster_mu_fg_Exc L2-3 LINC00507 FREM3', 'means_per_cluster_mu_fg_Exc L3-5 RORB COL22A1',
                        'means_per_cluster_mu_fg_Exc L3-5 RORB ESR1', 'means_per_cluster_mu_fg_Exc L3-5 RORB FILIP1L', 'means_per_cluster_mu_fg_Exc L3-5 RORB TWIST2'],
                  ncols=4, size=1.3,
                  img_key='hires',
                  vmin=0, vmax='p99.2', save='DH2a_first15yr_Exc_L3_5.pdf'
                 )
				 
with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):
    sc.pl.spatial(slide, cmap='magma',
                  color=['means_per_cluster_mu_fg_Exc L4-5 RORB DAPK2', 'means_per_cluster_mu_fg_Exc L4-5 RORB FOLH1B', 'means_per_cluster_mu_fg_Exc L4-6 FEZF2 IL26', 'means_per_cluster_mu_fg_Exc L4-6 RORB SEMA3E'], 
				  ncols=4, size=1.3,
                  img_key='hires',
                  vmin=0, vmax='p99.2', save='DH2a_first15yr_Exc_L4_6.pdf'
				 )
				 
with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):
    sc.pl.spatial(slide, cmap='magma',
                  color=['means_per_cluster_mu_fg_Exc L5-6 FEZF2 ABO', 'means_per_cluster_mu_fg_Exc L5-6 THEMIS C1QL3', 'means_per_cluster_mu_fg_Exc L5-6 THEMIS FGF10'
				        ],
                  ncols=4, size=1.3,
                  img_key='hires',
                  vmin=0, vmax='p99.2', save='DH2a_first15yr_Exc_L5_6.pdf'		
                  )


with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):
    sc.pl.spatial(slide, cmap='magma',
                  color=['means_per_cluster_mu_fg_Inh L1 SST NMBR', 'means_per_cluster_mu_fg_Inh L1-2 GAD1 MC4R', 'means_per_cluster_mu_fg_Inh L1-2 PAX6 CDH12', 
				  'means_per_cluster_mu_fg_Inh L1-2 SST BAGE2', 'means_per_cluster_mu_fg_Inh L1-2 VIP LBH', 'means_per_cluster_mu_fg_Inh L1-2 VIP PCDH20', 'means_per_cluster_mu_fg_Inh L1-2 VIP TSPAN12'],
                  ncols=4, size=1.3,
                  img_key='hires',
                  vmin=0, vmax='p99.2', save='DH2a_first15yr_Inh_L1_2.pdf'	
				  )
				  

with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):
    sc.pl.spatial(slide, cmap='magma',
                  color=['means_per_cluster_mu_fg_Inh L1-3 SST CALB1', 'means_per_cluster_mu_fg_Inh L1-3 VIP ADAMTSL1', 'means_per_cluster_mu_fg_Inh L1-3 VIP CCDC184', 'means_per_cluster_mu_fg_Inh L1-3 VIP CHRM2', 
				  'means_per_cluster_mu_fg_Inh L1-3 VIP GGH'],
                  ncols=4, size=1.3,
                  img_key='hires',
                  vmin=0, vmax='p99.2', save='DH2a_first15yr_Inh_L1_3.pdf'	
				  )
				  
				  
with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):
    sc.pl.spatial(slide, cmap='magma',
                  color=['means_per_cluster_mu_fg_Inh L1-4 LAMP5 LCP2', 'means_per_cluster_mu_fg_Inh L1-4 VIP CHRNA6', 'means_per_cluster_mu_fg_Inh L1-4 VIP OPRM1' 
				  ],
                  ncols=4, size=1.3,
                  img_key='hires',
                  vmin=0, vmax='p99.2', save='DH2a_first15yr_Inh_L1_4.pdf'  
                  )				  

with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):
    sc.pl.spatial(slide, cmap='magma',
                  color=['means_per_cluster_mu_fg_Inh L2-4 PVALB WFDC2', 'means_per_cluster_mu_fg_Inh L2-4 SST FRZB', 'means_per_cluster_mu_fg_Inh L2-4 VIP CBLN1', 'means_per_cluster_mu_fg_Inh L2-4 VIP SPAG17'],
                  ncols=4, size=1.3,
                  img_key='hires',
                  vmin=0, vmax='p99.2', save='DH2a_first15yr_Inh_L2_4.pdf'	
				  )
				  				 

with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):
    sc.pl.spatial(slide, cmap='magma',
                  color=['means_per_cluster_mu_fg_Inh L2-5 PVALB SCUBE3', 'means_per_cluster_mu_fg_Inh L2-5 VIP SERPINF1', 'means_per_cluster_mu_fg_Inh L2-5 VIP TYR', 'means_per_cluster_mu_fg_Inh L2-6 LAMP5 CA1', 
				  'means_per_cluster_mu_fg_Inh L2-6 VIP QPCT', 'means_per_cluster_mu_fg_Inh L3-5 SST ADGRG6', 'means_per_cluster_mu_fg_Inh L3-6 VIP HS3ST3A1'],
                  ncols=4, size=1.3,
                  img_key='hires',
                  vmin=0, vmax='p99.2', save='DH2a_first15yr_Inh_L2_6.pdf'	
				  )	  


with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):
    sc.pl.spatial(slide, cmap='magma',
                  color=['means_per_cluster_mu_fg_Inh L4-5 SST STK32A', 'means_per_cluster_mu_fg_Inh L4-6 PVALB SULF1', 'means_per_cluster_mu_fg_Inh L4-6 SST B3GAT2', 'means_per_cluster_mu_fg_Inh L4-6 SST GXYLT2', 
				  'means_per_cluster_mu_fg_Inh L5-6 PVALB LGR5', 'means_per_cluster_mu_fg_Inh L5-6 SST MIR548F2', 'means_per_cluster_mu_fg_Inh L5-6 SST NPM1P10', 'means_per_cluster_mu_fg_Micro L1-3 TYROBP'],
                  ncols=4, size=1.3,
                  img_key='hires',
                  vmin=0, vmax='p99.2', save='DH2a_first15yr_Inh_L5_6.pdf'  
				  )
			  

with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):
    sc.pl.spatial(slide, cmap='magma',
                  color=['means_per_cluster_mu_fg_Endo L2-6 NOSTRIN'],
                  ncols=4, size=1.3,
                  img_key='hires',
                  vmin=0, vmax='p99.2', save='DH2a_first15yr_Endo.pdf' 
				  )

with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):
    sc.pl.spatial(slide, cmap='magma',
                  color=['means_per_cluster_mu_fg_OPC L1-6 PDGFRA', 'means_per_cluster_mu_fg_Oligo L1-6 OPALIN', 'means_per_cluster_mu_fg_Inh L3-6 SST NPY'],
                  ncols=4, size=1.3,
                  img_key='hires',
                  vmin=0, vmax='p99.2', save='DH2a_first15yr_Oligo.pdf' 
				  )
				  
from cell2location.utils import select_slide
slide = select_slide(adata_vis, 'DH1a')

with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):
    sc.pl.spatial(slide, cmap='magma',
                  color=['means_per_cluster_mu_fg_Astro L1-2 FGFR3 GFAP', 'means_per_cluster_mu_fg_Endo L2-6 NOSTRIN', 'means_per_cluster_mu_fg_Micro L1-3 TYROBP', 'means_per_cluster_mu_fg_OPC L1-6 PDGFRA', 'means_per_cluster_mu_fg_Oligo L1-6 OPALIN','means_per_cluster_mu_fg_Exc L2 LAMP5 LTK', 'means_per_cluster_mu_fg_Exc L2-3 LINC00507 FREM3', 'means_per_cluster_mu_fg_Exc L4-6 RORB SEMA3E',
                        'means_per_cluster_mu_fg_Exc L5-6 FEZF2 ABO', 'means_per_cluster_mu_fg_Inh L4-5 SST STK32A', 'means_per_cluster_mu_fg_Inh L2-5 PVALB SCUBE3', 'means_per_cluster_mu_fg_Inh L2-5 VIP SERPINF1', 'means_per_cluster_mu_fg_Inh L2-6 LAMP5 CA1'],
                  ncols=6, size=1.3,
                  img_key='hires',
                  vmin=0, vmax='p99.2', save='DH1a_first15yr_celltypes_thesis.pdf'
                 )
				 
with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):
    sc.pl.spatial(slide, cmap='magma',
                  color=['means_per_cluster_mu_fg_Astro L1-2 FGFR3 GFAP', 'means_per_cluster_mu_fg_Endo L2-6 NOSTRIN', 'means_per_cluster_mu_fg_Micro L1-3 TYROBP', 'means_per_cluster_mu_fg_OPC L1-6 PDGFRA', 'means_per_cluster_mu_fg_Oligo L1-6 OPALIN','means_per_cluster_mu_fg_Exc L2 LAMP5 LTK',
                        'means_per_cluster_mu_fg_Inh L4-5 SST STK32A'],
                  ncols=6, size=1.3,
                  img_key='hires',
                  vmin=0, vmax='p99.2', save='DH1a_first15yr_celltypes_thesis2.pdf'
                 )		
				 
with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):
    sc.pl.spatial(slide, cmap='magma',
                  color=['means_per_cluster_mu_fg_Astro L1-6 FGFR3 SLC14A1', 'means_per_cluster_mu_fg_Exc L2 LAMP5 LTK', 'means_per_cluster_mu_fg_Exc L2-3 LINC00507 FREM3', 'means_per_cluster_mu_fg_Exc L3-5 RORB COL22A1',
                        'means_per_cluster_mu_fg_Exc L3-5 RORB ESR1', 'means_per_cluster_mu_fg_Exc L3-5 RORB FILIP1L', 'means_per_cluster_mu_fg_Exc L3-5 RORB TWIST2'],
                  ncols=4, size=1.3,
                  img_key='hires',
                  vmin=0, vmax='p99.2', save='DH1a_first15yr_Exc_L3_5.pdf'
                 )
				 
with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):
    sc.pl.spatial(slide, cmap='magma',
                  color=['means_per_cluster_mu_fg_Exc L4-5 RORB DAPK2', 'means_per_cluster_mu_fg_Exc L4-5 RORB FOLH1B', 'means_per_cluster_mu_fg_Exc L4-6 FEZF2 IL26', 'means_per_cluster_mu_fg_Exc L4-6 RORB SEMA3E'], 
				  ncols=4, size=1.3,
                  img_key='hires',
                  vmin=0, vmax='p99.2', save='DH1a_first15yr_Exc_L4_6.pdf'
				 )
				 
with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):
    sc.pl.spatial(slide, cmap='magma',
                  color=['means_per_cluster_mu_fg_Exc L5-6 FEZF2 ABO', 'means_per_cluster_mu_fg_Exc L5-6 THEMIS C1QL3', 'means_per_cluster_mu_fg_Exc L5-6 THEMIS FGF10'
				        ],
                  ncols=4, size=1.3,
                  img_key='hires',
                  vmin=0, vmax='p99.2', save='DH1a_first15yr_Exc_L5_6.pdf'		
                  )


with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):
    sc.pl.spatial(slide, cmap='magma',
                  color=['means_per_cluster_mu_fg_Inh L1 SST NMBR', 'means_per_cluster_mu_fg_Inh L1-2 GAD1 MC4R', 'means_per_cluster_mu_fg_Inh L1-2 PAX6 CDH12', 
				  'means_per_cluster_mu_fg_Inh L1-2 SST BAGE2', 'means_per_cluster_mu_fg_Inh L1-2 VIP LBH', 'means_per_cluster_mu_fg_Inh L1-2 VIP PCDH20', 'means_per_cluster_mu_fg_Inh L1-2 VIP TSPAN12'],
                  ncols=4, size=1.3,
                  img_key='hires',
                  vmin=0, vmax='p99.2', save='DH1a_first15yr_Inh_L1_2.pdf'	
				  )
				  

with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):
    sc.pl.spatial(slide, cmap='magma',
                  color=['means_per_cluster_mu_fg_Inh L1-3 SST CALB1', 'means_per_cluster_mu_fg_Inh L1-3 VIP ADAMTSL1', 'means_per_cluster_mu_fg_Inh L1-3 VIP CCDC184', 'means_per_cluster_mu_fg_Inh L1-3 VIP CHRM2', 
				  'means_per_cluster_mu_fg_Inh L1-3 VIP GGH'],
                  ncols=4, size=1.3,
                  img_key='hires',
                  vmin=0, vmax='p99.2', save='DH1a_first15yr_Inh_L1_3.pdf'	
				  )
				  
				  
with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):
    sc.pl.spatial(slide, cmap='magma',
                  color=['means_per_cluster_mu_fg_Inh L1-4 LAMP5 LCP2', 'means_per_cluster_mu_fg_Inh L1-4 VIP CHRNA6', 'means_per_cluster_mu_fg_Inh L1-4 VIP OPRM1' 
				  ],
                  ncols=4, size=1.3,
                  img_key='hires',
                  vmin=0, vmax='p99.2', save='DH1a_first15yr_Inh_L1_4.pdf'  
                  )				  

with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):
    sc.pl.spatial(slide, cmap='magma',
                  color=['means_per_cluster_mu_fg_Inh L2-4 PVALB WFDC2', 'means_per_cluster_mu_fg_Inh L2-4 SST FRZB', 'means_per_cluster_mu_fg_Inh L2-4 VIP CBLN1', 'means_per_cluster_mu_fg_Inh L2-4 VIP SPAG17'],
                  ncols=4, size=1.3,
                  img_key='hires',
                  vmin=0, vmax='p99.2', save='DH1a_first15yr_Inh_L2_4.pdf'	
				  )
				  
				  
				 

with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):
    sc.pl.spatial(slide, cmap='magma',
                  color=['means_per_cluster_mu_fg_Inh L2-5 PVALB SCUBE3', 'means_per_cluster_mu_fg_Inh L2-5 VIP SERPINF1', 'means_per_cluster_mu_fg_Inh L2-5 VIP TYR', 'means_per_cluster_mu_fg_Inh L2-6 LAMP5 CA1', 
				  'means_per_cluster_mu_fg_Inh L2-6 VIP QPCT', 'means_per_cluster_mu_fg_Inh L3-5 SST ADGRG6', 'means_per_cluster_mu_fg_Inh L3-6 VIP HS3ST3A1'],
                  ncols=4, size=1.3,
                  img_key='hires',
                  vmin=0, vmax='p99.2', save='DH1a_first15yr_Inh_L2_6.pdf'	
				  )			  


with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):
    sc.pl.spatial(slide, cmap='magma',
                  color=['means_per_cluster_mu_fg_Inh L4-5 SST STK32A', 'means_per_cluster_mu_fg_Inh L4-6 PVALB SULF1', 'means_per_cluster_mu_fg_Inh L4-6 SST B3GAT2', 'means_per_cluster_mu_fg_Inh L4-6 SST GXYLT2', 
				  'means_per_cluster_mu_fg_Inh L5-6 PVALB LGR5', 'means_per_cluster_mu_fg_Inh L5-6 SST MIR548F2', 'means_per_cluster_mu_fg_Inh L5-6 SST NPM1P10', 'means_per_cluster_mu_fg_Micro L1-3 TYROBP'],
                  ncols=4, size=1.3,
                  img_key='hires',
                  vmin=0, vmax='p99.2', save='DH1a_first15yr_Inh_L5_6.pdf'  
				  )
			  

with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):
    sc.pl.spatial(slide, cmap='magma',
                  color=['means_per_cluster_mu_fg_Endo L2-6 NOSTRIN'],
                  ncols=4, size=1.3,
                  img_key='hires',
                  vmin=0, vmax='p99.2', save='DH1a_first15yr_Endo.pdf' 
				  )

with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):
    sc.pl.spatial(slide, cmap='magma',
                  color=['means_per_cluster_mu_fg_OPC L1-6 PDGFRA', 'means_per_cluster_mu_fg_Oligo L1-6 OPALIN', 'means_per_cluster_mu_fg_Inh L3-6 SST NPY'],
                  ncols=4, size=1.3,
                  img_key='hires',
                  vmin=0, vmax='p99.2', save='DH1a_first15yr_Oligo.pdf' 
				  )
				  
from cell2location.utils import select_slide
slide = select_slide(adata_vis, 'DH3a')

with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):
    sc.pl.spatial(slide, cmap='magma',
                  color=['means_per_cluster_mu_fg_Astro L1-2 FGFR3 GFAP', 'means_per_cluster_mu_fg_Endo L2-6 NOSTRIN', 'means_per_cluster_mu_fg_Micro L1-3 TYROBP', 'means_per_cluster_mu_fg_OPC L1-6 PDGFRA', 'means_per_cluster_mu_fg_Oligo L1-6 OPALIN','means_per_cluster_mu_fg_Exc L2 LAMP5 LTK', 'means_per_cluster_mu_fg_Exc L2-3 LINC00507 FREM3', 'means_per_cluster_mu_fg_Exc L4-6 RORB SEMA3E',
                        'means_per_cluster_mu_fg_Exc L5-6 FEZF2 ABO', 'means_per_cluster_mu_fg_Inh L4-5 SST STK32A', 'means_per_cluster_mu_fg_Inh L2-5 PVALB SCUBE3', 'means_per_cluster_mu_fg_Inh L2-5 VIP SERPINF1', 'means_per_cluster_mu_fg_Inh L2-6 LAMP5 CA1'],
                  ncols=6, size=1.3,
                  img_key='hires',
                  vmin=0, vmax='p99.2', save='DH3a_31yr_celltypes_thesis.pdf'
                 )
				 
with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):
    sc.pl.spatial(slide, cmap='magma',
                  color=['means_per_cluster_mu_fg_Astro L1-2 FGFR3 GFAP', 'means_per_cluster_mu_fg_Endo L2-6 NOSTRIN', 'means_per_cluster_mu_fg_Micro L1-3 TYROBP', 'means_per_cluster_mu_fg_OPC L1-6 PDGFRA', 'means_per_cluster_mu_fg_Oligo L1-6 OPALIN','means_per_cluster_mu_fg_Exc L2 LAMP5 LTK',
                        'means_per_cluster_mu_fg_Inh L4-5 SST STK32A'],
                  ncols=6, size=1.3,
                  img_key='hires',
                  vmin=0, vmax='p99.2', save='DH3a_31yr_celltypes_thesis2.pdf'
                 )		

				 
with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):
    sc.pl.spatial(slide, cmap='magma',
                  color=['means_per_cluster_mu_fg_Astro L1-6 FGFR3 SLC14A1', 'means_per_cluster_mu_fg_Exc L2 LAMP5 LTK', 'means_per_cluster_mu_fg_Exc L2-3 LINC00507 FREM3', 'means_per_cluster_mu_fg_Exc L3-5 RORB COL22A1',
                        'means_per_cluster_mu_fg_Exc L3-5 RORB ESR1', 'means_per_cluster_mu_fg_Exc L3-5 RORB FILIP1L', 'means_per_cluster_mu_fg_Exc L3-5 RORB TWIST2'],
                  ncols=4, size=1.3,
                  img_key='hires',
                  vmin=0, vmax='p99.2', save='DH3a_31yr_Exc_L3_5.pdf'
                 )
				 
with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):
    sc.pl.spatial(slide, cmap='magma',
                  color=['means_per_cluster_mu_fg_Exc L4-5 RORB DAPK2', 'means_per_cluster_mu_fg_Exc L4-5 RORB FOLH1B', 'means_per_cluster_mu_fg_Exc L4-6 FEZF2 IL26', 'means_per_cluster_mu_fg_Exc L4-6 RORB SEMA3E'], 
				  ncols=4, size=1.3,
                  img_key='hires',
                  vmin=0, vmax='p99.2', save='DH3a_31yr_Exc_L4_6.pdf'
				 )
				 
with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):
    sc.pl.spatial(slide, cmap='magma',
                  color=['means_per_cluster_mu_fg_Exc L5-6 FEZF2 ABO', 'means_per_cluster_mu_fg_Exc L5-6 THEMIS C1QL3', 'means_per_cluster_mu_fg_Exc L5-6 THEMIS FGF10'
				        ],
                  ncols=4, size=1.3,
                  img_key='hires',
                  vmin=0, vmax='p99.2', save='DH3a_31yr_Exc_L5_6.pdf'		
                  )


with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):
    sc.pl.spatial(slide, cmap='magma',
                  color=['means_per_cluster_mu_fg_Inh L1 SST NMBR', 'means_per_cluster_mu_fg_Inh L1-2 GAD1 MC4R', 'means_per_cluster_mu_fg_Inh L1-2 PAX6 CDH12', 
				  'means_per_cluster_mu_fg_Inh L1-2 SST BAGE2', 'means_per_cluster_mu_fg_Inh L1-2 VIP LBH', 'means_per_cluster_mu_fg_Inh L1-2 VIP PCDH20', 'means_per_cluster_mu_fg_Inh L1-2 VIP TSPAN12'],
                  ncols=4, size=1.3,
                  img_key='hires',
                  vmin=0, vmax='p99.2', save='DH3a_31yr_Inh_L1_2.pdf'	
				  )
				  

with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):
    sc.pl.spatial(slide, cmap='magma',
                  color=['means_per_cluster_mu_fg_Inh L1-3 SST CALB1', 'means_per_cluster_mu_fg_Inh L1-3 VIP ADAMTSL1', 'means_per_cluster_mu_fg_Inh L1-3 VIP CCDC184', 'means_per_cluster_mu_fg_Inh L1-3 VIP CHRM2', 
				  'means_per_cluster_mu_fg_Inh L1-3 VIP GGH'],
                  ncols=4, size=1.3,
                  img_key='hires',
                  vmin=0, vmax='p99.2', save='DH3a_31yr_Inh_L1_3.pdf'	
				  )
				  
				  
with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):
    sc.pl.spatial(slide, cmap='magma',
                  color=['means_per_cluster_mu_fg_Inh L1-4 LAMP5 LCP2', 'means_per_cluster_mu_fg_Inh L1-4 VIP CHRNA6', 'means_per_cluster_mu_fg_Inh L1-4 VIP OPRM1' 
				  ],
                  ncols=4, size=1.3,
                  img_key='hires',
                  vmin=0, vmax='p99.2', save='DH3a_31yr_Inh_L1_4.pdf'  
                  )				  

with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):
    sc.pl.spatial(slide, cmap='magma',
                  color=['means_per_cluster_mu_fg_Inh L2-4 PVALB WFDC2', 'means_per_cluster_mu_fg_Inh L2-4 SST FRZB', 'means_per_cluster_mu_fg_Inh L2-4 VIP CBLN1', 'means_per_cluster_mu_fg_Inh L2-4 VIP SPAG17'],
                  ncols=4, size=1.3,
                  img_key='hires',
                  vmin=0, vmax='p99.2', save='DH3a_31yr_Inh_L2_4.pdf'	
				  )
				  
				  
				 

with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):
    sc.pl.spatial(slide, cmap='magma',
                  color=['means_per_cluster_mu_fg_Inh L2-5 PVALB SCUBE3', 'means_per_cluster_mu_fg_Inh L2-5 VIP SERPINF1', 'means_per_cluster_mu_fg_Inh L2-5 VIP TYR', 'means_per_cluster_mu_fg_Inh L2-6 LAMP5 CA1', 
				  'means_per_cluster_mu_fg_Inh L2-6 VIP QPCT', 'means_per_cluster_mu_fg_Inh L3-5 SST ADGRG6', 'means_per_cluster_mu_fg_Inh L3-6 VIP HS3ST3A1'],
                  ncols=4, size=1.3,
                  img_key='hires',
                  vmin=0, vmax='p99.2', save='DH3a_31yr_Inh_L2_6.pdf'	
				  )


with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):
    sc.pl.spatial(slide, cmap='magma',
                  color=['means_per_cluster_mu_fg_Inh L4-5 SST STK32A', 'means_per_cluster_mu_fg_Inh L4-6 PVALB SULF1', 'means_per_cluster_mu_fg_Inh L4-6 SST B3GAT2', 'means_per_cluster_mu_fg_Inh L4-6 SST GXYLT2', 
				  'means_per_cluster_mu_fg_Inh L5-6 PVALB LGR5', 'means_per_cluster_mu_fg_Inh L5-6 SST MIR548F2', 'means_per_cluster_mu_fg_Inh L5-6 SST NPM1P10', 'means_per_cluster_mu_fg_Micro L1-3 TYROBP'],
                  ncols=4, size=1.3,
                  img_key='hires',
                  vmin=0, vmax='p99.2', save='DH3a_31yr_Inh_L5_6.pdf'  
				  )
			  

with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):
    sc.pl.spatial(slide, cmap='magma',
                  color=['means_per_cluster_mu_fg_Endo L2-6 NOSTRIN'],
                  ncols=4, size=1.3,
                  img_key='hires',
                  vmin=0, vmax='p99.2', save='DH3a_31yr_Endo.pdf' 
				  )

with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):
    sc.pl.spatial(slide, cmap='magma',
                  color=['means_per_cluster_mu_fg_OPC L1-6 PDGFRA', 'means_per_cluster_mu_fg_Oligo L1-6 OPALIN', 'means_per_cluster_mu_fg_Inh L3-6 SST NPY'],
                  ncols=4, size=1.3,
                  img_key='hires',
                  vmin=0, vmax='p99.2', save='DH3a_31yr_Oligo.pdf' 
				  )
				  
from cell2location.utils import select_slide
slide = select_slide(adata_vis, 'DH4a')

with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):
    sc.pl.spatial(slide, cmap='magma',
                  color=['means_per_cluster_mu_fg_Astro L1-2 FGFR3 GFAP', 'means_per_cluster_mu_fg_Endo L2-6 NOSTRIN', 'means_per_cluster_mu_fg_Micro L1-3 TYROBP', 'means_per_cluster_mu_fg_OPC L1-6 PDGFRA', 'means_per_cluster_mu_fg_Oligo L1-6 OPALIN','means_per_cluster_mu_fg_Exc L2 LAMP5 LTK', 'means_per_cluster_mu_fg_Exc L2-3 LINC00507 FREM3', 'means_per_cluster_mu_fg_Exc L4-6 RORB SEMA3E',
                        'means_per_cluster_mu_fg_Exc L5-6 FEZF2 ABO', 'means_per_cluster_mu_fg_Inh L4-5 SST STK32A', 'means_per_cluster_mu_fg_Inh L2-5 PVALB SCUBE3', 'means_per_cluster_mu_fg_Inh L2-5 VIP SERPINF1', 'means_per_cluster_mu_fg_Inh L2-6 LAMP5 CA1'],
                  ncols=6, size=1.3,
                  img_key='hires',
                  vmin=0, vmax='p99.2', save='DH4a_31yr_celltypes_thesis.pdf'
                 )
				 
with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):
    sc.pl.spatial(slide, cmap='magma',
                  color=['means_per_cluster_mu_fg_Astro L1-2 FGFR3 GFAP', 'means_per_cluster_mu_fg_Endo L2-6 NOSTRIN', 'means_per_cluster_mu_fg_Micro L1-3 TYROBP', 'means_per_cluster_mu_fg_OPC L1-6 PDGFRA', 'means_per_cluster_mu_fg_Oligo L1-6 OPALIN','means_per_cluster_mu_fg_Exc L2 LAMP5 LTK',
                        'means_per_cluster_mu_fg_Inh L4-5 SST STK32A'],
                  ncols=6, size=1.3,
                  img_key='hires',
                  vmin=0, vmax='p99.2', save='DH4a_31yr_celltypes_thesis2.pdf'
                 )		

				 
with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):
    sc.pl.spatial(slide, cmap='magma',
                  color=['means_per_cluster_mu_fg_Astro L1-6 FGFR3 SLC14A1', 'means_per_cluster_mu_fg_Exc L2 LAMP5 LTK', 'means_per_cluster_mu_fg_Exc L2-3 LINC00507 FREM3', 'means_per_cluster_mu_fg_Exc L3-5 RORB COL22A1',
                        'means_per_cluster_mu_fg_Exc L3-5 RORB ESR1', 'means_per_cluster_mu_fg_Exc L3-5 RORB FILIP1L', 'means_per_cluster_mu_fg_Exc L3-5 RORB TWIST2'],
                  ncols=4, size=1.3,
                  img_key='hires',
                  vmin=0, vmax='p99.2', save='DH4a_31yr_Exc_L3_5.pdf'
                 )
				 
with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):
    sc.pl.spatial(slide, cmap='magma',
                  color=['means_per_cluster_mu_fg_Exc L4-5 RORB DAPK2', 'means_per_cluster_mu_fg_Exc L4-5 RORB FOLH1B', 'means_per_cluster_mu_fg_Exc L4-6 FEZF2 IL26', 'means_per_cluster_mu_fg_Exc L4-6 RORB SEMA3E'], 
				  ncols=4, size=1.3,
                  img_key='hires',
                  vmin=0, vmax='p99.2', save='DH4a_31yr_Exc_L4_6.pdf'
				 )
				 
with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):
    sc.pl.spatial(slide, cmap='magma',
                  color=['means_per_cluster_mu_fg_Exc L5-6 FEZF2 ABO', 'means_per_cluster_mu_fg_Exc L5-6 THEMIS C1QL3', 'means_per_cluster_mu_fg_Exc L5-6 THEMIS FGF10'
				        ],
                  ncols=4, size=1.3,
                  img_key='hires',
                  vmin=0, vmax='p99.2', save='DH4a_31yr_Exc_L5_6.pdf'		
                  )


with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):
    sc.pl.spatial(slide, cmap='magma',
                  color=['means_per_cluster_mu_fg_Inh L1 SST NMBR', 'means_per_cluster_mu_fg_Inh L1-2 GAD1 MC4R', 'means_per_cluster_mu_fg_Inh L1-2 PAX6 CDH12', 
				  'means_per_cluster_mu_fg_Inh L1-2 SST BAGE2', 'means_per_cluster_mu_fg_Inh L1-2 VIP LBH', 'means_per_cluster_mu_fg_Inh L1-2 VIP PCDH20', 'means_per_cluster_mu_fg_Inh L1-2 VIP TSPAN12'],
                  ncols=4, size=1.3,
                  img_key='hires',
                  vmin=0, vmax='p99.2', save='DH4a_31yr_Inh_L1_2.pdf'	
				  )
				  

with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):
    sc.pl.spatial(slide, cmap='magma',
                  color=['means_per_cluster_mu_fg_Inh L1-3 SST CALB1', 'means_per_cluster_mu_fg_Inh L1-3 VIP ADAMTSL1', 'means_per_cluster_mu_fg_Inh L1-3 VIP CCDC184', 'means_per_cluster_mu_fg_Inh L1-3 VIP CHRM2', 
				  'means_per_cluster_mu_fg_Inh L1-3 VIP GGH'],
                  ncols=4, size=1.3,
                  img_key='hires',
                  vmin=0, vmax='p99.2', save='DH4a_31yr_Inh_L1_3.pdf'	
				  )
				  
				  
with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):
    sc.pl.spatial(slide, cmap='magma',
                  color=['means_per_cluster_mu_fg_Inh L1-4 LAMP5 LCP2', 'means_per_cluster_mu_fg_Inh L1-4 VIP CHRNA6', 'means_per_cluster_mu_fg_Inh L1-4 VIP OPRM1' 
				  ],
                  ncols=4, size=1.3,
                  img_key='hires',
                  vmin=0, vmax='p99.2', save='DH4a_31yr_Inh_L1_4.pdf'  
                  )				  

with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):
    sc.pl.spatial(slide, cmap='magma',
                  color=['means_per_cluster_mu_fg_Inh L2-4 PVALB WFDC2', 'means_per_cluster_mu_fg_Inh L2-4 SST FRZB', 'means_per_cluster_mu_fg_Inh L2-4 VIP CBLN1', 'means_per_cluster_mu_fg_Inh L2-4 VIP SPAG17'],
                  ncols=4, size=1.3,
                  img_key='hires',
                  vmin=0, vmax='p99.2', save='DH4a_31yr_Inh_L2_4.pdf'	
				  )
				  
				  
				 

with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):
    sc.pl.spatial(slide, cmap='magma',
                  color=['means_per_cluster_mu_fg_Inh L2-5 PVALB SCUBE3', 'means_per_cluster_mu_fg_Inh L2-5 VIP SERPINF1', 'means_per_cluster_mu_fg_Inh L2-5 VIP TYR', 'means_per_cluster_mu_fg_Inh L2-6 LAMP5 CA1', 
				  'means_per_cluster_mu_fg_Inh L2-6 VIP QPCT', 'means_per_cluster_mu_fg_Inh L3-5 SST ADGRG6', 'means_per_cluster_mu_fg_Inh L3-6 VIP HS3ST3A1'],
                  ncols=4, size=1.3,
                  img_key='hires',
                  vmin=0, vmax='p99.2', save='DH4a_31yr_Inh_L2_6.pdf'	
				  )		  


with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):
    sc.pl.spatial(slide, cmap='magma',
                  color=['means_per_cluster_mu_fg_Inh L4-5 SST STK32A', 'means_per_cluster_mu_fg_Inh L4-6 PVALB SULF1', 'means_per_cluster_mu_fg_Inh L4-6 SST B3GAT2', 'means_per_cluster_mu_fg_Inh L4-6 SST GXYLT2', 
				  'means_per_cluster_mu_fg_Inh L5-6 PVALB LGR5', 'means_per_cluster_mu_fg_Inh L5-6 SST MIR548F2', 'means_per_cluster_mu_fg_Inh L5-6 SST NPM1P10', 'means_per_cluster_mu_fg_Micro L1-3 TYROBP'],
                  ncols=4, size=1.3,
                  img_key='hires',
                  vmin=0, vmax='p99.2', save='DH4a_31yr_Inh_L5_6.pdf'  
				  )
			  

with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):
    sc.pl.spatial(slide, cmap='magma',
                  color=['means_per_cluster_mu_fg_Endo L2-6 NOSTRIN'],
                  ncols=4, size=1.3,
                  img_key='hires',
                  vmin=0, vmax='p99.2', save='DH4a_31yr_Endo.pdf' 
				  )
				  

with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):
    sc.pl.spatial(slide, cmap='magma',
                  color=['means_per_cluster_mu_fg_OPC L1-6 PDGFRA', 'means_per_cluster_mu_fg_Oligo L1-6 OPALIN', 'means_per_cluster_mu_fg_Inh L3-6 SST NPY'],
                  ncols=4, size=1.3,
                  img_key='hires',
                  vmin=0, vmax='p99.2', save='DH3a_31yr_Oligo.pdf' 
				  )
				  
from cell2location.plt import plot_spatial

# select up to 6 clusters
clust_labels = ['means_per_cluster_mu_fg_Astro L1-2 FGFR3 GFAP', 'means_per_cluster_mu_fg_Astro L1-6 FGFR3 SLC14A1'
                        ]
clust_col = ['' + str(i) for i in clust_labels] # in case column names differ from labels

slide = select_slide(adata_vis, 'DH1')

with mpl.rc_context({'figure.figsize': (20, 20)}):
    fig = plot_spatial(
        adata=slide,
        # labels to show on a plot
        color=clust_col, labels=clust_labels,
		alpha_scaling=1.0,
        show_img=True,
		img_alpha=0.0,
		adjust_text=True,
        # 'fast' (white background) or 'dark_background'
        style='dark_background',
        # limit color scale at 99.2% quantile of cell abundance
        max_color_quantile=0.992,
        # size of locations (adjust depending on figure size)
        circle_diameter=9,
        colorbar_position='right'
    ), plt.savefig("DH1_4yr_Astro_clusters.pdf")

clust_labels = ['means_per_cluster_mu_fg_Exc L2 LAMP5 LTK', 'means_per_cluster_mu_fg_Exc L2-3 LINC00507 FREM3', 'means_per_cluster_mu_fg_Exc L3-5 RORB COL22A1',
                        'means_per_cluster_mu_fg_Exc L3-5 RORB ESR1', 'means_per_cluster_mu_fg_Exc L3-5 RORB FILIP1L', 'means_per_cluster_mu_fg_Exc L3-5 RORB TWIST2'
                        ]
clust_col = ['' + str(i) for i in clust_labels] # in case column names differ from labels

slide = select_slide(adata_vis, 'DH1')

with mpl.rc_context({'figure.figsize': (20, 20)}):
    fig = plot_spatial(
        adata=slide,
        # labels to show on a plot
        color=clust_col, labels=clust_labels,
		alpha_scaling=1.0,
        show_img=True,
		img_alpha=0.0,
		adjust_text=True,
        # 'fast' (white background) or 'dark_background'
        style='dark_background',
        # limit color scale at 99.2% quantile of cell abundance
        max_color_quantile=0.992,
        # size of locations (adjust depending on figure size)
        circle_diameter=9,
        colorbar_position='right'
    ), plt.savefig("DH1_4yr_Exc_L3_5_clusters.pdf")
	
clust_labels = ['means_per_cluster_mu_fg_Exc L4-5 RORB DAPK2', 'means_per_cluster_mu_fg_Exc L4-5 RORB FOLH1B', 'means_per_cluster_mu_fg_Exc L4-6 FEZF2 IL26', 'means_per_cluster_mu_fg_Exc L4-6 RORB SEMA3E'
                        ]
clust_col = ['' + str(i) for i in clust_labels] # in case column names differ from labels

slide = select_slide(adata_vis, 'DH1')

with mpl.rc_context({'figure.figsize': (20, 20)}):
    fig = plot_spatial(
        adata=slide,
        # labels to show on a plot
        color=clust_col, labels=clust_labels,
		alpha_scaling=1.0,
        show_img=True,
		img_alpha=0.0,
		adjust_text=True,
        # 'fast' (white background) or 'dark_background'
        style='dark_background',
        # limit color scale at 99.2% quantile of cell abundance
        max_color_quantile=0.992,
        # size of locations (adjust depending on figure size)
        circle_diameter=9,
        colorbar_position='right'
    ), plt.savefig("DH1_4yr_Exc_L4_6_clusters.pdf")
	
clust_labels = ['means_per_cluster_mu_fg_Exc L5-6 FEZF2 ABO', 'means_per_cluster_mu_fg_Exc L5-6 THEMIS C1QL3', 'means_per_cluster_mu_fg_Exc L5-6 THEMIS FGF10'
                        ]
clust_col = ['' + str(i) for i in clust_labels] # in case column names differ from labels

slide = select_slide(adata_vis, 'DH1')

with mpl.rc_context({'figure.figsize': (20, 20)}):
    fig = plot_spatial(
        adata=slide,
        # labels to show on a plot
        color=clust_col, labels=clust_labels,
		alpha_scaling=1.0,
        show_img=True,
		img_alpha=0.0,
		adjust_text=True,
        # 'fast' (white background) or 'dark_background'
        style='dark_background',
        # limit color scale at 99.2% quantile of cell abundance
        max_color_quantile=0.992,
        # size of locations (adjust depending on figure size)
        circle_diameter=9,
        colorbar_position='right'
    ), plt.savefig("DH1_4yr_Exc_L5_6_clusters.pdf")
	
clust_labels = ['means_per_cluster_mu_fg_Inh L1 SST NMBR', 'means_per_cluster_mu_fg_Inh L1-2 GAD1 MC4R', 'means_per_cluster_mu_fg_Inh L1-2 PAX6 CDH12', 
				  'means_per_cluster_mu_fg_Inh L1-2 SST BAGE2', 'means_per_cluster_mu_fg_Inh L1-2 VIP LBH', 'means_per_cluster_mu_fg_Inh L1-2 VIP PCDH20' 
                        ]
clust_col = ['' + str(i) for i in clust_labels] # in case column names differ from labels

slide = select_slide(adata_vis, 'DH1')

with mpl.rc_context({'figure.figsize': (20, 20)}):
    fig = plot_spatial(
        adata=slide,
        # labels to show on a plot
        color=clust_col, labels=clust_labels,
		alpha_scaling=1.0,
        show_img=True,
		img_alpha=0.0,
		adjust_text=True,
        # 'fast' (white background) or 'dark_background'
        style='dark_background',
        # limit color scale at 99.2% quantile of cell abundance
        max_color_quantile=0.992,
        # size of locations (adjust depending on figure size)
        circle_diameter=9,
        colorbar_position='right'
    ), plt.savefig("DH1_4yr_Inh_L1_2_clusters.pdf")
	
clust_labels = ['means_per_cluster_mu_fg_Inh L1-2 VIP TSPAN12', 'means_per_cluster_mu_fg_Inh L1-3 SST CALB1', 'means_per_cluster_mu_fg_Inh L1-3 VIP ADAMTSL1', 'means_per_cluster_mu_fg_Inh L1-3 VIP CCDC184', 'means_per_cluster_mu_fg_Inh L1-3 VIP CHRM2', 
				  'means_per_cluster_mu_fg_Inh L1-3 VIP GGH'
                        ]
clust_col = ['' + str(i) for i in clust_labels] # in case column names differ from labels

slide = select_slide(adata_vis, 'DH1')

with mpl.rc_context({'figure.figsize': (20, 20)}):
    fig = plot_spatial(
        adata=slide,
        # labels to show on a plot
        color=clust_col, labels=clust_labels,
		alpha_scaling=1.0,
        show_img=True,
		img_alpha=0.0,
		adjust_text=True,
        # 'fast' (white background) or 'dark_background'
        style='dark_background',
        # limit color scale at 99.2% quantile of cell abundance
        max_color_quantile=0.992,
        # size of locations (adjust depending on figure size)
        circle_diameter=9,
        colorbar_position='right'
    ), plt.savefig("DH1_4yr_Inh_L1_3_clusters.pdf")
	
clust_labels = ['means_per_cluster_mu_fg_Inh L1-4 LAMP5 LCP2', 'means_per_cluster_mu_fg_Inh L1-4 VIP CHRNA6', 'means_per_cluster_mu_fg_Inh L1-4 VIP OPRM1'
                        ]
clust_col = ['' + str(i) for i in clust_labels] # in case column names differ from labels

slide = select_slide(adata_vis, 'DH1')

with mpl.rc_context({'figure.figsize': (20, 20)}):
    fig = plot_spatial(
        adata=slide,
        # labels to show on a plot
        color=clust_col, labels=clust_labels,
		alpha_scaling=1.0,
        show_img=True,
		img_alpha=0.0,
		adjust_text=True,
        # 'fast' (white background) or 'dark_background'
        style='dark_background',
        # limit color scale at 99.2% quantile of cell abundance
        max_color_quantile=0.992,
        # size of locations (adjust depending on figure size)
        circle_diameter=9,
        colorbar_position='right'
    ), plt.savefig("DH1_4yr_Inh_L1_4_clusters.pdf")
	
clust_labels = ['means_per_cluster_mu_fg_Inh L2-4 PVALB WFDC2', 'means_per_cluster_mu_fg_Inh L2-4 SST FRZB', 'means_per_cluster_mu_fg_Inh L2-4 VIP CBLN1', 'means_per_cluster_mu_fg_Inh L2-4 VIP SPAG17'
                        ]
clust_col = ['' + str(i) for i in clust_labels] # in case column names differ from labels

slide = select_slide(adata_vis, 'DH1')

with mpl.rc_context({'figure.figsize': (20, 20)}):
    fig = plot_spatial(
        adata=slide,
        # labels to show on a plot
        color=clust_col, labels=clust_labels,
		alpha_scaling=1.0,
        show_img=True,
		img_alpha=0.0,
		adjust_text=True,
        # 'fast' (white background) or 'dark_background'
        style='dark_background',
        # limit color scale at 99.2% quantile of cell abundance
        max_color_quantile=0.992,
        # size of locations (adjust depending on figure size)
        circle_diameter=9,
        colorbar_position='right'
    ), plt.savefig("DH1_4yr_Inh_L2_4_clusters.pdf")
	
clust_labels = ['means_per_cluster_mu_fg_Inh L2-5 PVALB SCUBE3', 'means_per_cluster_mu_fg_Inh L2-5 VIP SERPINF1', 'means_per_cluster_mu_fg_Inh L2-5 VIP TYR', 'means_per_cluster_mu_fg_Inh L2-6 LAMP5 CA1', 
				  'means_per_cluster_mu_fg_Inh L2-6 VIP QPCT', 'means_per_cluster_mu_fg_Inh L3-5 SST ADGRG6', 'means_per_cluster_mu_fg_Inh L3-6 VIP HS3ST3A1'
                        ]
clust_col = ['' + str(i) for i in clust_labels] # in case column names differ from labels

slide = select_slide(adata_vis, 'DH1')

with mpl.rc_context({'figure.figsize': (20, 20)}):
    fig = plot_spatial(
        adata=slide,
        # labels to show on a plot
        color=clust_col, labels=clust_labels,
		alpha_scaling=1.0,
        show_img=True,
		img_alpha=0.0,
		adjust_text=True,
        # 'fast' (white background) or 'dark_background'
        style='dark_background',
        # limit color scale at 99.2% quantile of cell abundance
        max_color_quantile=0.992,
        # size of locations (adjust depending on figure size)
        circle_diameter=9,
        colorbar_position='right'
    ), plt.savefig("DH1_4yr_Inh_L2_6_clusters.pdf")
	

	
clust_labels = ['means_per_cluster_mu_fg_Inh L4-5 SST STK32A', 'means_per_cluster_mu_fg_Inh L4-6 PVALB SULF1', 'means_per_cluster_mu_fg_Inh L4-6 SST B3GAT2', 'means_per_cluster_mu_fg_Inh L4-6 SST GXYLT2', 
				  'means_per_cluster_mu_fg_Inh L5-6 PVALB LGR5', 'means_per_cluster_mu_fg_Inh L5-6 SST MIR548F2', 'means_per_cluster_mu_fg_Inh L5-6 SST NPM1P10'
                        ]
clust_col = ['' + str(i) for i in clust_labels] # in case column names differ from labels

slide = select_slide(adata_vis, 'DH1')

with mpl.rc_context({'figure.figsize': (20, 20)}):
    fig = plot_spatial(
        adata=slide,
        # labels to show on a plot
        color=clust_col, labels=clust_labels,
		alpha_scaling=1.0,
        show_img=True,
		img_alpha=0.0,
		adjust_text=True,
        # 'fast' (white background) or 'dark_background'
        style='dark_background',
        # limit color scale at 99.2% quantile of cell abundance
        max_color_quantile=0.992,
        # size of locations (adjust depending on figure size)
        circle_diameter=9,
        colorbar_position='right'
    ), plt.savefig("DH1_4yr_Inh_L5_6_clusters.pdf")
	
	
clust_labels = ['means_per_cluster_mu_fg_Micro L1-3 TYROBP', 'means_per_cluster_mu_fg_Endo L2-6 NOSTRIN'
               ]
clust_col = ['' + str(i) for i in clust_labels] # in case column names differ from labels

slide = select_slide(adata_vis, 'DH1')

with mpl.rc_context({'figure.figsize': (20, 20)}):
    fig = plot_spatial(
        adata=slide,
        # labels to show on a plot
        color=clust_col, labels=clust_labels,
		alpha_scaling=1.0,
        show_img=True,
		img_alpha=0.0,
		adjust_text=True,
        # 'fast' (white background) or 'dark_background'
        style='dark_background',
        # limit color scale at 99.2% quantile of cell abundance
        max_color_quantile=0.992,
        # size of locations (adjust depending on figure size)
        circle_diameter=9,
        colorbar_position='right'
    ), plt.savefig("DH1_4yr_Endo_clusters.pdf")
	
clust_labels = ['means_per_cluster_mu_fg_OPC L1-6 PDGFRA', 'means_per_cluster_mu_fg_Oligo L1-6 OPALIN', 'means_per_cluster_mu_fg_Inh L3-6 SST NPY'	
               ]
clust_col = ['' + str(i) for i in clust_labels] # in case column names differ from labels

slide = select_slide(adata_vis, 'DH1')

with mpl.rc_context({'figure.figsize': (20, 20)}):
    fig = plot_spatial(
        adata=slide,
        # labels to show on a plot
        color=clust_col, labels=clust_labels,
		alpha_scaling=1.0,
        show_img=True,
		img_alpha=0.0,
		adjust_text=True,
        # 'fast' (white background) or 'dark_background'
        style='dark_background',
        # limit color scale at 99.2% quantile of cell abundance
        max_color_quantile=0.992,
        # size of locations (adjust depending on figure size)
        circle_diameter=9,
        colorbar_position='right'
    ), plt.savefig("DH1_4yr_Oligo_clusters.pdf")

clust_labels = ['means_per_cluster_mu_fg_Micro L1-3 TYROBP', 'means_per_cluster_mu_fg_Endo L2-6 NOSTRIN'
               ] 			   
clust_col = ['' + str(i) for i in clust_labels] # in case column names differ from labels

slide = select_slide(adata_vis, 'DH1')

with mpl.rc_context({'figure.figsize': (20, 20)}):
    fig = plot_spatial(
        adata=slide,
        # labels to show on a plot
        color=clust_col, labels=clust_labels,
		alpha_scaling=1.0,
        show_img=True,
		img_alpha=0.0,
		adjust_text=True,
        # 'fast' (white background) or 'dark_background'
        style='dark_background',
        # limit color scale at 99.2% quantile of cell abundance
        max_color_quantile=0.992,
        # size of locations (adjust depending on figure size)
        circle_diameter=9,
        colorbar_position='right'
    ), plt.savefig("DH1_4yr_Endo_clusters.pdf")
	
clust_labels = ['means_per_cluster_mu_fg_Astro L1-2 FGFR3 GFAP', 'means_per_cluster_mu_fg_Astro L1-6 FGFR3 SLC14A1'
                        ]
clust_col = ['' + str(i) for i in clust_labels] # in case column names differ from labels

slide = select_slide(adata_vis, 'DH2')

with mpl.rc_context({'figure.figsize': (20, 20)}):
    fig = plot_spatial(
        adata=slide,
        # labels to show on a plot
        color=clust_col, labels=clust_labels,
		alpha_scaling=1.0,
        show_img=True,
		img_alpha=0.0,
		adjust_text=True,
        # 'fast' (white background) or 'dark_background'
        style='dark_background',
        # limit color scale at 99.2% quantile of cell abundance
        max_color_quantile=0.992,
        # size of locations (adjust depending on figure size)
        circle_diameter=9,
        colorbar_position='right'
    ), plt.savefig("DH2_4yr_Astro_clusters.pdf")

clust_labels = ['means_per_cluster_mu_fg_Exc L2 LAMP5 LTK', 'means_per_cluster_mu_fg_Exc L2-3 LINC00507 FREM3', 'means_per_cluster_mu_fg_Exc L3-5 RORB COL22A1',
                        'means_per_cluster_mu_fg_Exc L3-5 RORB ESR1', 'means_per_cluster_mu_fg_Exc L3-5 RORB FILIP1L', 'means_per_cluster_mu_fg_Exc L3-5 RORB TWIST2'
                        ]
clust_col = ['' + str(i) for i in clust_labels] # in case column names differ from labels

slide = select_slide(adata_vis, 'DH2')

with mpl.rc_context({'figure.figsize': (20, 20)}):
    fig = plot_spatial(
        adata=slide,
        # labels to show on a plot
        color=clust_col, labels=clust_labels,
		alpha_scaling=1.0,
        show_img=True,
		img_alpha=0.0,
		adjust_text=True,
        # 'fast' (white background) or 'dark_background'
        style='dark_background',
        # limit color scale at 99.2% quantile of cell abundance
        max_color_quantile=0.992,
        # size of locations (adjust depending on figure size)
        circle_diameter=9,
        colorbar_position='right'
    ), plt.savefig("DH2_4yr_Exc_L3_5_clusters.pdf")
	
clust_labels = ['means_per_cluster_mu_fg_Exc L4-5 RORB DAPK2', 'means_per_cluster_mu_fg_Exc L4-5 RORB FOLH1B', 'means_per_cluster_mu_fg_Exc L4-6 FEZF2 IL26', 'means_per_cluster_mu_fg_Exc L4-6 RORB SEMA3E'
                        ]
clust_col = ['' + str(i) for i in clust_labels] # in case column names differ from labels

slide = select_slide(adata_vis, 'DH2')

with mpl.rc_context({'figure.figsize': (20, 20)}):
    fig = plot_spatial(
        adata=slide,
        # labels to show on a plot
        color=clust_col, labels=clust_labels,
		alpha_scaling=1.0,
        show_img=True,
		img_alpha=0.0,
		adjust_text=True,
        # 'fast' (white background) or 'dark_background'
        style='dark_background',
        # limit color scale at 99.2% quantile of cell abundance
        max_color_quantile=0.992,
        # size of locations (adjust depending on figure size)
        circle_diameter=9,
        colorbar_position='right'
    ), plt.savefig("DH2_4yr_Exc_L4_6_clusters.pdf")
	
clust_labels = ['means_per_cluster_mu_fg_Exc L5-6 FEZF2 ABO', 'means_per_cluster_mu_fg_Exc L5-6 THEMIS C1QL3', 'means_per_cluster_mu_fg_Exc L5-6 THEMIS FGF10'
                        ]
clust_col = ['' + str(i) for i in clust_labels] # in case column names differ from labels

slide = select_slide(adata_vis, 'DH2')

with mpl.rc_context({'figure.figsize': (20, 20)}):
    fig = plot_spatial(
        adata=slide,
        # labels to show on a plot
        color=clust_col, labels=clust_labels,
		alpha_scaling=1.0,
        show_img=True,
		img_alpha=0.0,
		adjust_text=True,
        # 'fast' (white background) or 'dark_background'
        style='dark_background',
        # limit color scale at 99.2% quantile of cell abundance
        max_color_quantile=0.992,
        # size of locations (adjust depending on figure size)
        circle_diameter=9,
        colorbar_position='right'
    ), plt.savefig("DH2_4yr_Exc_L5_6_clusters.pdf")
	
clust_labels = ['means_per_cluster_mu_fg_Inh L1 SST NMBR', 'means_per_cluster_mu_fg_Inh L1-2 GAD1 MC4R', 'means_per_cluster_mu_fg_Inh L1-2 PAX6 CDH12', 
				  'means_per_cluster_mu_fg_Inh L1-2 SST BAGE2', 'means_per_cluster_mu_fg_Inh L1-2 VIP LBH', 'means_per_cluster_mu_fg_Inh L1-2 VIP PCDH20' 
                        ]
clust_col = ['' + str(i) for i in clust_labels] # in case column names differ from labels

slide = select_slide(adata_vis, 'DH2')

with mpl.rc_context({'figure.figsize': (20, 20)}):
    fig = plot_spatial(
        adata=slide,
        # labels to show on a plot
        color=clust_col, labels=clust_labels,
		alpha_scaling=1.0,
        show_img=True,
		img_alpha=0.0,
		adjust_text=True,
        # 'fast' (white background) or 'dark_background'
        style='dark_background',
        # limit color scale at 99.2% quantile of cell abundance
        max_color_quantile=0.992,
        # size of locations (adjust depending on figure size)
        circle_diameter=9,
        colorbar_position='right'
    ), plt.savefig("DH2_4yr_Inh_L1_2_clusters.pdf")
	
clust_labels = ['means_per_cluster_mu_fg_Inh L1-2 VIP TSPAN12', 'means_per_cluster_mu_fg_Inh L1-3 SST CALB1', 'means_per_cluster_mu_fg_Inh L1-3 VIP ADAMTSL1', 'means_per_cluster_mu_fg_Inh L1-3 VIP CCDC184', 'means_per_cluster_mu_fg_Inh L1-3 VIP CHRM2', 
				  'means_per_cluster_mu_fg_Inh L1-3 VIP GGH'
                        ]
clust_col = ['' + str(i) for i in clust_labels] # in case column names differ from labels

slide = select_slide(adata_vis, 'DH2')

with mpl.rc_context({'figure.figsize': (20, 20)}):
    fig = plot_spatial(
        adata=slide,
        # labels to show on a plot
        color=clust_col, labels=clust_labels,
		alpha_scaling=1.0,
        show_img=True,
		img_alpha=0.0,
		adjust_text=True,
        # 'fast' (white background) or 'dark_background'
        style='dark_background',
        # limit color scale at 99.2% quantile of cell abundance
        max_color_quantile=0.992,
        # size of locations (adjust depending on figure size)
        circle_diameter=9,
        colorbar_position='right'
    ), plt.savefig("DH2_4yr_Inh_L1_3_clusters.pdf")
	
clust_labels = ['means_per_cluster_mu_fg_Inh L1-4 LAMP5 LCP2', 'means_per_cluster_mu_fg_Inh L1-4 VIP CHRNA6', 'means_per_cluster_mu_fg_Inh L1-4 VIP OPRM1'
                        ]
clust_col = ['' + str(i) for i in clust_labels] # in case column names differ from labels

slide = select_slide(adata_vis, 'DH2')

with mpl.rc_context({'figure.figsize': (20, 20)}):
    fig = plot_spatial(
        adata=slide,
        # labels to show on a plot
        color=clust_col, labels=clust_labels,
		alpha_scaling=1.0,
        show_img=True,
		img_alpha=0.0,
		adjust_text=True,
        # 'fast' (white background) or 'dark_background'
        style='dark_background',
        # limit color scale at 99.2% quantile of cell abundance
        max_color_quantile=0.992,
        # size of locations (adjust depending on figure size)
        circle_diameter=9,
        colorbar_position='right'
    ), plt.savefig("DH2_4yr_Inh_L1_4_clusters.pdf")
	
clust_labels = ['means_per_cluster_mu_fg_Inh L2-4 PVALB WFDC2', 'means_per_cluster_mu_fg_Inh L2-4 SST FRZB', 'means_per_cluster_mu_fg_Inh L2-4 VIP CBLN1', 'means_per_cluster_mu_fg_Inh L2-4 VIP SPAG17'
                        ]
clust_col = ['' + str(i) for i in clust_labels] # in case column names differ from labels

slide = select_slide(adata_vis, 'DH2')

with mpl.rc_context({'figure.figsize': (20, 20)}):
    fig = plot_spatial(
        adata=slide,
        # labels to show on a plot
        color=clust_col, labels=clust_labels,
		alpha_scaling=1.0,
        show_img=True,
		img_alpha=0.0,
		adjust_text=True,
        # 'fast' (white background) or 'dark_background'
        style='dark_background',
        # limit color scale at 99.2% quantile of cell abundance
        max_color_quantile=0.992,
        # size of locations (adjust depending on figure size)
        circle_diameter=9,
        colorbar_position='right'
    ), plt.savefig("DH2_4yr_Inh_L2_4_clusters.pdf")
	
clust_labels = ['means_per_cluster_mu_fg_Inh L2-5 PVALB SCUBE3', 'means_per_cluster_mu_fg_Inh L2-5 VIP SERPINF1', 'means_per_cluster_mu_fg_Inh L2-5 VIP TYR', 'means_per_cluster_mu_fg_Inh L2-6 LAMP5 CA1', 
				  'means_per_cluster_mu_fg_Inh L2-6 VIP QPCT', 'means_per_cluster_mu_fg_Inh L3-5 SST ADGRG6', 'means_per_cluster_mu_fg_Inh L3-6 VIP HS3ST3A1'
                        ]
clust_col = ['' + str(i) for i in clust_labels] # in case column names differ from labels

slide = select_slide(adata_vis, 'DH2')

with mpl.rc_context({'figure.figsize': (20, 20)}):
    fig = plot_spatial(
        adata=slide,
        # labels to show on a plot
        color=clust_col, labels=clust_labels,
		alpha_scaling=1.0,
        show_img=True,
		img_alpha=0.0,
		adjust_text=True,
        # 'fast' (white background) or 'dark_background'
        style='dark_background',
        # limit color scale at 99.2% quantile of cell abundance
        max_color_quantile=0.992,
        # size of locations (adjust depending on figure size)
        circle_diameter=9,
        colorbar_position='right'
    ), plt.savefig("DH2_4yr_Inh_L2_6_clusters.pdf")
	
	
clust_labels = ['means_per_cluster_mu_fg_Inh L4-5 SST STK32A', 'means_per_cluster_mu_fg_Inh L4-6 PVALB SULF1', 'means_per_cluster_mu_fg_Inh L4-6 SST B3GAT2', 'means_per_cluster_mu_fg_Inh L4-6 SST GXYLT2', 
				  'means_per_cluster_mu_fg_Inh L5-6 PVALB LGR5', 'means_per_cluster_mu_fg_Inh L5-6 SST MIR548F2', 'means_per_cluster_mu_fg_Inh L5-6 SST NPM1P10'
                        ]
clust_col = ['' + str(i) for i in clust_labels] # in case column names differ from labels

slide = select_slide(adata_vis, 'DH2')

with mpl.rc_context({'figure.figsize': (20, 20)}):
    fig = plot_spatial(
        adata=slide,
        # labels to show on a plot
        color=clust_col, labels=clust_labels,
		alpha_scaling=1.0,
        show_img=True,
		img_alpha=0.0,
		adjust_text=True,
        # 'fast' (white background) or 'dark_background'
        style='dark_background',
        # limit color scale at 99.2% quantile of cell abundance
        max_color_quantile=0.992,
        # size of locations (adjust depending on figure size)
        circle_diameter=9,
        colorbar_position='right'
    ), plt.savefig("DH2_4yr_Inh_L5_6_clusters.pdf")
	
	
clust_labels = ['means_per_cluster_mu_fg_Micro L1-3 TYROBP', 'means_per_cluster_mu_fg_Endo L2-6 NOSTRIN'
               ]
clust_col = ['' + str(i) for i in clust_labels] # in case column names differ from labels

slide = select_slide(adata_vis, 'DH2')

with mpl.rc_context({'figure.figsize': (20, 20)}):
    fig = plot_spatial(
        adata=slide,
        # labels to show on a plot
        color=clust_col, labels=clust_labels,
		alpha_scaling=1.0,
        show_img=True,
		img_alpha=0.0,
		adjust_text=True,
        # 'fast' (white background) or 'dark_background'
        style='dark_background',
        # limit color scale at 99.2% quantile of cell abundance
        max_color_quantile=0.992,
        # size of locations (adjust depending on figure size)
        circle_diameter=9,
        colorbar_position='right'
    ), plt.savefig("DH2_4yr_Endo_clusters.pdf")
	
clust_labels = ['means_per_cluster_mu_fg_OPC L1-6 PDGFRA', 'means_per_cluster_mu_fg_Oligo L1-6 OPALIN', 'means_per_cluster_mu_fg_Inh L3-6 SST NPY'	
               ]
clust_col = ['' + str(i) for i in clust_labels] # in case column names differ from labels

slide = select_slide(adata_vis, 'DH2')

with mpl.rc_context({'figure.figsize': (20, 20)}):
    fig = plot_spatial(
        adata=slide,
        # labels to show on a plot
        color=clust_col, labels=clust_labels,
		alpha_scaling=1.0,
        show_img=True,
		img_alpha=0.0,
		adjust_text=True,
        # 'fast' (white background) or 'dark_background'
        style='dark_background',
        # limit color scale at 99.2% quantile of cell abundance
        max_color_quantile=0.992,
        # size of locations (adjust depending on figure size)
        circle_diameter=9,
        colorbar_position='right'
    ), plt.savefig("DH2_4yr_Oligo_clusters.pdf")	
	
clust_labels = ['means_per_cluster_mu_fg_Astro L1-2 FGFR3 GFAP', 'means_per_cluster_mu_fg_Astro L1-6 FGFR3 SLC14A1'
                        ]
clust_col = ['' + str(i) for i in clust_labels] # in case column names differ from labels

slide = select_slide(adata_vis, 'DH3')

with mpl.rc_context({'figure.figsize': (20, 20)}):
    fig = plot_spatial(
        adata=slide,
        # labels to show on a plot
        color=clust_col, labels=clust_labels,
		alpha_scaling=1.0,
        show_img=True,
		img_alpha=0.0,
		adjust_text=True,
        # 'fast' (white background) or 'dark_background'
        style='dark_background',
        # limit color scale at 99.2% quantile of cell abundance
        max_color_quantile=0.992,
        # size of locations (adjust depending on figure size)
        circle_diameter=9,
        colorbar_position='right'
    ), plt.savefig("DH3_2nd_15yr_Astro_clusters.pdf")

clust_labels = ['means_per_cluster_mu_fg_Exc L2 LAMP5 LTK', 'means_per_cluster_mu_fg_Exc L2-3 LINC00507 FREM3', 'means_per_cluster_mu_fg_Exc L3-5 RORB COL22A1',
                        'means_per_cluster_mu_fg_Exc L3-5 RORB ESR1', 'means_per_cluster_mu_fg_Exc L3-5 RORB FILIP1L', 'means_per_cluster_mu_fg_Exc L3-5 RORB TWIST2'
                        ]
clust_col = ['' + str(i) for i in clust_labels] # in case column names differ from labels

slide = select_slide(adata_vis, 'DH3')

with mpl.rc_context({'figure.figsize': (20, 20)}):
    fig = plot_spatial(
        adata=slide,
        # labels to show on a plot
        color=clust_col, labels=clust_labels,
		alpha_scaling=1.0,
        show_img=True,
		img_alpha=0.0,
		adjust_text=True,
        # 'fast' (white background) or 'dark_background'
        style='dark_background',
        # limit color scale at 99.2% quantile of cell abundance
        max_color_quantile=0.992,
        # size of locations (adjust depending on figure size)
        circle_diameter=9,
        colorbar_position='right'
    ), plt.savefig("DH3_2nd_15yr_Exc_L3_5_clusters.pdf")
	
clust_labels = ['means_per_cluster_mu_fg_Exc L4-5 RORB DAPK2', 'means_per_cluster_mu_fg_Exc L4-5 RORB FOLH1B', 'means_per_cluster_mu_fg_Exc L4-6 FEZF2 IL26', 'means_per_cluster_mu_fg_Exc L4-6 RORB SEMA3E'
                        ]
clust_col = ['' + str(i) for i in clust_labels] # in case column names differ from labels

slide = select_slide(adata_vis, 'DH3')

with mpl.rc_context({'figure.figsize': (20, 20)}):
    fig = plot_spatial(
        adata=slide,
        # labels to show on a plot
        color=clust_col, labels=clust_labels,
		alpha_scaling=1.0,
        show_img=True,
		img_alpha=0.0,
		adjust_text=True,
        # 'fast' (white background) or 'dark_background'
        style='dark_background',
        # limit color scale at 99.2% quantile of cell abundance
        max_color_quantile=0.992,
        # size of locations (adjust depending on figure size)
        circle_diameter=9,
        colorbar_position='right'
    ), plt.savefig("DH3_2nd_15yr_Exc_L4_6_clusters.pdf")
	
clust_labels = ['means_per_cluster_mu_fg_Exc L5-6 FEZF2 ABO', 'means_per_cluster_mu_fg_Exc L5-6 THEMIS C1QL3', 'means_per_cluster_mu_fg_Exc L5-6 THEMIS FGF10'
                        ]
clust_col = ['' + str(i) for i in clust_labels] # in case column names differ from labels

slide = select_slide(adata_vis, 'DH3')

with mpl.rc_context({'figure.figsize': (20, 20)}):
    fig = plot_spatial(
        adata=slide,
        # labels to show on a plot
        color=clust_col, labels=clust_labels,
		alpha_scaling=1.0,
        show_img=True,
		img_alpha=0.0,
		adjust_text=True,
        # 'fast' (white background) or 'dark_background'
        style='dark_background',
        # limit color scale at 99.2% quantile of cell abundance
        max_color_quantile=0.992,
        # size of locations (adjust depending on figure size)
        circle_diameter=9,
        colorbar_position='right'
    ), plt.savefig("DH3_2nd_15yr_Exc_L5_6_clusters.pdf")
	
clust_labels = ['means_per_cluster_mu_fg_Inh L1 SST NMBR', 'means_per_cluster_mu_fg_Inh L1-2 GAD1 MC4R', 'means_per_cluster_mu_fg_Inh L1-2 PAX6 CDH12', 
				  'means_per_cluster_mu_fg_Inh L1-2 SST BAGE2', 'means_per_cluster_mu_fg_Inh L1-2 VIP LBH', 'means_per_cluster_mu_fg_Inh L1-2 VIP PCDH20' 
                        ]
clust_col = ['' + str(i) for i in clust_labels] # in case column names differ from labels

slide = select_slide(adata_vis, 'DH3')

with mpl.rc_context({'figure.figsize': (20, 20)}):
    fig = plot_spatial(
        adata=slide,
        # labels to show on a plot
        color=clust_col, labels=clust_labels,
		alpha_scaling=1.0,
        show_img=True,
		img_alpha=0.0,
		adjust_text=True,
        # 'fast' (white background) or 'dark_background'
        style='dark_background',
        # limit color scale at 99.2% quantile of cell abundance
        max_color_quantile=0.992,
        # size of locations (adjust depending on figure size)
        circle_diameter=9,
        colorbar_position='right'
    ), plt.savefig("DH3_2nd_15yr_Inh_L1_2_clusters.pdf")
	
clust_labels = ['means_per_cluster_mu_fg_Inh L1-2 VIP TSPAN12', 'means_per_cluster_mu_fg_Inh L1-3 SST CALB1', 'means_per_cluster_mu_fg_Inh L1-3 VIP ADAMTSL1', 'means_per_cluster_mu_fg_Inh L1-3 VIP CCDC184', 'means_per_cluster_mu_fg_Inh L1-3 VIP CHRM2', 
				  'means_per_cluster_mu_fg_Inh L1-3 VIP GGH'
                        ]
clust_col = ['' + str(i) for i in clust_labels] # in case column names differ from labels

slide = select_slide(adata_vis, 'DH3')

with mpl.rc_context({'figure.figsize': (20, 20)}):
    fig = plot_spatial(
        adata=slide,
        # labels to show on a plot
        color=clust_col, labels=clust_labels,
		alpha_scaling=1.0,
        show_img=True,
		img_alpha=0.0,
		adjust_text=True,
        # 'fast' (white background) or 'dark_background'
        style='dark_background',
        # limit color scale at 99.2% quantile of cell abundance
        max_color_quantile=0.992,
        # size of locations (adjust depending on figure size)
        circle_diameter=9,
        colorbar_position='right'
    ), plt.savefig("DH3_2nd_15yr_Inh_L1_3_clusters.pdf")
	
clust_labels = ['means_per_cluster_mu_fg_Inh L1-4 LAMP5 LCP2', 'means_per_cluster_mu_fg_Inh L1-4 VIP CHRNA6', 'means_per_cluster_mu_fg_Inh L1-4 VIP OPRM1'
                        ]
clust_col = ['' + str(i) for i in clust_labels] # in case column names differ from labels

slide = select_slide(adata_vis, 'DH3')

with mpl.rc_context({'figure.figsize': (20, 20)}):
    fig = plot_spatial(
        adata=slide,
        # labels to show on a plot
        color=clust_col, labels=clust_labels,
		alpha_scaling=1.0,
        show_img=True,
		img_alpha=0.0,
		adjust_text=True,
        # 'fast' (white background) or 'dark_background'
        style='dark_background',
        # limit color scale at 99.2% quantile of cell abundance
        max_color_quantile=0.992,
        # size of locations (adjust depending on figure size)
        circle_diameter=9,
        colorbar_position='right'
    ), plt.savefig("DH3_2nd_15yr_Inh_L1_4_clusters.pdf")
	
clust_labels = ['means_per_cluster_mu_fg_Inh L2-4 PVALB WFDC2', 'means_per_cluster_mu_fg_Inh L2-4 SST FRZB', 'means_per_cluster_mu_fg_Inh L2-4 VIP CBLN1', 'means_per_cluster_mu_fg_Inh L2-4 VIP SPAG17'
                        ]
clust_col = ['' + str(i) for i in clust_labels] # in case column names differ from labels

slide = select_slide(adata_vis, 'DH3')

with mpl.rc_context({'figure.figsize': (20, 20)}):
    fig = plot_spatial(
        adata=slide,
        # labels to show on a plot
        color=clust_col, labels=clust_labels,
		alpha_scaling=1.0,
        show_img=True,
		img_alpha=0.0,
		adjust_text=True,
        # 'fast' (white background) or 'dark_background'
        style='dark_background',
        # limit color scale at 99.2% quantile of cell abundance
        max_color_quantile=0.992,
        # size of locations (adjust depending on figure size)
        circle_diameter=9,
        colorbar_position='right'
    ), plt.savefig("DH3_2nd_15yr_Inh_L2_4_clusters.pdf")
	
clust_labels = ['means_per_cluster_mu_fg_Inh L2-5 PVALB SCUBE3', 'means_per_cluster_mu_fg_Inh L2-5 VIP SERPINF1', 'means_per_cluster_mu_fg_Inh L2-5 VIP TYR', 'means_per_cluster_mu_fg_Inh L2-6 LAMP5 CA1', 
				  'means_per_cluster_mu_fg_Inh L2-6 VIP QPCT', 'means_per_cluster_mu_fg_Inh L3-5 SST ADGRG6', 'means_per_cluster_mu_fg_Inh L3-6 VIP HS3ST3A1'
                        ]
clust_col = ['' + str(i) for i in clust_labels] # in case column names differ from labels

slide = select_slide(adata_vis, 'DH3')

with mpl.rc_context({'figure.figsize': (20, 20)}):
    fig = plot_spatial(
        adata=slide,
        # labels to show on a plot
        color=clust_col, labels=clust_labels,
		alpha_scaling=1.0,
        show_img=True,
		img_alpha=0.0,
		adjust_text=True,
        # 'fast' (white background) or 'dark_background'
        style='dark_background',
        # limit color scale at 99.2% quantile of cell abundance
        max_color_quantile=0.992,
        # size of locations (adjust depending on figure size)
        circle_diameter=9,
        colorbar_position='right'
    ), plt.savefig("DH3_2nd_15yr_Inh_L2_6_clusters.pdf")
	

	
clust_labels = ['means_per_cluster_mu_fg_Inh L4-5 SST STK32A', 'means_per_cluster_mu_fg_Inh L4-6 PVALB SULF1', 'means_per_cluster_mu_fg_Inh L4-6 SST B3GAT2', 'means_per_cluster_mu_fg_Inh L4-6 SST GXYLT2', 
				  'means_per_cluster_mu_fg_Inh L5-6 PVALB LGR5', 'means_per_cluster_mu_fg_Inh L5-6 SST MIR548F2', 'means_per_cluster_mu_fg_Inh L5-6 SST NPM1P10'
                        ]
clust_col = ['' + str(i) for i in clust_labels] # in case column names differ from labels

slide = select_slide(adata_vis, 'DH3')

with mpl.rc_context({'figure.figsize': (20, 20)}):
    fig = plot_spatial(
        adata=slide,
        # labels to show on a plot
        color=clust_col, labels=clust_labels,
		alpha_scaling=1.0,
        show_img=True,
		img_alpha=0.0,
		adjust_text=True,
        # 'fast' (white background) or 'dark_background'
        style='dark_background',
        # limit color scale at 99.2% quantile of cell abundance
        max_color_quantile=0.992,
        # size of locations (adjust depending on figure size)
        circle_diameter=9,
        colorbar_position='right'
    ), plt.savefig("DH3_2nd_15yr_Inh_L5_6_clusters.pdf")
	
	
clust_labels = ['means_per_cluster_mu_fg_Micro L1-3 TYROBP', 'means_per_cluster_mu_fg_Endo L2-6 NOSTRIN'
               ]
clust_col = ['' + str(i) for i in clust_labels] # in case column names differ from labels

slide = select_slide(adata_vis, 'DH3')

with mpl.rc_context({'figure.figsize': (20, 20)}):
    fig = plot_spatial(
        adata=slide,
        # labels to show on a plot
        color=clust_col, labels=clust_labels,
		alpha_scaling=1.0,
        show_img=True,
		img_alpha=0.0,
		adjust_text=True,
        # 'fast' (white background) or 'dark_background'
        style='dark_background',
        # limit color scale at 99.2% quantile of cell abundance
        max_color_quantile=0.992,
        # size of locations (adjust depending on figure size)
        circle_diameter=9,
        colorbar_position='right'
    ), plt.savefig("DH3_2nd_15yr_Endo_clusters.pdf")
	
clust_labels = ['means_per_cluster_mu_fg_OPC L1-6 PDGFRA', 'means_per_cluster_mu_fg_Oligo L1-6 OPALIN', 'means_per_cluster_mu_fg_Inh L3-6 SST NPY'	
               ]
clust_col = ['' + str(i) for i in clust_labels] # in case column names differ from labels

slide = select_slide(adata_vis, 'DH3')

with mpl.rc_context({'figure.figsize': (20, 20)}):
    fig = plot_spatial(
        adata=slide,
        # labels to show on a plot
        color=clust_col, labels=clust_labels,
		alpha_scaling=1.0,
        show_img=True,
		img_alpha=0.0,
		adjust_text=True,
        # 'fast' (white background) or 'dark_background'
        style='dark_background',
        # limit color scale at 99.2% quantile of cell abundance
        max_color_quantile=0.992,
        # size of locations (adjust depending on figure size)
        circle_diameter=9,
        colorbar_position='right'
    ), plt.savefig("DH3_2nd_15yr_Oligo_clusters.pdf")
	
clust_labels = ['means_per_cluster_mu_fg_Astro L1-2 FGFR3 GFAP', 'means_per_cluster_mu_fg_Astro L1-6 FGFR3 SLC14A1'
                        ]
clust_col = ['' + str(i) for i in clust_labels] # in case column names differ from labels

slide = select_slide(adata_vis, 'DH4')

with mpl.rc_context({'figure.figsize': (20, 20)}):
    fig = plot_spatial(
        adata=slide,
        # labels to show on a plot
        color=clust_col, labels=clust_labels,
		alpha_scaling=1.0,
        show_img=True,
		img_alpha=0.0,
		adjust_text=True,
        # 'fast' (white background) or 'dark_background'
        style='dark_background',
        # limit color scale at 99.2% quantile of cell abundance
        max_color_quantile=0.992,
        # size of locations (adjust depending on figure size)
        circle_diameter=9,
        colorbar_position='right'
    ), plt.savefig("DH4_2nd_15yr_Astro_clusters.pdf")

clust_labels = ['means_per_cluster_mu_fg_Exc L2 LAMP5 LTK', 'means_per_cluster_mu_fg_Exc L2-3 LINC00507 FREM3', 'means_per_cluster_mu_fg_Exc L3-5 RORB COL22A1',
                        'means_per_cluster_mu_fg_Exc L3-5 RORB ESR1', 'means_per_cluster_mu_fg_Exc L3-5 RORB FILIP1L', 'means_per_cluster_mu_fg_Exc L3-5 RORB TWIST2'
                        ]
clust_col = ['' + str(i) for i in clust_labels] # in case column names differ from labels

slide = select_slide(adata_vis, 'DH4')

with mpl.rc_context({'figure.figsize': (20, 20)}):
    fig = plot_spatial(
        adata=slide,
        # labels to show on a plot
        color=clust_col, labels=clust_labels,
		alpha_scaling=1.0,
        show_img=True,
		img_alpha=0.0,
		adjust_text=True,
        # 'fast' (white background) or 'dark_background'
        style='dark_background',
        # limit color scale at 99.2% quantile of cell abundance
        max_color_quantile=0.992,
        # size of locations (adjust depending on figure size)
        circle_diameter=9,
        colorbar_position='right'
    ), plt.savefig("DH4_2nd_15yr_Exc_L3_5_clusters.pdf")
	
clust_labels = ['means_per_cluster_mu_fg_Exc L4-5 RORB DAPK2', 'means_per_cluster_mu_fg_Exc L4-5 RORB FOLH1B', 'means_per_cluster_mu_fg_Exc L4-6 FEZF2 IL26', 'means_per_cluster_mu_fg_Exc L4-6 RORB SEMA3E'
                        ]
clust_col = ['' + str(i) for i in clust_labels] # in case column names differ from labels

slide = select_slide(adata_vis, 'DH4')

with mpl.rc_context({'figure.figsize': (20, 20)}):
    fig = plot_spatial(
        adata=slide,
        # labels to show on a plot
        color=clust_col, labels=clust_labels,
		alpha_scaling=1.0,
        show_img=True,
		img_alpha=0.0,
		adjust_text=True,
        # 'fast' (white background) or 'dark_background'
        style='dark_background',
        # limit color scale at 99.2% quantile of cell abundance
        max_color_quantile=0.992,
        # size of locations (adjust depending on figure size)
        circle_diameter=9,
        colorbar_position='right'
    ), plt.savefig("DH4_2nd_15yr_Exc_L4_6_clusters.pdf")
	
clust_labels = ['means_per_cluster_mu_fg_Exc L5-6 FEZF2 ABO', 'means_per_cluster_mu_fg_Exc L5-6 THEMIS C1QL3', 'means_per_cluster_mu_fg_Exc L5-6 THEMIS FGF10'
                        ]
clust_col = ['' + str(i) for i in clust_labels] # in case column names differ from labels

slide = select_slide(adata_vis, 'DH4')

with mpl.rc_context({'figure.figsize': (20, 20)}):
    fig = plot_spatial(
        adata=slide,
        # labels to show on a plot
        color=clust_col, labels=clust_labels,
		alpha_scaling=1.0,
        show_img=True,
		img_alpha=0.0,
		adjust_text=True,
        # 'fast' (white background) or 'dark_background'
        style='dark_background',
        # limit color scale at 99.2% quantile of cell abundance
        max_color_quantile=0.992,
        # size of locations (adjust depending on figure size)
        circle_diameter=9,
        colorbar_position='right'
    ), plt.savefig("DH4_2nd_15yr_Exc_L5_6_clusters.pdf")
	
clust_labels = ['means_per_cluster_mu_fg_Inh L1 SST NMBR', 'means_per_cluster_mu_fg_Inh L1-2 GAD1 MC4R', 'means_per_cluster_mu_fg_Inh L1-2 PAX6 CDH12', 
				  'means_per_cluster_mu_fg_Inh L1-2 SST BAGE2', 'means_per_cluster_mu_fg_Inh L1-2 VIP LBH', 'means_per_cluster_mu_fg_Inh L1-2 VIP PCDH20' 
                        ]
clust_col = ['' + str(i) for i in clust_labels] # in case column names differ from labels

slide = select_slide(adata_vis, 'DH4')

with mpl.rc_context({'figure.figsize': (20, 20)}):
    fig = plot_spatial(
        adata=slide,
        # labels to show on a plot
        color=clust_col, labels=clust_labels,
		alpha_scaling=1.0,
        show_img=True,
		img_alpha=0.0,
		adjust_text=True,
        # 'fast' (white background) or 'dark_background'
        style='dark_background',
        # limit color scale at 99.2% quantile of cell abundance
        max_color_quantile=0.992,
        # size of locations (adjust depending on figure size)
        circle_diameter=9,
        colorbar_position='right'
    ), plt.savefig("DH4_2nd_15yr_Inh_L1_2_clusters.pdf")
	
clust_labels = ['means_per_cluster_mu_fg_Inh L1-2 VIP TSPAN12', 'means_per_cluster_mu_fg_Inh L1-3 SST CALB1', 'means_per_cluster_mu_fg_Inh L1-3 VIP ADAMTSL1', 'means_per_cluster_mu_fg_Inh L1-3 VIP CCDC184', 'means_per_cluster_mu_fg_Inh L1-3 VIP CHRM2', 
				  'means_per_cluster_mu_fg_Inh L1-3 VIP GGH'
                        ]
clust_col = ['' + str(i) for i in clust_labels] # in case column names differ from labels

slide = select_slide(adata_vis, 'DH4')

with mpl.rc_context({'figure.figsize': (20, 20)}):
    fig = plot_spatial(
        adata=slide,
        # labels to show on a plot
        color=clust_col, labels=clust_labels,
		alpha_scaling=1.0,
        show_img=True,
		img_alpha=0.0,
		adjust_text=True,
        # 'fast' (white background) or 'dark_background'
        style='dark_background',
        # limit color scale at 99.2% quantile of cell abundance
        max_color_quantile=0.992,
        # size of locations (adjust depending on figure size)
        circle_diameter=9,
        colorbar_position='right'
    ), plt.savefig("DH4_2nd_15yr_Inh_L1_3_clusters.pdf")
	
clust_labels = ['means_per_cluster_mu_fg_Inh L1-4 LAMP5 LCP2', 'means_per_cluster_mu_fg_Inh L1-4 VIP CHRNA6', 'means_per_cluster_mu_fg_Inh L1-4 VIP OPRM1'
                        ]
clust_col = ['' + str(i) for i in clust_labels] # in case column names differ from labels

slide = select_slide(adata_vis, 'DH4')

with mpl.rc_context({'figure.figsize': (20, 20)}):
    fig = plot_spatial(
        adata=slide,
        # labels to show on a plot
        color=clust_col, labels=clust_labels,
		alpha_scaling=1.0,
        show_img=True,
		img_alpha=0.0,
		adjust_text=True,
        # 'fast' (white background) or 'dark_background'
        style='dark_background',
        # limit color scale at 99.2% quantile of cell abundance
        max_color_quantile=0.992,
        # size of locations (adjust depending on figure size)
        circle_diameter=9,
        colorbar_position='right'
    ), plt.savefig("DH4_2nd_15yr_Inh_L1_4_clusters.pdf")
	
clust_labels = ['means_per_cluster_mu_fg_Inh L2-4 PVALB WFDC2', 'means_per_cluster_mu_fg_Inh L2-4 SST FRZB', 'means_per_cluster_mu_fg_Inh L2-4 VIP CBLN1', 'means_per_cluster_mu_fg_Inh L2-4 VIP SPAG17'
                        ]
clust_col = ['' + str(i) for i in clust_labels] # in case column names differ from labels

slide = select_slide(adata_vis, 'DH4')

with mpl.rc_context({'figure.figsize': (20, 20)}):
    fig = plot_spatial(
        adata=slide,
        # labels to show on a plot
        color=clust_col, labels=clust_labels,
		alpha_scaling=1.0,
        show_img=True,
		img_alpha=0.0,
		adjust_text=True,
        # 'fast' (white background) or 'dark_background'
        style='dark_background',
        # limit color scale at 99.2% quantile of cell abundance
        max_color_quantile=0.992,
        # size of locations (adjust depending on figure size)
        circle_diameter=9,
        colorbar_position='right'
    ), plt.savefig("DH4_2nd_15yr_Inh_L2_4_clusters.pdf")
	
clust_labels = ['means_per_cluster_mu_fg_Inh L2-5 PVALB SCUBE3', 'means_per_cluster_mu_fg_Inh L2-5 VIP SERPINF1', 'means_per_cluster_mu_fg_Inh L2-5 VIP TYR', 'means_per_cluster_mu_fg_Inh L2-6 LAMP5 CA1', 
				  'means_per_cluster_mu_fg_Inh L2-6 VIP QPCT', 'means_per_cluster_mu_fg_Inh L3-5 SST ADGRG6', 'means_per_cluster_mu_fg_Inh L3-6 VIP HS3ST3A1'
                        ]
clust_col = ['' + str(i) for i in clust_labels] # in case column names differ from labels

slide = select_slide(adata_vis, 'DH4')

with mpl.rc_context({'figure.figsize': (20, 20)}):
    fig = plot_spatial(
        adata=slide,
        # labels to show on a plot
        color=clust_col, labels=clust_labels,
		alpha_scaling=1.0,
        show_img=True,
		img_alpha=0.0,
		adjust_text=True,
        # 'fast' (white background) or 'dark_background'
        style='dark_background',
        # limit color scale at 99.2% quantile of cell abundance
        max_color_quantile=0.992,
        # size of locations (adjust depending on figure size)
        circle_diameter=9,
        colorbar_position='right'
    ), plt.savefig("DH4_2nd_15yr_Inh_L2_6_clusters.pdf")
	

	
clust_labels = ['means_per_cluster_mu_fg_Inh L4-5 SST STK32A', 'means_per_cluster_mu_fg_Inh L4-6 PVALB SULF1', 'means_per_cluster_mu_fg_Inh L4-6 SST B3GAT2', 'means_per_cluster_mu_fg_Inh L4-6 SST GXYLT2', 
				  'means_per_cluster_mu_fg_Inh L5-6 PVALB LGR5', 'means_per_cluster_mu_fg_Inh L5-6 SST MIR548F2', 'means_per_cluster_mu_fg_Inh L5-6 SST NPM1P10'
                        ]
clust_col = ['' + str(i) for i in clust_labels] # in case column names differ from labels

slide = select_slide(adata_vis, 'DH4')

with mpl.rc_context({'figure.figsize': (20, 20)}):
    fig = plot_spatial(
        adata=slide,
        # labels to show on a plot
        color=clust_col, labels=clust_labels,
		alpha_scaling=1.0,
        show_img=True,
		img_alpha=0.0,
		adjust_text=True,
        # 'fast' (white background) or 'dark_background'
        style='dark_background',
        # limit color scale at 99.2% quantile of cell abundance
        max_color_quantile=0.992,
        # size of locations (adjust depending on figure size)
        circle_diameter=9,
        colorbar_position='right'
    ), plt.savefig("DH4_2nd_15yr_Inh_L5_6_clusters.pdf")
	
	
clust_labels = ['means_per_cluster_mu_fg_Micro L1-3 TYROBP', 'means_per_cluster_mu_fg_Endo L2-6 NOSTRIN'
               ]
clust_col = ['' + str(i) for i in clust_labels] # in case column names differ from labels

slide = select_slide(adata_vis, 'DH4')

with mpl.rc_context({'figure.figsize': (20, 20)}):
    fig = plot_spatial(
        adata=slide,
        # labels to show on a plot
        color=clust_col, labels=clust_labels,
		alpha_scaling=1.0,
        show_img=True,
		img_alpha=0.0,
		adjust_text=True,
        # 'fast' (white background) or 'dark_background'
        style='dark_background',
        # limit color scale at 99.2% quantile of cell abundance
        max_color_quantile=0.992,
        # size of locations (adjust depending on figure size)
        circle_diameter=9,
        colorbar_position='right'
    ), plt.savefig("DH4_2nd_15yr_Endo_clusters.pdf")
	
clust_labels = ['means_per_cluster_mu_fg_OPC L1-6 PDGFRA', 'means_per_cluster_mu_fg_Oligo L1-6 OPALIN', 'means_per_cluster_mu_fg_Inh L3-6 SST NPY'	
               ]
clust_col = ['' + str(i) for i in clust_labels] # in case column names differ from labels

slide = select_slide(adata_vis, 'DH4')

with mpl.rc_context({'figure.figsize': (20, 20)}):
    fig = plot_spatial(
        adata=slide,
        # labels to show on a plot
        color=clust_col, labels=clust_labels,
		alpha_scaling=1.0,
        show_img=True,
		img_alpha=0.0,
		adjust_text=True,
        # 'fast' (white background) or 'dark_background'
        style='dark_background',
        # limit color scale at 99.2% quantile of cell abundance
        max_color_quantile=0.992,
        # size of locations (adjust depending on figure size)
        circle_diameter=9,
        colorbar_position='right'
    ), plt.savefig("DH4_2nd_15yr_Oligo_clusters.pdf")
	
clust_labels = ['means_per_cluster_mu_fg_Astro L1-2 FGFR3 GFAP', 'means_per_cluster_mu_fg_Astro L1-6 FGFR3 SLC14A1'
                        ]
clust_col = ['' + str(i) for i in clust_labels] # in case column names differ from labels

slide = select_slide(adata_vis, 'DH1a')

with mpl.rc_context({'figure.figsize': (20, 20)}):
    fig = plot_spatial(
        adata=slide,
        # labels to show on a plot
        color=clust_col, labels=clust_labels,
		alpha_scaling=1.0,
        show_img=True,
		img_alpha=0.0,
		adjust_text=True,
        # 'fast' (white background) or 'dark_background'
        style='dark_background',
        # limit color scale at 99.2% quantile of cell abundance
        max_color_quantile=0.992,
        # size of locations (adjust depending on figure size)
        circle_diameter=9,
        colorbar_position='right'
    ), plt.savefig("DH1a_first15yr_Astro_clusters.pdf")

clust_labels = ['means_per_cluster_mu_fg_Exc L2 LAMP5 LTK', 'means_per_cluster_mu_fg_Exc L2-3 LINC00507 FREM3', 'means_per_cluster_mu_fg_Exc L3-5 RORB COL22A1',
                        'means_per_cluster_mu_fg_Exc L3-5 RORB ESR1', 'means_per_cluster_mu_fg_Exc L3-5 RORB FILIP1L', 'means_per_cluster_mu_fg_Exc L3-5 RORB TWIST2'
                        ]
clust_col = ['' + str(i) for i in clust_labels] # in case column names differ from labels

slide = select_slide(adata_vis, 'DH1a')

with mpl.rc_context({'figure.figsize': (20, 20)}):
    fig = plot_spatial(
        adata=slide,
        # labels to show on a plot
        color=clust_col, labels=clust_labels,
		alpha_scaling=1.0,
        show_img=True,
		img_alpha=0.0,
		adjust_text=True,
        # 'fast' (white background) or 'dark_background'
        style='dark_background',
        # limit color scale at 99.2% quantile of cell abundance
        max_color_quantile=0.992,
        # size of locations (adjust depending on figure size)
        circle_diameter=9,
        colorbar_position='right'
    ), plt.savefig("DH1a_first15yr_Exc_L3_5_clusters.pdf")
	
clust_labels = ['means_per_cluster_mu_fg_Exc L4-5 RORB DAPK2', 'means_per_cluster_mu_fg_Exc L4-5 RORB FOLH1B', 'means_per_cluster_mu_fg_Exc L4-6 FEZF2 IL26', 'means_per_cluster_mu_fg_Exc L4-6 RORB SEMA3E'
                        ]
clust_col = ['' + str(i) for i in clust_labels] # in case column names differ from labels

slide = select_slide(adata_vis, 'DH1a')

with mpl.rc_context({'figure.figsize': (20, 20)}):
    fig = plot_spatial(
        adata=slide,
        # labels to show on a plot
        color=clust_col, labels=clust_labels,
		alpha_scaling=1.0,
        show_img=True,
		img_alpha=0.0,
		adjust_text=True,
        # 'fast' (white background) or 'dark_background'
        style='dark_background',
        # limit color scale at 99.2% quantile of cell abundance
        max_color_quantile=0.992,
        # size of locations (adjust depending on figure size)
        circle_diameter=9,
        colorbar_position='right'
    ), plt.savefig("DH1a_first15yr_Exc_L4_6_clusters.pdf")
	
clust_labels = ['means_per_cluster_mu_fg_Exc L5-6 FEZF2 ABO', 'means_per_cluster_mu_fg_Exc L5-6 THEMIS C1QL3', 'means_per_cluster_mu_fg_Exc L5-6 THEMIS FGF10'
                        ]
clust_col = ['' + str(i) for i in clust_labels] # in case column names differ from labels

slide = select_slide(adata_vis, 'DH1a')

with mpl.rc_context({'figure.figsize': (20, 20)}):
    fig = plot_spatial(
        adata=slide,
        # labels to show on a plot
        color=clust_col, labels=clust_labels,
		alpha_scaling=1.0,
        show_img=True,
		img_alpha=0.0,
		adjust_text=True,
        # 'fast' (white background) or 'dark_background'
        style='dark_background',
        # limit color scale at 99.2% quantile of cell abundance
        max_color_quantile=0.992,
        # size of locations (adjust depending on figure size)
        circle_diameter=9,
        colorbar_position='right'
    ), plt.savefig("DH1a_first15yr_Exc_L5_6_clusters.pdf")
	
clust_labels = ['means_per_cluster_mu_fg_Inh L1 SST NMBR', 'means_per_cluster_mu_fg_Inh L1-2 GAD1 MC4R', 'means_per_cluster_mu_fg_Inh L1-2 PAX6 CDH12', 
				  'means_per_cluster_mu_fg_Inh L1-2 SST BAGE2', 'means_per_cluster_mu_fg_Inh L1-2 VIP LBH', 'means_per_cluster_mu_fg_Inh L1-2 VIP PCDH20' 
                        ]
clust_col = ['' + str(i) for i in clust_labels] # in case column names differ from labels

slide = select_slide(adata_vis, 'DH1a')

with mpl.rc_context({'figure.figsize': (20, 20)}):
    fig = plot_spatial(
        adata=slide,
        # labels to show on a plot
        color=clust_col, labels=clust_labels,
		alpha_scaling=1.0,
        show_img=True,
		img_alpha=0.0,
		adjust_text=True,
        # 'fast' (white background) or 'dark_background'
        style='dark_background',
        # limit color scale at 99.2% quantile of cell abundance
        max_color_quantile=0.992,
        # size of locations (adjust depending on figure size)
        circle_diameter=9,
        colorbar_position='right'
    ), plt.savefig("DH1a_first15yr_Inh_L1_2_clusters.pdf")
	
clust_labels = ['means_per_cluster_mu_fg_Inh L1-2 VIP TSPAN12', 'means_per_cluster_mu_fg_Inh L1-3 SST CALB1', 'means_per_cluster_mu_fg_Inh L1-3 VIP ADAMTSL1', 'means_per_cluster_mu_fg_Inh L1-3 VIP CCDC184', 'means_per_cluster_mu_fg_Inh L1-3 VIP CHRM2', 
				  'means_per_cluster_mu_fg_Inh L1-3 VIP GGH'
                        ]
clust_col = ['' + str(i) for i in clust_labels] # in case column names differ from labels

slide = select_slide(adata_vis, 'DH1a')

with mpl.rc_context({'figure.figsize': (20, 20)}):
    fig = plot_spatial(
        adata=slide,
        # labels to show on a plot
        color=clust_col, labels=clust_labels,
		alpha_scaling=1.0,
        show_img=True,
		img_alpha=0.0,
		adjust_text=True,
        # 'fast' (white background) or 'dark_background'
        style='dark_background',
        # limit color scale at 99.2% quantile of cell abundance
        max_color_quantile=0.992,
        # size of locations (adjust depending on figure size)
        circle_diameter=9,
        colorbar_position='right'
    ), plt.savefig("DH1a_first15yr_Inh_L1_3_clusters.pdf")
	
clust_labels = ['means_per_cluster_mu_fg_Inh L1-4 LAMP5 LCP2', 'means_per_cluster_mu_fg_Inh L1-4 VIP CHRNA6', 'means_per_cluster_mu_fg_Inh L1-4 VIP OPRM1'
                        ]
clust_col = ['' + str(i) for i in clust_labels] # in case column names differ from labels

slide = select_slide(adata_vis, 'DH1a')

with mpl.rc_context({'figure.figsize': (20, 20)}):
    fig = plot_spatial(
        adata=slide,
        # labels to show on a plot
        color=clust_col, labels=clust_labels,
		alpha_scaling=1.0,
        show_img=True,
		img_alpha=0.0,
		adjust_text=True,
        # 'fast' (white background) or 'dark_background'
        style='dark_background',
        # limit color scale at 99.2% quantile of cell abundance
        max_color_quantile=0.992,
        # size of locations (adjust depending on figure size)
        circle_diameter=9,
        colorbar_position='right'
    ), plt.savefig("DH1a_first15yr_Inh_L1_4_clusters.pdf")
	
clust_labels = ['means_per_cluster_mu_fg_Inh L2-4 PVALB WFDC2', 'means_per_cluster_mu_fg_Inh L2-4 SST FRZB', 'means_per_cluster_mu_fg_Inh L2-4 VIP CBLN1', 'means_per_cluster_mu_fg_Inh L2-4 VIP SPAG17'
                        ]
clust_col = ['' + str(i) for i in clust_labels] # in case column names differ from labels

slide = select_slide(adata_vis, 'DH1a')

with mpl.rc_context({'figure.figsize': (20, 20)}):
    fig = plot_spatial(
        adata=slide,
        # labels to show on a plot
        color=clust_col, labels=clust_labels,
		alpha_scaling=1.0,
        show_img=True,
		img_alpha=0.0,
		adjust_text=True,
        # 'fast' (white background) or 'dark_background'
        style='dark_background',
        # limit color scale at 99.2% quantile of cell abundance
        max_color_quantile=0.992,
        # size of locations (adjust depending on figure size)
        circle_diameter=9,
        colorbar_position='right'
    ), plt.savefig("DH1a_first15yr_Inh_L2_4_clusters.pdf")
	
clust_labels = ['means_per_cluster_mu_fg_Inh L2-5 PVALB SCUBE3', 'means_per_cluster_mu_fg_Inh L2-5 VIP SERPINF1', 'means_per_cluster_mu_fg_Inh L2-5 VIP TYR', 'means_per_cluster_mu_fg_Inh L2-6 LAMP5 CA1', 
				  'means_per_cluster_mu_fg_Inh L2-6 VIP QPCT', 'means_per_cluster_mu_fg_Inh L3-5 SST ADGRG6', 'means_per_cluster_mu_fg_Inh L3-6 VIP HS3ST3A1'
                        ]
clust_col = ['' + str(i) for i in clust_labels] # in case column names differ from labels

slide = select_slide(adata_vis, 'DH1a')

with mpl.rc_context({'figure.figsize': (20, 20)}):
    fig = plot_spatial(
        adata=slide,
        # labels to show on a plot
        color=clust_col, labels=clust_labels,
		alpha_scaling=1.0,
        show_img=True,
		img_alpha=0.0,
		adjust_text=True,
        # 'fast' (white background) or 'dark_background'
        style='dark_background',
        # limit color scale at 99.2% quantile of cell abundance
        max_color_quantile=0.992,
        # size of locations (adjust depending on figure size)
        circle_diameter=9,
        colorbar_position='right'
    ), plt.savefig("DH1a_first15yr_Inh_L2_6_clusters.pdf")
	

	
clust_labels = ['means_per_cluster_mu_fg_Inh L4-5 SST STK32A', 'means_per_cluster_mu_fg_Inh L4-6 PVALB SULF1', 'means_per_cluster_mu_fg_Inh L4-6 SST B3GAT2', 'means_per_cluster_mu_fg_Inh L4-6 SST GXYLT2', 
				  'means_per_cluster_mu_fg_Inh L5-6 PVALB LGR5', 'means_per_cluster_mu_fg_Inh L5-6 SST MIR548F2', 'means_per_cluster_mu_fg_Inh L5-6 SST NPM1P10'
                        ]
clust_col = ['' + str(i) for i in clust_labels] # in case column names differ from labels

slide = select_slide(adata_vis, 'DH1a')

with mpl.rc_context({'figure.figsize': (20, 20)}):
    fig = plot_spatial(
        adata=slide,
        # labels to show on a plot
        color=clust_col, labels=clust_labels,
		alpha_scaling=1.0,
        show_img=True,
		img_alpha=0.0,
		adjust_text=True,
        # 'fast' (white background) or 'dark_background'
        style='dark_background',
        # limit color scale at 99.2% quantile of cell abundance
        max_color_quantile=0.992,
        # size of locations (adjust depending on figure size)
        circle_diameter=9,
        colorbar_position='right'
    ), plt.savefig("DH1a_first15yr_Inh_L5_6_clusters.pdf")
	
	
clust_labels = ['means_per_cluster_mu_fg_Micro L1-3 TYROBP', 'means_per_cluster_mu_fg_Endo L2-6 NOSTRIN'
               ]
clust_col = ['' + str(i) for i in clust_labels] # in case column names differ from labels

slide = select_slide(adata_vis, 'DH1a')

with mpl.rc_context({'figure.figsize': (20, 20)}):
    fig = plot_spatial(
        adata=slide,
        # labels to show on a plot
        color=clust_col, labels=clust_labels,
		alpha_scaling=1.0,
        show_img=True,
		img_alpha=0.0,
		adjust_text=True,
        # 'fast' (white background) or 'dark_background'
        style='dark_background',
        # limit color scale at 99.2% quantile of cell abundance
        max_color_quantile=0.992,
        # size of locations (adjust depending on figure size)
        circle_diameter=9,
        colorbar_position='right'
    ), plt.savefig("DH1a_first15yr_Endo_clusters.pdf")
	
clust_labels = ['means_per_cluster_mu_fg_OPC L1-6 PDGFRA', 'means_per_cluster_mu_fg_Oligo L1-6 OPALIN', 'means_per_cluster_mu_fg_Inh L3-6 SST NPY'	
               ]
clust_col = ['' + str(i) for i in clust_labels] # in case column names differ from labels

slide = select_slide(adata_vis, 'DH1a')

with mpl.rc_context({'figure.figsize': (20, 20)}):
    fig = plot_spatial(
        adata=slide,
        # labels to show on a plot
        color=clust_col, labels=clust_labels,
		alpha_scaling=1.0,
        show_img=True,
		img_alpha=0.0,
		adjust_text=True,
        # 'fast' (white background) or 'dark_background'
        style='dark_background',
        # limit color scale at 99.2% quantile of cell abundance
        max_color_quantile=0.992,
        # size of locations (adjust depending on figure size)
        circle_diameter=9,
        colorbar_position='right'
    ), plt.savefig("DH1a_first15yr_Oligo_clusters.pdf")
	
clust_labels = ['means_per_cluster_mu_fg_Astro L1-2 FGFR3 GFAP', 'means_per_cluster_mu_fg_Astro L1-6 FGFR3 SLC14A1'
                        ]
clust_col = ['' + str(i) for i in clust_labels] # in case column names differ from labels

slide = select_slide(adata_vis, 'DH2a')

with mpl.rc_context({'figure.figsize': (20, 20)}):
    fig = plot_spatial(
        adata=slide,
        # labels to show on a plot
        color=clust_col, labels=clust_labels,
		alpha_scaling=1.0,
        show_img=True,
		img_alpha=0.0,
		adjust_text=True,
        # 'fast' (white background) or 'dark_background'
        style='dark_background',
        # limit color scale at 99.2% quantile of cell abundance
        max_color_quantile=0.992,
        # size of locations (adjust depending on figure size)
        circle_diameter=9,
        colorbar_position='right'
    ), plt.savefig("DH2a_first15yr_Astro_clusters.pdf")

clust_labels = ['means_per_cluster_mu_fg_Exc L2 LAMP5 LTK', 'means_per_cluster_mu_fg_Exc L2-3 LINC00507 FREM3', 'means_per_cluster_mu_fg_Exc L3-5 RORB COL22A1',
                        'means_per_cluster_mu_fg_Exc L3-5 RORB ESR1', 'means_per_cluster_mu_fg_Exc L3-5 RORB FILIP1L', 'means_per_cluster_mu_fg_Exc L3-5 RORB TWIST2'
                        ]
clust_col = ['' + str(i) for i in clust_labels] # in case column names differ from labels

slide = select_slide(adata_vis, 'DH2a')

with mpl.rc_context({'figure.figsize': (20, 20)}):
    fig = plot_spatial(
        adata=slide,
        # labels to show on a plot
        color=clust_col, labels=clust_labels,
		alpha_scaling=1.0,
        show_img=True,
		img_alpha=0.0,
		adjust_text=True,
        # 'fast' (white background) or 'dark_background'
        style='dark_background',
        # limit color scale at 99.2% quantile of cell abundance
        max_color_quantile=0.992,
        # size of locations (adjust depending on figure size)
        circle_diameter=9,
        colorbar_position='right'
    ), plt.savefig("DH2a_first15yr_Exc_L3_5_clusters.pdf")
	
clust_labels = ['means_per_cluster_mu_fg_Exc L4-5 RORB DAPK2', 'means_per_cluster_mu_fg_Exc L4-5 RORB FOLH1B', 'means_per_cluster_mu_fg_Exc L4-6 FEZF2 IL26', 'means_per_cluster_mu_fg_Exc L4-6 RORB SEMA3E'
                        ]
clust_col = ['' + str(i) for i in clust_labels] # in case column names differ from labels

slide = select_slide(adata_vis, 'DH2a')

with mpl.rc_context({'figure.figsize': (20, 20)}):
    fig = plot_spatial(
        adata=slide,
        # labels to show on a plot
        color=clust_col, labels=clust_labels,
		alpha_scaling=1.0,
        show_img=True,
		img_alpha=0.0,
		adjust_text=True,
        # 'fast' (white background) or 'dark_background'
        style='dark_background',
        # limit color scale at 99.2% quantile of cell abundance
        max_color_quantile=0.992,
        # size of locations (adjust depending on figure size)
        circle_diameter=9,
        colorbar_position='right'
    ), plt.savefig("DH2a_first15yr_Exc_L4_6_clusters.pdf")
	
clust_labels = ['means_per_cluster_mu_fg_Exc L5-6 FEZF2 ABO', 'means_per_cluster_mu_fg_Exc L5-6 THEMIS C1QL3', 'means_per_cluster_mu_fg_Exc L5-6 THEMIS FGF10'
                        ]
clust_col = ['' + str(i) for i in clust_labels] # in case column names differ from labels

slide = select_slide(adata_vis, 'DH2a')

with mpl.rc_context({'figure.figsize': (20, 20)}):
    fig = plot_spatial(
        adata=slide,
        # labels to show on a plot
        color=clust_col, labels=clust_labels,
		alpha_scaling=1.0,
        show_img=True,
		img_alpha=0.0,
		adjust_text=True,
        # 'fast' (white background) or 'dark_background'
        style='dark_background',
        # limit color scale at 99.2% quantile of cell abundance
        max_color_quantile=0.992,
        # size of locations (adjust depending on figure size)
        circle_diameter=9,
        colorbar_position='right'
    ), plt.savefig("DH2a_first15yr_Exc_L5_6_clusters.pdf")
	
clust_labels = ['means_per_cluster_mu_fg_Inh L1 SST NMBR', 'means_per_cluster_mu_fg_Inh L1-2 GAD1 MC4R', 'means_per_cluster_mu_fg_Inh L1-2 PAX6 CDH12', 
				  'means_per_cluster_mu_fg_Inh L1-2 SST BAGE2', 'means_per_cluster_mu_fg_Inh L1-2 VIP LBH', 'means_per_cluster_mu_fg_Inh L1-2 VIP PCDH20' 
                        ]
clust_col = ['' + str(i) for i in clust_labels] # in case column names differ from labels

slide = select_slide(adata_vis, 'DH2a')

with mpl.rc_context({'figure.figsize': (20, 20)}):
    fig = plot_spatial(
        adata=slide,
        # labels to show on a plot
        color=clust_col, labels=clust_labels,
		alpha_scaling=1.0,
        show_img=True,
		img_alpha=0.0,
		adjust_text=True,
        # 'fast' (white background) or 'dark_background'
        style='dark_background',
        # limit color scale at 99.2% quantile of cell abundance
        max_color_quantile=0.992,
        # size of locations (adjust depending on figure size)
        circle_diameter=9,
        colorbar_position='right'
    ), plt.savefig("DH2a_first15yr_Inh_L1_2_clusters.pdf")
	
clust_labels = ['means_per_cluster_mu_fg_Inh L1-2 VIP TSPAN12', 'means_per_cluster_mu_fg_Inh L1-3 SST CALB1', 'means_per_cluster_mu_fg_Inh L1-3 VIP ADAMTSL1', 'means_per_cluster_mu_fg_Inh L1-3 VIP CCDC184', 'means_per_cluster_mu_fg_Inh L1-3 VIP CHRM2', 
				  'means_per_cluster_mu_fg_Inh L1-3 VIP GGH'
                        ]
clust_col = ['' + str(i) for i in clust_labels] # in case column names differ from labels

slide = select_slide(adata_vis, 'DH2a')

with mpl.rc_context({'figure.figsize': (20, 20)}):
    fig = plot_spatial(
        adata=slide,
        # labels to show on a plot
        color=clust_col, labels=clust_labels,
		alpha_scaling=1.0,
        show_img=True,
		img_alpha=0.0,
		adjust_text=True,
        # 'fast' (white background) or 'dark_background'
        style='dark_background',
        # limit color scale at 99.2% quantile of cell abundance
        max_color_quantile=0.992,
        # size of locations (adjust depending on figure size)
        circle_diameter=9,
        colorbar_position='right'
    ), plt.savefig("DH2a_first15yr_Inh_L1_3_clusters.pdf")
	
clust_labels = ['means_per_cluster_mu_fg_Inh L1-4 LAMP5 LCP2', 'means_per_cluster_mu_fg_Inh L1-4 VIP CHRNA6', 'means_per_cluster_mu_fg_Inh L1-4 VIP OPRM1'
                        ]
clust_col = ['' + str(i) for i in clust_labels] # in case column names differ from labels

slide = select_slide(adata_vis, 'DH2a')

with mpl.rc_context({'figure.figsize': (20, 20)}):
    fig = plot_spatial(
        adata=slide,
        # labels to show on a plot
        color=clust_col, labels=clust_labels,
		alpha_scaling=1.0,
        show_img=True,
		img_alpha=0.0,
		adjust_text=True,
        # 'fast' (white background) or 'dark_background'
        style='dark_background',
        # limit color scale at 99.2% quantile of cell abundance
        max_color_quantile=0.992,
        # size of locations (adjust depending on figure size)
        circle_diameter=9,
        colorbar_position='right'
    ), plt.savefig("DH2a_first15yr_Inh_L1_4_clusters.pdf")
	
clust_labels = ['means_per_cluster_mu_fg_Inh L2-4 PVALB WFDC2', 'means_per_cluster_mu_fg_Inh L2-4 SST FRZB', 'means_per_cluster_mu_fg_Inh L2-4 VIP CBLN1', 'means_per_cluster_mu_fg_Inh L2-4 VIP SPAG17'
                        ]
clust_col = ['' + str(i) for i in clust_labels] # in case column names differ from labels

slide = select_slide(adata_vis, 'DH2a')

with mpl.rc_context({'figure.figsize': (20, 20)}):
    fig = plot_spatial(
        adata=slide,
        # labels to show on a plot
        color=clust_col, labels=clust_labels,
		alpha_scaling=1.0,
        show_img=True,
		img_alpha=0.0,
		adjust_text=True,
        # 'fast' (white background) or 'dark_background'
        style='dark_background',
        # limit color scale at 99.2% quantile of cell abundance
        max_color_quantile=0.992,
        # size of locations (adjust depending on figure size)
        circle_diameter=9,
        colorbar_position='right'
    ), plt.savefig("DH2a_first15yr_Inh_L2_4_clusters.pdf")
	
clust_labels = ['means_per_cluster_mu_fg_Inh L2-5 PVALB SCUBE3', 'means_per_cluster_mu_fg_Inh L2-5 VIP SERPINF1', 'means_per_cluster_mu_fg_Inh L2-5 VIP TYR', 'means_per_cluster_mu_fg_Inh L2-6 LAMP5 CA1', 
				  'means_per_cluster_mu_fg_Inh L2-6 VIP QPCT', 'means_per_cluster_mu_fg_Inh L3-5 SST ADGRG6', 'means_per_cluster_mu_fg_Inh L3-6 VIP HS3ST3A1'
                        ]
clust_col = ['' + str(i) for i in clust_labels] # in case column names differ from labels

slide = select_slide(adata_vis, 'DH2a')

with mpl.rc_context({'figure.figsize': (20, 20)}):
    fig = plot_spatial(
        adata=slide,
        # labels to show on a plot
        color=clust_col, labels=clust_labels,
		alpha_scaling=1.0,
        show_img=True,
		img_alpha=0.0,
		adjust_text=True,
        # 'fast' (white background) or 'dark_background'
        style='dark_background',
        # limit color scale at 99.2% quantile of cell abundance
        max_color_quantile=0.992,
        # size of locations (adjust depending on figure size)
        circle_diameter=9,
        colorbar_position='right'
    ), plt.savefig("DH2a_first15yr_Inh_L2_6_clusters.pdf")
	

	
clust_labels = ['means_per_cluster_mu_fg_Inh L4-5 SST STK32A', 'means_per_cluster_mu_fg_Inh L4-6 PVALB SULF1', 'means_per_cluster_mu_fg_Inh L4-6 SST B3GAT2', 'means_per_cluster_mu_fg_Inh L4-6 SST GXYLT2', 
				  'means_per_cluster_mu_fg_Inh L5-6 PVALB LGR5', 'means_per_cluster_mu_fg_Inh L5-6 SST MIR548F2', 'means_per_cluster_mu_fg_Inh L5-6 SST NPM1P10'
                        ]
clust_col = ['' + str(i) for i in clust_labels] # in case column names differ from labels

slide = select_slide(adata_vis, 'DH2a')

with mpl.rc_context({'figure.figsize': (20, 20)}):
    fig = plot_spatial(
        adata=slide,
        # labels to show on a plot
        color=clust_col, labels=clust_labels,
		alpha_scaling=1.0,
        show_img=True,
		img_alpha=0.0,
		adjust_text=True,
        # 'fast' (white background) or 'dark_background'
        style='dark_background',
        # limit color scale at 99.2% quantile of cell abundance
        max_color_quantile=0.992,
        # size of locations (adjust depending on figure size)
        circle_diameter=9,
        colorbar_position='right'
    ), plt.savefig("DH2a_first15yr_Inh_L5_6_clusters.pdf")
	
	
clust_labels = ['means_per_cluster_mu_fg_Micro L1-3 TYROBP', 'means_per_cluster_mu_fg_Endo L2-6 NOSTRIN'
               ]
clust_col = ['' + str(i) for i in clust_labels] # in case column names differ from labels

slide = select_slide(adata_vis, 'DH2a')

with mpl.rc_context({'figure.figsize': (20, 20)}):
    fig = plot_spatial(
        adata=slide,
        # labels to show on a plot
        color=clust_col, labels=clust_labels,
		alpha_scaling=1.0,
        show_img=True,
		img_alpha=0.0,
		adjust_text=True,
        # 'fast' (white background) or 'dark_background'
        style='dark_background',
        # limit color scale at 99.2% quantile of cell abundance
        max_color_quantile=0.992,
        # size of locations (adjust depending on figure size)
        circle_diameter=9,
        colorbar_position='right'
    ), plt.savefig("DH2a_first15yr_Endo_clusters.pdf")
	
clust_labels = ['means_per_cluster_mu_fg_OPC L1-6 PDGFRA', 'means_per_cluster_mu_fg_Oligo L1-6 OPALIN', 'means_per_cluster_mu_fg_Inh L3-6 SST NPY'	
               ]
clust_col = ['' + str(i) for i in clust_labels] # in case column names differ from labels

slide = select_slide(adata_vis, 'DH2a')

with mpl.rc_context({'figure.figsize': (20, 20)}):
    fig = plot_spatial(
        adata=slide,
        # labels to show on a plot
        color=clust_col, labels=clust_labels,
		alpha_scaling=1.0,
        show_img=True,
		img_alpha=0.0,
		adjust_text=True,
        # 'fast' (white background) or 'dark_background'
        style='dark_background',
        # limit color scale at 99.2% quantile of cell abundance
        max_color_quantile=0.992,
        # size of locations (adjust depending on figure size)
        circle_diameter=9,
        colorbar_position='right'
    ), plt.savefig("DH2a_first15yr_Oligo_clusters.pdf")
	
	
clust_labels = ['means_per_cluster_mu_fg_Astro L1-2 FGFR3 GFAP', 'means_per_cluster_mu_fg_Astro L1-6 FGFR3 SLC14A1'
                        ]
clust_col = ['' + str(i) for i in clust_labels] # in case column names differ from labels

slide = select_slide(adata_vis, 'DH3a')

with mpl.rc_context({'figure.figsize': (20, 20)}):
    fig = plot_spatial(
        adata=slide,
        # labels to show on a plot
        color=clust_col, labels=clust_labels,
		alpha_scaling=1.0,
        show_img=True,
		img_alpha=0.0,
		adjust_text=True,
        # 'fast' (white background) or 'dark_background'
        style='dark_background',
        # limit color scale at 99.2% quantile of cell abundance
        max_color_quantile=0.992,
        # size of locations (adjust depending on figure size)
        circle_diameter=9,
        colorbar_position='right'
    ), plt.savefig("DH3a_31yr_Astro_clusters.pdf")

clust_labels = ['means_per_cluster_mu_fg_Exc L2 LAMP5 LTK', 'means_per_cluster_mu_fg_Exc L2-3 LINC00507 FREM3', 'means_per_cluster_mu_fg_Exc L3-5 RORB COL22A1',
                        'means_per_cluster_mu_fg_Exc L3-5 RORB ESR1', 'means_per_cluster_mu_fg_Exc L3-5 RORB FILIP1L', 'means_per_cluster_mu_fg_Exc L3-5 RORB TWIST2'
                        ]
clust_col = ['' + str(i) for i in clust_labels] # in case column names differ from labels

slide = select_slide(adata_vis, 'DH3a')

with mpl.rc_context({'figure.figsize': (20, 20)}):
    fig = plot_spatial(
        adata=slide,
        # labels to show on a plot
        color=clust_col, labels=clust_labels,
		alpha_scaling=1.0,
        show_img=True,
		img_alpha=0.0,
		adjust_text=True,
        # 'fast' (white background) or 'dark_background'
        style='dark_background',
        # limit color scale at 99.2% quantile of cell abundance
        max_color_quantile=0.992,
        # size of locations (adjust depending on figure size)
        circle_diameter=9,
        colorbar_position='right'
    ), plt.savefig("DH3a_31yr_Exc_L3_5_clusters.pdf")
	
clust_labels = ['means_per_cluster_mu_fg_Exc L4-5 RORB DAPK2', 'means_per_cluster_mu_fg_Exc L4-5 RORB FOLH1B', 'means_per_cluster_mu_fg_Exc L4-6 FEZF2 IL26', 'means_per_cluster_mu_fg_Exc L4-6 RORB SEMA3E'
                        ]
clust_col = ['' + str(i) for i in clust_labels] # in case column names differ from labels

slide = select_slide(adata_vis, 'DH3a')

with mpl.rc_context({'figure.figsize': (20, 20)}):
    fig = plot_spatial(
        adata=slide,
        # labels to show on a plot
        color=clust_col, labels=clust_labels,
		alpha_scaling=1.0,
        show_img=True,
		img_alpha=0.0,
		adjust_text=True,
        # 'fast' (white background) or 'dark_background'
        style='dark_background',
        # limit color scale at 99.2% quantile of cell abundance
        max_color_quantile=0.992,
        # size of locations (adjust depending on figure size)
        circle_diameter=9,
        colorbar_position='right'
    ), plt.savefig("DH3a_31yr_Exc_L4_6_clusters.pdf")
	
clust_labels = ['means_per_cluster_mu_fg_Exc L5-6 FEZF2 ABO', 'means_per_cluster_mu_fg_Exc L5-6 THEMIS C1QL3', 'means_per_cluster_mu_fg_Exc L5-6 THEMIS FGF10'
                        ]
clust_col = ['' + str(i) for i in clust_labels] # in case column names differ from labels

slide = select_slide(adata_vis, 'DH3a')

with mpl.rc_context({'figure.figsize': (20, 20)}):
    fig = plot_spatial(
        adata=slide,
        # labels to show on a plot
        color=clust_col, labels=clust_labels,
		alpha_scaling=1.0,
        show_img=True,
		img_alpha=0.0,
		adjust_text=True,
        # 'fast' (white background) or 'dark_background'
        style='dark_background',
        # limit color scale at 99.2% quantile of cell abundance
        max_color_quantile=0.992,
        # size of locations (adjust depending on figure size)
        circle_diameter=9,
        colorbar_position='right'
    ), plt.savefig("DH3a_31yr_Exc_L5_6_clusters.pdf")
	
clust_labels = ['means_per_cluster_mu_fg_Inh L1 SST NMBR', 'means_per_cluster_mu_fg_Inh L1-2 GAD1 MC4R', 'means_per_cluster_mu_fg_Inh L1-2 PAX6 CDH12', 
				  'means_per_cluster_mu_fg_Inh L1-2 SST BAGE2', 'means_per_cluster_mu_fg_Inh L1-2 VIP LBH', 'means_per_cluster_mu_fg_Inh L1-2 VIP PCDH20' 
                        ]
clust_col = ['' + str(i) for i in clust_labels] # in case column names differ from labels

slide = select_slide(adata_vis, 'DH3a')

with mpl.rc_context({'figure.figsize': (20, 20)}):
    fig = plot_spatial(
        adata=slide,
        # labels to show on a plot
        color=clust_col, labels=clust_labels,
		alpha_scaling=1.0,
        show_img=True,
		img_alpha=0.0,
		adjust_text=True,
        # 'fast' (white background) or 'dark_background'
        style='dark_background',
        # limit color scale at 99.2% quantile of cell abundance
        max_color_quantile=0.992,
        # size of locations (adjust depending on figure size)
        circle_diameter=9,
        colorbar_position='right'
    ), plt.savefig("DH3a_31yr_Inh_L1_2_clusters.pdf")
	
clust_labels = ['means_per_cluster_mu_fg_Inh L1-2 VIP TSPAN12', 'means_per_cluster_mu_fg_Inh L1-3 SST CALB1', 'means_per_cluster_mu_fg_Inh L1-3 VIP ADAMTSL1', 'means_per_cluster_mu_fg_Inh L1-3 VIP CCDC184', 'means_per_cluster_mu_fg_Inh L1-3 VIP CHRM2', 
				  'means_per_cluster_mu_fg_Inh L1-3 VIP GGH'
                        ]
clust_col = ['' + str(i) for i in clust_labels] # in case column names differ from labels

slide = select_slide(adata_vis, 'DH3a')

with mpl.rc_context({'figure.figsize': (20, 20)}):
    fig = plot_spatial(
        adata=slide,
        # labels to show on a plot
        color=clust_col, labels=clust_labels,
		alpha_scaling=1.0,
        show_img=True,
		img_alpha=0.0,
		adjust_text=True,
        # 'fast' (white background) or 'dark_background'
        style='dark_background',
        # limit color scale at 99.2% quantile of cell abundance
        max_color_quantile=0.992,
        # size of locations (adjust depending on figure size)
        circle_diameter=9,
        colorbar_position='right'
    ), plt.savefig("DH3a_31yr_Inh_L1_3_clusters.pdf")
	
clust_labels = ['means_per_cluster_mu_fg_Inh L1-4 LAMP5 LCP2', 'means_per_cluster_mu_fg_Inh L1-4 VIP CHRNA6', 'means_per_cluster_mu_fg_Inh L1-4 VIP OPRM1'
                        ]
clust_col = ['' + str(i) for i in clust_labels] # in case column names differ from labels

slide = select_slide(adata_vis, 'DH3a')

with mpl.rc_context({'figure.figsize': (20, 20)}):
    fig = plot_spatial(
        adata=slide,
        # labels to show on a plot
        color=clust_col, labels=clust_labels,
		alpha_scaling=1.0,
        show_img=True,
		img_alpha=0.0,
		adjust_text=True,
        # 'fast' (white background) or 'dark_background'
        style='dark_background',
        # limit color scale at 99.2% quantile of cell abundance
        max_color_quantile=0.992,
        # size of locations (adjust depending on figure size)
        circle_diameter=9,
        colorbar_position='right'
    ), plt.savefig("DH3a_31yr_Inh_L1_4_clusters.pdf")
	
clust_labels = ['means_per_cluster_mu_fg_Inh L2-4 PVALB WFDC2', 'means_per_cluster_mu_fg_Inh L2-4 SST FRZB', 'means_per_cluster_mu_fg_Inh L2-4 VIP CBLN1', 'means_per_cluster_mu_fg_Inh L2-4 VIP SPAG17'
                        ]
clust_col = ['' + str(i) for i in clust_labels] # in case column names differ from labels

slide = select_slide(adata_vis, 'DH3a')

with mpl.rc_context({'figure.figsize': (20, 20)}):
    fig = plot_spatial(
        adata=slide,
        # labels to show on a plot
        color=clust_col, labels=clust_labels,
		alpha_scaling=1.0,
        show_img=True,
		img_alpha=0.0,
		adjust_text=True,
        # 'fast' (white background) or 'dark_background'
        style='dark_background',
        # limit color scale at 99.2% quantile of cell abundance
        max_color_quantile=0.992,
        # size of locations (adjust depending on figure size)
        circle_diameter=9,
        colorbar_position='right'
    ), plt.savefig("DH3a_31yr_Inh_L2_4_clusters.pdf")
	
clust_labels = ['means_per_cluster_mu_fg_Inh L2-5 PVALB SCUBE3', 'means_per_cluster_mu_fg_Inh L2-5 VIP SERPINF1', 'means_per_cluster_mu_fg_Inh L2-5 VIP TYR', 'means_per_cluster_mu_fg_Inh L2-6 LAMP5 CA1', 
				  'means_per_cluster_mu_fg_Inh L2-6 VIP QPCT', 'means_per_cluster_mu_fg_Inh L3-5 SST ADGRG6', 'means_per_cluster_mu_fg_Inh L3-6 VIP HS3ST3A1'
                        ]
clust_col = ['' + str(i) for i in clust_labels] # in case column names differ from labels

slide = select_slide(adata_vis, 'DH3a')

with mpl.rc_context({'figure.figsize': (20, 20)}):
    fig = plot_spatial(
        adata=slide,
        # labels to show on a plot
        color=clust_col, labels=clust_labels,
		alpha_scaling=1.0,
        show_img=True,
		img_alpha=0.0,
		adjust_text=True,
        # 'fast' (white background) or 'dark_background'
        style='dark_background',
        # limit color scale at 99.2% quantile of cell abundance
        max_color_quantile=0.992,
        # size of locations (adjust depending on figure size)
        circle_diameter=9,
        colorbar_position='right'
    ), plt.savefig("DH3a_31yr_Inh_L2_6_clusters.pdf")
	
	
clust_labels = ['means_per_cluster_mu_fg_Inh L4-5 SST STK32A', 'means_per_cluster_mu_fg_Inh L4-6 PVALB SULF1', 'means_per_cluster_mu_fg_Inh L4-6 SST B3GAT2', 'means_per_cluster_mu_fg_Inh L4-6 SST GXYLT2', 
				  'means_per_cluster_mu_fg_Inh L5-6 PVALB LGR5', 'means_per_cluster_mu_fg_Inh L5-6 SST MIR548F2', 'means_per_cluster_mu_fg_Inh L5-6 SST NPM1P10'
                        ]
clust_col = ['' + str(i) for i in clust_labels] # in case column names differ from labels

slide = select_slide(adata_vis, 'DH3a')

with mpl.rc_context({'figure.figsize': (20, 20)}):
    fig = plot_spatial(
        adata=slide,
        # labels to show on a plot
        color=clust_col, labels=clust_labels,
		alpha_scaling=1.0,
        show_img=True,
		img_alpha=0.0,
		adjust_text=True,
        # 'fast' (white background) or 'dark_background'
        style='dark_background',
        # limit color scale at 99.2% quantile of cell abundance
        max_color_quantile=0.992,
        # size of locations (adjust depending on figure size)
        circle_diameter=9,
        colorbar_position='right'
    ), plt.savefig("DH3a_31yr_Inh_L5_6_clusters.pdf")
	
	
clust_labels = ['means_per_cluster_mu_fg_Micro L1-3 TYROBP', 'means_per_cluster_mu_fg_Endo L2-6 NOSTRIN'
               ]
clust_col = ['' + str(i) for i in clust_labels] # in case column names differ from labels

slide = select_slide(adata_vis, 'DH3a')

with mpl.rc_context({'figure.figsize': (20, 20)}):
    fig = plot_spatial(
        adata=slide,
        # labels to show on a plot
        color=clust_col, labels=clust_labels,
		alpha_scaling=1.0,
        show_img=True,
		img_alpha=0.0,
		adjust_text=True,
        # 'fast' (white background) or 'dark_background'
        style='dark_background',
        # limit color scale at 99.2% quantile of cell abundance
        max_color_quantile=0.992,
        # size of locations (adjust depending on figure size)
        circle_diameter=9,
        colorbar_position='right'
    ), plt.savefig("DH3a_31yr_Endo_clusters.pdf")
	
clust_labels = ['means_per_cluster_mu_fg_OPC L1-6 PDGFRA', 'means_per_cluster_mu_fg_Oligo L1-6 OPALIN', 'means_per_cluster_mu_fg_Inh L3-6 SST NPY'	
               ]
clust_col = ['' + str(i) for i in clust_labels] # in case column names differ from labels

slide = select_slide(adata_vis, 'DH3a')

with mpl.rc_context({'figure.figsize': (20, 20)}):
    fig = plot_spatial(
        adata=slide,
        # labels to show on a plot
        color=clust_col, labels=clust_labels,
		alpha_scaling=1.0,
        show_img=True,
		img_alpha=0.0,
		adjust_text=True,
        # 'fast' (white background) or 'dark_background'
        style='dark_background',
        # limit color scale at 99.2% quantile of cell abundance
        max_color_quantile=0.992,
        # size of locations (adjust depending on figure size)
        circle_diameter=9,
        colorbar_position='right'
    ), plt.savefig("DH3a_31yr_Oligo_clusters.pdf")
	

clust_labels = ['means_per_cluster_mu_fg_Astro L1-2 FGFR3 GFAP', 'means_per_cluster_mu_fg_Astro L1-6 FGFR3 SLC14A1'
                        ]
clust_col = ['' + str(i) for i in clust_labels] # in case column names differ from labels

slide = select_slide(adata_vis, 'DH4a')

with mpl.rc_context({'figure.figsize': (20, 20)}):
    fig = plot_spatial(
        adata=slide,
        # labels to show on a plot
        color=clust_col, labels=clust_labels,
		alpha_scaling=1.0,
        show_img=True,
		img_alpha=0.0,
		adjust_text=True,
        # 'fast' (white background) or 'dark_background'
        style='dark_background',
        # limit color scale at 99.2% quantile of cell abundance
        max_color_quantile=0.992,
        # size of locations (adjust depending on figure size)
        circle_diameter=9,
        colorbar_position='right'
    ), plt.savefig("DH4a_31yr_Astro_clusters.pdf")

clust_labels = ['means_per_cluster_mu_fg_Exc L2 LAMP5 LTK', 'means_per_cluster_mu_fg_Exc L2-3 LINC00507 FREM3', 'means_per_cluster_mu_fg_Exc L3-5 RORB COL22A1',
                        'means_per_cluster_mu_fg_Exc L3-5 RORB ESR1', 'means_per_cluster_mu_fg_Exc L3-5 RORB FILIP1L', 'means_per_cluster_mu_fg_Exc L3-5 RORB TWIST2'
                        ]
clust_col = ['' + str(i) for i in clust_labels] # in case column names differ from labels

slide = select_slide(adata_vis, 'DH4a')

with mpl.rc_context({'figure.figsize': (20, 20)}):
    fig = plot_spatial(
        adata=slide,
        # labels to show on a plot
        color=clust_col, labels=clust_labels,
		alpha_scaling=1.0,
        show_img=True,
		img_alpha=0.0,
		adjust_text=True,
        # 'fast' (white background) or 'dark_background'
        style='dark_background',
        # limit color scale at 99.2% quantile of cell abundance
        max_color_quantile=0.992,
        # size of locations (adjust depending on figure size)
        circle_diameter=9,
        colorbar_position='right'
    ), plt.savefig("DH4a_31yr_Exc_L3_5_clusters.pdf")
	
clust_labels = ['means_per_cluster_mu_fg_Exc L4-5 RORB DAPK2', 'means_per_cluster_mu_fg_Exc L4-5 RORB FOLH1B', 'means_per_cluster_mu_fg_Exc L4-6 FEZF2 IL26', 'means_per_cluster_mu_fg_Exc L4-6 RORB SEMA3E'
                        ]
clust_col = ['' + str(i) for i in clust_labels] # in case column names differ from labels

slide = select_slide(adata_vis, 'DH4a')

with mpl.rc_context({'figure.figsize': (20, 20)}):
    fig = plot_spatial(
        adata=slide,
        # labels to show on a plot
        color=clust_col, labels=clust_labels,
		alpha_scaling=1.0,
        show_img=True,
		img_alpha=0.0,
		adjust_text=True,
        # 'fast' (white background) or 'dark_background'
        style='dark_background',
        # limit color scale at 99.2% quantile of cell abundance
        max_color_quantile=0.992,
        # size of locations (adjust depending on figure size)
        circle_diameter=9,
        colorbar_position='right'
    ), plt.savefig("DH4a_31yr_Exc_L4_6_clusters.pdf")
	
clust_labels = ['means_per_cluster_mu_fg_Exc L5-6 FEZF2 ABO', 'means_per_cluster_mu_fg_Exc L5-6 THEMIS C1QL3', 'means_per_cluster_mu_fg_Exc L5-6 THEMIS FGF10'
                        ]
clust_col = ['' + str(i) for i in clust_labels] # in case column names differ from labels

slide = select_slide(adata_vis, 'DH4a')

with mpl.rc_context({'figure.figsize': (20, 20)}):
    fig = plot_spatial(
        adata=slide,
        # labels to show on a plot
        color=clust_col, labels=clust_labels,
		alpha_scaling=1.0,
        show_img=True,
		img_alpha=0.0,
		adjust_text=True,
        # 'fast' (white background) or 'dark_background'
        style='dark_background',
        # limit color scale at 99.2% quantile of cell abundance
        max_color_quantile=0.992,
        # size of locations (adjust depending on figure size)
        circle_diameter=9,
        colorbar_position='right'
    ), plt.savefig("DH4a_31yr_Exc_L5_6_clusters.pdf")
	
clust_labels = ['means_per_cluster_mu_fg_Inh L1 SST NMBR', 'means_per_cluster_mu_fg_Inh L1-2 GAD1 MC4R', 'means_per_cluster_mu_fg_Inh L1-2 PAX6 CDH12', 
				  'means_per_cluster_mu_fg_Inh L1-2 SST BAGE2', 'means_per_cluster_mu_fg_Inh L1-2 VIP LBH', 'means_per_cluster_mu_fg_Inh L1-2 VIP PCDH20' 
                        ]
clust_col = ['' + str(i) for i in clust_labels] # in case column names differ from labels

slide = select_slide(adata_vis, 'DH4a')

with mpl.rc_context({'figure.figsize': (20, 20)}):
    fig = plot_spatial(
        adata=slide,
        # labels to show on a plot
        color=clust_col, labels=clust_labels,
		alpha_scaling=1.0,
        show_img=True,
		img_alpha=0.0,
		adjust_text=True,
        # 'fast' (white background) or 'dark_background'
        style='dark_background',
        # limit color scale at 99.2% quantile of cell abundance
        max_color_quantile=0.992,
        # size of locations (adjust depending on figure size)
        circle_diameter=9,
        colorbar_position='right'
    ), plt.savefig("DH4a_31yr_Inh_L1_2_clusters.pdf")
	
clust_labels = ['means_per_cluster_mu_fg_Inh L1-2 VIP TSPAN12', 'means_per_cluster_mu_fg_Inh L1-3 SST CALB1', 'means_per_cluster_mu_fg_Inh L1-3 VIP ADAMTSL1', 'means_per_cluster_mu_fg_Inh L1-3 VIP CCDC184', 'means_per_cluster_mu_fg_Inh L1-3 VIP CHRM2', 
				  'means_per_cluster_mu_fg_Inh L1-3 VIP GGH'
                        ]
clust_col = ['' + str(i) for i in clust_labels] # in case column names differ from labels

slide = select_slide(adata_vis, 'DH4a')

with mpl.rc_context({'figure.figsize': (20, 20)}):
    fig = plot_spatial(
        adata=slide,
        # labels to show on a plot
        color=clust_col, labels=clust_labels,
		alpha_scaling=1.0,
        show_img=True,
		img_alpha=0.0,
		adjust_text=True,
        # 'fast' (white background) or 'dark_background'
        style='dark_background',
        # limit color scale at 99.2% quantile of cell abundance
        max_color_quantile=0.992,
        # size of locations (adjust depending on figure size)
        circle_diameter=9,
        colorbar_position='right'
    ), plt.savefig("DH4a_31yr_Inh_L1_3_clusters.pdf")
	
clust_labels = ['means_per_cluster_mu_fg_Inh L1-4 LAMP5 LCP2', 'means_per_cluster_mu_fg_Inh L1-4 VIP CHRNA6', 'means_per_cluster_mu_fg_Inh L1-4 VIP OPRM1'
                        ]
clust_col = ['' + str(i) for i in clust_labels] # in case column names differ from labels

slide = select_slide(adata_vis, 'DH4a')

with mpl.rc_context({'figure.figsize': (20, 20)}):
    fig = plot_spatial(
        adata=slide,
        # labels to show on a plot
        color=clust_col, labels=clust_labels,
		alpha_scaling=1.0,
        show_img=True,
		img_alpha=0.0,
		adjust_text=True,
        # 'fast' (white background) or 'dark_background'
        style='dark_background',
        # limit color scale at 99.2% quantile of cell abundance
        max_color_quantile=0.992,
        # size of locations (adjust depending on figure size)
        circle_diameter=9,
        colorbar_position='right'
    ), plt.savefig("DH4a_31yr_Inh_L1_4_clusters.pdf")
	
clust_labels = ['means_per_cluster_mu_fg_Inh L2-4 PVALB WFDC2', 'means_per_cluster_mu_fg_Inh L2-4 SST FRZB', 'means_per_cluster_mu_fg_Inh L2-4 VIP CBLN1', 'means_per_cluster_mu_fg_Inh L2-4 VIP SPAG17'
                        ]
clust_col = ['' + str(i) for i in clust_labels] # in case column names differ from labels

slide = select_slide(adata_vis, 'DH4a')

with mpl.rc_context({'figure.figsize': (20, 20)}):
    fig = plot_spatial(
        adata=slide,
        # labels to show on a plot
        color=clust_col, labels=clust_labels,
		alpha_scaling=1.0,
        show_img=True,
		img_alpha=0.0,
		adjust_text=True,
        # 'fast' (white background) or 'dark_background'
        style='dark_background',
        # limit color scale at 99.2% quantile of cell abundance
        max_color_quantile=0.992,
        # size of locations (adjust depending on figure size)
        circle_diameter=9,
        colorbar_position='right'
    ), plt.savefig("DH4a_31yr_Inh_L2_4_clusters.pdf")
	
clust_labels = ['means_per_cluster_mu_fg_Inh L2-5 PVALB SCUBE3', 'means_per_cluster_mu_fg_Inh L2-5 VIP SERPINF1', 'means_per_cluster_mu_fg_Inh L2-5 VIP TYR', 'means_per_cluster_mu_fg_Inh L2-6 LAMP5 CA1', 
				  'means_per_cluster_mu_fg_Inh L2-6 VIP QPCT', 'means_per_cluster_mu_fg_Inh L3-5 SST ADGRG6', 'means_per_cluster_mu_fg_Inh L3-6 VIP HS3ST3A1'
                        ]
clust_col = ['' + str(i) for i in clust_labels] # in case column names differ from labels

slide = select_slide(adata_vis, 'DH4a')

with mpl.rc_context({'figure.figsize': (20, 20)}):
    fig = plot_spatial(
        adata=slide,
        # labels to show on a plot
        color=clust_col, labels=clust_labels,
		alpha_scaling=1.0,
        show_img=True,
		img_alpha=0.0,
		adjust_text=True,
        # 'fast' (white background) or 'dark_background'
        style='dark_background',
        # limit color scale at 99.2% quantile of cell abundance
        max_color_quantile=0.992,
        # size of locations (adjust depending on figure size)
        circle_diameter=9,
        colorbar_position='right'
    ), plt.savefig("DH4a_31yr_Inh_L2_6_clusters.pdf")
	
	
clust_labels = ['means_per_cluster_mu_fg_Inh L4-5 SST STK32A', 'means_per_cluster_mu_fg_Inh L4-6 PVALB SULF1', 'means_per_cluster_mu_fg_Inh L4-6 SST B3GAT2', 'means_per_cluster_mu_fg_Inh L4-6 SST GXYLT2', 
				  'means_per_cluster_mu_fg_Inh L5-6 PVALB LGR5', 'means_per_cluster_mu_fg_Inh L5-6 SST MIR548F2', 'means_per_cluster_mu_fg_Inh L5-6 SST NPM1P10'
                        ]
clust_col = ['' + str(i) for i in clust_labels] # in case column names differ from labels

slide = select_slide(adata_vis, 'DH4a')

with mpl.rc_context({'figure.figsize': (20, 20)}):
    fig = plot_spatial(
        adata=slide,
        # labels to show on a plot
        color=clust_col, labels=clust_labels,
		alpha_scaling=1.0,
        show_img=True,
		img_alpha=0.0,
		adjust_text=True,
        # 'fast' (white background) or 'dark_background'
        style='dark_background',
        # limit color scale at 99.2% quantile of cell abundance
        max_color_quantile=0.992,
        # size of locations (adjust depending on figure size)
        circle_diameter=9,
        colorbar_position='right'
    ), plt.savefig("DH4a_31yr_Inh_L5_6_clusters.pdf")
	
	
clust_labels = ['means_per_cluster_mu_fg_Micro L1-3 TYROBP', 'means_per_cluster_mu_fg_Endo L2-6 NOSTRIN'
               ]
clust_col = ['' + str(i) for i in clust_labels] # in case column names differ from labels

slide = select_slide(adata_vis, 'DH4a')

with mpl.rc_context({'figure.figsize': (20, 20)}):
    fig = plot_spatial(
        adata=slide,
        # labels to show on a plot
        color=clust_col, labels=clust_labels,
		alpha_scaling=1.0,
        show_img=True,
		img_alpha=0.0,
		adjust_text=True,
        # 'fast' (white background) or 'dark_background'
        style='dark_background',
        # limit color scale at 99.2% quantile of cell abundance
        max_color_quantile=0.992,
        # size of locations (adjust depending on figure size)
        circle_diameter=9,
        colorbar_position='right'
    ), plt.savefig("DH4a_31yr_Endo_clusters.pdf")
	
clust_labels = ['means_per_cluster_mu_fg_OPC L1-6 PDGFRA', 'means_per_cluster_mu_fg_Oligo L1-6 OPALIN', 'means_per_cluster_mu_fg_Inh L3-6 SST NPY'	
               ]
clust_col = ['' + str(i) for i in clust_labels] # in case column names differ from labels

slide = select_slide(adata_vis, 'DH4a')

with mpl.rc_context({'figure.figsize': (20, 20)}):
    fig = plot_spatial(
        adata=slide,
        # labels to show on a plot
        color=clust_col, labels=clust_labels,
		alpha_scaling=1.0,
        show_img=True,
		img_alpha=0.0,
		adjust_text=True,
        # 'fast' (white background) or 'dark_background'
        style='dark_background',
        # limit color scale at 99.2% quantile of cell abundance
        max_color_quantile=0.992,
        # size of locations (adjust depending on figure size)
        circle_diameter=9,
        colorbar_position='right'
    ), plt.savefig("DH4a_31yr_Oligo_clusters.pdf")
	
# compute KNN using the cell2location output stored in adata.obsm
sc.pp.neighbors(adata_vis, use_rep='q05_cell_abundance_w_sf',
                n_neighbors = 21)

# Cluster spots into regions using scanpy
sc.tl.leiden(adata_vis, resolution=0.8)

# add region as categorical variable
adata_vis.obs["region_cluster"] = adata_vis.obs["leiden"].astype("category")
# compute UMAP using KNN graph based on the cell2location output
sc.tl.umap(adata_vis, min_dist = 0.3, spread = 1)	
	
with mpl.rc_context({'axes.facecolor':  'white',
                     'figure.figsize': [8, 8]}):
    sc.pl.umap(adata_vis, color=['region_cluster'], size=30,
               color_map = 'RdPu', ncols = 2, legend_loc='on data',
               legend_fontsize=20, save="umap_region_cluster_n21_res0.8.pdf")
    sc.pl.umap(adata_vis, color=['sample'], size=30,
               color_map = 'RdPu', ncols = 2,
               legend_fontsize=20, save="umap_all_the_samples_res0.8.pdf")	

slide = select_slide(adata_vis, 'DH1')
			   
with mpl.rc_context({'axes.facecolor':  'white',
                     'figure.figsize': [8, 8]}):
    sc.pl.umap(slide, color=['region_cluster'], size=30,
               color_map = 'RdPu', ncols = 2, legend_loc='on data',
               legend_fontsize=20, save="DH1_4yr_cluster_n21_res0.8.pdf")
    sc.pl.umap(slide, color=['sample'], size=30,
               color_map = 'RdPu', ncols = 2,
               legend_fontsize=20, save="DH1_4yr_umap_sample_res_n21_0.8.pdf")
			   
with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):
    sc.pl.spatial(slide, color=['region_cluster'],
                  size=1.3, img_key='hires', alpha=0.5, save="DH1_spatial_cluster_region_4yr_n21_res0.8.pdf")

slide = select_slide(adata_vis, 'DH2')

with mpl.rc_context({'axes.facecolor':  'white',
                     'figure.figsize': [8, 8]}):
    sc.pl.umap(slide, color=['region_cluster'], size=30,
               color_map = 'RdPu', ncols = 2, legend_loc='on data',
               legend_fontsize=20, save="DH2_4yr_umap_region_cluster_n21_res0.8.pdf")
    sc.pl.umap(slide, color=['sample'], size=30,
               color_map = 'RdPu', ncols = 2,
               legend_fontsize=20, save="DH2_umap_sample_n21_res0.8.pdf")
			   
with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):
    sc.pl.spatial(slide, color=['region_cluster'],
                  size=1.3, img_key='hires', alpha=0.5, save="DH2_spatial_cluster_region_4yr_n21_res0.8.pdf")

slide = select_slide(adata_vis, 'DH3')

with mpl.rc_context({'axes.facecolor':  'white',
                     'figure.figsize': [8, 8]}):
    sc.pl.umap(slide, color=['region_cluster'], size=30,
               color_map = 'RdPu', ncols = 2, legend_loc='on data',
               legend_fontsize=20, save="DH3_2nd_15yr_umap_region_cluster_n21_res0.8.pdf")
    sc.pl.umap(slide, color=['sample'], size=30,
               color_map = 'RdPu', ncols = 2,
               legend_fontsize=20, save="DH3_2nd_15yr_umap_sample_n21_res0.8.pdf")
			   
with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):
    sc.pl.spatial(slide, color=['region_cluster'],
                  size=1.3, img_key='hires', alpha=0.5, save="DH3_spatial_cluster_region_2nd_15yr_n21_res0.8.pdf")

slide = select_slide(adata_vis, 'DH4')

with mpl.rc_context({'axes.facecolor':  'white',
                     'figure.figsize': [8, 8]}):
    sc.pl.umap(slide, color=['region_cluster'], size=30,
               color_map = 'RdPu', ncols = 2, legend_loc='on data',
               legend_fontsize=20, save="DH4_2nd_15yr_umap_region_cluster_n21_res0.8.pdf")
    sc.pl.umap(slide, color=['sample'], size=30,
               color_map = 'RdPu', ncols = 2,
               legend_fontsize=20, save="DH4_2nd_15yr_umap_sample_n21_res0.8.pdf")
			   
with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):
    sc.pl.spatial(slide, color=['region_cluster'],
                  size=1.3, img_key='hires', alpha=0.5, save="DH4_spatial_cluster_region_2nd_15yr_n21_res0.8.pdf")
				  
slide = select_slide(adata_vis, 'DH1a')

with mpl.rc_context({'axes.facecolor':  'white',
                     'figure.figsize': [8, 8]}):
    sc.pl.umap(slide, color=['region_cluster'], size=30,
               color_map = 'RdPu', ncols = 2, legend_loc='on data',
               legend_fontsize=20, save="DH1a_first15yr_umap_region_cluster_n21_res0.8.pdf")
    sc.pl.umap(slide, color=['sample'], size=30,
               color_map = 'RdPu', ncols = 2,
               legend_fontsize=20, save="DH1a_first15yr_umap_sample_n21_res0.8.pdf")
			   
with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):
    sc.pl.spatial(slide, color=['region_cluster'],
                  size=1.3, img_key='hires', alpha=0.5, save="DH1a_spatial_cluster_region_first15yr_n21_res0.8.pdf")

slide = select_slide(adata_vis, 'DH2a')

with mpl.rc_context({'axes.facecolor':  'white',
                     'figure.figsize': [8, 8]}):
    sc.pl.umap(slide, color=['region_cluster'], size=30,
               color_map = 'RdPu', ncols = 2, legend_loc='on data',
               legend_fontsize=20, save="DH2a_first15yr_umap_region_cluster_n21_res0.8.pdf")
    sc.pl.umap(slide, color=['sample'], size=30,
               color_map = 'RdPu', ncols = 2,
               legend_fontsize=20, save="DH2a_first15yr_umap_sample_n21_res0.8.pdf")
			   
with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):
    sc.pl.spatial(slide, color=['region_cluster'],
                  size=1.3, img_key='hires', alpha=0.5, save="DH2_spatial_cluster_region_first15yr_n21_res0.8.pdf")

slide = select_slide(adata_vis, 'DH3a')

with mpl.rc_context({'axes.facecolor':  'white',
                     'figure.figsize': [8, 8]}):
    sc.pl.umap(slide, color=['region_cluster'], size=30,
               color_map = 'RdPu', ncols = 2, legend_loc='on data',
               legend_fontsize=20, save="DH3a_31yr_umap_region_cluster_n21_res0.8.pdf")
    sc.pl.umap(slide, color=['sample'], size=30,
               color_map = 'RdPu', ncols = 2,
               legend_fontsize=20, save="DH3a_31yr_umap_sample_n21_res0.8.pdf")
			   
with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):
    sc.pl.spatial(slide, color=['region_cluster'],
                  size=1.3, img_key='hires', alpha=0.5, save="DH3_spatial_cluster_region_31yr_n21_res0.8.pdf")

slide = select_slide(adata_vis, 'DH4a')

with mpl.rc_context({'axes.facecolor':  'white',
                     'figure.figsize': [8, 8]}):
    sc.pl.umap(slide, color=['region_cluster'], size=30,
               color_map = 'RdPu', ncols = 2, legend_loc='on data',
               legend_fontsize=20, save="DH4a_31yr_umap_region_cluster_n21_res0.8.pdf")
    sc.pl.umap(slide, color=['sample'], size=30,
               color_map = 'RdPu', ncols = 2,
               legend_fontsize=20, save="DH4a_31yr_umap_sample_n21_res0.8.pdf")
			   
with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):
    sc.pl.spatial(slide, color=['region_cluster'],
                  size=1.3, img_key='hires', alpha=0.5, save="DH4_spatial_cluster_region_31yr_n21_res0.8.pdf")


from cell2location import run_colocation
res_dict, adata_vis = run_colocation(
    adata_vis,
    model_name='CoLocatedGroupsSklearnNMF',
    train_args={
      'n_fact': np.arange(5, 30), # IMPORTANT: use a wider range of the number of factors (5-30)
      'sample_name_col': 'sample', # columns in adata_vis.obs that identifies sample
      'n_restarts': 3 # number of training restarts
    },
    export_args={'path': f'{run_name}/CoLocatedComb/'}
)

sam = np.array(adata_vis.obs['sample'])
for i in np.unique(sam):
    
    s1 = adata_vis.obs[['region_cluster']]
    s1 = s1.loc[sam == i]
    s1.index = [x[10:] for x in s1.index]
    s1.index.name = 'Barcode'
    
    s1.to_csv(f'{run_name}/region_cluster29_ + i + .csv')
	
adata_file = f"{run_name}/sp_with_clusters.h5ad"
adata_vis.write(adata_file)
adata_file
