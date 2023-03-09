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
results_folder = '/scratch/mshruv003/cell2loc_integrated_with_Allen_MTG/'
sp_results_folder = f'{results_folder}'
sc_results_folder = f'{results_folder}reference_signature_cluster_label/'

run_name = '/scratch/mshruv003/cell2loc_integrated_with_Allen_MTG/cell2location_map'

sp_data_file = '/scratch/mshruv003/cell2loc_integrated_with_Allen_MTG/cell2location_map/sp_with_clusters_n21_res08.h5ad'
adata_vis2 = anndata.read(sp_data_file)


def unpickle_model(path, mod_name):
    r""" Unpickle model
    """
    file = path + 'model_' + mod_name + ".p"
    
    mod1_ann = pickle.load(file = open(file, "rb"))
    return mod1_ann['mod']

n_fact = 15
mod_path = f'{results_folder}/cell2location_map/CoLocatedComb/CoLocatedGroupsSklearnNMF_27703locations_57factors/models/'
adata_file = f'{results_folder}/cell2location_map/CoLocatedComb/CoLocatedGroupsSklearnNMF_27703locations_57factors/anndata/sp.h5ad'


mod_sk = unpickle_model(mod_path, f'n_fact{n_fact}')

adata_vis_sk = anndata.read(adata_file)

adata_snrna_raw = sc.read_h5ad('/scratch/mshruv003/cell2loc_integrated_with_Allen_MTG/reference_signature_cluster_label_batch/sc.h5ad')

adata_vis_sk.uns['mod']['factor_names']

b_dev_sel = [
       'means_per_cluster_mu_fg_Astro L1-2 FGFR3 GFAP',
       'means_per_cluster_mu_fg_Astro L1-6 FGFR3 SLC14A1',
       'means_per_cluster_mu_fg_Endo L2-6 NOSTRIN',
       'means_per_cluster_mu_fg_Exc L2 LAMP5 LTK',
       'means_per_cluster_mu_fg_Exc L2-3 LINC00507 FREM3',
       'means_per_cluster_mu_fg_Exc L3-5 RORB COL22A1',
       'means_per_cluster_mu_fg_Exc L3-5 RORB ESR1',
       'means_per_cluster_mu_fg_Exc L3-5 RORB FILIP1L',
       'means_per_cluster_mu_fg_Exc L3-5 RORB TWIST2',
       'means_per_cluster_mu_fg_Exc L4-5 RORB DAPK2',
       'means_per_cluster_mu_fg_Exc L4-5 RORB FOLH1B',
       'means_per_cluster_mu_fg_Exc L4-6 FEZF2 IL26',
       'means_per_cluster_mu_fg_Exc L4-6 RORB C1R',
       'means_per_cluster_mu_fg_Exc L4-6 RORB SEMA3E',
       'means_per_cluster_mu_fg_Exc L5-6 FEZF2 ABO',
       'means_per_cluster_mu_fg_Exc L5-6 FEZF2 EFTUD1P1',
       'means_per_cluster_mu_fg_Exc L5-6 THEMIS C1QL3',
       'means_per_cluster_mu_fg_Exc L5-6 THEMIS FGF10',
       'means_per_cluster_mu_fg_Inh L1 SST CHRNA4',
       'means_per_cluster_mu_fg_Inh L1 SST NMBR',
       'means_per_cluster_mu_fg_Inh L1-2 GAD1 MC4R',
       'means_per_cluster_mu_fg_Inh L1-2 PAX6 CDH12',
       'means_per_cluster_mu_fg_Inh L1-2 SST BAGE2',
       'means_per_cluster_mu_fg_Inh L1-2 VIP LBH',
       'means_per_cluster_mu_fg_Inh L1-2 VIP PCDH20',
       'means_per_cluster_mu_fg_Inh L1-2 VIP TSPAN12',
       'means_per_cluster_mu_fg_Inh L1-3 SST CALB1',
       'means_per_cluster_mu_fg_Inh L1-3 VIP ADAMTSL1',
       'means_per_cluster_mu_fg_Inh L1-3 VIP CCDC184',
       'means_per_cluster_mu_fg_Inh L1-3 VIP CHRM2',
       'means_per_cluster_mu_fg_Inh L1-3 VIP GGH',
       'means_per_cluster_mu_fg_Inh L1-4 LAMP5 LCP2',
       'means_per_cluster_mu_fg_Inh L1-4 VIP CHRNA6',
       'means_per_cluster_mu_fg_Inh L1-4 VIP OPRM1',
       'means_per_cluster_mu_fg_Inh L2-3 VIP CASC6',
       'means_per_cluster_mu_fg_Inh L2-4 PVALB WFDC2',
       'means_per_cluster_mu_fg_Inh L2-4 SST FRZB',
       'means_per_cluster_mu_fg_Inh L2-4 VIP CBLN1',
       'means_per_cluster_mu_fg_Inh L2-4 VIP SPAG17',
       'means_per_cluster_mu_fg_Inh L2-5 PVALB SCUBE3',
       'means_per_cluster_mu_fg_Inh L2-5 VIP SERPINF1',
       'means_per_cluster_mu_fg_Inh L2-5 VIP TYR',
       'means_per_cluster_mu_fg_Inh L2-6 LAMP5 CA1',
       'means_per_cluster_mu_fg_Inh L2-6 VIP QPCT',
       'means_per_cluster_mu_fg_Inh L3-5 SST ADGRG6',
       'means_per_cluster_mu_fg_Inh L3-6 VIP HS3ST3A1',
       'means_per_cluster_mu_fg_Inh L4-5 SST STK32A',
       'means_per_cluster_mu_fg_Inh L4-6 PVALB SULF1',
       'means_per_cluster_mu_fg_Inh L4-6 SST B3GAT2',
       'means_per_cluster_mu_fg_Inh L4-6 SST GXYLT2',
       'means_per_cluster_mu_fg_Inh L5-6 PVALB LGR5',
       'means_per_cluster_mu_fg_Inh L5-6 SST MIR548F2',
       'means_per_cluster_mu_fg_Inh L5-6 SST NPM1P10',
       'means_per_cluster_mu_fg_Micro L1-3 TYROBP',
       'means_per_cluster_mu_fg_NA',
       'means_per_cluster_mu_fg_OPC L1-6 PDGFRA',
       'means_per_cluster_mu_fg_Oligo L1-6 OPALIN']



from cell2location.cluster_averages.cluster_averages import get_cluster_averages_df
from cell2location.plt.plot_heatmap import clustermap
#plt.rcParams['xtick.bottom'] = plt.rcParams['xtick.labelbottom'] = True
#plt.rcParams['xtick.top'] = plt.rcParams['xtick.labeltop'] = False

ct_aver = get_cluster_averages_df(X=adata_vis.obs[['' + i 
                                                      for i in adata_vis.uns['mod']['factor_names']]],
                        cluster_col=adata_vis.obs["region_cluster"])
ct_aver.index = adata_vis.uns['mod']['factor_names']
ct_aver.columns = ['region_' + c for c in ct_aver.columns]

# normalise to get 10% of each cell type in each location
ct_aver = (ct_aver.T / ct_aver.max(1)).T

rcParams['pdf.fonttype'] = 42 # enables correct plotting of text
clustermap(ct_aver.loc[b_dev_sel, :],
           cluster_rows=False, cluster_cols=True, 
           figure_size=[5.9 + 0.12 * mod_sk.n_fact, 5.9 + 0.1 * mod_sk.n_var],
           fun_type='dotplot', array_size=None)

plt.savefig(Fig4D_suppl_cluster_dotplot.pdf', bbox_inches='tight')
plt.show()