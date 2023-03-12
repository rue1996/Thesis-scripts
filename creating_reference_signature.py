#4th-step
#Creating reference signature using the Hockman_lab snRNAseq datasets
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

fpath = "/scratch/mshruv003/cell2loc_integrated_with_Allen_MTG"
ref_run_name = f'{fpath}/reference_signature_cluster_label_batch'


adata_ref = sc.read("/scratch/mshruv003/cell2loc_integrated_with_Allen_MTG/Allen_MTG_seurat.h5ad")

from cell2location.utils.filtering import filter_genes
selected = filter_genes(adata_ref, cell_count_cutoff=5, cell_percentage_cutoff2=0.03, nonz_mean_cutoff=1.12)

# filter the object
adata_ref = adata_ref[:, selected].copy()

scvi.data.setup_anndata(adata=adata_ref,
                        # 10X reaction / sample / batch
                        batch_key='sample',
                        # cell type, covariate used for constructing signatures
                        labels_key='Allen_high_res',
                        # multiplicative technical effects (platform, 3' vs 5', donor effect)
                        categorical_covariate_keys=['single_cell_chemistry']
                       )
scvi.data.view_anndata_setup(adata_ref)

from cell2location.models import RegressionModel
mod = RegressionModel(adata_ref)

mod.train(max_epochs=250, batch_size=2500, train_size=1, lr=0.002, use_gpu=False)
mod.plot_history(20)

adata_ref = mod.export_posterior(
    adata_ref, sample_kwargs={'num_samples': 1000, 'batch_size': 2500, 'use_gpu': False}
)
mod.save(f"{ref_run_name}", overwrite=True)
adata_file = f"{ref_run_name}/sc.h5ad"
adata_ref.__dict__['_raw'].__dict__['_var'] = adata_ref.__dict__['_raw'].__dict__['_var'].rename(columns={'_index': 'features'})
adata_ref.write(adata_file)
adata_file



