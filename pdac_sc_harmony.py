#!/usr/bin/env python3

import numpy as np 
import pandas as pd
import scanpy as sc
import scanpy.external as sce
import pdac_utils
import bbknn
import sys
sys.path.append('/lustre/home/mcchen/PDAC/process/PDAC/script')


sc.settings.verbosity = 1
#saving path
sc.settings.figdir = '/lustre/home/mcchen/PDAC/process/PDAC/fig1/'
#figure paras
sc.settings.set_figure_params(dpi=100, fontsize=10, dpi_save=400,facecolor = 'white', figsize=(8,8), format='png')
#data read in
pdac_pids =["PDAC01_Tumor_GEX","PDAC02_Tumor_GEX","PDAC03_Tumor_GEX","PDAC04_Normal_GEX","PDAC04_Tumor_GEX","PDAC05_Tumor_GEX","PDAC06_Normal_GEX","PDAC06_Tumor_GEX",
          "PDAC07_Tumor_GEX","PDAC08_Normal_GEX","PDAC08_Tumor_GEX","PDAC09_Normal_GEX","PDAC09_Tumor_GEX","PDAC10_Normal_GEX","PDAC10_Tumor_GEX","PDAC11_Tumor_GEX"]
Normal_pids =["PDAC04_Normal_GEX","PDAC06_Normal_GEX","PDAC08_Normal_GEX","PDAC09_Normal_GEX","PDAC10_Normal_GEX"]
Tumor_pids =["PDAC01_Tumor_GEX","PDAC02_Tumor_GEX","PDAC03_Tumor_GEX","PDAC04_Tumor_GEX","PDAC05_Tumor_GEX","PDAC06_Tumor_GEX","PDAC07_Tumor_GEX","PDAC08_Tumor_GEX",
             "PDAC09_Tumor_GEX","PDAC10_Tumor_GEX","PDAC11_Tumor_GEX"]
pdac_adatas = []
for pid in pdac_pids:
    adata = sc.read_10x_h5('/lustre/home/whhou/00.datasets/2021_Data/process/mapping/PDAC/RNA/%s/%s/outs/filtered_feature_bc_matrix.h5'%(pid,pid))
    ncells, ngenes = adata.shape 
    adata.obs['n_counts'] = adata.X.sum(1)
    adata.obs['n_genes'] = (adata.X > 0).sum(1)
    adata.obs['pid'] = [pid]*adata.shape[0]
    if pid in Tumor_pids:
        adata.obs['type'] = ['pdac']*adata.shape[0]
    else:
        adata.obs['type'] = ['normal']*adata.shape[0]
    adata.var_names_make_unique()
    pdac_utils.calculate_qc(adata)
    pdac_utils.preprocess(adata)
    pdac_adatas.append(adata)    
#adding pdac12
adata = sc.read_10x_h5('/lustre/home/mcchen/data/PDAC/PDAC12_Tumor_GEX/PDAC12-Tumor-GEX/outs/filtered_feature_bc_matrix.h5')
ncells, ngenes = adata.shape
adata.obs['n_counts'] = adata.X.sum(1)
adata.obs['n_genes'] = (adata.X > 0).sum(1)
adata.obs['pid'] = ["PDAC12_Tumor_GEX"]*adata.shape[0]
adata.obs['type'] = ['pdac']*adata.shape[0]
adata.var_names_make_unique()
pdac_utils.calculate_qc(adata)
pdac_utils.preprocess(adata)
pdac_adatas.append(adata)
#combined and save
pdac_combined = pdac_adatas[0].concatenate(pdac_adatas[1:]) 
pdac_utils.basic_filter(pdac_combined)
pdac_utils.normalization(pdac_combined)
sc.pp.highly_variable_genes(pdac_combined,flavor = 'seurat', min_mean=0.0125, max_mean=3, min_disp=0.25,batch_key='pid')
pdac_combined.raw = pdac_combined
sc.pp.scale(pdac_combined, max_value=10)
sc.tl.pca(pdac_combined, svd_solver='arpack')
sce.pp.harmony_integrate(pdac_combined, 'pid')
pdac_combined.obsm['X_pca'] = pdac_combined.obsm['X_pca_harmony']
pdac_combined.write('/lustre/home/mcchen/PDAC/process/PDAC/obj/pdac_combined_harmony.h5ad')
#neighbor/umap
pdac_utils.analysis_data(pdac_combined,1.6)
pdac_utils.analysis_plots(pdac_combined)
pdac_combined.write('/lustre/home/mcchen/PDAC/process/PDAC/obj/pdac_combined_harmony.h5ad')
#marker_gene
genemarkers = {}
filename="/lustre/home/mcchen/ref/pdac_label1.txt"
f = open(filename)
for line in f:
    tokens = line.strip().split('\t')
    if len(tokens) > 1:
        markers = [gene.upper() for gene in tokens[1:]]
        markers = [x for x in markers if x in  pdac_combined.raw.var_names]
        genemarkers[tokens[0]] = markers
sc.pl.dotplot(pdac_combined, genemarkers, groupby ='leiden_harmony',dendrogram=True,save = 'label1_umap')
immune = {}
filename="/lustre/home/mcchen/ref/immune_label.txt"
f = open(filename)
for line in f:
    tokens = line.strip().split(',')
    if len(tokens) > 1:
        markers = [gene.upper() for gene in tokens[1:]]
        markers = [x for x in markers if x in  pdac_combined.raw.var_names]
        immune[tokens[0]] = markers
sc.pl.dotplot(pdac_combined, immune, groupby ='leiden_harmony', dendrogram=True,save = 'label_immune_umap')
                 