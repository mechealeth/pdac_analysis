#!/usr/bin/env python3


import numpy as np
import pandas as pd
import scanpy as sc
import scrublet as scr
from anndata import AnnData
from collections import Counter
from matplotlib import pyplot as plt
import scanpy.external as sce


def load_genemarkers():
    genemarkers = {}
    filename="/lustre/home/mcchen/ref/pdac_marker_gene.txt"
    f = open(filename)
    for line in f:
        tokens = line.strip().split('\t')
        if len(tokens) > 1:
            markers = [gene.upper() for gene in tokens[1:]]
            genemarkers[tokens[0]] = markers
    return genemarkers


def cell_cycle_score(adata):
    cell_cycle_genes = [x.strip() for x in open('/lustre/home/mcchen/ref/regev_lab_cell_cycle_genes.txt')]
    s_genes = cell_cycle_genes[:43]
    g2m_genes = cell_cycle_genes[43:]
    cell_cycle_genes = [x for x in cell_cycle_genes if x in adata.var_names]
    sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)
    return adata


def calculate_qc(adata):
    scrub = scr.Scrublet(adata.X)
    doublet_scores, predicted_doublets = scrub.scrub_doublets()
    adata.obs['scrublet_scores'] = doublet_scores
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    return adata


def qc_plot(adata):
    sc.pl.scatter(adata, 'n_counts', 'n_genes', color='pct_counts_mt')
    sc.pl.scatter(adata, 'n_counts', 'n_genes', color='scrublet_scores')
    
    

def basic_filter(adata):
    sc.pp.filter_genes(adata, min_cells=3)
    sc.pp.filter_cells(adata, min_genes=200)
    adata = adata[adata.obs['scrublet_scores'] <= 0.4,:].copy() 
    adata = adata[adata.obs.n_genes_by_counts < 5000,:].copy()
    adata = adata[adata.obs.pct_counts_mt < 20,:].copy
    return adata

def normalization(adata):
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    return adata

def hvg_scale(adata):
    sc.pp.highly_variable_genes(adata, flavor = 'seurat', min_mean=0.0125, max_mean=3, min_disp=0.25,batch_key='pid')
    sc.pp.scale(adata, max_value=10)
    return adata

def preprocess(adata):
    qc_plot(adata)
    basic_filter(adata)
    normalization(adata)
    adata.raw = adata
    hvg_scale(adata)
    return adata

def analysis_data(adata,resolution = None):
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
    sc.tl.umap(adata, min_dist=0.5)
    if resolution:
        sc.tl.leiden(adata, resolution=resolution,key_added='leiden')
    else:
        sc.tl.leiden(adata, resolution,key_added='leiden')
    sc.tl.rank_genes_groups(adata,groupby = 'leiden',key_added='rank_genes_groups', method='t-test', use_raw=True)
    return adata
                
def analysis_plots(adata):
    sc.pl.pca(adata, color=['pct_counts_mt', 'n_genes'], cmap="jet",save = 'pca_qc')
    sc.pl.pca_variance_ratio(adata, n_pcs = 50,log=True, save = 'pca_variance_ratio')
    sc.pl.umap(adata, color=['n_genes', 'scrublet_scores'], save = '_gene_doublet')
    sc.pl.umap(adata, color=['pid'],  save = '_pid_leiden_umap')
    sc.pl.umap(adata, color=['leiden'],  save = '_pid_leiden')
    sc.pl.rank_genes_groups(adata,key='rank_genes_groups', n_genes=25, sharey=False, save = '_rank_gene_group')
    return adata


        
    







   
    

           
            
    