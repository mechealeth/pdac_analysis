#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 27 17:06:44 2022

@author: mecheal
"""

myloid.obs['label'] = celltype_labels
TNK_subset.obs['label'] = TNK_subset.obs['label'].astype(str) 
pdac_combined.obs['label1'] = pdac_combined.obs['label1'].astype(str)
pdac_combined.obs['label1'].update(TNK_subset.obs['label'])

import scanpy as sc
import infercnvpy as cnv
import matplotlib.pyplot as plt

gtf_path = /lustre/home/whhou/00.datasets/10x/refdata-cellranger-GRCh38-3.0.0/genes/genes.gtf

cnv.io.genomic_position_from_gtf(gtf_file='/lustre/home/whhou/00.datasets/10x/refdata-cellranger-GRCh38-3.0.0/genes/genes.gtf', adata=pdac_combined)
pdac_combined.var.loc[:, [ "chromosome", "start", "end"]].head()

exclude_chromosomes = ['chrnan','chrMT','chrKI270734.1','chrKI270726.1','chrGL000194.1','chrGL000219.1','chrGL000195.1','chrGL000219.1','chrGL000009.2',
                      'chrGL000218.1','chrKI270721.1','chrKI270728.1','chrKI270711.1' ,'chrKI270713.1','chrY','chrX']
pdac_combined.var['chromosome'] = ['chr'+str(i) for i in pdac_combined.var['chromosome']]

cnv.tl.infercnv(
    pdac_combined,
    reference_key="label1",
    reference_cat=[
        "CD4_Tn",
        "MAIT",
        "Bcell",
        "CD8_Tn",
        "CD8_Tem_TNFSF9",
        "Endothelial",
        "CD4_Tem_CREM",
        "Mast",
        "CD4_Tfh",
        "Fib_COL11A1",
        "Plasma",
        "CD8_Tex",
        "Fib_COL11A1",
        "SMC_CCL19"
    ],
    window_size=150,
    exclude_chromosomes=exclude_chromosomes
)

cnv.pl.chromosome_heatmap(DA, groupby="leiden",save = 'duc_cnV')
cnv.tl.pca(DA)
cnv.pp.neighbors(DA,n_neighbors=8, n_pcs=20)
cnv.tl.umap(DA, min_dist=0.4)
cnv.tl.leiden(DA,resolution = 0.2)
cnv.tl.cnv_score(DA)
cnv.pl.umap(DA,color= ['cnv_score','cnv_leiden'],save ='umap_cnv')
cnv.pl.chromosome_heatmap(DA, groupby="leiden", dendrogram=True,save ='DA_leiden')
sc.pl.umap(Ductal_Acinar,color= ['cnv_score','cnv_leiden','leiden'],save ='ad')

sc.pl.violin(Ductal_Acinar, ['cnv_score'], groupby='cnv_leiden',save = 'vio' )

Ductal_Acinar.uns['cnv'] = DA.uns['cnv']

cancer = Ductal_Acinar[Ductal_Acinar.obs.cnv_leiden.isin(['1','5']), :].copy()
DA.obs['leiden'] = Ductal_Acinar.obs['leiden']