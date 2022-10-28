#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 20 09:25:51 2022

@author: mecheal
"""

import scanpy as sc
import scanpy.external as sce
import numpy as np

import matplotlib as mpl
from copy import copy
Greens = copy(mpl.cm.Greens)
Greens.set_under("lightgray")

sc.settings.verbosity = 1
#saving path
sc.settings.figdir = '/lustre/home/mcchen/PDAC/process/PDAC/fig2/'
#figure paras
sc.settings.set_figure_params(dpi=100, fontsize=10, dpi_save=200,facecolor = 'white', figsize=(8,8), format='png')


pdac_combined = sc.read('/lustre/home/mcchen/PDAC/process/PDAC/obj/pdac_combined1010.h5ad')

Myeloid = pdac_combined[pdac_combined.obs.label1.isin(['Myeloid']), :].copy()
sc.pp.normalize_total(Myeloid, target_sum=1e4)
sc.pp.log1p(Myeloid)
Myeloid.raw = Myeloid.copy()
My =Myeloid.raw.to_adata()


#calculate new hvg
My.uns['log1p']["base"] = None 
sc.pp.highly_variable_genes(My,  min_mean=0.0125, max_mean=3, min_disp=0.25,batch_key='pid')
remove =  Myeloid.var_names[Myeloid.var_names.str.startswith(("RPS","RPL","MT-",'MALAT1'))]
hvg = [x for x in My.var['highly_variable'][My.var.highly_variable].index.values if x not in remove]
My.var['highly_variable'] = False
My.var.loc[My.var_names.isin(hvg), 'highly_variable'] = True

#
sc.pp.regress_out(My, ['total_counts', 'pct_counts_mt'])
sc.pp.scale(My, max_value=10)
sc.tl.pca(My, svd_solver='arpack' )
sce.pp.harmony_integrate(My, 'pid')
My.obsm['X_pca'] =My.obsm['X_pca_harmony']
sc.pp.neighbors(My, n_neighbors=10, n_pcs=30)
sc.tl.umap(My, min_dist=0.4)
sc.tl.leiden(My,0.7,key_added='leiden')
sc.pl.umap(My, color=['leiden'],legend_loc = 'on data',save ='umap_My')

#transfer
Myeloid.uns = My.uns
Myeloid.obsm = My.obsm
Myeloid.varm = My.varm
Myeloid.obsp = My.obsp
Myeloid.obs['leiden'] =My.obs['leiden']


#check
sc.pl.umap(Myeloid, color=['CD14','ITGAM','MNDA','MPEG1','ITGAX'],frameon=False,color_map=Greens,vmin=0.00001,save ='umap_leidenc1')
sc.pl.umap(Myeloid, color=['leiden'],legend_loc = 'on data',save ='umap_my')
#Mono
sc.pl.umap(Myeloid, color=['S100A9', 'LYZ', 'FCN1'],frameon=False,color_map=Greens,vmin=0.00001,save ='mono')
#Macro
sc.pl.umap(Myeloid, color=['C1QA', 'CD68', 'TREM2'],frameon=False,color_map=Greens,vmin=0.00001,save ='macro')
#DC
sc.pl.umap(Myeloid, color=['IRF7', 'HLA-DRA', 'LYZ', 'CST3'],frameon=False,color_map=Greens,vmin=0.00001,save ='dc')

sc.pl.umap(My, color=['CTSD'],frameon=False,color_map=Greens,vmin=0.00001,save ='ctsd')

myloid_label = {'0':'Myeloid',
              '1':'Myeloid',
              '2':'Myeloid',
              '3':'Myeloid',
              '4':'Myeloid',
              '5':'Myeloid',
              '6':'Myeloid',
              '7':'Myeloid',
              '8':'Myeloid',
              '9':'Myeloid',
              '10':'Myeloid',
              '11':'Myeloid' ,
              '12':'uncharacterize',
              '13':'Myeloid',
              '14':'Myeloid',
              '15':'Myeloid', 
     }

celltype_labels = []
for leiden in Myeloid.obs.leiden:
    celltype_labels.append(myloid_label.get(leiden))
Myeloid.obs['label'] = celltype_labels
Myeloid.obs['label'] = Myeloid.obs['label'].astype(str)  
pdac_combined.obs['label1'] = pdac_combined.obs['label1'].astype(str)
pdac_combined.obs['label1'].update(Myeloid.obs['label'])    

sc.tl.leiden(Myeloid, restrict_to=('leiden', ['10']), resolution=0.4, key_added='leiden')

genemarkers = {}
filename="/lustre/home/mcchen/ref/My_marker_gene.txt"
f = open(filename)
for line in f:
     tokens = line.strip().split('\t')
     if len(tokens) > 1:
         markers = [gene.upper() for gene in tokens[1:]]
         markers = [x for x in markers if x in Myeloid.var_names]
         genemarkers[tokens[0]] = markers
sc.tl.dendrogram(Myeloid, groupby ='leiden')        
sc.pl.dotplot(Myeloid, genemarkers, groupby ='leiden',dendrogram=True,save = 'my_marker')


genemarkers = {}
filename="/lustre/home/mcchen/ref/mye_marker.txt"
f = open(filename)
for line in f:
     tokens = line.strip().split(',')
     if len(tokens) > 1:
         markers = [gene.upper() for gene in tokens[1:]]
         markers = [x for x in markers if x in Myeloid.var_names]
         genemarkers[tokens[0]] = markers
sc.tl.dendrogram(Myeloid, groupby ='leiden')        
sc.pl.dotplot(Myeloid, genemarkers, groupby ='leiden',dendrogram=True,save = 'myc_marker')



#plot
Myeloid.uns['log1p']["base"] = None 
sc.tl.rank_genes_groups(Myeloid,groupby = 'leiden',key_added='rank_genes_groups1',pts = True, method='t-test', use_raw=True)
sc.pl.rank_genes_groups(Myeloid,key='rank_genes_groups1', n_genes=30, sharey=False, save = '_rank_gene_group')
sc.pl.umap(Myeloid, color=['n_genes', 'pct_counts_mt','scrublet_scores'], save = '_gene_doublet')
sc.pl.umap(Myeloid, color=['pid'],  save = '_pid_leiden_umap')
sc.pl.umap(Myeloid, color=['leiden'],save ='umap_leiden')
sc.pl.umap(Myeloid, color=['leiden'], legend_loc ='on data',save ='umap_leiden_num')

Myeloid.uns['log1p']["base"] = None 
sc.tl.rank_genes_groups(Myeloid,groupby = 'label',key_added='rank_genes_groups',pts = True, method='t-test', use_raw=True)
sc.pl.rank_genes_groups(Myeloid,key='rank_genes_groups', n_genes=30, sharey=False, save = '_rank_gene_group')

de12= sc.get.rank_genes_groups_df(Myeloid, group='12', key='rank_genes_groups1', pval_cutoff=0.05, log2fc_min=0.25)
defc= sc.get.rank_genes_groups_df(Myeloid, group='Macro_FCN1', key='rank_genes_groups', pval_cutoff=0.05, log2fc_min=0.25)
#compare two cluster
sc.tl.rank_genes_groups(Myeloid, 'label', groups=['Macro'], reference='Macro_C1QC',key_added='1v6', pts = True,method='t-test')
de16 = sc.get.rank_genes_groups_df(Myeloid, group='Macro', key='1v6', pval_cutoff=0.05, log2fc_min=0.25)
sc.pl.rank_genes_groups(Myeloid, groups=['1'], n_genes=40, save = '1v6')
sc.pl.umap(Myeloid, color=['CCL5'],save ='umap_leiden')
#Mono_CD14
sc.pl.umap(Myeloid, color=['FCN1','S100A8','VCAN','S100A12'],frameon=False,color_map=Greens,vmin=0.00001,save ='Mono_CD14')
genemarkers = ['FCN1','S100A8','VCAN','S100A12','THBS1']
sc.pl.dotplot(Myeloid, genemarkers, groupby ='leiden',save = 'Mono_CD14')
#cDC1_CLEC9A
sc.pl.umap(Myeloid, color=['CLEC9A','C1orf54','CPNE3','IDO1'],frameon=False,color_map=Greens,vmin=0.00001,save ='cDC1_CLEC9A')
genemarkers = ['CLEC9A','C1orf54','CPNE3','IDO1','IRF8']
sc.pl.dotplot(Myeloid, genemarkers, groupby ='leiden',save = 'cDC1_CLEC9A')
#cDC2_CD1C
sc.pl.umap(Myeloid, color=['CD1C','CLEC10A','FCER1A','JAML'],frameon=False,color_map=Greens,vmin=0.00001,save ='cDC2_CD1C')
genemarkers = ['CD1C','CLEC10A','FCER1A','JAML','CNN2','PLD4']
sc.pl.dotplot(Myeloid, genemarkers, groupby ='leiden',save = 'cDC2_CD1C')
#cDC3_LAMP3
sc.pl.umap(Myeloid, color=['LAMP3','IL7R','CCR7','FSCN1'],frameon=False,color_map=Greens,vmin=0.00001,save ='cDC3_LAMP3')
genemarkers = ['LAMP3','IL7R','CCR7','FSCN1','RFTN1']
sc.pl.dotplot(Myeloid, genemarkers, groupby ='leiden',save = 'cDC3_LAMP3')
#Macro_C1QC
sc.pl.umap(Myeloid, color=['C1QC','C1QA','C1QB','APOE'],frameon=False,color_map=Greens,vmin=0.00001,save ='Macro_C1QC')
genemarkers = ['C1QC','C1QA','C1QB','APOE','PLTP']
sc.pl.dotplot(Myeloid, genemarkers, groupby ='leiden',save = 'Macro_C1QC')
#Macro_ISG15
sc.pl.umap(Myeloid, color=['ISG15','IFI6','IFITM1','MX1'],frameon=False,color_map=Greens,vmin=0.00001,save ='Macro_ISG15')
genemarkers = ['ISG15','IFI6','IFITM1','MX1','ISG20','PLTP','APOBEC3A']
sc.pl.dotplot(Myeloid, genemarkers, groupby ='leiden',save = 'Macro_ISG15')
#Macro_SPP1
sc.pl.umap(Myeloid, color=['SPP1','MARCO','C15orf48','AQP9'],frameon=False,color_map=Greens,vmin=0.00001,save ='Macro_SPP1')
genemarkers = ['SPP1','MARCO','C15orf48','AQP9','IL1RN','FN1','RETN']
sc.pl.dotplot(Myeloid, genemarkers, groupby ='leiden',save = 'Macro_SPP1')
#c13
sc.pl.umap(Myeloid, color=['LILRB2','LILRA5','MTSS1','HES4'],frameon=False,color_map=Greens,vmin=0.00001,save ='c13')
genemarkers = ['LILRB2','LILRA5','MTSS1','HES4','IL1RN','TCF7L2','PLAC8']
sc.pl.dotplot(Myeloid, genemarkers, groupby ='leiden',save = 'c13')
#Macro_C1QC_SPP1
sc.pl.umap(Myeloid, color=['SPP1','C1QC','NUPR1','RNASE1'],frameon=False,color_map=Greens,vmin=0.00001,save ='Macro_C1QC_SPP1')
genemarkers = ['C1QA','C1QC','SPP1','NUPR1','RNASE1','SELENOP']
sc.pl.dotplot(Myeloid, genemarkers, groupby ='leiden',save = 'Macro_C1QC_SPP1')



genemarkers = {}
filename="/lustre/home/mcchen/ref/pdac_label.txt"
f = open(filename)
for line in f:
    tokens = line.strip().split('\t')
    if len(tokens) > 1:
        markers = [gene.upper() for gene in tokens[1:]]
        markers = [x for x in markers if x in Myeloid.raw.var_names]
        genemarkers[tokens[0]] = markers
sc.pl.dotplot(Myeloid, genemarkers, groupby ='leiden',save = 'label')


genemarkers = {}
filename="/lustre/home/mcchen/ref/pdac_label1.txt"
f = open(filename)
for line in f:
    tokens = line.strip().split('\t')
    if len(tokens) > 1:
        markers = [gene.upper() for gene in tokens[1:]]
        markers = [x for x in markers if x in Myeloid.raw.var_names]
        genemarkers[tokens[0]] = markers
sc.pl.dotplot(Myeloid, genemarkers, groupby ='leiden',save = 'label1')


sc.pl.dotplot(Fibroblast, genemarkers, groupby ='leiden',save = 'c1')
#Myloid_pro
sc.pl.umap(Myeloid, color=['MKI67', 'TOP2A'],frameon=False,color_map=Greens,vmin=0.00001,save ='proferating')
genemarkers = ['UBE2C','MKI67', 'TOP2A']
sc.pl.dotplot(Myeloid, genemarkers, groupby ='leiden',save = 'proferating')

#annotattion
myloid_label = {'0':'Macro_FCN1',
              '1':'Macro_C1QC',
              '2':'Macro_C1QC',
              '3':'cDC2_CD1C',
              '4':'Macro_SPP1',
              '5':'Macro_SPP1',
              '6':'Macro_FCN1',
              '7':'Macro_ISG15',
              '8':'Macro_C1QC',
              '9':'Macro_C1QC',
              '10':'Macro_C1QC',
              '11':'Myeloid_Proferating' ,
              '12':'Mono_Macro',
              '13':'Macro_LIRB2',
              '14':'cDC1_CLEC9A',
              '15':'cDC3_LAMP3', 
     }

celltype_labels = []
for leiden in Myeloid.obs.leiden:
    celltype_labels.append(myloid_label.get(leiden))
Myeloid.obs['label'] = celltype_labels
sc.pl.umap(Myeloid, color=['leiden','label'], legend_loc ='on data',save ='umap_leiden_num')


Myeloid.obs['label'] = Myeloid.obs['label'].astype(str)  
pdac_combined.obs['label1'] = pdac_combined.obs['label1'].astype(str)
pdac_combined.obs['label1'].update(Myeloid.obs['label'])    
sc.pl.umap(pdac_combined, color=['label1'],save ='umap_leiden')

de= sc.get.rank_genes_groups_df(Myeloid, group='Macro', key='rank_genes_groups', pval_cutoff=0.05, log2fc_min=0.25)


tmp = pd.crosstab(adata.obs['leiden_0.6'],adata.obs['type'], normalize='index')
tmp.plot.bar(stacked=True).legend(loc='upper right')