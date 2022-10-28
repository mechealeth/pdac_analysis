#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 24 11:26:28 2022

@author: mecheal
"""

label1 = {
"PDAC04_Normal_GEX" :'N',
"PDAC11_Tumor_GEX":'High_NI',
"PDAC06_Normal_GEX" :'N',
"PDAC08_Normal_GEX" :'N',
"PDAC08_Tumor_GEX " :'Low_NI',
"PDAC04_Tumor_GEX":'High_NI',
"PDAC09_Normal_GEX":'N',
"PDAC12_Tumor_GEX ":'No_NI',
"PDAC10_Normal_GEX":'N',
"PDAC10_Tumor_GEX":'High_NI',
"PDAC06_Tumor_GEX":'High_NI',
"PDAC09_Tumor_GEX":'High_NI',
"PDAC01_Tumor_GEX":'Low_NI',
"PDAC07_Tumor_GEX":'High_NI',
"PDAC02_Tumor_GEX":'High_NI',
"PDAC03_Tumor_GEX":'High_NI',
"PDAC05_Tumor_GEX":'High_NI'
}
celltype_labels = []
for leiden in pdac_combined.obs.pid:
    celltype_labels.append(label1.get(leiden))
pdac_combined.obs['NI'] = celltype_labels

sc.pl.umap(pdac_combined, color='NI', groups=['Low_NI'],na_in_legend=False,frameon=False, save = 'Low_NI')

genemarkers = {}
filename="/lustre/home/mcchen/ref/pdac_heatmap.txt"
f = open(filename)
for line in f:
    tokens = line.strip().split(',')
    if len(tokens) > 1:
        markers = [gene.upper() for gene in tokens[1:]]
        markers = [x for x in markers if x in pdac_combined.raw.var_names]
        genemarkers[tokens[0]] = markers
sc.pl.heatmap(pdac_combined,genemarkers, groupby='label', cmap='viridis', save ='heatmap')

markers = ['PRSS1', 'CTRB2', 'CTRB1', 'REG1B','CD79A', 'MS4A1', 'CD79B','UBE2C', 'TOP2A', 'MKI67','SPP1', 'MMP7', 'KRT7', 'CFTR',
          'PLVAP', 'ENG', 'PECAM1', 'SERPINE1', 'CLDN5', 'DCN', 'LUM', 'COL3A1','KIT', 'MS4A2', 'GATA2','CD14', 'ITGAM', 'MNDA','JCHAIN', 'CD79A', 'MZB1',
        'NGFR', 'CDH19', 'SOX10', 'PLP1','CD3E', 'CD3G', 'CD3D']
sc.pl.matrixplot(pdac_combined, markers, groupby='label',save ='heatmap')