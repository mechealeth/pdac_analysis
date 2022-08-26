#!/usr/bin/Rscript
#load library
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(harmony)
#load data
setwd("/lustre/home/whhou/00.datasets/2021_Data/process/mapping/PDAC/RNA")
out_data_dir <- "/lustre/home/mcchen/PDAC/process/PDAC/obj/"
dir() -> filename
patterns <- c("Tumor_GEX", "Normal_GEX")
f_PDAC_GEX <- filename[grep(paste(patterns, collapse="|"),filename)]
PDAC_GEX <- vector("list",length = length(f_PDAC_GEX))
names(PDAC_GEX) <- f_PDAC_GEX
for (i in f_PDAC_GEX) {
  PDAC_GEX[[i]] <- Read10X(data.dir = paste0("./",i,"/",i,"/outs","/filtered_feature_bc_matrix/"))
}
#create seurat object  
for (i in f_PDAC_GEX) {
  PDAC_GEX[[i]] <- CreateSeuratObject(counts = PDAC_GEX[[i]], project = names(PDAC_GEX[i]), min.cell = 3, min.features = 200)
}
#pdac meta data editing
PDAC_GEX[[1]]$ID <- "PDAC01_Tumor"
PDAC_GEX[[2]]$ID <- "PDAC02_Tumor"
PDAC_GEX[[3]]$ID <- "PDAC03_Tumor"
PDAC_GEX[[4]]$ID <- "PDAC04_Normal"
PDAC_GEX[[5]]$ID <- "PDAC04_Tumor"
PDAC_GEX[[6]]$ID <- "PDAC05_Tumor"
PDAC_GEX[[7]]$ID <- "PDAC06_Normal"
PDAC_GEX[[8]]$ID <- "PDAC06_Tumor"
PDAC_GEX[[9]]$ID <- "PDAC07_Tumor"
PDAC_GEX[[10]]$ID <- "PDAC08_Normal"
PDAC_GEX[[11]]$ID <- "PDAC08_Tumor"
PDAC_GEX[[12]]$ID <- "PDAC09_Normal"
PDAC_GEX[[13]]$ID <- "PDAC09_Tumor"
PDAC_GEX[[14]]$ID <- "PDAC10_Normal"
PDAC_GEX[[15]]$ID <- "PDAC10_Tumor"
PDAC_GEX[[16]]$ID <- "PDAC11_Tumor"
#disease state
PDAC_GEX[[1]]$DiseaseState <- "PDAC"
PDAC_GEX[[2]]$DiseaseState <- "PDAC"
PDAC_GEX[[3]]$DiseaseState <- "PDAC"
PDAC_GEX[[4]]$DiseaseState <- "Adj_Normal"
PDAC_GEX[[5]]$DiseaseState <- "PDAC"
PDAC_GEX[[6]]$DiseaseState <- "PDAC"
PDAC_GEX[[7]]$DiseaseState <- "Adj_Normal"
PDAC_GEX[[8]]$DiseaseState <- "PDAC"
PDAC_GEX[[9]]$DiseaseState <- "PDAC"
PDAC_GEX[[10]]$DiseaseState <- "Adj_Normal"
PDAC_GEX[[11]]$DiseaseState <- "PDAC"
PDAC_GEX[[12]]$DiseaseState <- "Adj_Normal"
PDAC_GEX[[13]]$DiseaseState <- "PDAC"
PDAC_GEX[[14]]$DiseaseState <- "Adj_Normal"
PDAC_GEX[[15]]$DiseaseState <- "PDAC"
PDAC_GEX[[16]]$DiseaseState <- "PDAC"
#PDAC12
PDAC_GEX <- Read10X('/lustre/home/mcchen/data/PDAC/PDAC12_Tumor_GEX/PDAC12-Tumor-GEX/outs/filtered_feature_bc_matrix/')
PDAC_GEX_merge <- CreateSeuratObject(PDAC_GEX, min.cell = 3, min.features = 200)
PDAC_GEX_merge$ID <- "PDAC12_Tumor"
PDAC_GEX_merge$DiseaseState <- "PDAC"
#PDAC merge data and filter
for (i in 1:length(f_PDAC_GEX)) {
  PDAC_GEX_merge <- merge(x = PDAC_GEX_merge, y = PDAC_GEX[[i]])
}
PDAC_GEX_merge [["percent.mt"]] <- PercentageFeatureSet(PDAC_GEX_merge , pattern = "^MT-")
plot1 <- FeatureScatter(PDAC_GEX_merge, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(PDAC_GEX_merge, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
p1 <- plot1 + plot2
ggsave(p1, file=paste0(out_data_dir, 'qc_plot.pdf'))
PDAC_GEX_merge <- subset(PDAC_GEX_merge, subset = Percent_mito < 20)
PDAC_GEX_merge <- subset(PDAC_GEX_merge,subset = nCount_RNA < 60000)
PDAC_GEX_merge %>% NormalizeData() %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA() -> PDAC_GEX_merge
saveRDS(PDAC_GEX_merge, file=paste0(out_data_dir, 'PDAC_GEX_merge.rds'))
#clean
rm(PDAC_GEX)
gc()
#with harmony batch correction
PDAC_GEX_merge <- RunHarmony(PDAC_GEX_merge, group.by.vars = "ID")
p <- ElbowPlot(PDAC_GEX_merge)
ggsave(p, file=paste0(out_data_dir, 'pdac1_harmony_elbow_1.2.pdf'))
PDAC_GEX_merge_harmony <- RunUMAP(PDAC_GEX_merge,reduction = "harmony", dims = 1:50,verbos = F)
PDAC_GEX_merge_harmony <- FindNeighbors(PDAC_GEX_merge_harmony,reduction = "harmony", dims = 1:50)
PDAC_GEX_merge_harmony <- FindClusters(PDAC_GEX_merge_harmony,resolution = 1.2)
p2 <- DimPlot(PDAC_GEX_merge_harmony,reduction = "umap",label = T) + NoLegend()
p3 <- DimPlot(PDAC_GEX_merge_harmony,reduction = "umap", group.by = "ID") + ggtitle("patient_ID")
ggsave(p2, file=paste0(out_data_dir, 'pdac1_harmony_umap_1.2.pdf'))
ggsave(p3, file=paste0(out_data_dir, 'pdac_harmony_umap_sample.pdf'))
saveRDS(PDAC_GEX_merge_harmony, file=paste0(out_data_dir, 'PDAC_GEX_merge_harmony.rds'))
# cluster cell
clustermarker <- FindAllMarkers(PDAC_GEX_merge_harmony,logfc.threshold = 0.25)
write.table(clustermarker,file = sprintf("%s/pdac1_harmony_cluster_maker.xls",out_data_dir),quote = F,sep = "\t",row.names = F)
top5 <- clustermarker %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
p4<- DoHeatmap(PDAC_GEX_merge_harmony, features = top5$gene,slot = "scale.data", disp.min = -2, disp.max = 2, group.by = "ident" , group.bar = T) + scale_fill_gradientn(colors = c("blue", "white", "red")) + theme(axis.text.y = element_text(size = 4))
ggsave(p4, file=paste0(out_data_dir, 'pdac1_seurat_cluster_1.2_marker_heatmap.png'), width=18, height=16, dpi = 300, units = "in", device='png')
saveRDS(PDAC_GEX_merge_harmony, "PDAC_GEX_merge_harmony_cluster.rds")
#marker gene plot
Cell_Types <- c("Epi","T_Cell","Myeloid","B_Cell","Fibroblast","Macro","NK", "Endo","Acinar",'DC','Neural','Mono','Plasma')
Epi_Markers <- c("KRT7","KRT8","KRT18","KRT19","EPCAM","CDH1")
T_Cell_Markers <- c("CD3E","CD3G","CD3D","CD4","IL7R","CD8A","LEF1")
Myeloid_Markers <- c("CD14","ITGAM","MNDA","MPEG1","ITGAX")
B_Cell_Markers <- c("CD79A","MS4A1","CD19")
Fibroblast_Markers <- c("CDH11","PDGFRA","PDGFRB","ACTA2")
NK_Markers <- c("NCR3","FCGR3A","NCAM1","KLRF1","KLRC1","CD38","KLRC1")
Endo_Markers <- c("CDH5","PECAM1")
Acinar_Markers <- c("TRY4","SPINK1","AMY2A")
Macro_Markers <- c('C1QA','CD68','TREM2')
DC_Markers <- c('IRF7','HLA-DRA','LYZ','CST3')
Neural_Markers <- c('S100B','PLP1','SOX10')
Mono_Markers <- c('S100A9','LYZ','FCN1')
Plasma_Markers <- c('JCHAIN','CD79A')
All_Markers <- list(Epi_Markers,T_Cell_Markers,Myeloid_Markers,B_Cell_Markers,Fibroblast_Markers,Macro_Markers,NK_Markers,Endo_Markers,Acinar_Markers,DC_Markers,Neural_Markers,Mono_Markers,Plasma_Markers)
Markers = c(Epi_Markers,T_Cell_Markers,Myeloid_Markers,B_Cell_Markers,Fibroblast_Markers,Macro_Markers,NK_Markers,Endo_Markers,Acinar_Markers,DC_Markers,Neural_Markers,Mono_Markers,Plasma_Markers)
#plot marker
for(i in 1: length(All_Markers)){
  marker_to_test <- All_Markers[[i]]
  PDAC_GEX_merge_harmony <- AddModuleScore(PDAC_GEX_merge_harmony, features = marker_to_test, name = Cell_Types[i])
  i <- FeaturePlot(PDAC_GEX_merge_harmony, label = T,features = paste0(Cell_Types[i], "1"))
  ggsave(i, file=paste0(out_data_dir, filename = paste0(Cell_Types[i], "_pdf")))
  }
p5 <- DoHeatmap(PDAC_GEX_merge_harmony,features = Markers,slot = data)
ggsave(p2, file=paste0(out_data_dir, 'pdac1_harmony_celltype_marker_heatmap.png'), width=18, height=14, dpi = 300, units = "in", device='png')
