#!/usr/bin/Rscript
#load library
library(SeuratDisk)
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(RColorBrewer)
library(tidyr)
library(dplyr)
#setting path
setwd("/lustre/home/whhou/00.datasets/2021_Data/process/mapping/PDAC/RNA")
out_data_dir <- "/lustre/home/mcchen/PDAC/process/PBMC/result/"
#load data
dir() -> filename
f_PBMC_GEX <- filename[grep("PBMC_GEX",filename)]
PBMC_GEX <- vector("list",length = length(f_PBMC_GEX))
names(PBMC_GEX) <- f_PBMC_GEX
for (i in f_PBMC_GEX) {
  PBMC_GEX[[i]]<- Read10X(data.dir = paste0("./",i,"/",i,"/outs","/filtered_feature_bc_matrix/"))
}
PBMC_GEX$PDAC08_PBMC_GEX <- Read10X('/lustre/home/mcchen/data/PDAC/PDAC08_PBMC_GEX/PDAC08_PBMC_GEX/outs/filtered_feature_bc_matrix/')
PBMC_GEX$PDAC12_PBMC_GEX <- Read10X('/lustre/home/mcchen/data/PDAC/PDAC12_PBMC_GEX/PDAC12-PBMC-GEX/outs/filtered_feature_bc_matrix/')
#create seurat object  
for (i in names(PBMC_GEX)) {
  PBMC_GEX[[i]] <- CreateSeuratObject(counts = PBMC_GEX[[i]], project = names(PBMC_GEX[i]), min.cell = 3, min.features = 200)
}
#pbmc meta data editing
PBMC_GEX[[1]]$ID <- "PDAC05_PBMC"
PBMC_GEX[[2]]$ID <- "PDAC06_PBMC"
PBMC_GEX[[3]]$ID <- "PDAC07_PBMC"
PBMC_GEX[[4]]$ID <- "PDAC09_PBMC"
PBMC_GEX[[5]]$ID <- "PDAC10_PBMC"
PBMC_GEX[[6]]$ID <- "PDAC11_PBMC"
PBMC_GEX[[7]]$ID <- "PDAC08_PBMC"
PBMC_GEX[[8]]$ID <- "PDAC12_PBMC"
#do filter
pbmc_combined <- merge(PBMC_GEX[[1]], PBMC_GEX[2:length(PBMC_GEX)])
pbmc_combined <- PercentageFeatureSet(pbmc_combined, pattern = "^MT-",col.name = "Percent_mito")
pbmc_combined <- subset(pbmc_combined,subset = Percent_mito < 20 )
pbmc_combined <- subset(pbmc_combined,subset = nFeature_RNA < 2500 )
pbmc_combined <- SCTransform(pbmc_combined,verbose = F)
#load data 
ref_pbmc <- LoadH5Seurat("/lustre/home/mcchen/ref/pbmc_multimodal.h5seurat")
anchors <- FindTransferAnchors(
  reference = ref_pbmc,
  query = pbmc_combined,
  normalization.method = "SCT",
  reference.reduction = "spca",
  dims = 1:50
)
pbmc_combined <- MapQuery(
  anchorset = anchors,
  query = pbmc_combined,
  reference = ref_pbmc,
  refdata = list(
    celltype.l1 = "celltype.l1",
    celltype.l2 = "celltype.l2",
    celltype.l3 = "celltype.l3"
  ),
  reference.reduction = "spca", 
  reduction.model = "wnn.umap"
)
# Merge the batches 
plot1 <- DimPlot(pbmc_combined, reduction = "ref.umap", group.by = "celltype.l1", label = TRUE, repel = TRUE, label.size = 3)
plot2 <- DimPlot(pbmc_combined, reduction = "ref.umap", group.by = "celltype.l2", label = TRUE, repel = TRUE, label.size = 3) 
p2 <- plot1 + plot2
ggsave(p2, file=paste0(out_data_dir, 'pbmc_mapping_plot.pdf'))
p3 <- DimPlot(pbmc_combined, reduction = "ref.umap", group.by = "ID", label = TRUE, repel = TRUE, label.size = 3) 
ggsave(p3, file=paste0(out_data_dir, 'pbmc_pid_umap_plot.pdf'))
saveRDS(pbmc_combined, file=paste0(out_data_dir, 'pbmc_combined_mapping.rds'))


