library(Seurat)
library(dplyr)
library(tidyverse)
library(patchwork)
library(devtools)
library(harmony)
Bas.combined <- readRDS('../Bas.combined.rds')
Bas_normal <- Bas.combined[,Bas.combined$state=='Normal']
Bas_BPH <- Bas.combined[,Bas.combined$state=='BPH']
scRNAlist <- list()
scRNAlist[[1]] <- CreateSeuratObject(Bas_normal@assays$RNA@counts, project='Normal', min.cells=1, min.features = 100,meta.data = Bas_normal@meta.data)
scRNAlist[[1]][["percent.mt"]] <- PercentageFeatureSet(scRNAlist[[1]], pattern = "^MT-")
scRNAlist[[1]] <- subset(scRNAlist[[1]], subset =percent.mt < 10 ) 
scRNAlist[[2]] <- CreateSeuratObject(Bas_BPH@assays$RNA@counts, project='BPH', min.cells=1, min.features = 100,meta.data = Bas_BPH@meta.data)
scRNAlist[[2]][["percent.mt"]] <- PercentageFeatureSet(scRNAlist[[2]], pattern = "^MT-")
scRNAlist[[2]] <- subset(scRNAlist[[2]], subset = percent.mt < 10)   
scRNA_harmony <- merge(scRNAlist[[1]], y=c(scRNAlist[[2]]))
scRNA_harmony <- NormalizeData(scRNA_harmony) %>% FindVariableFeatures(nfeatures=1200) %>% ScaleData() %>% RunPCA(verbose=FALSE)
scRNA_harmony <- RunHarmony(scRNA_harmony,group.by.vars = c('state'),plot_convergence	=T)
scRNA_harmony <- RunUMAP(scRNA_harmony, reduction = "harmony", dims = 1:10)
scRNA_harmony <- FindNeighbors(scRNA_harmony, reduction = "harmony", dims = 1:10) %>% FindClusters(resolution = 0.03)
marker <- FindAllMarkers(scRNA_harmony,only.pos = T)
marker %>%
  group_by(cluster) %>%
  top_n(n = 50, wt = avg_log2FC) -> top50
DoHeatmap(scRNA_harmony, features = top50$gene) + NoLegend()

cluster = as.character(scRNA_harmony$seurat_clusters)
cluster[cluster==0]='Basal-B'
cluster[cluster==1]='Basal-A'
scRNA_harmony$cluster <- factor(cluster,levels = c('Basal-A','Basal-B'))
Idents(scRNA_harmony) <- scRNA_harmony$cluster
plot1 = DimPlot(scRNA_harmony, reduction = "umap", label=T) 
plot2 = DimPlot(scRNA_harmony, reduction = "umap", group.by='state') 
plot1+plot2
plot3=VlnPlot(scRNA_harmony,features = c('KRT5','KRT15','TP63','NKX3-1','C1R'))
plot3
plot4=FeaturePlot(scRNA_harmony,features = c('KRT5','KRT15','TP63','NKX3-1','C1R'),reduction = 'umap')
plot5 = DimPlot(scRNA_harmony, reduction = "umap", split.by = 'state') 
plot1+plot2+plot4
plot4
table(scRNA_harmony$seurat_clusters,scRNA_harmony$state)
#       BPH Normal
#Bas-B 160    158
#Bas-A  18    121
#0.8988764 0.5663082
#0.1011236 0.4336918
saveRDS(scRNA_harmony,file = 'scRNA_harmony.rds')
scRNA_harmony <- readRDS('scRNA_harmony.rds')
scRNA_harmony
FeaturePlot(scRNA_harmony,features = c('KRT13', 'OLFM4', 'LY6D', 'S100P'),reduction = 'umap')
