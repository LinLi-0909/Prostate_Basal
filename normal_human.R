library(Seurat)
human <- readRDS('human_prostate.rds')
human
table(human$seurat_clusters)
DimPlot(human,label = T)
human_ep <- human[,human$seurat_clusters%in%c(0,1,2)]
human_ep
saveRDS(human_ep,'human_ep.rds')
human_ep <- readRDS('human_ep.rds')
human_ep
human_ep$cluster <- human_ep$seurat_clusters
human_ep[["percent.mt"]] <- PercentageFeatureSet(human_ep, pattern = "^MT-")
VlnPlot(human_ep, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
human_ep#18497 features across 9825 samples within 1 assay
human_ep <- subset(human_ep, subset =  percent.mt < 15&nFeature_RNA >500)
human_ep#18497 features across 9730 samples within 1 assay
human_ep <- NormalizeData(human_ep, normalization.method = "LogNormalize", scale.factor = 10000)
human_ep <- FindVariableFeatures(human_ep, selection.method = "vst",nfeatures = 1000)
all.genes <- rownames(human_ep)
human_ep <- ScaleData(human_ep, features = all.genes)
human_ep <- RunPCA(human_ep, features = VariableFeatures(object = human_ep))
ElbowPlot(human_ep)
human_ep <- FindNeighbors(human_ep, dims = 1:6)
human_ep <- FindClusters(human_ep, resolution = 0.08)
human_ep <- RunUMAP(human_ep, dims = 1:6)
DimPlot(human_ep, reduction = "umap", label = TRUE,pt.size = 1.5)#3
FeaturePlot(human_ep,features = c('KRT8','KRT18','PIGR','PSCA','KRT4','NKX3-1','HOXB13','KRT15','KRT5','TP63'),cols = c('grey','red'))
VlnPlot(human_ep,features =c('KRT8','KRT18','PIGR','PSCA','KRT4','NKX3-1','HOXB13','KRT15','KRT5','TP63'),pt.size = 0 )
cluster <- as.character(human_ep$seurat_clusters)
cluster[cluster==0]='Basal'
cluster[cluster==1]='Luminal-A/B'
cluster[cluster==2]='Luminal-C'
human_ep$cluster <- factor(cluster,levels = c('Luminal-C','Luminal-A/B','Basal'))
Idents(human_ep) <- human_ep$cluster
DimPlot(human_ep, reduction = "umap", label = TRUE,pt.size = 1.5)#3
DotPlot(human_ep, features = c('KRT5','TP63','KRT15','ACPP','NKX3-1','HOXB13','PIGR','PSCA','KRT4'),group.by = 'cluster',dot.scale = 4,dot.min=0)+RotatedAxis()+
  scale_x_discrete("")+scale_y_discrete()+scale_color_gradientn(colours = rev(c("chocolate","white","slateblue")))+
  scale_size_continuous(range=c(1, 12))+theme(panel.grid = element_line(color="lightgrey"))
saveRDS(human_ep,file = 'human_ep.rds')

Basal_normal <- human_ep[,human_ep$cluster=='Basal']
saveRDS(Basal_normal,file = 'Basal_normal.rds')
Basal_normal <- readRDS('Basal_normal.rds')
