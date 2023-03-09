data1 <- read.table('GSM4556601_MM033_human_clean_raw_counts_matrix.txt',header = T,row.names = 1)
head(data1[,1:5])
dim(data1)
data2 <- read.table('GSM4556602_MM037_clean_raw_counts_matrix.txt',header = T,row.names = 1)
head(data1[,1:5])
dim(data1)
dim(data2)
h1 <- CreateSeuratObject(counts = data1,min.cells = 3,min.features = 200)
h1$batch <- rep('h1',2291)
h2 <- CreateSeuratObject(counts = data2,min.cells = 3,min.features = 200)
h2$batch <- rep('h2',1597)
BPH <- merge(h1,h2)
BPH[["percent.mt"]] <- PercentageFeatureSet(BPH, pattern = "^MT-")
VlnPlot(BPH, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
BPH <- subset(BPH, subset =  percent.mt < 20)
BPH#19575 features across 3837 samples within 1 assay
BPH <- NormalizeData(BPH, normalization.method = "LogNormalize", scale.factor = 1000000)
BPH <- FindVariableFeatures(BPH, selection.method = "vst",nfeatures = 1000)
all.genes <- rownames(BPH)
BPH <- ScaleData(BPH, features = all.genes)
BPH <- RunPCA(BPH, features = VariableFeatures(object = BPH))
ElbowPlot(BPH)
BPH <- FindNeighbors(BPH, dims = 1:7)
BPH <- FindClusters(BPH, resolution = 0.01)
BPH <- RunUMAP(BPH, dims = 1:7)
DimPlot(BPH, reduction = "umap", label = TRUE,pt.size = 1.5)#3
cluster <- as.character(BPH$seurat_clusters)
cluster[cluster==0]='EP'
cluster[cluster==1]='T cell'
cluster[cluster==2]='Macrophage'
BPH$cluster <- factor(cluster,levels = c('EP','Macrophage','T cell'))
Idents(BPH) <- BPH$cluster
gene= c('EPCAM','KRT8','KRT18','KRT14','TP63','KRT5','KRT15','NKX3-1','CCL3','CD8A','LTB')
DotPlot(BPH, features = gene,dot.scale = 4,dot.min=0)+RotatedAxis()+
  scale_x_discrete("")+scale_y_discrete()+scale_color_gradientn(colours = rev(c("chocolate","white","slateblue")))+
  scale_size_continuous(range=c(1, 12))+theme(panel.grid = element_line(color="lightgrey"))
VlnPlot(BPH,features =gene,pt.size = 0)
FeaturePlot(BPH,features =gene,reduction = 'umap',cols = c('grey','red'))
FeaturePlot(BPH,features ='KRT5',reduction = 'umap',cols = c('grey','red'))
saveRDS(BPH,file = 'BPH.rds')##3== basal cells
BPH <- readRDS('Integrated/BPH.rds')
BPH_EP <- BPH[,BPH$seurat_clusters==0]
BPH_EP#19575 features across 2967 samples within 1 assa
BPH_EP <- NormalizeData(BPH_EP, normalization.method = "LogNormalize", scale.factor = 1000000)
BPH_EP <- FindVariableFeatures(BPH_EP, selection.method = "vst",nfeatures = 1000)
all.genes <- rownames(BPH_EP)
BPH_EP <- ScaleData(BPH_EP, features = all.genes)
BPH_EP <- RunPCA(BPH_EP, features = VariableFeatures(object = BPH_EP))
ElbowPlot(BPH_EP)
BPH_EP <- FindNeighbors(BPH_EP, dims = 1:10)
BPH_EP <- FindClusters(BPH_EP, resolution = 0.04)
BPH_EP <- RunUMAP(BPH_EP, dims = 1:10)
DimPlot(BPH_EP, reduction = "umap", label = TRUE,pt.size = 1.5)#3
VlnPlot(BPH_EP,features = c('KRT8','KRT18','HOXB13','NKX3-1','TACSTD2','PSCA','PIGR','KRT5','KRT15','TP63'))
gene= c('KRT5','TP63','KRT15','ACPP','NKX3-1','HOXB13','PIGR','PSCA','KRT4')
DotPlot(BPH_EP, features = gene,dot.scale = 4,dot.min=0)+RotatedAxis()+
  scale_x_discrete("")+scale_y_discrete()+scale_color_gradientn(colours = rev(c("chocolate","white","slateblue")))+
  scale_size_continuous(range=c(1, 12))+theme(panel.grid = element_line(color="lightgrey"))
FeaturePlot(BPH_EP,features =gene,cols = c('grey','red'))
cluster = as.character(BPH_EP$seurat_clusters)
cluster[cluster==0]='Luminal-A/B'
cluster[cluster==1]='Luminal-C'
cluster[cluster==2]='Basal'
BPH_EP$cluster <- factor(cluster,levels = c('Luminal-C','Luminal-A/B','Basal'))
Idents(BPH_EP)<- BPH_EP$cluster
saveRDS(BPH_EP,file = 'BPH_EP.rds')
BPH_EP <- readRDS('../elife_BPH/BPH_EP.rds')
saveRDS(BPH_EP,file = '../elife_BPH/BPH_EP.rds')
##Integration
library(patchwork)
library(Seurat)
BPH <- readRDS('human1.rds')
h2 <- readRDS('human2.rds')
Bas1 <- BPH[,BPH$seurat_clusters==4]
Bas2 <- h2[,h2$seurat_clusters==6]
human_ep <- readRDS('H:/Basal/human_basal/human_ep.rds')##2,3 basal
Bas3 <- human_ep[,human_ep$seurat_clusters%in%c(0,3)]
Bas1 <- readRDS('Bas1.rds')
Bas2 <- readRDS('Bas2.rds')
Bas3 <- readRDS('Bas3.rds')
Bas1#18428 features across 89 samples within 1 assay 
Bas2#18321 features across 77 samples within 1 assay 
Bas3#17482 features across 3680 samples within 1 assay 
Bas1$batch <- rep('BatcBPH',89)
Bas2$batch <- rep('Batch2',77)
Bas3$batch <- rep('Batch3',3680)
Bas1$state <-  rep('BPH',89)
Bas2$state <- rep('BPH',77)
Bas3$state <- rep('Normal',3680)
Bas3$cluster <- as.character(Bas3$seurat_clusters)
table(Bas3$seurat_clusters)#3077    603 
##250 50
##sampling
cell_id <- colnames(Bas3)
cell_id1 <- cell_id[Bas3$seurat_clusters==0]
cell_id2 <- cell_id[Bas3$seurat_clusters==3]
set.seed(10000)
cell_id1s <- sample(cell_id1,250)
set.seed(10000)
cell_id2s <- sample(cell_id2,50)
Bas31 <- Bas3[,colnames(Bas3)%in%c(cell_id1s,cell_id2s)]
Bas31
Bas12 <- merge(Bas1,Bas2)
Bas <- merge(Bas12,Bas31)
Bas#22531 features across 466 samples within 1 assay 
table(Bas$batch)
Bas <- CreateSeuratObject(counts = Bas@assays$RNA@counts,min.cells = 1,min.features = 200,meta.data = Bas@meta.data)
Bas
Bas.list <- SplitObject(Bas, split.by = "batch")
Bas.list <- lapply(X = Bas.list, FUN = function(x) {
  x <- NormalizeData(x, scale.factor = 1000000)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)
})
features <- SelectIntegrationFeatures(object.list = Bas.list)
features
Bas.anchors <- FindIntegrationAnchors(object.list = Bas.list, anchor.features = features,k.anchor = 50)
Bas.combined <- IntegrateData(anchorset = Bas.anchors,k.weight = 70)
DefaultAssay(Bas.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
Bas.combined <- ScaleData(Bas.combined, verbose = FALSE)
Bas.combined <- RunPCA(Bas.combined, npcs = 30, verbose = FALSE)
Bas.combined <- RunUMAP(Bas.combined, reduction = "pca", dims = 1:12)
Bas.combined <- RunTSNE(Bas.combined, reduction = "pca", dims = 1:12)
Bas.combined <- FindNeighbors(Bas.combined, reduction = "pca", dims = 1:28)
Bas.combined <- FindClusters(Bas.combined, resolution = 0.2,algorithm = 1)
table(Bas.combined$seurat_clusters,Bas.combined$state)
p1 <- DimPlot(Bas.combined, reduction = "tsne", group.by = "batch",pt.size = 3)
p2 <- DimPlot(Bas.combined, reduction = "tsne", label = TRUE, repel = TRUE,pt.size = 3)
#p3 <- DimPlot(Bas.combined, reduction = "tsne", label = TRUE, repel = TRUE,pt.size = 2,group.by = 'cluster')
p1 + p2
table(Bas.combined$seurat_clusters,Bas.combined$batch)
table(Bas.combined$batch)
DefaultAssay(Bas.combined)='RNA'
FeaturePlot(Bas.combined,features = c('KRT5','KRT14','TP63','NKX3-1','C1R'),reduction = 'tsne',cols = c('grey','red'))
DimPlot(Bas.combined, reduction = "umap", split.by = "batch",pt.size = 4)
DimPlot(Bas.combined, reduction = "tsne", split.by = "batch",pt.size = 4)
DimPlot(Bas.combined, reduction = "tsne", split.by = "state",pt.size = 4)
meta <- Bas.combined@meta.data
write.csv(meta,'meta_combined.csv')
meta<- read.csv('meta_combined.csv',header = T,row.names = 1)
Bas.combined$celltype <- meta$cell_type
DimPlot(Bas.combined, reduction = "tsne", group.by = "celltype",pt.size = 4)
VlnPlot(Bas.combined,features = c('NKX3-1','C1R'),pt.size = 0)
cell_type <- as.character(Bas.combined$seurat_clusters)
cell_type[cell_type==0]='Bas-B'
cell_type[cell_type==1]='Bas-A'
Bas.combined$celltype <- as.factor(cell_type)
VlnPlot(Bas.combined,features = c('NKX3-1','C1R'),pt.size = 0,group.by = 'celltype')
per <- table(Bas.combined$seurat_clusters,Bas.combined$state)
per

#       BPH Normal
#Bas-B  162  113
#Bas-A   4  187

table(Bas.combined$state)#166    300 
t(per)/c(166,300)
#         Bas-B       Bas-A
#BPH    0.95783133 0.04216867
#Normal 0.31666667 0.68333333
