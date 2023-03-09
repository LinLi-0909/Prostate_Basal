rm(list=ls())
library(ggplot2)
library(Seurat)
library(dplyr)
library(Matrix)
data <- Read10X(data.dir = "outs+YFP/outs/filtered_gene_bc_matrices_mex/output1_genome/")
dim(data)
cluster <- read.csv('cluster_new.csv',header = T)
head(cluster)
basal <- data[,cluster$cluster==3]
dim(basal)
write.csv(basal,'basal.csv')

basal <- read.csv('basal.csv',header = T,row.names = 1)
basal <- CreateSeuratObject(counts  = basal, min.cells = 3, min.features = 200, 
                           project = "10X_basal")
basal
# 13876 features across 747 samples within 1 assay 
##QC and selecting cells for further analysis
mito.features <- grep(pattern = "^mt-", x = rownames(x = basal), value = TRUE)
percent.mito <- Matrix::colSums(x = GetAssayData(object = basal, slot = "counts")[mito.features,
                                                                                 ])/Matrix::colSums(x = GetAssayData(object =basal, slot = "counts"))

# The [[ operator can add columns to object metadata, and is a great place
# to stash QC stats
basal[["percent.mito"]] <- percent.mito
VlnPlot(object =basal, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"),
        ncol = 3)
FeatureScatter(object = basal, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(object = basal, feature1 = "nCount_RNA", feature2 = "percent.mito")
basal <- subset(x = basal, subset = nFeature_RNA > 200  & nFeature_RNA < 5000)
basal #16401 features across 8541 samples within 1 assay 
#%16401 features across 8537 samples within 1 assay 

basal <- NormalizeData(object = basal, normalization.method = "LogNormalize",
                             scale.factor = 100000)
##Detection of variable genes across the single cells
basal <- FindVariableFeatures(object = basal, selection.method = "mean.var.plot",
                             mean.cutoff = c(0.25, 8), dispersion.cutoff = c(0.5, 10))
VariableFeaturePlot(basal)

basal <- ScaleData(object = basal, features = rownames(x = basal),vars.to.regress = c("nCount_RNA", "percent.mito"))
#Scaling the data and removing unwanted sources of variation
basal <- RunPCA(object = basal, features = VariableFeatures(object =basal))
DimPlot(object = basal)
basal <- ProjectDim(object = basal)
basal <- JackStraw(object = basal, num.replicate = 100)
basal <- ScoreJackStraw(object = basal, dims = 1:10)
JackStrawPlot(object = basal, dims = 1:10)
ElbowPlot(object = basal)

basal <- FindNeighbors(object = basal, dims = 1:10)
basal <- FindClusters(object = basal, resolution = 0.16)
basal <- RunTSNE(object = basal, dims = 1:10)
DimPlot(object = basal, reduction = "tsne", label = TRUE, pt.size = 2)
ggsave('basal_tsne.pdf',dpi = 600,width = 8,height = 6)
FeaturePlot(object = basal, cols = c("lightgrey", "red"), features ='Krt5' , pt.size = 2)
FeaturePlot(object = basal, cols = c("lightgrey", "red"), features ='Krt14' , pt.size = 2)
FeaturePlot(object = basal, cols = c("lightgrey", "red"), features ='Trp63' , pt.size = 2)
basal_cluster <- basal$RNA_snn_res.0.16
write.csv(basal_cluster,'basal_cluster.csv')
save(basal,file ='basal.RDS')
FeaturePlot(object = basal, cols = c("lightgrey", "red"), features ='Asrgl1' , pt.size = 2)
FeaturePlot(object = basal, cols = c("lightgrey", "red"), features ='Gsdma' , pt.size = 2)
FeaturePlot(object = basal, cols = c("lightgrey", "red"), features =c('Krt17') , pt.size = 2)
FeaturePlot(object = basal, cols = c("lightgrey", "red"), features =c('Cpe') , pt.size = 2)
FeaturePlot(object = basal, cols = c("lightgrey", "red"), features =c('Cnn2') , pt.size = 2)
FeaturePlot(object = basal, cols = c("lightgrey", "red"), features =c('Prelp') , pt.size = 2)
FeaturePlot(object = basal, cols = c("lightgrey", "red"), features =c('Dlk2') , pt.size = 2)

load('basal.RDS')
cluster <- as.character(basal$seurat_clusters)
cluster[cluster==0]= 'Bas-A'
cluster[cluster==1]= 'Bas-B'
basal$cluster <- cluster
Idents(basal)<-basal$cluster
gene =c('Krt17','Dlk2','Prelp','Nkx3-1','Pbsn','C1rb')
DotPlot(basal, features = gene,dot.scale = 4,dot.min=0)+RotatedAxis()+
  scale_x_discrete("")+scale_y_discrete()+scale_color_gradientn(colours = rev(c("chocolate","white","slateblue")))+
  scale_size_continuous(range=c(1, 12))+theme(panel.grid = element_line(color="lightgrey"))

markers <- FindAllMarkers(basal,only.pos = T,logfc.threshold = 0.25)
markers %>%
  group_by(cluster) %>%
  top_n(n = 30, wt = avg_log2FC) -> top30
DoHeatmap(basal,features = top30$gene,group.colors = brewer.pal(8,"Set2")[c(1,2)])+scale_fill_gradientn(colors = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100))
write.csv(markers,'Basal_RNA_marker.csv')
basal_norm <- GetAssayData(basal,slot = 'data')
head(basal_norm[,1:5])
write.table(basal_norm,'basal_norm_data.txt')
markers <- read.csv('Basal_RNA_marker.csv',header=T,row.names=1)
tf <- read.table('mouse_tf_trrust.txt',header = F)
head(tf)
Basal_B <- markers[markers$cluster==1,]
Basal_B
Basal_B_TF <- Basal_B[Basal_B$gene%in%tf$V1,]
dim(Basal_B_TF)
write.csv(Basal_B_TF,'Basal_B_TF.csv')

Basal_A <- markers[markers$cluster==0,]
Basal_A
Basal_A_TF <- Basal_A[Basal_A$gene%in%tf$V1,]
dim(Basal_A_TF)
write.csv(Basal_A_TF,'Basal_A_TF.csv')


library(edgeR)
head(cluster)
#TGACAACGTTTGTTGG.1
basal <- GetAssayData(basal,slot = 'counts')
head(basal[,1:5])
basal <- as.matrix(basal)
head(basal[,1:4])
dim(basal) #13876   746
a <- apply(basal,1,function(x){sum(x!=0)})
length(a)
data_filter <- basal[a>=10,]
dim(data_filter) #11679   746
cluster <- read.csv('basal_cluster.csv',header = T)
table(cluster$x)
group <- cluster$x
group[group==0] <-'c0'
group[group==1] <-'c1'

group <- as.factor(group)
table(group)

y <- DGEList(counts=data_filter,group=group,genes=rownames(data_filter))
dge <- calcNormFactors(y)
design <- model.matrix(~0+group)
colnames(design) <- levels(group)
dge <- estimateDisp(dge, design,robust=TRUE)
fit <- glmQLFit(dge,design,robust=TRUE)
sig_cut=0.01


AB<-makeContrasts(A=c0-c1,levels = design)
lrtAB <- glmQLFTest(fit, contrast=AB)
diffAB <- topTags(lrtAB,n=40000,adjust.method="BH", sort.by="PValue", p.value=sig_cut)$table
dim(diffAB) #3272    6
c0_sig <- diffAB[diffAB$logFC>=0.5&diffAB$FDR<=0.01,]
dim(c0_sig) #1167    6
c1_sig <- diffAB[diffAB$logFC<=-0.5&diffAB$FDR<=0.01,]
dim(c1_sig)
write.csv(c0_sig,'basal_c0_signatures.csv')
write.csv(c1_sig,'basal_c1_signatures.csv')
c0_sig <- c0_sig[order(c0_sig$logFC,decreasing = T),]
head(c0_sig)
c1_sig$logFC <- -c1_sig$logFC
c1_sig <- c1_sig[order(c1_sig$logFC,decreasing = T),]
head(c1_sig)
cluster <- read.csv('basal_cluster.csv',header = T,row.names = 1)
g0 <- data.frame(gene= c0_sig$genes)
g1 <- data.frame(gene= c1_sig$genes)
d1 <- basal_norm[rownames(basal_norm)%in%c0_sig$genes,]
d2 <- basal_norm[rownames(basal_norm)%in%c1_sig$genes,]
dim(d1)
dim(d2)
d1 <- as.data.frame(d1)
d2 <- as.data.frame(d2)
d1 <- cbind(rownames(d1),d1)
d2 <- cbind(rownames(d2),d2)
head(d1[,1:4])
head(d2[,1:5])
colnames(d1)[1] <- 'gene'
colnames(d2)[1] <- 'gene'
t1 <- merge(g0,d1,all = F,sort = F)
dim(t1)
t2 <- merge(g1,d2,all = F,sort = F)
dim(t2)
t1 <- t1[1:15,]
t2<- t2[1:15,]
dim(t1)
dim(t2)
m <- rbind(t1,t2)
dim(m)
rownames(m) <- m$gene
m <- m[,-1]
m0 <- m[,cluster$x==0]
m1 <- m[,cluster$x==1]
dim(m0)
dim(m1)
head(m0[,1:4])
head(m1[,1:5])
m2 <- cbind(m0,m1)

library(gplots)

d2 <- log2(m2+1)
m2 <- t(scale(t(d2)))
colors <- colorRampPalette(c("deepskyblue4",'white', "brown3"))(400)
col <- colorRampPalette(c('darkblue','gold'))(400)
breaks <- c(seq(-3,-0.5,length.out = 200),seq(-0.4,4,length.out = 201))
heatmap.2(m2, col=colors,key=T,keysize = 1, symkey=FALSE,  density.info="none", trace="none", cexRow=0.8,labCol = '',labRow = '',srtRow =1,
          Colv=F ,Rowv=F,dendrogram="none", adjCol = c(0.5,1),symbreaks=F,breaks = breaks)
heatmap.2(m2, col=colors,key=T,keysize = 1, symkey=FALSE,  density.info="none", trace="none", cexRow=0.8,labCol = '',srtRow =1,
          Colv=F ,Rowv=F,dendrogram="none", adjCol = c(0.5,1),symbreaks=F,breaks = breaks)
