# ライブラリの読み込み----------------------------------------------------------
library(dplyr)
library(Seurat)
library(patchwork)

# split the dataset into a list of two seurat objects (stim and CTRL)----
DCA.list<-list()
data.sap <- Read10X(data.dir = "/home/cardiovascular/Desktop/DCA_human/CellRanger210615/1SA/outs/filtered_feature_bc_matrix/")
data.uapmi<- Read10X(data.dir = "/home/cardiovascular/Desktop/DCA_human/CellRanger210615/2ACS/outs/filtered_feature_bc_matrix/")

DCA.list$sap<- CreateSeuratObject(counts = data.sap, project = "pbmc3k",
                                   min.cells = 3, min.features = 200)
DCA.list$sap<- AddMetaData(DCA.list$sap , "sa", col.name = "stim")
DCA.list$uapmi <- CreateSeuratObject(counts = data.uapmi, project = "pbmc3k",
                                   min.cells = 3, min.features = 200)
DCA.list$uapmi <- AddMetaData(DCA.list$uapmi , "uapmi", col.name = "stim")

# low-quality cellの確認①------------------------------------------------------
# MT-から始まるミトコンドリアRNAを"percent.mt"列としてデータに追加

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
DCA.list$sap[["percent.mt"]] <- PercentageFeatureSet(DCA.list$sap, pattern = "^MT-")
DCA.list$uapmi[["percent.mt"]] <- PercentageFeatureSet(DCA.list$uapmi, pattern = "^MT-")
VlnPlot(DCA.list$sap, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(DCA.list$uapmi, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# low-quality cellの確認②------------------------------------------------------
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

# percent.mtが極端に大きい（10%以上）、またはnFeature_RNAが小さすぎる（200以下）・大きすぎる(2500以上)細胞は、
# 死細胞やdoubletの可能性が高い
plot1 <- FeatureScatter(DCA.list$sap, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(DCA.list$sap, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

plot1 <- FeatureScatter(DCA.list$uapmi, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(DCA.list$uapmi, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2


# low-quality cellのフィルタリング---------------------------------------------
# ミトコンドリアRNAが７以下、遺伝子数が200~2500の間の細胞のみ抽出
DCA.list$sap <- subset(DCA.list$sap, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 8)
DCA.list$uapmi <- subset(DCA.list$uapmi, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 8)


# 先程と比べlow-qualityデータがフィルターされていることを確認
plot1 <- FeatureScatter(DCA.list$sap, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(DCA.list$uapmi, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

plot1 <- FeatureScatter(DCA.list$sap, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(DCA.list$uapmi, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2


# データの正規化（normalization）-----------------------------------------------
DCA.list$sap<- NormalizeData(DCA.list$sap, normalization.method = "LogNormalize", scale.factor = 10000)
DCA.list$uapmi <- NormalizeData(DCA.list$uapmi, normalization.method = "LogNormalize", scale.factor = 10000)

all.genes <- rownames(DCA.list$sap)
DCA.list$sap <- ScaleData(DCA.list$sap, features = all.genes)

all.genes <- rownames(DCA.list$uapmi)
DCA.list$uapmi <- ScaleData(DCA.list$uapmi, features = all.genes)

DCA.list$sap
DCA.list$uapmi 

# UMAPによるクラスタリング　2つ合わせて------------------------------------------------------------------------------

# normalize and identify variable features for each dataset independently
all.genes <- rownames(DCA.list["sap"]) +  rownames(DCA.list["uapmi"])

DCA.list <- lapply(X = DCA.list, FUN = function(x) {
  x <- NormalizeData(x, normalization.method = "LogNormalize", scale.factor = 10000)
  x <- ScaleData(x, features = all.genes)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})


# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = DCA.list)

# Integrationの実行
DCA.anchors <- FindIntegrationAnchors(object.list = DCA.list, anchor.features = features)

# this command creates an 'integrated' data assay
DCA.combined <- IntegrateData(anchorset =DCA.anchors)

# specify that we will perform downstream analysis on the corrected data note that the original
# unmodified data still resides in the 'RNA' assay
DefaultAssay(DCA.combined) <- "integrated"

DCA.combined <- ScaleData(DCA.combined, verbose = FALSE)
DCA.combined <- RunPCA(DCA.combined, npcs = 30, verbose = FALSE)
DCA.combined <- RunUMAP(DCA.combined, reduction = "pca", dims = 1:30)
DCA.combined <- FindNeighbors(DCA.combined, reduction = "pca", dims = 1:30)
DCA.combined <- FindClusters(DCA.combined, resolution = 0.5)

p1 <- DimPlot(DCA.combined, reduction = "umap", group.by = "stim")
p2 <- DimPlot(DCA.combined, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2

# Visualization 別々------
p1 <- DimPlot(DCA.combined, reduction = "umap", split.by = "stim")
p1
p2 <- DimPlot(DCA.combined, reduction = "umap", repel = TRUE)
p2

# cluster change name----
library(ggplot2)

levels(DCA.combined) # [1] "0" "1" "2" "3" "4" "5" "6" "7 "8" "9"
levels(DCA.combined) <- c("0", "1", "2", "5", "9", "3", "4", "6", "7", "8")

new.cluster.ids <- c("0","1","2","3","4", "5", "6","7","8","9")
names(new.cluster.ids) <- levels(DCA.combined)
DCA.combined<- RenameIdents(DCA.combined, new.cluster.ids)
DimPlot(DCA.combined, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

saveRDS(DCA.combined, "/home/cardiovascular/Desktop/DCA_human/RDS/DCA_allcells.rds")
DCA.combined <- readRDS ("/home/cardiovascular/Desktop/DCA_human/RDS/DCA_allcells.rds")


# Myeloid----
VlnPlot(DCA.combined, features = c("CD68","CD14","CSF1R","IL1B","S100A8", "S100A9") )

# Bcells-----
VlnPlot(DCA.combined, features = c("CD79A","CD79B","FCER2","CD22","MYH11","CD34") )

# Tcells-----
VlnPlot(DCA.combined, features = c("CD3E","CD4","CD8A","IL7R","LEF1","GZMK") )

VlnPlot(DCA.combined, features = c("GLS"), split.by = "stim" , cols = c("sa" ="deepskyblue", "uapmi"="red"))

# ヒートマップの作成------------------------------------------------------------------------------
DCA.combined.markers <- FindAllMarkers(DCA.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
DCA.combined.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

top10 <- DCA.combined.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(DCA.combined, features = top10$gene) + NoLegend()
write.table(top10, "table.txt", quote=F, col.names=F, append=T)

#Vinplot----
VlnPlot(DCA.combined, features = c("SLC25A44"), split.by = "stim" )
VlnPlot(DCA.combined, features = c("FCGR3A") )
VlnPlot(DCA.combined, features = c("Th17A") )
VlnPlot(DCA.combined, features = c("CD68") , levels = new_order)

#myeloid-------
DCA.combined_1 <- subset(DCA.combined, idents = c(5, 6, 7))
DCA.combined_1 <- ScaleData(DCA.combined_1, verbose = FALSE)
DCA.combined_1<- RunPCA(DCA.combined_1, npcs = 30, verbose = FALSE)
DCA.combined_1 <- RunUMAP(DCA.combined_1 , reduction = "pca", dims = 1:30)
DCA.combined_1 <- FindNeighbors(DCA.combined_1 , reduction = "pca", dims = 1:30)
DCA.combined_1 <- FindClusters(DCA.combined_1, resolution = 0.8)
p1 <- DimPlot(DCA.combined_1, reduction = "umap", split.by = "stim")
p1
p2 <- DimPlot(DCA.combined_1, reduction = "umap", label = FALSE, repel = TRUE) 
p2

levels(DCA.combined_1) # [1] "0" "1" "2" "3" "4" "5" 
levels(DCA.combined_1) <- c("1", "4", "3", "0", "2", "5" )

new.cluster.ids <- c("0","1","2","3","4", "5")
names(new.cluster.ids) <- levels(DCA.combined_1)
DCA.combined_1<- RenameIdents(DCA.combined_1, new.cluster.ids)
DimPlot(DCA.combined_1, reduction = "umap", label = TRUE, pt.size = 1.0) + NoLegend()

ident.colors <- (scales::hue_pal())(n = length(x = levels(x = DCA.combined_1)))
names(x = ident.colors) <- levels(x = DCA.combined_1)
cell.colors <- ident.colors[Idents(object = DCA.combined_1)]
names(x = cell.colors) <- colnames(x = DCA.combined_1)
cell.colors <- c("#F8766D","#ABA300","#0CB702","#00A9FF","00C19A","FF61CC")

DCA.combined_1.markers <- FindAllMarkers(DCA.combined_1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
DCA.combined_1.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

top10 <- DCA.combined_1.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(DCA.combined_1, features = top10$gene) + NoLegend()
write.table(top10, "table.txt", quote=F, col.names=F, append=T)


VlnPlot(DCA.combined_1, features = c("FCGR3A","FCGR2A","ITGB2","CD16","CEACAM8","CD55","CD44","ELANE") )
VlnPlot(DCA.combined_1, features = c("PTPRC","CD14","MPO","LYZ","ITGB2","CD16","FCGR3A","FCGR3B","CD55","CD44","ELANE","VCAN", "CD52","FCN1","S100A12","IL1B","CSF1R","S100A8","S100A9","CAMP","DEFB4A","CTSS","ITGAM","CD68","CYBB") )
VlnPlot(DCA.combined_1, features = c("DEFB1","DEFA1","DEFB103B") )
VlnPlot(DCA.combined_1, features = c("S100A6","FOS","S100A12","S100A9","HMGB2","KLF6"))
VlnPlot(DCA.combined_1, features = c("HLA-DPA1","HLA-DPB1","HLA-DQA1","CLEC10A","FCER1A","CD1C") )
VlnPlot(DCA.combined_1, features = c("TLR4","TLR2","TLR7","TLR3","TLR12", "TLR9") )
VlnPlot(DCA.combined_1, features = c("HDC","KIT","TPSAB1","CCL5") )
VlnPlot(DCA.combined_1, features = c("IL1B","TNF","CXCL8","NR1H3", "CEBPA","CEBPB") )
VlnPlot(DCA.combined_1, features = c("CD68","CSF1R","ITGAM","S100A8","S100A9","S100A12") )
VlnPlot(DCA.combined_1, features = c("CD14","FCGR3A","IL1B","LYZ","VCAN","FCN1") )
VlnPlot(DCA.combined_1, features = c("CASP1","CASP4","DEFB4A","CTSS","CYBB","IL1B","VEGFA") )
VlnPlot(DCA.combined_1, features = c("CCL2","CCL20","CCL4","HMGB2","KLF6","VEGFA") )
VlnPlot(DCA.combined_1, features = c("IL1B","TNF","CXCL8","CXCL3", "ABCG1","ABCA1") )
VlnPlot(DCA.combined_1, features = c("CD9","TREM2","MRC1","FCER2","C1QA", "C1QB") )
VlnPlot(DCA.combined_1, features = c("C1QC","CD59","APOE","APOC1","NUPR1", "C1QB") )
VlnPlot(DCA.combined_1, features = c("CEACAM8","MPO","IL6","CD63","S100A8","S100A9") )
VlnPlot(DCA.combined_1, features = c("MKI67","ASP175","CSF1R","IL1B","S100A8","S100A9") )
VlnPlot(DCA.combined_1, features = c("CD80","CD86","CD40","CD209","IL10","CD83") )
VlnPlot(DCA.combined_1, features = c("MT2A","MALAT1","HNRNPU","CD3E","KLF4","CCL5") )
VlnPlot(DCA.combined_1, features = c("CD80","CD86","CD40","CD209","IL10","CD83") )
VlnPlot(DCA.combined_1, features = c("CD36","CCR2","CD14","CD61","CX3CR1","CCL5","PTPRC") )
VlnPlot(DCA.combined_1, features = c("HNRNPU","MALAT1","CD14","ICAM1","CX3CR1","CCL5","PTPRC") )
VlnPlot(DCA.combined_1, features = c("HDC","KIT","TPSAB1","GZMB","CCL5","TPSB2") )
VlnPlot(DCA.combined_1, features = c("ENPP3","MITF","CD63","ITGA2","CD33","ITGAM") )
VlnPlot(DCA.combined_1, features = c("TNF","TNFSF13","CXCL8","CXCL3","CCL3") )
VlnPlot(DCA.combined_1, features = c("CD9","MRC1","FCER2","CCL22","AREG","EREG") )
VlnPlot(DCA.combined_1, features = c("MKI67","TUBB","STMN1","TYMS","NFKB1","STAT4") )
VlnPlot(DCA.combined_1, features = c("C1QA","C1QB","C1QC","TYMS","NFKB1","STAT4") )
VlnPlot(DCA.combined_1, features = c("TET2","DNMT3A","ASXL1","PPM1D","NLRP3","STING1","CGAS" ))
VlnPlot(DCA.combined_1, features = c("OLR1","SCARB1","SCARB2","MSR1","MARCO","SRSF2") )
VlnPlot(DCA.combined_1, features = c("NR1H3","CTSD","CTSL","SPP1","MARCO","FABP4") )
VlnPlot(DCA.combined_1, features = c("TET2","TET1","TET3","SPP1","MARCO","SLC25A44") )
VlnPlot(DCA.combined_1, features = c("TNF","TNFSF13"),split.by = "stim",  cols = c("sa" ="deepskyblue", "uapmi"="red"))
VlnPlot(DCA.combined_1, features = c("CD163","F13A1","MS4A4A","SPP1","MARCO","SLC25A44") )


#DC----
DCA.combined_1_DC_only <- subset(DCA.combined_1, idents = c(5))
saveRDS(DCA.combined_1_DC_only , "/home/cardiovascular/Desktop/DCA_human/RDS/DCA_DC.rds")


#Tcells--------
DCA.combined_2 <- subset(DCA.combined, idents = c(0,1,2,3,4))
DCA.combined_2 <- ScaleData(DCA.combined_2, verbose = FALSE)
DCA.combined_2<- RunPCA(DCA.combined_2, npcs = 30, verbose = FALSE)
DCA.combined_2 <- RunUMAP(DCA.combined_2 , reduction = "pca", dims = 1:30)
DCA.combined_2 <- FindNeighbors(DCA.combined_2, reduction = "pca", dims = 1:30)
DCA.combined_2 <- FindClusters(DCA.combined_2, resolution = 0.5)
p1 <- DimPlot(DCA.combined_2, reduction = "umap", split.by = "stim")
p1
p2 <- DimPlot(DCA.combined_2, reduction = "umap", label = FALSE, repel = TRUE)
p2

levels(DCA.combined_2) # [1] "0" "1" "2" "3" "4" 
levels(DCA.combined_2) <- c("0", "2", "1", "3", "4")

new.cluster.ids <- c("0","1","2","3","4")
names(new.cluster.ids) <- levels(DCA.combined_2)
DCA.combined_2<- RenameIdents(DCA.combined_2, new.cluster.ids)


VlnPlot(DCA.combined_2, features = c("CD3E","CD4","CD8A","GZMB","TYROBP","NKG7") )
VlnPlot(DCA.combined_2, features = c("GZMA","GZMK","PRF1","LEF1","LEF1","SELL") )
VlnPlot(DCA.combined_2, features = c("RORA","GATA3","CD40LG","IL10","IL4","IFNG") )
VlnPlot(DCA.combined_2, features = c("IL10","IL4","IL6","IFNG","IL17A") )
VlnPlot(DCA.combined_2, features = c("CCR7","GATA3","RORC","FOXP3") )
VlnPlot(DCA.combined_2, features = c("GZMB","TBX21","NKG7","GNLY","CX3CR1","CD69") )

VlnPlot(DCA.combined_2, features = c("RORA","GATA3","PD1","CD40LG","IL10","IL4","LILRB1R") ,split.by = "stim")

DCA.combined.markers_2 <- FindAllMarkers(DCA.combined_2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
DCA.combined.markers_2 %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

top10 <- DCA.combined.markers_2 %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(DCA.combined_2, features = top10$gene) 
write.table(top10, "table.txt", quote=F, col.names=F, append=T)

#CD8----

DCA.combined_2_1 <- subset(DCA.combined_2, idents = c(0,1))
DCA.combined_2_1 <- ScaleData(DCA.combined_2_1, verbose = FALSE)
DCA.combined_2_1<- RunPCA(DCA.combined_2_1, npcs = 30, verbose = FALSE)
DCA.combined_2_1 <- RunUMAP(DCA.combined_2_1, reduction = "pca", dims = 1:30)
DCA.combined_2_1 <- FindNeighbors(DCA.combined_2_1, reduction = "pca", dims = 1:30)
DCA.combined_2_1 <- FindClusters(DCA.combined_2_1, resolution = 0.5)
p1 <- DimPlot(DCA.combined_2_1, reduction = "umap", split.by = "stim")
p1
p2 <- DimPlot(DCA.combined_2_1, reduction = "umap", label = FALSE, repel = TRUE)
p2

VlnPlot(DCA.combined_2_1, features = c("GZMA","GZMK","GZMB","TBX21","CX3CR1","NKG7") )
VlnPlot(DCA.combined_2_1, features = c("CD69","ITGAE","IL7R","SELL","CD27","CD44") )
VlnPlot(DCA.combined_2_1, features = c("IFNG","EOMES","SELL","CD27","LTB","CCL4") )
VlnPlot(DCA.combined_2_1, features = c("TNF","CD74","IFNG","CD74","CD74","CXCR6") )
VlnPlot(DCA.combined_2_1, features = c("IL7R","CD27","PRF1","CXCR6","SELL","CCR7") )

DCA.combined.markers_2_1 <- FindAllMarkers(DCA.combined_2_1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
DCA.combined.markers_2_1 %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

top10 <- DCA.combined.markers_2_1 %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(DCA.combined_2_1, features = top10$gene) 
write.table(top10, "table.txt", quote=F, col.names=F, append=T)

prop.table(table(Idents(DCA.combined_2_1),DCA.combined_2_1$stim))

#CD4----

DCA.combined_2_2 <- subset(DCA.combined_2, idents = c(2,3))
DCA.combined_2_2 <- ScaleData(DCA.combined_2_2, verbose = FALSE)
DCA.combined_2_2<- RunPCA(DCA.combined_2_2, npcs = 30, verbose = FALSE)
DCA.combined_2_2 <- RunUMAP(DCA.combined_2_2, reduction = "pca", dims = 1:30)
DCA.combined_2_2 <- FindNeighbors(DCA.combined_2_2, reduction = "pca", dims = 1:30)
DCA.combined_2_2 <- FindClusters(DCA.combined_2_2, resolution = 0.7)

DimPlot(DCA.combined_2_2, reduction = "umap", split.by = "stim")
DimPlot(DCA.combined_2_2, reduction = "umap", label = FALSE, repel = TRUE)


VlnPlot(DCA.combined_2_2, features = c("GZMA","GZMK","PRF1","LEF1","IL7R","SELL") )
VlnPlot(DCA.combined_2_2, features = c("GZMB","PRF1","CD28","RORC","FOXP3","CTLA4") )
VlnPlot(DCA.combined_2_2, features = c("IL17A","IL23","IL2","IL10","IL4","IFNG") )
VlnPlot(DCA.combined_2_2, features = c("GZMB","CD28","PRF1","CD69","IL7R","SELL") )
VlnPlot(DCA.combined_2_2, features = c("CXCR3","CCR4","CCR6","TER1"), split.by ="stim",cols = c("sa" ="deepskyblue", "uapmi"="red") )
VlnPlot(DCA.combined_2_2, features = c("CCL4","IFNG","TNFAIP3","SELL","LEF1","CCR7") )
VlnPlot(DCA.combined_2_2, features = c("ZAP70"), split.by ="stim",cols = c("sa" ="deepskyblue", "uapmi"="red") )
VlnPlot(DCA.combined_2_2, features = c("GZMA","GZMK","CD69","CCL5","IL7R","SELL") )
VlnPlot(DCA.combined_2_2, features = c("IL4","IL10","IL17RA","TBX21","GATA3","RORC") )
VlnPlot(DCA.combined_2_2, features = c("LCK","ZAP70","CD40LG","TBX21","GATA3","RORC") )


FeaturePlot(DCA.combined_2_2, features=c('CD69'), min.cutoff=2.0, max.cutoff='q90',split.by = "stim" )
FeaturePlot(DCA.combined_2_2, features=c('CD69'), min.cutoff=2.0, max.cutoff='q90')

levels(DCA.combined_2_2) # [1] "0" "1" "2" "3" "4" 
levels(DCA.combined_2_2) <- c("0", "1", "3", "2", "4")

DCA.combined_2_2_prop <- prop.table(table(Idents(DCA.combined_2_2),DCA.combined_2_2 $stim))
write.csv(DCA.combined_2_2_prop, "prop.DCA.combined_2_2_prop.csv")


new.cluster.ids <- c("0","1","2","3","4")
names(new.cluster.ids) <- levels(DCA.combined_2_2)
DCA.combined_2_2<- RenameIdents(DCA.combined_2_2, new.cluster.ids)


DCA.combined.markers_2_2 <- FindAllMarkers(DCA.combined_2_2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
DCA.combined.markers_2_2 %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

top10 <- DCA.combined.markers_2_2 %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(DCA.combined_2_2, features = top10$gene)
write.table(top10, "table.txt", quote=F, col.names=F, append=T)

DCA.combined_2_2 <- readRDS ("/home/cardiovascular/Desktop/DCA_human/RDS/DCA_CD4Tcells.rds")


#Bcells---------
DCA.combined_3 <- subset(DCA.combined, idents = c(8,9))
DCA.combined_3 <- ScaleData(DCA.combined_3, verbose = FALSE)
DCA.combined_3 <- RunPCA(DCA.combined_3, npcs = 30, verbose = FALSE)
DCA.combined_3 <- RunUMAP(DCA.combined_3 , reduction = "pca", dims = 1:30)
DCA.combined_3 <- FindNeighbors(DCA.combined_3 , reduction = "pca", dims = 1:30)
DCA.combined_3 <- FindClusters(DCA.combined_3, resolution = 0.5)
p1 <- DimPlot(DCA.combined_3, reduction = "umap", split.by = "stim")
p1
p2 <- DimPlot(DCA.combined_3, reduction = "umap", label = FALSE, repel = TRUE)
p2 

VlnPlot(DCA.combined_3, features = c("CD79A","CD79B","CD27","FCER2","CD22") )
VlnPlot(DCA.combined_3, features = c("CD38","PTPRC","CD40","CD24","CD79B","FCER2") )
VlnPlot(DCA.combined_3, features = c("IGHG1","IGHG3","IGHG","IGHM","IGHD","IGHE"))

VlnPlot(DCA.combined_3, features = c("IGHG1","IGHG3","IGHG","IGHM","IGHD","IGHE"),split.by = "stim",  cols = c("sa" ="deepskyblue", "uapmi"="red"))
VlnPlot(DCA.combined_3, features = c("CD40","IL4R","IL21R","PAX5","BACH2","ICOSL","CD80","CD86","HLA-DQB1","HLA-DQB1","HLA-DPA1"), split.by = "stim",  cols = c("sa" ="deepskyblue", "uapmi"="red"))
VlnPlot(DCA.combined_3, features = c("IGHG1","IGHG2","IGHG3","IGLG1"),split.by = "stim",  cols = c("sa" ="deepskyblue", "uapmi"="red"))


DCA.combined_3_1 <- subset(DCA.combined_3, idents = c(1))
DCA.combined_3_1 <- ScaleData(DCA.combined_3_1, verbose = FALSE)
DCA.combined_3_1<- RunPCA(DCA.combined_3_1, npcs = 30, verbose = FALSE)
DCA.combined_3_1 <- RunUMAP(DCA.combined_3_1 , reduction = "pca", dims = 1:30)
DCA.combined_3_1 <- FindNeighbors(DCA.combined_3_1 , reduction = "pca", dims = 1:30)
DCA.combined_3_1 <- FindClusters(DCA.combined_3_1, resolution = 0.5)
p1 <- DimPlot(DCA.combined_3_1, reduction = "umap", split.by = "stim")
p1

# ヒートマップの作成------------------------------------------------------------------------------
DCA.combined_1.markers <- FindAllMarkers(DCA.combined_1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
DCA.combined_1.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

top10 <- DCA.combined_1.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(DCA.combined_1, features = top10$gene) + NoLegend()
write.table(top10, "table.txt", quote=F, col.names=F, append=T)

#RDSファイルの読み込み
DCA.combined_1 <- readRDS (file =  "/home/cardiovascular/Desktop/DCA_human/RDS/DCA_Myeoloid_cluster.rds")



# VInplot----
VlnPlot(DCA.combined_1, features = c("CXCL10","S100A8","HMGB2","HLA-DPA1","HLA-DQA1","HLA-DPB1") ,  split.by = "stim" , cols = c("sa" ="deepskyblue", "uapmi"="red"))
VlnPlot(DCA.combined_1, features = c("APOE") ,  split.by = "stim" , cols = c("sa" ="deepskyblue", "uapmi"="red"))
VlnPlot(DCA.combined_1, features = c("CXCL10") ,  split.by = "stim" , cols = c("sa" ="deepskyblue", "uapmi"="red"))

DCA.combined_1_DC <- subset(DCA.combined_1, idents = c("5_sa", "5_uapmi"))
saveRDS(DCA.combined_1_DC , "/home/cardiovascular/Desktop/DCA_human/RDS/DCA.combined_1_DC.rds")
DCA.combined_1_DC_SA <- subset(DCA.combined_1, idents = c("5_sa"))
saveRDS(DCA.combined_1_DC_SA , "/home/cardiovascular/Desktop/DCA_human/RDS/DCA.combined_1_DC_SA.rds")

DCA.combined_1_DC_ACS <- subset(DCA.combined_1, idents = c("5_uapmi"))
saveRDS(DCA.combined_1_DC_ACS , "/home/cardiovascular/Desktop/DCA_human/RDS/DCA.combined_1_DC_ACS.rds")




DCA.combined <- readRDS ("/home/cardiovascular/Desktop/DCA_human/RDS/DCA_allcells.rds")


VlnPlot(DCA.combined_1_DC, features = c("CD40","CD80","CD86","HLA-DPA1","HLA-DQA1","HLA-DPB1") ,  split.by = "stim" , cols = c("sa" ="deepskyblue", "uapmi"="red"))
VlnPlot(DCA.combined_1_DC, features = c("CLEC10A","FCER1A","CD1C","IL3RA","RELB") ,  split.by = "stim" , cols = c("sa" ="deepskyblue", "uapmi"="red"))
VlnPlot(DCA.combined_1_DC, features = c("NECTIN2","PVR","TNFSF4","IL3RA","RELB") ,  split.by = "stim" , cols = c("sa" ="deepskyblue", "uapmi"="red"))


DCA.combined_2_2_cmCD4 <- subset(DCA.combined_2_2, idents = c("3_sa", "3_uapmi"))
VlnPlot(DCA.combined_2_2_cmCD4, features = c("CD28","CD40LG","CTLA4","TIGHT","FOXP3","TNFRSF4") ,  split.by = "stim" , cols = c("sa" ="deepskyblue", "uapmi"="red"))
VlnPlot(DCA.combined_2_2_cmCD4, features = c("TIGIT","CD226" ,"IL17A","IL4","IL10") ,  split.by = "stim" , cols = c("sa" ="deepskyblue", "uapmi"="red"))

DCA.combined_2_2_eCD4 <- subset(DCA.combined_2_2, idents = c("2_sa", "2_uapmi"))
VlnPlot(DCA.combined_2_eCD4, features = c("CD28","CD40LG","PDCD1","LAG3","CTLA4","TNFRSF4") ,  split.by = "stim" , cols = c("sa" ="deepskyblue", "uapmi"="red"))
VlnPlot(DCA.combined_2_eCD4, features = c("TIGIT","CD226","HAVCR2","PDCD1") ,  split.by = "stim" , cols = c("sa" ="deepskyblue", "uapmi"="red"))

DCA.combined_2_2_nCD4 <- subset(DCA.combined_2_2, idents = c("0_sa","3_sa","1_sa", "0_uapmi", "1_uapmi", "3_uapmi"))
VlnPlot(DCA.combined_2_2_nCD4, features = c("CD28","CD40LG","CTLA4","LAG3","CTLA4","TNFRSF4") ,  split.by = "stim" , cols = c("sa" ="deepskyblue", "uapmi"="red"))
VlnPlot(DCA.combined_2_2_nCD4, features = c("TIGIT","CD226","HAVCR2","PDCD1") ,  split.by = "stim" , cols = c("sa" ="deepskyblue", "uapmi"="red"))

DCA.combined_2_2_e1CD4 <- subset(DCA.combined_2_2, idents = c("0_sa", "0_uapmi"))
VlnPlot(DCA.combined_2_2_e1CD4, features = c("CD28","CD40LG","PDCD1","LAG3","CTLA4","TNFRSF4") ,  split.by = "stim" , cols = c("sa" ="deepskyblue", "uapmi"="red"))
VlnPlot(DCA.combined_2_2_e1CD4, features = c("TIGIT","CD226","HAVCR2","PDCD1") ,  split.by = "stim" , cols = c("sa" ="deepskyblue", "uapmi"="red"))

DCA.combined_2_2_regCD4 <- subset(DCA.combined_2_2, idents = c("4_sa", "4_uapmi"))
VlnPlot(DCA.combined_2_2_regCD4, features = c("CD28","CD40LG","PDCD1","LAG3","CTLA4","TNFRSF4") ,  split.by = "stim" , cols = c("sa" ="deepskyblue", "uapmi"="red"))
VlnPlot(DCA.combined_2_2_regCD4, features = c("TIGIT","CD226","HAVCR2","PDCD1") ,  split.by = "stim" , cols = c("sa" ="deepskyblue", "uapmi"="red"))


DCA.combined_2.markers <- FindAllMarkers(DCA.combined_2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
DCA.combined_2.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

top10 <- DCA.combined.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(DCA.combined_2, features = top10$gene) + NoLegend()
write.table(top10, "table.txt", quote=F, col.names=F, append=T)
        
        
DCA.combined_3.markers <- FindAllMarkers(DCA.combined_3, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
DCA.combined_3.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

top10 <- DCA.combined.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(DCA.combined_3, features = top10$gene) + NoLegend()
write.table(top10, "table.txt", quote=F, col.names=F, append=T)

# volcano all myeloid------

DCA.combined$celltype <- Idents(DCA.combined)
DCA.combined$celltype.stim <- paste(Idents(DCA.combined), DCA.combined$stim, sep="_")
Idents(DCA.combined) <- "celltype.stim" # Idents()関数で変更
levels(DCA.combined) # 変更できているか確認_1


DCA.combined_Myeloid_sa<- subset(DCA.combined, idents = c("5_sa","6_sa","7_sa"))
DCA.combined_Myeloid_uapmi<- subset(DCA.combined, idents = c("5_uapmi","6_uapmi","7_uapmi"))

DCA.combined_Tcells_sa<- subset(DCA.combined, idents = c("0_sa","1_sa","2_sa","3_sa","4_sa"))
DCA.combined_Tcells_uapmi<- subset(DCA.combined, idents = c("0_uapmi","1_uapmi","2_uapmi","3_uapmi","4_uapmi"))

DCA.combined_Bcells_sa<- subset(DCA.combined, idents = c("8_sa","9_sa"))
DCA.combined_Bcells_uapmi<- subset(DCA.combined, idents = c("8_uapmi","9_uapmi"))



DCA.table_myeloid <- FindMarkers(DCA.combined, ident.1 = c("5_sa","6_sa","7_sa"), ident.2 =c("5_uapmi","6_uapmi","7_uapmi"), verbose = FALSE, logfc.threshold = 0)

DCA.table_myeloid$logp <- -log10(DCA.table_myeloid$p_val)

DCA.table_myeloid_filtered_left = subset(DCA.table_myeloid, logp>=5 & avg_log2FC <= -1.0)
DCA.table_myeloid_filtered_right = subset(DCA.table_myeloid, logp>=5 & avg_log2FC >= 1.0)

genes.to.label.left <- rownames(DCA.table_myeloid_filtered_left)
genes.to.label.right <- rownames(DCA.table_myeloid_filtered_right)

p1 <- ggplot(DCA.table_myeloid, aes(avg_log2FC, logp, label)) + geom_point() 
p1 <- LabelPoints(plot = p1, points = genes.to.label.right,color="blue", repel = TRUE, xnudge=0)
p1 <- LabelPoints(plot = p1, points = genes.to.label.left,color="red", repel = TRUE, xnudge=0)
p1

# VOLCANO Trem2 macs----



# volcano all CD4------

DCA.combined_2$celltype <- Idents(DCA.combined_2)
DCA.combined_2$celltype.stim <- paste(Idents(DCA.combined_2), DCA.combined_2$stim, sep="_")
Idents(DCA.combined_2) <- "celltype.stim" # Idents()関数で変更
levels(DCA.combined_2) # 変更できているか確認_1

DCA.table_CD4 <- FindMarkers(DCA.combined_2, ident.1 = c("2_sa"), ident.2 =c("2_uapmi"), verbose = FALSE, logfc.threshold = 0)

DCA.table_CD4$logp <- -log10(DCA.table_CD4$p_val)

DCA.table_CD4_filtered_left = subset(DCA.table_CD4, logp>=2 & avg_log2FC <= -0.5)
DCA.table_CD4_filtered_right = subset(DCA.table_CD4, logp>=2 & avg_log2FC >= 0.5)

genes.to.label.left <- rownames(DCA.table_CD4_filtered_left)
genes.to.label.right <- rownames(DCA.table_CD4_filtered_right)

p1 <- ggplot(DCA.table_CD4, aes(avg_log2FC, logp, label)) + geom_point() 
p1 <- LabelPoints(plot = p1, points = genes.to.label.right,color="blue", repel = TRUE, xnudge=0)
p1 <- LabelPoints(plot = p1, points = genes.to.label.left,color="red", repel = TRUE, xnudge=0)
p1

write.table(DCA.table_CD4, "effecCD4.txt", quote=F, col.names=F, append=T)

# central memory CD4

DCA.combined_2_2$celltype <- Idents(DCA.combined_2_2)
DCA.combined_2_2$celltype.stim <- paste(Idents(DCA.combined_2_2), DCA.combined_2_2$stim, sep="_")
Idents(DCA.combined_2_2) <- "celltype.stim" # Idents()関数で変更
levels(DCA.combined_2_2) # 変更できているか確認

DCA.table_cCD4 <- FindMarkers(DCA.combined_2_2, ident.1 = c("3_sa"), ident.2 =c("3_uapmi"), verbose = FALSE, logfc.threshold = 0)

DCA.table_cCD4$logp <- -log10(DCA.table_cCD4$p_val)

DCA.table_cCD4_filtered_left = subset(DCA.table_cCD4, logp>=2 & avg_log2FC <= -0.4)
DCA.table_cCD4_filtered_right = subset(DCA.table_cCD4, logp>=2 & avg_log2FC >= 0.4)

genes.to.label.left <- rownames(DCA.table_cCD4_filtered_left)
genes.to.label.right <- rownames(DCA.table_cCD4_filtered_right)

p1 <- ggplot(DCA.table_cCD4, aes(avg_log2FC, logp, label)) + geom_point() 
p1 <- LabelPoints(plot = p1, points = genes.to.label.right,color="blue", repel = TRUE, xnudge=-0.1)
p1 <- LabelPoints(plot = p1, points = genes.to.label.left,color="red", repel = TRUE, xnudge=-0.1)
p1

DCA.combined_2_2_SA <- subset(DCA.combined_2_2, idents = c("0_sa","1_sa","2_sa","3_sa","4_sa"))
saveRDS(DCA.combined_2_2_SA, "/home/cardiovascular/Desktop/DCA_human/RDS/DCA_CD4Tcells_SA.rds")
DCA.combined_2_2_ACS <- subset(DCA.combined_2_2, idents = c("0_uapmi","1_uapmi","2_uapmi","3_uapmi","4_uapmi"))
saveRDS(DCA.combined_2_2_ACS , "/home/cardiovascular/Desktop/DCA_human/RDS/DCA_CD4Tcells_ACS.rds")



write.table(DCA.table1, "table.txt", quote=F, col.names=F, append=T)

# volcano myeloid subset
DCA.combined_1$celltype <- Idents(DCA.combined_1)
DCA.combined_1$celltype.stim <- paste(Idents(DCA.combined_1), DCA.combined_1$stim, sep="_")
Idents(DCA.combined_1) <- "celltype.stim" # Idents()関数で変更
levels(DCA.combined_1)

DCA.table1 <- FindMarkers(DCA.combined_1, ident.1 = "5_uapmi", ident.2 = "5_sa", verbose = FALSE, logfc.threshold = 0)
DCA.table2 <- FindMarkers(DCA.combined_1, ident.1 = "0_uapmi", ident.2 = "0_sa", verbose = FALSE, logfc.threshold = 0)
DCA.table3 <- FindMarkers(DCA.combined_1, ident.1 = "1_uapmi", ident.2 = "1_sa", verbose = FALSE, logfc.threshold = 0)

DCA.table1$logp <- -log10(DCA.table1$p_val)
DCA.table2$logp <- -log10(DCA.table2$p_val)
DCA.table3$logp <- -log10(DCA.table3$p_val)
DCA.table4$logp <- -log10(DCA.table4$p_val)

# cluster1----
DCA.table1_filtered_left = subset(DCA.table1, logp>=1.30103& avg_log2FC <= -0.8)
DCA.table1_filtered_right = subset(DCA.table1, logp>=1.30103 & avg_log2FC >= 0.8)

genes.to.label.left <- rownames(DCA.table1_filtered_left)
genes.to.label.right <- rownames(DCA.table1_filtered_right)

p1 <- ggplot(DCA.table1, aes(avg_log2FC, logp, label)) + geom_point() 
p1 <- LabelPoints(plot = p1, points = genes.to.label.right,color="red", repel = TRUE, xnudge=0)
p1 <- LabelPoints(plot = p1, points = genes.to.label.left,color="blue", repel = TRUE, xnudge=0)
p1

write.table(DCA.table1, "table.txt", quote=F, col.names=F, append=T)

VlnPlot(DCA.combined_1, ident.1 = "5_sa", ident.2 = "5_uapmi", features = c("CXCL10") ,  split.by = "stim" , cols = c("sa" ="deepskyblue", "uapmi"="red"))

# cluster2----
DCA.table2_filtered_left = subset(DCA.table2, logp>=2 & avg_log2FC <= -0.8)
DCA.table2_filtered_right = subset(DCA.table2, logp>=2 & avg_log2FC >= 0.8)

genes.to.label.left <- rownames(DCA.table2_filtered_left)
genes.to.label.right <- rownames(DCA.table2_filtered_right)

p1 <- ggplot(DCA.table2, aes(avg_log2FC, logp, label)) + geom_point() 
p1 <- LabelPoints(plot = p1, points = genes.to.label.right,color="red", repel = TRUE, xnudge=0)
p1 <- LabelPoints(plot = p1, points = genes.to.label.left,color="blue", repel = TRUE, xnudge=0)
p1

write.table(DCA.table2, "table.txt", quote=F, col.names=F, append=T)

# cluster3----

DCA.table3_filtered_left = subset(DCA.table3, logp>=2 & avg_log2FC <= -0.8)
DCA.table3_filtered_right = subset(DCA.table3, logp>=2 & avg_log2FC >= 0.8)

genes.to.label.left <- rownames(DCA.table3_filtered_left)
genes.to.label.right <- rownames(DCA.table3_filtered_right)

p1 <- ggplot(DCA.table3, aes(avg_log2FC, logp, label)) + geom_point() 
p1 <- LabelPoints(plot = p1, points = genes.to.label.right,color="red", repel = TRUE, xnudge=0)
p1 <- LabelPoints(plot = p1, points = genes.to.label.left,color="blue", repel = TRUE, xnudge=0)
p1

write.table(DCA.table2, "table.txt", quote=F, col.names=F, append=T)


plot1 <- FeatureScatter(DCA.combined_3, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(DCA.combined_3, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

DCA.combined_3



