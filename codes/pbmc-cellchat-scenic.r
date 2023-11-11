#### pbmc ####

#Erythrocyte
FeaturePlot(pbmc, features = "HBA1", order = TRUE, reduction = "umap")
FeaturePlot(pbmc, features = "HBA2", order = TRUE, reduction = "umap")
#Megakaryocyte
FeaturePlot(pbmc, features = "PPBP", order = TRUE, reduction = "umap")
FeaturePlot(pbmc, features = "TUBB1", order = TRUE, reduction = "umap")
FeaturePlot(pbmc, features = "PF4", order = TRUE, reduction = "umap")
FeaturePlot(pbmc, features = "SPARC", order = TRUE, reduction = "umap")
#B cell & Plasm cell
FeaturePlot(pbmc, features = "MS4A1", order = TRUE, reduction = "umap")
FeaturePlot(pbmc, features = "MZB1", order = TRUE, reduction = "umap")
FeaturePlot(pbmc, features = "CD79A", order = TRUE, reduction = "umap")
FeaturePlot(pbmc, features = "IL4R", order = TRUE, reduction = "umap")
FeaturePlot(pbmc, features = "CD27", order = TRUE, reduction = "umap")
FeaturePlot(pbmc, features = "JCHAIN", order = TRUE, reduction = "umap") 
#B Monocyte
FeaturePlot(pbmc, features = "CD14", order = TRUE, reduction = "umap")
FeaturePlot(pbmc, features = "LYZ", order = TRUE, reduction = "umap")
FeaturePlot(pbmc, features = "FCGR3A", order = TRUE, reduction = "umap")
FeaturePlot(pbmc, features = "MS4A7", order = TRUE, reduction = "umap")
# pDC
FeaturePlot(pbmc, features = "CLEC4C", order = TRUE, reduction = "umap")
FeaturePlot(pbmc, features = "LILRA4", order = TRUE, reduction = "umap")
FeaturePlot(pbmc, features = "SERPINF1", order = TRUE, reduction = "umap")
FeaturePlot(pbmc, features = "GZMB", order = TRUE, reduction = "umap")
FeaturePlot(pbmc, features = "MZB1", order = TRUE, reduction = "umap")
# cDC
FeaturePlot(pbmc, features = "CD1C", order = TRUE, reduction = "umap")
FeaturePlot(pbmc, features = "CLEC10A", order = TRUE, reduction = "umap")
# NK
FeaturePlot(pbmc, features = "NCAM1", order = TRUE, reduction = "umap")
FeaturePlot(pbmc, features = "FCGR3A", order = TRUE, reduction = "umap")
FeaturePlot(pbmc, features = "KLRF1", order = TRUE, reduction = "umap")
FeaturePlot(pbmc, features = "GNLY", order = TRUE, reduction = "umap")
#CD8 T
FeaturePlot(pbmc, features = "CD3D", order = TRUE, reduction = "umap")
FeaturePlot(pbmc, features = "CD8A", order = TRUE, reduction = "umap")
FeaturePlot(pbmc, features = "GZMK", order = TRUE, reduction = "umap")
#CD4 T
FeaturePlot(pbmc, features = "CD3D", order = TRUE, reduction = "umap")
FeaturePlot(pbmc, features = "CD4", order = TRUE, reduction = "umap")


# 未注释UMAP图
umap_unannotation <- DimPlot(pbmc, cols = celltype.color,
                             reduction = "umap",
                             label = "T", 
                             pt.size = 0.15,
                             label.size = 4)
ggsave(filename = "./output/umap_unannotation.pdf", umap_unannotation, wi = 9, he = 7)
# 细胞注释
pbmc <- RenameIdents(object = pbmc, 
                     "0" = "CD4 T", #CD4
                     "1" = "CD4 T", #CD4
                     "2" = "NK", #NCAM1-, FCGR3A, KLRF1, GNLY
                     "3" = "NK", #NCAM1, FCGR3A, KLRF1, GNLY
                     "4" = "Naive B", #MS4A1, MZB1, CD79A, CD27-
                     "5" = "CD14+ Monocyte", #CD14, LYZ
                     "6" = "CD8 T", #CD8A
                     "7" = "CD8 T", #CD8A
                     "8" = "Memory B", #CD27, MS4A1, MZB1, CD79A
                     "9" = "FCGR3A+ Monocyte", #FCGR3A, MS4A7, CD14-
                     "10" = "CD8 T", #CD8A
                     "11" = "NK", #NCAM1-, FCGR3A, KLRF1, GNLY
                     "12" = "Megakaryocyte",#PPBP, TUBB1, PF4, SPARC
                     "13" = "NK", #NCAM1, FCGR3A-, KLRF1, GNLY
                     "14" = "Erythrocyte", #HBA1, HBA2
                     "15" = "cDC", #CD1C, CLEC10A
                     "16" = "NK", #NCAM1, FCGR3A, KLRF1, GNLY
                     "17" = "Plasm cell", #JCHAIN, MS4A1-, MZB1, CD79A
                     "18" = "CD4 T", #CD4
                     "19" = "Naive B",#MS4A1, MZB1, CD79A, CD27-
                     "20" = "Naive B")#MS4A1, MZB1, CD79A, CD27-
# 注释后UMAP图
umap_celltype <- DimPlot(pbmc, cols = celltype.color,
                         reduction = "umap",
                         label = "T", 
                         pt.size = 0.15,
                         label.size = 4)
umap_celltype
ggsave(filename = "./output/umap_celltype.pdf", umap_celltype, wi = 9, he = 7)

umap_celltype_splitbygroup <- DimPlot(pbmc, cols = celltype.color,
                         reduction = "umap",
                         label = "T", 
                         pt.size = 0.15,
                         label.size = 4,
                         split.by = "group")
umap_celltype_splitbygroup
ggsave(filename = "./output/umap_celltype_splitbygroup.pdf", umap_celltype_splitbygroup, wi = 16, he = 7)
# FeaturePlot图
FeaturePlotgene <- c("CD3D",  #CD4 T
                     "NCAM1", "KLRF1", "GNLY",#NK
                     "MS4A1",  "CD79A", #Naive B
                     "CD14", "MS4A7", #CD14+ Monocyte
                     "CD8A", "GZMK", #CD8 T
                     "FCGR3A", #FCGR3A+ Monocyte
                     "PPBP",  #Megakaryocyte
                     "HBA1", #Erythrocyte
                     "CLEC10A", #cDC
                     "MZB1", "LILRA4") #Plasma cell
p1 <- FeaturePlot_scCustom(pbmc, features = FeaturePlotgene, reduction = "umap")
ggsave(filename = './output/FeaturePlot.pdf',plot = p1,he=9,wi=12)
# FeaturePlot图2
FeaturePlotgene2 <- c("CD3D",  #CD4 T
                     "MS4A1",   #Naive B
                     "CD14",  #CD14+ Monocyte
                     "CD8A",  #CD8 T
                     "FCGR3A", #FCGR3A+ Monocyte
                     "MZB1") #Plasma cell 
p3 <- FeaturePlot_scCustom(pbmc, features = "CD3D", reduction = "umap")
ggsave(filename = './output/FeaturePlot CD3D.pdf',plot = p3,he=6,wi=6)
p4 <- FeaturePlot_scCustom(pbmc, features = "MS4A1", reduction = "umap")
ggsave(filename = './output/FeaturePlot MS4A1.pdf',plot = p4,he=6,wi=6)
p5 <- FeaturePlot_scCustom(pbmc, features = "CD14", reduction = "umap")
ggsave(filename = './output/FeaturePlot CD14.pdf',plot = p5,he=6,wi=6)
p6 <- FeaturePlot_scCustom(pbmc, features = "FCGR3A", reduction = "umap")
ggsave(filename = './output/FeaturePlot FCGR3A.pdf',plot = p6,he=6,wi=6)

# Dotplot图
DotPlotgene <- c("CD3D", #CD4 T
                 "NCAM1", "KLRF1", "GNLY",#NK
                 "MS4A1",  "CD79A", "IL4R", #Naive B
                 "CD14", "LYZ", "MS4A7", #CD14+ Monocyte
                 "CD8A", "GZMK", #CD8 T
                 "CD27",  #Memory B
                 "FCGR3A", #FCGR3A+ Monocyte
                 "PPBP", "TUBB1", #Megakaryocyte
                 "HBA1", "HBA2", #Erythrocyte
                 "CD1C", "CLEC10A", #cDC
                 "MZB1", "JCHAIN", "CLEC4C", "LILRA4") #Plasma cell
p2 <- DotPlot(pbmc, features = DotPlotgene,cols = c("lightgrey", "darkblue"),
              dot.min = 0, dot.scale = 6,  col.min = 0,
              col.max = 5) + RotatedAxis()+theme(axis.text.y = element_text(size = 8), axis.text.x = element_text(size = 8))
p2
ggsave(filename = './output/DotPlot.pdf',plot = p2,he = 3,wi = 9)

#### pbmc.CD4 ####
pbmc.CD4 <- subset(pbmc, ident = c("CD4 T"))
pbmc.CD4 <- NormalizeData(pbmc.CD4, verbose = FALSE)
pbmc.CD4 <- FindVariableFeatures(pbmc.CD4, selection.method = "vst", nfeatures = 2000)
pbmc.CD4 <- ScaleData(pbmc.CD4, verbose = FALSE)
pbmc.CD4 <- RunPCA(pbmc.CD4, pc.genes = pbmc@var.genes, npcs = 20, verbose = FALSE, reduction.name = "pca")
pbmc.CD4 <- RunHarmony(pbmc.CD4 , reduction = "pca",group.by.vars = "Samples", reduction.save = "harmony")

  # UMAP analysis
pbmc.CD4 <- RunUMAP(pbmc.CD4, reduction = "harmony", dims = 1:20, reduction.name = "umap")
pbmc.CD4 <- FindNeighbors(pbmc.CD4, reduction = "harmony", dims = 1:20)
pbmc.CD4 <- FindClusters(pbmc.CD4, resolution = 1.1) 

# Featureplot图
FeaturePlot_scCustom(pbmc.CD4, features = "CXCR5", order = TRUE, reduction = "umap")
# 细胞注释
pbmc.CD4 <- RenameIdents(object = pbmc.CD4, 
                     "0" = "Naive T", #CCR7, SELL
                     "1" = "Naive T", #CCR7, SELL
                     "2" = "Tfh", #CXCR5+ ICOS CCR10-
                     "3" = "Th17", #CCR6+ CXCR3- KLRB1+
                     "4" = "Tfh", #CXCR5+ ICOS CCR10-
                     "5" = "Tcm", #CD27, MYC, NOSIP 
                     "6" = "Tcm", #CD27, MYC, NOSIP
                     "7" = "Tem", #GZMK, GZMA, IFNG
                     "8" = "Tem", #GZMK, GZMA, IFNG
                     "9" = "Treg", #FOXP3
                     "10" = "Tem", #GZMK, GZMA, IFNG
                     "11" = "Treg", #FOXP3
                     "12" = "Tem", #GZMK, GZMA, IFNG
                     "13" = "Naive T") #CCR7, SELL
# 注释后UMAP图
umap_celltype_pbmc.CD4 <- DimPlot(pbmc.CD4, cols = celltype.color,
                         reduction = "umap",
                         label = "T", 
                         pt.size = 0.25,
                         label.size = 4)
umap_celltype_pbmc.CD4
ggsave(filename = "./output/umap_celltype_pbmc.CD4.pdf", umap_celltype_pbmc.CD4, wi = 6, he = 5)

umap_celltype_pbmc.CD4_splitbygroup <- DimPlot(pbmc.CD4, cols = celltype.color,
                                      reduction = "umap",
                                      label = "T", 
                                      pt.size = 0.25,
                                      label.size = 4,
                                      split.by = "group")
umap_celltype_pbmc.CD4_splitbygroup
ggsave(filename = "./output/umap_celltype_pbmc.CD4_splitbygroup.pdf", umap_celltype_pbmc.CD4_splitbygroup, wi = 9, he = 5)                     

# Dotplot图
DotPlotgene <- c( "CCR7", "SELL", #Naive CD4+ T CCR7+ SELL+
                  "CXCR5", "ICOS",  #Tfh CXCR5+ ICOS CCR10-
                  "CCR6", "KLRB1",   #Th17 CCR6+ CXCR3- KLRB1+ #Th1 CCR6- CXCR3+ CCR5+ IFNG+ #Th2 CCR6- CXCR3- 
                  "CD27", "NOSIP", #Central memory CD4+ T(Tcm),  CD27+ MYC+ NOSIP+ 
                  "GZMK", "GZMA", "CXCR3","IFNG", #Effector memory CD4+ T (Tem), GZMK+ GZMA+ IFNG+                        
                  "FOXP3","CCR10") #Treg  FOX?P3+
p5 <- DotPlot(pbmc.CD4, features = DotPlotgene,cols = c("lightgrey", "darkblue"),
              dot.min = 0, dot.scale = 7,  col.min = 0,
              col.max = 5) + RotatedAxis()+theme(axis.text.y = element_text(size = 9), axis.text.x = element_text(size = 9))
ggsave(filename = './output/DotPlot_pbmc.CD4.pdf',plot = p5,he = 5,wi = 9)
# 堆积小提琴图
p6 <- Stacked_VlnPlot(pbmc.CD4 , features = DotPlotgene , colors_use = celltype.color)+
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))
ggsave(filename = "./output/Stacked_VlnPlot_pbmc.CD4.pdf", p6, width = 6, height = 10) 
# FeaturePlot图
FeaturePlotgene <- c("CCR6", "KLRB1",   #Th17 CCR6+ CXCR3- KLRB1+ 
                  "GZMA", "CXCR3", #Effector memory CD4+ T (Tem), GZMK+ GZMA+ IFNG+                        
                  "FOXP3","CCR10") #Treg  FOXP3+
p7 <- FeaturePlot_scCustom(pbmc.CD4, features = FeaturePlotgene, reduction = "umap")
p7
ggsave(filename = './output/FeaturePlot_pbmc.CD4.pdf',plot = p7,he=9, wi=8)
# 热图
pbmc.CD4 <- ScaleData(pbmc.CD4, features = rownames(pbmc.CD4))
DoHeatmap(pbmc.CD4, features = DotPlotgene,size = 3, angle = 90)+
  theme(axis.text.y = element_text(size = 8))
AverageHeatmap(object = pbmc.CD4, 
               markerGene = DotPlotgene,
               clusterAnnoName = F,
               htCol = c("#0099CC", "white", "#CC0033"),
               htRange = c(-2, 0, 2),
               annoCol = T,
               myanCol = c("#DC143C","#20B2AA","#FFA500","#9370DB","#F08080","#1E90FF"))
# markGenes = "FOXP3"

# 细胞比例
head(Idents(pbmc.CD4))
pbmc.CD4@meta.data$celltype <- Idents(pbmc.CD4) #metadata里加一列celltype
metadata <- pbmc.CD4@meta.data
data <- as.data.frame(table(pbmc.CD4$group,pbmc.CD4$celltype))
colnames(data) <- c("group","CellType","Freq")
# library(dplyr)
df <- data %>% 
   group_by(group) %>% 
   mutate(Total = sum(Freq)) %>% 
   ungroup() %>% 
   mutate(Percent = Freq/Total)
 df$CellType  <- factor(df$CellType,levels = unique(df$CellType))
write.csv(df, file = "./output/pbmc.CD4_cell_percent.csv",row.names = F,quote = F)
write.table(df, file = "./output/pbmc.CD4_cell_percent.txt", row.names = T,sep = "\t")
# library(ggplot2)
# library(ggpubr)
p <- ggplot(df, aes(x = group, y = Percent, fill = CellType)) +
   geom_bar(position = "fill",stat="identity",  alpha = 1, width = 0.95) + #color = 'white',
   scale_y_continuous(expand = c(0,0)) +
   theme_classic()+
   theme(axis.text.x = element_text(angle = 45,vjust = 0.5, hjust=0.5)) 
p1 <- p + scale_fill_manual(values = celltype.color) 
p1
ggsave(filename = './output/pbmc.CD4_cell_percent_barplot.pdf', plot = p1, he=7.5, wi=5)

#### 基于分组找差异基因 ####
pbmc.CD4@meta.data$celltype <- Idents(pbmc.CD4) #加一列celltype
pbmc.CD4@meta.data$celltype_condition=paste(pbmc.CD4@meta.data$celltype, pbmc.CD4@meta.data$group, sep = "_")
marker_condition=data.frame()
Idents(pbmc.CD4)="celltype_condition"
for ( ci in sort(as.character(unique(pbmc.CD4@meta.data$celltype))) ) {
  tmp.marker <- FindMarkers(
    pbmc.CD4, logfc.threshold = 0, min.pct = 0.1, #logfc不筛选
    only.pos = F, test.use = "wilcox", #wilcox
    ident.1=paste0(ci,"_PBMC_BD"), ident.2=paste0(ci,"_PBMC_HC") #ident.1比ident.2
  )
  
  tmp.marker$gene=rownames(tmp.marker)
  tmp.marker$condition=ifelse(tmp.marker$avg_log2FC > 0, paste0(ci,"_PBMC_BD"),paste0(ci,"_PBMC_HC"))
  tmp.marker$cluster=ci
  
  #tmp.marker=tmp.marker%>%filter(p_val_adj < 0.01) #p_val_adj值不筛选
  tmp.marker=as.data.frame(tmp.marker)
  tmp.marker=tmp.marker%>%arrange(desc(avg_log2FC))
  
  marker_condition=marker_condition%>%rbind(tmp.marker)
}
  write.table(marker_condition,file = "./output/pbmc.CD4_markers.BasedOncondition.txt",quote = F,sep = "\t",row.names = F,col.names = T)
  marker_condition=read.table("./output/pbmc.CD4_markers.BasedOncondition.txt",header = T,sep = "\t",stringsAsFactors = F)
marker_condition$sig=""
marker_condition$sig[abs(marker_condition$avg_log2FC) > 0.25 & marker_condition$p_val_adj < 0.01] = "sig"
marker_condition$sig2=paste(marker_condition$cluster,marker_condition$sig,sep = "_")
marker_condition$sig2[str_detect(marker_condition$sig2,"_$")]="not_sig"
marker_condition$sig2=str_replace(marker_condition$sig2,"_sig","")
#控制顺序
marker_condition$sig2=factor(marker_condition$sig2,levels = c("not",sort(unique(marker_condition$cluster))))
marker_condition$cluster=factor(marker_condition$cluster,levels = sort(unique(marker_condition$cluster)))
marker_condition=marker_condition%>%arrange(cluster,sig2)
#控制范围
marker_condition$avg_log2FC[marker_condition$avg_log2FC > 3]=3
marker_condition$avg_log2FC[marker_condition$avg_log2FC < c(-3)]= -3
#配色
library(RColorBrewer)
library(scales)
color_ct=c(brewer.pal(12, "Set3")[-c(2,3,9,12)],
           brewer.pal(5, "Set1")[2],
           brewer.pal(3, "Dark2")[1])
color_ct=celltype.color
names(color_ct)=sort(unique(as.character(marker_condition$cluster)))
#画图
marker_condition %>% ggplot(aes(x=cluster,y=avg_log2FC,color=sig2))+geom_jitter(width = 0.25,size=0.5)+
  scale_color_manual(values = c(color_ct,"not"="#dee1e6"))+
  scale_y_continuous("HC VS BD, average log2FC",expand = c(0.02,0))+
  theme_bw()+
  theme(
    panel.grid = element_blank(),
    legend.position = "none",
    axis.text.x.bottom = element_text(angle = 45,hjust = 1,size = 14,color = "black"),
    axis.text.y.left = element_text(size = 14,color = "black"),
    axis.title.x.bottom = element_blank(),
    axis.title.y.left = element_text(size = 16)
  )
ggsave(filename = "./output/pbmc.CD4_DEG_jitter1.pdf", width = 4, height = 3)

#### 差异基因GO分析 ####
library(clusterProfiler)
library(org.Hs.eg.db)
deg <- subset(marker_condition, p_val_adj < 0.05 & abs(avg_log2FC) > 0.25) #筛选所有的差异基因
ids=bitr(deg$gene,'SYMBOL','ENTREZID','org.Hs.eg.db') ## 将SYMBOL转成ENTREZID
sce.markers = merge(deg, ids, by.x='gene',by.y='SYMBOL')
sce.markers_up = subset(sce.markers, avg_log2FC > 0) ## BD组上调的差异基因
sce.markers_down = subset(sce.markers, avg_log2FC < 0) ## BD组下调的差异基因

gcSample = split(sce.markers$ENTREZID, sce.markers$cluster) 
# 函数 split() 可以按照分组因子，把向量，矩阵和数据框进行适当的分组
# 它的返回值是一个列表，代表分组变量每个水平的观测
## GO
xx <- compareCluster(gcSample,
                     fun = "enrichGO",
                     OrgDb = "org.Hs.eg.db",
                     ont = "ALL", #BP, CC, MF
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.01,
                     qvalueCutoff = 0.05)
p <- dotplot(xx)
p1 <- p + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5))
ggsave(filename = "./output/pbmc.CD4_GO.pdf", plot = p1, width = 9, height = 9)





table(sce.markers$cluster)
table(sce.markers_up$cluster)
table(sce.markers_down$cluster)

#### Cellchat-CD4 T ####
# 细胞通讯 
table(pbmc.CD4@active.ident)
Idents(pbmc.CD4) <- pbmc.CD4$group
Idents(pbmc.CD4) <- pbmc.CD4$celltype
# 提取两个分组的Seurat对象
pbmc.CD4_BD <- subset(x = pbmc.CD4, idents = "PBMC_BD")
pbmc.CD4_HC <- subset(x = pbmc.CD4, idents = "PBMC_HC")
# 创建Cellchat对象
pbmc.CD4_BD <- createCellChat(pbmc.CD4_BD@assays$RNA@data, meta = pbmc.CD4_BD@meta.data, group.by = "celltype")
pbmc.CD4_HC <- createCellChat(pbmc.CD4_HC@assays$RNA@data, meta = pbmc.CD4_HC@meta.data, group.by = "celltype")
save(pbmc.CD4_BD, pbmc.CD4_HC, file = "./output/pbmc.CD4_Cellchat")
# 分析pbmc.CD4_BD的细胞通讯网络
cellchat <- pbmc.CD4_BD
cellchat@DB <- CellChatDB.human
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
# cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat, raw.use = T, population.size = T)
# cellchat <- filterCommunication(cellchat, min.cells = 5)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
# cellchat <- computeNetSimilarity(cellchat, type = "functional")
# cellchat <- netEmbedding(cellchat, type = "functional")
# cellchat <- netClustering(cellchat, type = "functional")
# cellchat <- computeNetSimilarity(cellchat, type = "structural")
# cellchat <- netEmbedding(cellchat, type = "structural")
# cellchat <- netClustering(cellchat, type = "structural")
# 保存rds文件 
pbmc.CD4_BD_Cellchat <- cellchat
saveRDS(pbmc.CD4_BD_Cellchat, file = "./output/pbmc.CD4_BD_Cellchat.rds")
# 分析pbmc.CD4_HC的细胞通讯网络
cellchat <- pbmc.CD4_HC
cellchat@DB <- CellChatDB.human
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
# cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat, raw.use = T, population.size = T)
# cellchat <- filterCommunication(cellchat, min.cells = 5)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
# cellchat <- computeNetSimilarity(cellchat, type = "functional")
# cellchat <- netEmbedding(cellchat, type = "functional")
# cellchat <- netClustering(cellchat, type = "functional")
# cellchat <- computeNetSimilarity(cellchat, type = "structural")
# cellchat <- netEmbedding(cellchat, type = "structural")
# cellchat <- netClustering(cellchat, type = "structural")
# 保存rds文件
pbmc.CD4_HC_Cellchat <- cellchat
saveRDS(pbmc.CD4_HC_Cellchat, file = "./output/pbmc.CD4_HC_Cellchat.rds")
#合并Cellchat对象 保存rds文件 (readRDS下面3个)
pbmc.CD4_BD_Cellchat <- readRDS(file = "./output/pbmc.CD4_BD_Cellchat.rds")
pbmc.CD4_HC_Cellchat <- readRDS(file = "./output/pbmc.CD4_HC_Cellchat.rds")
pbmc.CD4_Cellchat <- list(pbmc.CD4_HC = pbmc.CD4_HC_Cellchat, pbmc.CD4_BD = pbmc.CD4_BD_Cellchat)
# cellchat <- mergeCellChat(pbmc.CD4_Cellchat, add.names = names(pbmc.CD4_Cellchat), cell.prefix = T)
# saveRDS(cellchat, file = "./output/pbmc.CD4_Cellchat.rds")
cellchat <- readRDS(file = "./output/pbmc.CD4_Cellchat.rds")

# 3.1.所有细胞群总体观：通讯数量与强度对比
# 所有细胞群总体观：通讯数量与强度对比
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "count")
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
p1 <- gg1 + gg2
ggsave(filename = "./output/pbmc.CD4_Cellchat_Overview_number_strength.pdf", plot = p1, width = 6, height = 4)
# 数量与强度差异网络图
par(mfrow = c(1,2))
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")
# pbmc.CD4_Cellchat_netVisual_diffInteraction_circle
# 数量与强度差异热图
par(mfrow = c(1,1))
h1 <- netVisual_heatmap(cellchat)
h2 <- netVisual_heatmap(cellchat, measure = "weight")
h1 + h2
# pbmc.CD4_Cellchat_netVisual_diffInteraction_heatmap

# 细胞互作数量对比网络图
par(mfrow = c(1,2))
weight.max <- getMaxWeight(pbmc.CD4_Cellchat, attribute = c("idents", "count"))
for(i in 1:length(pbmc.CD4_Cellchat)){
  netVisual_circle(pbmc.CD4_Cellchat[[i]]@net$count, weight.scale = T, label.edge = F, edge.weight.max = weight.max[2],
                   edge.width.max = 12, title.name = paste0("Number of interactions - ", names(pbmc.CD4_Cellchat)[i]))
}
# pbmc.CD4_Cellchat_netVisual_circle_number
# 细胞互作强度对比网络图
par(mfrow = c(1,2))
weight.max <- getMaxWeight(pbmc.CD4_Cellchat, attribute = c("idents", "count"))
for(i in 1:length(pbmc.CD4_Cellchat)){
  netVisual_circle(pbmc.CD4_Cellchat[[i]]@net$weigh, weight.scale = T, label.edge = F, edge.weight.max = weight.max[2],
                   edge.width.max = 12, title.name = paste0("Number of strength - ", names(pbmc.CD4_Cellchat)[i]))
}
# pbmc.CD4_Cellchat_netVisual_circle_strength

# 3.2.制定特定细胞互作数量对比网络图
par(mfrow = c(1,2))
s.cell <- c("Treg", "Th17", "Tem","Naive T")
count1 <- pbmc.CD4_Cellchat[[1]]@net$count[s.cell, s.cell]
count2 <- pbmc.CD4_Cellchat[[2]]@net$count[s.cell, s.cell]
weight.max <- max(max(count1), max(count2))
netVisual_circle(count1, weight.scale = T, label.edge = T, edge.weight.max = weight.max, edge.width.max = 12,
                 title.name = paste0("Number of interactions-", names(pbmc.CD4_Cellchat)[1]))
netVisual_circle(count2, weight.scale = T, label.edge = T, edge.weight.max = weight.max, edge.width.max = 12,
                 title.name = paste0("Number of interactions-", names(pbmc.CD4_Cellchat)[2]))
# pbmc.CD4_Cellchat_netVisual_circle_specific cell

# 3.3.比较主要的source和target
num.link <- sapply(pbmc.CD4_Cellchat, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets 
gg <- list()
for (i in 1:length(pbmc.CD4_Cellchat)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(pbmc.CD4_Cellchat[[i]], title = names(pbmc.CD4_Cellchat)[i],
                                               weight.MinMax = weight.MinMax)
}
patchwork::wrap_plots(plots = gg)
# pbmc.CD4_Cellchat_netAnalysis_signalingRole_scatter_source&target
# 3.4.鉴定差异来源
# 鉴定指定细胞亚群的信号变化
gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Treg",
                                            signaling.exclude = "MIF")
gg2 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Th17",
                                            signaling.exclude = "MIF")
gg3 <- gg1 + gg2
patchwork::wrap_plots(plots = gg3)
# pbmc.CD4_Cellchat_netAnalysis_signalingRole_scatter_source&target_specific cell

# 3.5.保守和特异性信号通路的识别与可视化
# 通路信号强度对比分析
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = T)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = T)
p <- gg1 + gg2
ggsave(filename = "./output/pbmc.CD4_Cellchat_Compare_pathway_strength.pdf",plot = p, width = 10,height = 6)
# 3.4.流行学习识别差异信号通路(没跑出来)
cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional")
cellchat <- netClustering(cellchat, type = "functional")
# netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5)
# netVisual_embeddingPairwiseZoomIn(cellchat, type = "functional", nCol = 2)
cellchat <- computeNetSimilarityPairwise(cellchat, type = "structural")
cellchat <- netEmbedding(cellchat, type = "structural")
cellchat <- netClustering(cellchat, type = "structural")
netVisual_embeddingPairwise(cellchat, type = "structural", label.size = 3.5)
netVisual_embeddingPairwiseZoomIn(cellchat, type = "structural", nCol = 2)
p <- rankSimilarity(cellchat, type = "structural") + ggtitle("structural similarity of pathway")
ggsave(filename = "./output/pbmc.CD4_Cellchat_Pathway_Similarity.pdf", p, width = 8, height = 5)

## 3.6.细胞信号模式对比
library(ComplexHeatmap)
# 总体信号模式对比
pathway.union <- union(pbmc.CD4_Cellchat[[1]]@netP$pathways, pbmc.CD4_Cellchat[[2]]@netP$pathways)
ht1 <- netAnalysis_signalingRole_heatmap(pbmc.CD4_Cellchat[[1]], pattern = "all", signaling = pathway.union,
                                         title = names(pbmc.CD4_Cellchat)[1], width = 8, height = 10)
ht2 <- netAnalysis_signalingRole_heatmap(pbmc.CD4_Cellchat[[2]], pattern = "all", signaling = pathway.union,
                                         title = names(pbmc.CD4_Cellchat)[2], width = 8, height = 10)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
# pbmc.CD4_Cellchat_netAnalysis_signalingRole_heatmap
# 输出信号模式对比
pathway.union <- union(pbmc.CD4_Cellchat[[1]]@netP$pathways, pbmc.CD4_Cellchat[[2]]@netP$pathways)
ht1 <- netAnalysis_signalingRole_heatmap(pbmc.CD4_Cellchat[[1]], pattern = "outgoing", signaling = pathway.union,
                                         title = names(pbmc.CD4_Cellchat)[1], width = 8, height = 10)
ht2 <- netAnalysis_signalingRole_heatmap(pbmc.CD4_Cellchat[[2]], pattern = "outgoing", signaling = pathway.union,
                                         title = names(pbmc.CD4_Cellchat)[2], width = 8, height = 10)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
# pbmc.CD4_Cellchat_netAnalysis_signalingRole_heatmap_outgoing
# 输入信号模式对比
pathway.union <- union(pbmc.CD4_Cellchat[[1]]@netP$pathways, pbmc.CD4_Cellchat[[2]]@netP$pathways)
ht1 <- netAnalysis_signalingRole_heatmap(pbmc.CD4_Cellchat[[1]], pattern = "incoming", signaling = pathway.union,
                                         title = names(pbmc.CD4_Cellchat)[1], width = 8, height = 10) # color.heatmap = "GnBu"
ht2 <- netAnalysis_signalingRole_heatmap(pbmc.CD4_Cellchat[[2]], pattern = "incoming", signaling = pathway.union,
                                         title = names(pbmc.CD4_Cellchat)[2], width = 8, height = 10)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
# pbmc.CD4_Cellchat_netAnalysis_signalingRole_heatmap_incoming

## 3.7.特定信号通路的对比 (3种形式图是一个意思)
pathway.union
pathways.show <- c("ITGB2")
weight.max <- getMaxWeight(pbmc.CD4_Cellchat, slot.name = c("netP"), attribute = pathways.show)
par(mfrow = c(1,2), xpd = TRUE)
for(i in 1:length(pbmc.CD4_Cellchat)){
  netVisual_aggregate(pbmc.CD4_Cellchat[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1],
                      edge.width.max = 10, signaling.name = paste(pathways.show, names(pbmc.CD4_Cellchat)[i]))
}
# pbmc.CD4_Cellchat_netVisual_aggregate_ITGB2_circle
# 热图
par(mfrow = c(1,2), xpd = TRUE)
ht <- list()
for(i in 1:length(pbmc.CD4_Cellchat)){
  ht[[i]] <- netVisual_heatmap(pbmc.CD4_Cellchat[[i]], signaling = pathways.show, color.heatmap = "Reds",
                               title.name = paste(pathways.show, "signaling ", names(pbmc.CD4_Cellchat)[i]))
}
ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))
# pbmc.CD4_Cellchat_netVisual_aggregate_ITGB2_heatmap
# 和弦图
par(mfrow = c(1,2), xpd = T)
for(i in 1:length(pbmc.CD4_Cellchat)){
  netVisual_aggregate(pbmc.CD4_Cellchat[[i]], signaling = pathways.show, layout = "chord", pt.title = 3,
                      title.space = 0.05, vertex.label.cex = 0.6, signaling.name = paste(pathways.show, names(pbmc.CD4_Cellchat)[i]))
}
# pbmc.CD4_Cellchat_netVisual_aggregate_ITGB2_chord

## 3.8.配体-受体对比分析
# 气泡图展示所有配体受体对的差异(基于概率)
levels(cellchat@idents$joint)
p4 <- netVisual_bubble(cellchat, sources.use = c(1,2,3,4,5,6), targets.use = c(3,6), comparison = c(1,2), angle.x = 45)
ggsave(filename = "./output/pbmc.CD4_Cellchat_Compare_LR_bubble 123456-36.pdf", plot = p4, width = 12, height = 6)
# 气泡图展示上调或下调的配体受体对
p5 <- netVisual_bubble(cellchat, sources.use = c(1,2,3,4,5,6), targets.use = c(3,6), comparison = c(1,2),
                       max.dataset = 2, title.name = "Increased signaling in TIL", angle.x = 45, remove.isolate = T)
p6 <- netVisual_bubble(cellchat, sources.use = c(1,2,3,4,5,6), targets.use = c(3,6), comparison = c(1,2),
                       max.dataset = 1, title.name = "Decreased signaling in TIL", angle.x = 45, remove.isolate = T)
pc <- p5 + p6
ggsave(filename = "./output/pbmc.CD4_Cellchat_Compare_LR_regulated 123456-36.pdf", plot = pc, width = 12, height = 5.5)

# 交换source和target
p7 <- netVisual_bubble(cellchat, sources.use = c(3,6), targets.use = c(1,2,3,4,5,6), comparison = c(1,2), angle.x = 45)
ggsave(filename = "./output/pbmc.CD4_Cellchat_Compare_LR_bubble 36-123456.pdf", plot = p7, width = 12, height = 6)
# 气泡图展示上调或下调的配体受体对
p8 <- netVisual_bubble(cellchat, sources.use = c(3,6), targets.use = c(1,2,3,4,5,6), comparison = c(1,2),
                       max.dataset = 2, title.name = "Increased signaling in TIL", angle.x = 45, remove.isolate = T)
p9 <- netVisual_bubble(cellchat, sources.use = c(3,6), targets.use = c(1,2,3,4,5,6), comparison = c(1,2),
                       max.dataset = 1, title.name = "Decreased signaling in TIL", angle.x = 45, remove.isolate = T)
pc1 <- p8 + p9
ggsave(filename = "./output/pbmc.CD4_Cellchat_Compare_LR_regulated 36-123456.pdf", plot = pc1, width = 12, height = 5.5)

## 3.9.气泡图展示所有配体受体对的差异(基于信号分子表达)
# 识别上调/下调的配体-受体对（基于基因表达）
pos.dataset = "pbmc.CD4_BD"
features.name = pos.dataset
cellchat1 <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets",
                                        pos.dataset = pos.dataset, features.name = features.name,
                                        only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1,
                                        thresh.p = 1)
net <- netMappingDEG(cellchat1, features.name = features.name)
net.up <- subsetCommunication(cellchat1, net = net, datasets = "pbmc.CD4_BD", ligand.logFC = 0.2,
                              receptor.logFC = NULL)
net.down <- subsetCommunication(cellchat1, net = net, datasets = "pbmc.CD4_HC", ligand.logFC = -0.1,
                              receptor.logFC = -0.1)
# 反卷积从复杂多亚基信号中解单个基因表达
gene.up <- extractGeneSubsetFromPair(net.up, cellchat1)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat1)

pairLR.use.up = net.up[, "interaction_name", drop = F]
gg10 <- netVisual_bubble(cellchat1, pairLR.use = pairLR.use.up, sources.use = c(1,2,6),
                         targets.use = c(3,4,5), comparison = c(1,2), angle.x = 45,
                         remove.isolate = T, title.name = paste0("Upregulated signaling in ",
                                                                 names(pbmc.CD4_Cellchat)[1]))
pairLR.use.down = net.down[, "interaction_name", drop = F]
gg11 <- netVisual_bubble(cellchat1, pairLR.use = pairLR.use.down, sources.use = c(1,2,6),
                         targets.use = c(3,4,5), comparison = c(1,2), angle.x = 45,
                         remove.isolate = T, title.name = paste0("Downregulated signaling in ",
                                                                 names(pbmc.CD4_Cellchat)[1]))
gg12 = gg10 + gg11

#### Cellchat CD4T-Monocyte ####
# 读取Seurat对象
pbmc.CD4_Monocyte_BD <- readRDS(file = "./output/pbmc_CD4T_mon_BD.rds")
pbmc.CD4_Monocyte_HC <- readRDS(file = "./output/pbmc_CD4T_mon_HC.rds")
# 设置celltype为facotr
pbmc.CD4_Monocyte_BD@meta.data$celltype2 <- as.factor(pbmc.CD4_Monocyte_BD@meta.data$subclass)
pbmc.CD4_Monocyte_HC@meta.data$celltype2 <- as.factor(pbmc.CD4_Monocyte_HC@meta.data$subclass)
table(pbmc.CD4_Monocyte_BD$celltype2)
# 创建Cellchat对象
pbmc.CD4_Monocyte_BD <- createCellChat(pbmc.CD4_Monocyte_BD@assays$RNA@data, meta = pbmc.CD4_Monocyte_BD@meta.data, group.by = "celltype2")
pbmc.CD4_Monocyte_HC <- createCellChat(pbmc.CD4_Monocyte_HC@assays$RNA@data, meta = pbmc.CD4_Monocyte_HC@meta.data, group.by = "celltype2")
save(pbmc.CD4_Monocyte_BD, pbmc.CD4_Monocyte_HC, file = "./output/Cellchat_CD4_Monocyte/pbmc.CD4_Monocyte_Cellchat")
# 分析pbmc.CD4_Monocyte_BD的细胞通讯网络
cellchat <- pbmc.CD4_Monocyte_BD
cellchat@DB <- CellChatDB.human
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
# cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat, raw.use = T, population.size = T)
# cellchat <- filterCommunication(cellchat, min.cells = 5)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
# cellchat <- computeNetSimilarity(cellchat, type = "functional")
# cellchat <- netEmbedding(cellchat, type = "functional")
# cellchat <- netClustering(cellchat, type = "functional")
# cellchat <- computeNetSimilarity(cellchat, type = "structural")
# cellchat <- netEmbedding(cellchat, type = "structural")
# cellchat <- netClustering(cellchat, type = "structural")
# 保存rds文件 
pbmc.CD4_Monocyte_BD_Cellchat <- cellchat
saveRDS(pbmc.CD4_Monocyte_BD_Cellchat, file = "./output/Cellchat_CD4_Monocyte/pbmc.CD4_Monocyte_BD_Cellchat.rds")
# 分析pbmc.CD4_Monocyte_HC的细胞通讯网络
cellchat <- pbmc.CD4_Monocyte_HC
cellchat@DB <- CellChatDB.human
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
# cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat, raw.use = T, population.size = T)
# cellchat <- filterCommunication(cellchat, min.cells = 5)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
# cellchat <- computeNetSimilarity(cellchat, type = "functional")
# cellchat <- netEmbedding(cellchat, type = "functional")
# cellchat <- netClustering(cellchat, type = "functional")
# cellchat <- computeNetSimilarity(cellchat, type = "structural")
# cellchat <- netEmbedding(cellchat, type = "structural")
# cellchat <- netClustering(cellchat, type = "structural")
# 保存rds文件 
pbmc.CD4_Monocyte_HC_Cellchat <- cellchat
saveRDS(pbmc.CD4_Monocyte_HC_Cellchat, file = "./output/Cellchat_CD4_Monocyte/pbmc.CD4_Monocyte_HC_Cellchat.rds")
# 合并Cellchat对象 保存rds文件(readRDS下面3个)
pbmc.CD4_Monocyte_BD_Cellchat <- readRDS(file = "./output/Cellchat_CD4_Monocyte/pbmc.CD4_Monocyte_BD_Cellchat.rds")
pbmc.CD4_Monocyte_HC_Cellchat <- readRDS(file = "./output/Cellchat_CD4_Monocyte/pbmc.CD4_Monocyte_HC_Cellchat.rds")
pbmc.CD4_Monocyte_Cellchat <- list(pbmc.CD4_Monocyte_HC = pbmc.CD4_Monocyte_HC_Cellchat, pbmc.CD4_Monocyte_BD = pbmc.CD4_Monocyte_BD_Cellchat)
# cellchat <- mergeCellChat(pbmc.CD4_Monocyte_Cellchat, add.names = names(pbmc.CD4_Monocyte_Cellchat), cell.prefix = T)
# saveRDS(cellchat, file = "./output/Cellchat_CD4_Monocyte/pbmc.CD4_Monocyte_Cellchat.rds")
cellchat <- readRDS(file = "./output/Cellchat_CD4_Monocyte/pbmc.CD4_Monocyte_Cellchat.rds")

# 3.1.所有细胞群总体观：通讯数量与强度对比
# 所有细胞群总体观：通讯数量与强度对比
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "count")
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
p1 <- gg1 + gg2
ggsave(filename = "./output/Cellchat_CD4_Monocyte/pbmc.CD4_Monocyte_Cellchat_Overview_number_strength.pdf", plot = p1, width = 6, height = 4)
# 数量与强度差异网络图
par(mfrow = c(1,2))
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")
# pbmc.CD4_Monocyte_Cellchat_netVisual_diffInteraction_circle
# 数量与强度差异热图
par(mfrow = c(1,1))
h1 <- netVisual_heatmap(cellchat)
h2 <- netVisual_heatmap(cellchat, measure = "weight")
h1 + h2
# pbmc.CD4_Monocyte_Cellchat_netVisual_diffInteraction_heatmap

# 细胞互作数量对比网络图
par(mfrow = c(1,2))
weight.max <- getMaxWeight(pbmc.CD4_Monocyte_Cellchat, attribute = c("idents", "count"))
for(i in 1:length(pbmc.CD4_Monocyte_Cellchat)){
  netVisual_circle(pbmc.CD4_Monocyte_Cellchat[[i]]@net$count, weight.scale = T, label.edge = F, edge.weight.max = weight.max[2],
                   edge.width.max = 12, title.name = paste0("Number of interactions - ", names(pbmc.CD4_Monocyte_Cellchat)[i]))
}
# pbmc.CD4_Monocyte_Cellchat_netVisual_circle_number
# 细胞互作强度对比网络图
par(mfrow = c(1,2))
weight.max <- getMaxWeight(pbmc.CD4_Monocyte_Cellchat, attribute = c("idents", "count"))
for(i in 1:length(pbmc.CD4_Monocyte_Cellchat)){
  netVisual_circle(pbmc.CD4_Monocyte_Cellchat[[i]]@net$weigh, weight.scale = T, label.edge = F, edge.weight.max = weight.max[2],
                   edge.width.max = 12, title.name = paste0("Number of strength - ", names(pbmc.CD4_Monocyte_Cellchat)[i]))
}
# pbmc.CD4_Monocyte_Cellchat_netVisual_circle_strength

# 3.3.比较主要的source和target
num.link <- sapply(pbmc.CD4_Monocyte_Cellchat, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets 
gg <- list()
for (i in 1:length(pbmc.CD4_Monocyte_Cellchat)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(pbmc.CD4_Monocyte_Cellchat[[i]], title = names(pbmc.CD4_Monocyte_Cellchat)[i],
                                               weight.MinMax = weight.MinMax)
}
patchwork::wrap_plots(plots = gg)
# pbmc.CD4_Monocyte_Cellchat_netAnalysis_signalingRole_scatter_source&target

# 3.4.鉴定差异来源
# 鉴定指定细胞亚群的信号变化
gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Treg",
                                            signaling.exclude = "MIF")
gg2 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Th17",
                                            signaling.exclude = "MIF")
gg3 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "FCGR3A+ Monocyte",
                                            signaling.exclude = "MIF")
gg4 <- gg1 + gg2 + gg3
patchwork::wrap_plots(plots = gg4)
# pbmc.CD4_Monocyte_Cellchat_netAnalysis_signalingRole_scatter_source&target_specific cell

# 3.5.保守和特异性信号通路的识别与可视化
# 通路信号强度对比分析
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = T)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = T)
p <- gg1 + gg2
ggsave(filename = "./output/Cellchat_CD4_Monocyte/pbmc.CD4_Monocyte_Cellchat_Compare_pathway_strength.pdf",plot = p, width = 10,height = 6)

## 3.6.细胞信号模式对比
library(ComplexHeatmap)
# 总体信号模式对比
pathway.union <- union(pbmc.CD4_Monocyte_Cellchat[[1]]@netP$pathways, pbmc.CD4_Monocyte_Cellchat[[2]]@netP$pathways)
ht1 <- netAnalysis_signalingRole_heatmap(pbmc.CD4_Monocyte_Cellchat[[1]], pattern = "all", signaling = pathway.union,
                                         title = names(pbmc.CD4_Monocyte_Cellchat)[1], width = 8, height = 10)
ht2 <- netAnalysis_signalingRole_heatmap(pbmc.CD4_Monocyte_Cellchat[[2]], pattern = "all", signaling = pathway.union,
                                         title = names(pbmc.CD4_Monocyte_Cellchat)[2], width = 8, height = 10)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
# pbmc.CD4_Monocyte_Cellchat_netAnalysis_signalingRole_heatmap
# 输出信号模式对比
pathway.union <- union(pbmc.CD4_Monocyte_Cellchat[[1]]@netP$pathways, pbmc.CD4_Monocyte_Cellchat[[2]]@netP$pathways)
ht1 <- netAnalysis_signalingRole_heatmap(pbmc.CD4_Monocyte_Cellchat[[1]], pattern = "outgoing", signaling = pathway.union,
                                         title = names(pbmc.CD4_Monocyte_Cellchat)[1], width = 8, height = 10)
ht2 <- netAnalysis_signalingRole_heatmap(pbmc.CD4_Monocyte_Cellchat[[2]], pattern = "outgoing", signaling = pathway.union,
                                         title = names(pbmc.CD4_Monocyte_Cellchat)[2], width = 8, height = 10)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
# pbmc.CD4_Monocyte_Cellchat_netAnalysis_signalingRole_heatmap_outgoing
# 输入信号模式对比
pathway.union <- union(pbmc.CD4_Monocyte_Cellchat[[1]]@netP$pathways, pbmc.CD4_Monocyte_Cellchat[[2]]@netP$pathways)
ht1 <- netAnalysis_signalingRole_heatmap(pbmc.CD4_Monocyte_Cellchat[[1]], pattern = "incoming", signaling = pathway.union,
                                         title = names(pbmc.CD4_Monocyte_Cellchat)[1], width = 8, height = 10) # color.heatmap = "GnBu"
ht2 <- netAnalysis_signalingRole_heatmap(pbmc.CD4_Monocyte_Cellchat[[2]], pattern = "incoming", signaling = pathway.union,
                                         title = names(pbmc.CD4_Monocyte_Cellchat)[2], width = 8, height = 10)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
# pbmc.CD4_Monocyte_Cellchat_netAnalysis_signalingRole_heatmap_incoming

## 3.7.特定信号通路的对比 (3种形式图是一个意思)
pathway.union
pathways.show <- c("SELPLG")
weight.max <- getMaxWeight(pbmc.CD4_Monocyte_Cellchat, slot.name = c("netP"), attribute = pathways.show)
par(mfrow = c(1,2), xpd = TRUE)
for(i in 1:length(pbmc.CD4_Monocyte_Cellchat)){
  netVisual_aggregate(pbmc.CD4_Monocyte_Cellchat[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1],
                      edge.width.max = 10, signaling.name = paste(pathways.show, names(pbmc.CD4_Monocyte_Cellchat)[i]))
}
# pbmc.CD4_Monocyte_Cellchat_netVisual_aggregate_ITGB2_circle
# 热图
par(mfrow = c(1,2), xpd = TRUE)
ht <- list()
for(i in 1:length(pbmc.CD4_Monocyte_Cellchat)){
  ht[[i]] <- netVisual_heatmap(pbmc.CD4_Monocyte_Cellchat[[i]], signaling = pathways.show, color.heatmap = "Reds",
                               title.name = paste(pathways.show, "signaling ", names(pbmc.CD4_Monocyte_Cellchat)[i]))
}
ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))
# pbmc.CD4_Monocyte_Cellchat_netVisual_aggregate_ITGB2_heatmap
# 和弦图
par(mfrow = c(1,2), xpd = T)
for(i in 1:length(pbmc.CD4_Monocyte_Cellchat)){
  netVisual_aggregate(pbmc.CD4_Monocyte_Cellchat[[i]], signaling = pathways.show, layout = "chord", pt.title = 3,
                      title.space = 0.05, vertex.label.cex = 0.6, signaling.name = paste(pathways.show, names(pbmc.CD4_Monocyte_Cellchat)[i]))
}
# pbmc.CD4_Monocyte_Cellchat_netVisual_aggregate_ITGB2_chord

## 3.8.配体-受体对比分析
# 气泡图展示所有配体受体对的差异(基于概率)
levels(cellchat@idents$joint)
p4 <- netVisual_bubble(cellchat, sources.use = c(2,3,4,5,6,7), targets.use = c(1), comparison = c(1,2), angle.x = 45)
ggsave(filename = "./output/Cellchat_CD4_Monocyte/pbmc.CD4_Monocyte_Cellchat_Compare_LR_bubble 234567-1.pdf", plot = p4, width = 12, height = 6)
# 气泡图展示上调或下调的配体受体对
p5 <- netVisual_bubble(cellchat, sources.use = c(2,3,4,5,6,7), targets.use = c(1), comparison = c(1,2),
                       max.dataset = 2, title.name = "Increased signaling in BD", angle.x = 45, remove.isolate = T)
p6 <- netVisual_bubble(cellchat, sources.use = c(2,3,4,5,6,7), targets.use = c(1), comparison = c(1,2),
                       max.dataset = 1, title.name = "Decreased signaling in BD", angle.x = 45, remove.isolate = T)
pc <- p5 + p6
ggsave(filename = "./output/Cellchat_CD4_Monocyte/pbmc.CD4_Monocyte_Cellchat_Compare_LR_regulated 234567-1.pdf", plot = pc, width = 12, height = 5.5)

# 交换source和target
p7 <- netVisual_bubble(cellchat, sources.use = c(1), targets.use = c(2,3,4,5,6,7), comparison = c(1,2), angle.x = 45)
ggsave(filename = "./output/Cellchat_CD4_Monocyte/pbmc.CD4_Monocyte_Cellchat_Compare_LR_bubble 1-234567.pdf", plot = p7, width = 12, height = 6)
# 气泡图展示上调或下调的配体受体对
p8 <- netVisual_bubble(cellchat, sources.use = c(1), targets.use = c(2,3,4,5,6,7), comparison = c(1,2),
                       max.dataset = 2, title.name = "Increased signaling in BD", angle.x = 45, remove.isolate = T)
p9 <- netVisual_bubble(cellchat, sources.use = c(1), targets.use = c(2,3,4,5,6,7), comparison = c(1,2),
                       max.dataset = 1, title.name = "Decreased signaling in BD", angle.x = 45, remove.isolate = T)
pc1 <- p8 + p9
ggsave(filename = "./output/Cellchat_CD4_Monocyte/pbmc.CD4_Monocyte_Cellchat_Compare_LR_regulated 1-234567.pdf", plot = pc1, width = 12, height = 5.5)







