## Author: Yu Shuhuan
## Date: 2023-10-16
## Brief description: preprocessing of PBMC single cell data using removeDoublets, seurat, singleR etc.

```{r}
library(SingleR)
library(Seurat)
options(Seurat.object.assay.version = "v5")
library(EnsDb.Mmulatta.v109)
library(Signac)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(viridis)
library(DoubletFinder)
library(patchwork)
library(harmony)
library(cowplot)
library(SeuratWrappers)
library(celldex)
library(scCustomize)

### pre-processing data  

for(i in 1:8){
  samples=c("PBMC_BD_1","PBMC_BD_2","PBMC_BD_3","PBMC_BD_4","PBMC_HC_1","PBMC_HC_2","PBMC_HC_3","PBMC_HC_4")
  sample <- samples[i]
  print(paste0("Now processing: ",sample))

  pbmc.data <- Read10X(data.dir = paste0("./raw-data/",sample,"/"))
  pbmc <- CreateSeuratObject(counts = pbmc.data)
  pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
  # Visualize QC metrics as a violin plot
  p <- VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  ggsave(paste0("./result/Plots/",sample,"-QC-violin.pdf"),p)

  raw_cell_number <- dim(pbmc)[2]

  #quality control
  pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 25000 & percent.mt < 5)
  
  cell_numbers_after_qc <- dim(pbmc)[2]

  pbmc <- NormalizeData(pbmc)
  pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
  pbmc <- ScaleData(pbmc)
  pbmc <- RunPCA(pbmc)
  # ElbowPlot(pbmc)
  pbmc <- RunUMAP(pbmc, dims = 1:20)

  # cluster cells
  pbmc <- FindNeighbors(pbmc, dims = 1:20)
  pbmc <- FindClusters(pbmc,resolution=0.3)
  
  # pK Identification (no ground-truth) 
  sweep.res.list <- paramSweep_v3(pbmc, PCs = 1:20, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  mpK<-as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)]))

  # remofve doublets by DoubletFinder
  pbmc$celltype <- Idents(pbmc)
  DoubletRate=ncol(pbmc)*8*1e-6
  homotypic.prop=modelHomotypic(pbmc@meta.data$celltype)
  nExp_poi=round(DoubletRate*length(pbmc$celltype))
  nExp_poi.adj=round(nExp_poi*(1-homotypic.prop))
  pbmc=doubletFinder_v3(pbmc, PCs = 1:20, pN = 0.25, pK = mpK, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)
  DF.clarrsications=names(pbmc@meta.data)[grep("DF.classifications",names(pbmc@meta.data))]
  
  # Plot results
  p1 <- DimPlot(pbmc, reduction = "umap", group.by = DF.clarrsications) + ggtitle("Doublets")
  ggsave(paste0("./result/Plots/",sample,"_DoubletFinder.pdf"),p1,width = 10,height = 10)

  # filter doublets
  Idents(pbmc) <- DF.clarrsications
  pbmc<-subset(x = pbmc, idents="Singlet")   #过滤细胞，只保留经过检测的Singlet
  
  filtered_cell_number <- dim(pbmc)[2]

  print(paste0("raw cell numbers: ",raw_cell_number))
  print(paste0("cell numbers after quality control: ", cell_numbers_after_qc))
  print(paste0("filtered cell numbers: ",filtered_cell_number))

  # save filtered rna_counts 
  filtered_barcodes <- colnames(pbmc[["RNA"]]@counts)
  rna_counts <- pbmc[["RNA"]]@counts
  filtered.rna.counts = rna_counts[,filtered_barcodes]
  colnames(filtered.rna.counts) = paste0(sample,".",colnames(filtered.rna.counts))
  write.table(filtered.rna.counts, paste0("./result/",sample,"_filtered.rna.counts.txt"),row.names = T,sep = "\t",quote = F)

  rm(list = ls())
  gc()
}


### Integration BD samples

BD1.counts <- read.table("./result/PBMC_BD_1_filtered.rna.counts.txt", header = T,sep = "\t",row.names = 1)
BD2.counts <- read.table("./result/PBMC_BD_2_filtered.rna.counts.txt", header = T,sep = "\t",row.names = 1)
BD3.counts <- read.table("./result/PBMC_BD_3_filtered.rna.counts.txt", header = T,sep = "\t",row.names = 1)
BD4.counts <- read.table("./result/PBMC_BD_4_filtered.rna.counts.txt", header = T,sep = "\t",row.names = 1)
BD4.counts <- BD4.counts[rownames(BD1.counts),]
HC1.counts <- read.table("./result/PBMC_HC_1_filtered.rna.counts.txt", header = T,sep = "\t",row.names = 1)
HC1.counts <- HC1.counts[rownames(BD1.counts),]
HC2.counts <- read.table("./result/PBMC_HC_2_filtered.rna.counts.txt", header = T,sep = "\t",row.names = 1)
HC3.counts <- read.table("./result/PBMC_HC_3_filtered.rna.counts.txt", header = T,sep = "\t",row.names = 1)
HC4.counts <- read.table("./result/PBMC_HC_4_filtered.rna.counts.txt", header = T,sep = "\t",row.names = 1)

all.samples.counts <- cbind(BD1.counts,BD2.counts,BD3.counts,BD4.counts,HC1.counts,HC2.counts,HC3.counts,HC4.counts)
rm("BD1.counts", "BD2.counts", "BD3.counts", "BD4.counts","HC1.counts","HC2.counts","HC3.counts","HC4.counts")

pbmc <- CreateSeuratObject(counts = all.samples.counts, project = "pbmc", min.cells = 5) %>%
    Seurat::NormalizeData(verbose = FALSE) %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
    ScaleData(verbose = FALSE) %>% 
    RunPCA(pc.genes = pbmc@var.genes, npcs = 20, verbose = FALSE)

pbmc@meta.data$Samples <- c(rep("PBMC_BD_1", 4409), rep("PBMC_BD_2", 4145), rep("PBMC_BD_3", 2809), rep("PBMC_BD_4", 5759),
                            rep("PBMC_HC_1", 2202), rep("PBMC_HC_2", 5204), rep("PBMC_HC_3", 5864), rep("PBMC_HC_4", 4624)
                            )

p1 <- DimPlot(object = pbmc, reduction = "pca", pt.size = .1, group.by = "Samples")
p2 <- VlnPlot(object = pbmc, features = "PC_1", group.by = "Samples", pt.size = .1)
ggsave("./result/Plots/plot_uncorrected_PCs_samples.pdf", p1+p2, width = 12, height = 6)

# run Harmony
options(repr.plot.height = 2.5, repr.plot.width = 6)
pbmc <- pbmc %>% 
    RunHarmony("Samples", plot_convergence = TRUE)

harmony_embeddings <- Embeddings(pbmc, 'harmony')
harmony_embeddings[1:5, 1:5]

options(repr.plot.height = 5, repr.plot.width = 12)
p1 <- DimPlot(object = pbmc, reduction = "harmony", pt.size = .1, group.by = "Samples")
p2 <- VlnPlot(object = pbmc, features = "harmony_1", group.by = "Samples", pt.size = .1)
ggsave("./result/Plots/plot_corrected_PCs_samples.pdf", plot_grid(p1,p2), width = 12, height = 6)

# UMAP analysis
pbmc <- RunUMAP(pbmc, reduction = "harmony", dims = 1:20)
pbmc <- FindNeighbors(pbmc, reduction = "harmony", dims = 1:20)
pbmc <- FindClusters(pbmc, resolution = 0.8) 

p1 <- DimPlot(pbmc, reduction = "umap", group.by = "Samples", pt.size = .1, split.by = 'Samples')
ggsave("./result/Plots/plot_unlabel_Harmony_dimplot.pdf", p1, width = 20,height = 7)

p2 <- DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = .1)
ggsave("./result/Plots/plot_unlabel_Harmony_umap.pdf", p2, width = 7,height = 7)

# annotation clusters by SingleR
hpca.se <- HumanPrimaryCellAtlasData()

pbmc_for_SingleR <- GetAssayData(pbmc, slot="data") ##获取标准化矩阵
pbmc.hesc <- SingleR(test = pbmc_for_SingleR, ref = hpca.se, labels = hpca.se$label.main) #
pbmc.hesc

# annotate manually
marker.genes <- c("CD3D","CD27","IL7R", # CD4 T
                 "CD14","LYZ",  # CD16- Monocyte
                 "CD79A","MS4A1","IL4R",  # Naive B
                 "CD8A",  # CD8 T
                 "KLRF1", # Activate NK
                 "SPON2","GZMK","GZMB",  # Memory B
                 "MYOM2",  # Resting NK
                 "PPBP",   # Megakaryocyte
                 "FCGR3A","MS4A7",   # CD16+ Monocyte
                 "CLEC10A","CD1C",  # cDC
                 "HBA1",  # Erythrocyte
                 "LILRA4","CLEC4C", #pDC
                 "MZB1" # Plasma B: "IGHA1" not expressed
                 )
cd_genes=unique(marker.genes)

cell.type.cols <- c("#DC143C","#20B2AA","#FFA500","#9370DB","#F08080","#1E90FF",
            "#808000","#CCCCFF","#A0522D","#800080","#D2691E",
            "#87CEEB","#40E0D0","#008B8B","#228B22","#E9967A",
            "#4682B4","#F0E68C","#EE82EE","#FF6347","#8B4513"
            )
p1 <- Stacked_VlnPlot(pbmc,features = cd_genes,colors_use = cell.type.cols)+
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))
pbmc@meta.data$labels <-pbmc.hesc$labels
p2 <- DimPlot(pbmc, group.by = c("seurat_clusters", "labels"),reduction = "umap",cols = cell.type.cols, label = TRUE)
ggsave("./result/Plots/singleR.anno_marker.genes.vinplot.pdf", p2+p1, width = 30, height = 10)

#seurat 和 singleR的table表
meta=pbmc@meta.data
table(pbmc.hesc$labels,meta$seurat_clusters)

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

pbmc$celltype <- Idents(pbmc)

p <- DimPlot(pbmc, reduction = "umap",label=TRUE, cols = cell.type.cols[1:11])
ggsave("./result/Plots/Plots-label-dimplot.pdf",p)

saveRDS(pbmc,"./result/pbmc_scRNA_harmony.rds")
# save.image("./result/PBMC_Tcells_part1.Rdata")
```
