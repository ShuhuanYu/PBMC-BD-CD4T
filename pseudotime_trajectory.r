## Author: Yu Shuhuan
## Date: 2023-10-17
## Brief description: inference the differiation trajectory of CD4T subtypes, especiallty Treg and Th17

```{r}
### monocle3

library(monocle3)
library(Seurat)
library(SeuratData)
library(SeuratWrappers)
library(ggplot2)
library(patchwork)
library(magrittr)

# load seurat object
CD4T <- readRDS("./result/annotaion_pbmc.CD4_scRNA_hamony.rds")
CD4T$subtype <- Idents(CD4T)

# Asign partions
cds <- as.cell_data_set(CD4T)

recreate.partitions <- c(rep(1, length(cds@colData@rownames)))
names(recreate.partitions) <- cds@colData@rownames
recreate.partitions <- as.factor(recreate.partitions)
recreate.partitions

cds@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partitions

# Assign cluster information
list.cluster <- CD4T@active.ident
cds@clusters@listData[["UMAP"]][["clusters"]] <- list.cluster
# Assign UMAP coordinates
cds@int_colData@listData[["reducedDims"]]@listData[["UMAP"]] <- CD4T@reductions$umap@cell.embeddings

# visualization
cluster.before.traj <-plot_cells(cds, color_cells_by = "cluster", label_groups_by_cluster = F, 
           group_label_size = 5) + theme(legend.position = "right")

cds <- learn_graph(cds, use_partition = F)
p1 <- plot_cells(cds, color_cells_by = "cluster", label_groups_by_cluster = F,
           label_branch_points = T, label_roots = T, label_leaves = F,
           group_label_size = 5)
# Order cells in Pseudotime
cds <- order_cells(cds, reduction_method = "UMAP", root_cells = colnames(cds[, clusters(cds) == "Naive T"]))
p2 <- plot_cells(cds, color_cells_by = "pseudotime", label_groups_by_cluster = T,
           label_branch_points = T, label_roots = F, label_leaves = F)

# Cells ordered by Monocle3 Pseudotime
cds$monocle3_pseudotime <- pseudotime(cds)
data.pseudo <- as.data.frame(colData(cds))

ggplot(data.pseudo, aes(monocle3_pseudotime, subtype, fill = subtype)) + geom_boxplot()
p3 <- ggplot(data.pseudo, aes(monocle3_pseudotime, reorder(subtype, monocle3_pseudotime), fill = subtype)) + geom_boxplot()
ggsave("./result/Plots/monocle3-CD4T_subtype_pseudoime-trajectory.pdf",p1+p2+p3,width=15,height=6)

saveRDS(cds, "./result/monocle3-PBMC_CD4T.cds")
```

##
To explore the differentiation of CD4T subtypes, we conducted pseudotime trajectory inference using
Monocle3.Naive T cells were selected as root node manually based on biological knowledge.The developmental
trajectory of CD4 T subtypes suggested a binary branched structure with two end states:fate1,Treg,and 
fate2,Th17(Figure).Inference from Slingshot also supported that(Supplementary Fig.).Reording the pseudotime
revealed Th17 cells are the most mature subtype in CD4 T cells.To explore further the differences between two 
development branches, each branch were choosen using "choose_graph_segments" function.

beanplot: the p-value of BD and HC pseudotime test in fate1 and fate2 are 0.2326 and 0.005242 respectively. wilcox.test function

pseudotime heatmap for each fate: monocle3 algothrithm and monocle2 visualization.
##
