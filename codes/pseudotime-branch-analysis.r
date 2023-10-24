## Author: Yu Shuhuan
## Date: 2023-10-23

``{r}


cds <- readRDS("./result/monocle3-PBMC_CD4T.rds")

# compare the pseudotime between BD and HC samples
cds_fate1 <- choose_graph_segments(cds) # Treg
cds_fate2 <- choose_graph_segments(cds) # Th17

fate1.pseudo <- as.data.frame(colData(cds_fate1))
fate2.pseudo <- as.data.frame(colData(cds_fate2))

data <- data.frame(rbind(fate1.pseudo[,c("monocle3_pseudotime","group")],fate2.pseudo[,c("monocle3_pseudotime","group")]))
data$fate <- rep(c("fate1","fate2"),c(dim(fate1.pseudo)[1],dim(fate2.pseudo)[1]))
data$group <- ifelse(data$group=="PBMC_BD","BD","HC")
data$x <- paste(data$fate,data$group,sep = " ")

# visualization(violin plot)
library(vioplot)
library(beanplot)

pdf("./result/Plots/trajectory-part4.pdf")
beanplot(monocle3_pseudotime ~ x, data = data, ll = 0.04,
         main = "beanplot", ylab = "pseudotime", side= "both",
         border = NA,horizontal = F,
         col = list(c("#006ebc","#006ebc"),c("#f5633e", "#f5633e")))
legend("bottomleft", fill =c("#006ebc", "#f5633e"),
       legend = c("HC", "BD"))
dev.off()

# differential test
wilcox.test(x = data[which(data$x == "fate1 BD"),1], y = data[which(data$x == "fate1 HC"),1])   # p=0.2326
wilcox.test(x = data[which(data$x == "fate2 BD"),1], y = data[which(data$x == "fate2 HC"),1])   # p=0.005242

saveRDS(cds_fate1,"./result/fate1.rds")
saveRDS(cds_fate2,"./result/fate2.rds")

### Differential expression analysis ### fate1

fate1 <- readRDS("./result/fate1.rds")

expr_matrix <- exprs(fate1)
pheno_feature <- as.data.frame(colData(fate1))
pheno_feature$Cluster <- pheno_feature$subtype
pheno_feature$Cluster <- gsub("Naive T",1,pheno_feature$Cluster)
pheno_feature$Cluster <- gsub("Tfh",2,pheno_feature$Cluster)
pheno_feature$Cluster <- gsub("Th17",3,pheno_feature$Cluster)
pheno_feature$Cluster <- gsub("Tcm",4,pheno_feature$Cluster)
pheno_feature$Cluster <- gsub("Tem",5,pheno_feature$Cluster)
pheno_feature$Cluster <- gsub("Treg",6,pheno_feature$Cluster)
colnames(pheno_feature)[19] <- "Pseudotime"
pd <- new("AnnotatedDataFrame", data = pheno_feature)
gene_feature <- data.frame(gene_short_name=rownames(expr_matrix))
rownames(gene_feature) <- rownames(expr_matrix)
fd <- new("AnnotatedDataFrame", data = gene_feature)

detach("package:monocle3")
library(monocle)
monocle2.fate1.object <- newCellDataSet(as.matrix(expr_matrix),phenoData = pd,featureData = fd)

monocle2.fate1.object <- estimateSizeFactors(monocle2.fate1.object)
monocle2.fate1.object <- estimateDispersions(monocle2.fate1.object)

diff_test_res <- differentialGeneTest(monocle2.fate1.object,fullModelFormulaStr = "~sm.ns(Pseudotime)")
sig_genes_all <- subset(diff_test_res, qval < 0.1)
sig_genes_all <- sig_genes_all[order(sig_genes_all$qval,decreasing = F),]
sig_gene_names <- row.names(sig_genes_all[1:20,])
# sig_gene_names <- row.names(subset(diff_test_res, qval < 0.1))
enrichment_genes_fate1 <- data.frame(symbol=row.names(sig_genes_all[1:50,]))

write.table(enrichment_genes_fate1,"./result/enrichment_genes_fate1.txt",row.names = F,quote = F)
saveRDS(monocle2.fate1.object,"./result/monocle2.fate1.object.rds")

library(pheatmap)
library(colorRamps)

x_monocle <- monocle2.fate1.object
newdata<-data.frame(Pseudotime = seq(min(pData(x_monocle)$Pseudotime), 
                                     max(pData(x_monocle)$Pseudotime), length.out = 100))
m <- genSmoothCurves(x_monocle[sig_gene_names,],
                                     new_data = newdata,
                                     cores = 10,
                                     trend_formula = "~sm.ns(Pseudotime, df=3)")
m = m[!apply(m, 1, sum) == 0, ]
m = log10(m + 1)
m = m[!apply(m, 1, sd) == 0, ]
m = Matrix::t(scale(Matrix::t(m), center = TRUE))
m = m[is.na(row.names(m)) == FALSE, ]
m[is.nan(m)] = 0
m[m > 3] = 3
m[m < -3] = -3
heatmap_matrix <- m
row_dist <- as.dist((1 - cor(Matrix::t(heatmap_matrix)))/2)
row_dist[is.na(row_dist)] <- 1
bks <- seq(-3.1, 3.1, by = 0.1)
hmcols <- blue2green2red(length(bks) - 1)
ph <- pheatmap(heatmap_matrix, useRaster = T, cluster_cols = FALSE, 
        cluster_rows = TRUE, show_rownames = F, show_colnames = F, 
        clustering_distance_rows = row_dist, clustering_method = "ward.D2", 
        cutree_rows = 6, silent = TRUE, filename = NA, 
        breaks = bks, border_color = NA, color = hmcols)
annotation_row <- data.frame(Cluster = factor(cutree(ph$tree_row, 
            6)))
annotation_col <- NA
feature_label <- as.character(fData(x_monocle[sig_gene_names,])[row.names(heatmap_matrix), 
       "gene_short_name"])
feature_label[is.na(feature_label)] <- row.names(heatmap_matrix)
row_ann_labels <- as.character(fData(x_monocle[sig_gene_names,])[row.names(annotation_row), 
       "gene_short_name"])
row_ann_labels[is.na(row_ann_labels)] <- row.names(annotation_row)
row_ann_labels <- row.names(annotation_row)
row.names(heatmap_matrix) <- feature_label
colnames(heatmap_matrix) <- c(1:ncol(heatmap_matrix))
ph_res <- pheatmap(heatmap_matrix[, ], useRaster = T, cluster_cols = FALSE, 
        cluster_rows = TRUE, show_rownames = TRUE, 
        show_colnames = F, clustering_distance_rows = row_dist, 
        clustering_method = "ward.D2", cutree_rows = 6, 
        annotation_row = annotation_row, annotation_col = annotation_col, 
        treeheight_row = 20, breaks = bks, fontsize = 6, color = hmcols, 
        border_color = NA, silent = TRUE, filename = NA)
# grid::grid.rect(gp = grid::gpar("fill", col = NA))
# grid::grid.draw(ph_res$gtable)

ggsave("./result/Plots/fate1-pseudotime-heatmap.pdf",ph_res)

### fate2

fate2 <- readRDS("./result/fate2.rds")

expr_matrix <- exprs(fate2)
pheno_feature <- as.data.frame(colData(fate2))
pheno_feature$Cluster <- pheno_feature$subtype
pheno_feature$Cluster <- gsub("Naive T",1,pheno_feature$Cluster)
pheno_feature$Cluster <- gsub("Tfh",2,pheno_feature$Cluster)
pheno_feature$Cluster <- gsub("Th17",3,pheno_feature$Cluster)
pheno_feature$Cluster <- gsub("Tcm",4,pheno_feature$Cluster)
pheno_feature$Cluster <- gsub("Tem",5,pheno_feature$Cluster)
pheno_feature$Cluster <- gsub("Treg",6,pheno_feature$Cluster)
colnames(pheno_feature)[19] <- "Pseudotime"
pd <- new("AnnotatedDataFrame", data = pheno_feature)
gene_feature <- data.frame(gene_short_name=rownames(expr_matrix))
rownames(gene_feature) <- rownames(expr_matrix)
fd <- new("AnnotatedDataFrame", data = gene_feature)

detach("package:monocle3")
library(monocle)
monocle2.fate2.object <- newCellDataSet(as.matrix(expr_matrix),phenoData = pd,featureData = fd)

monocle2.fate2.object <- estimateSizeFactors(monocle2.fate2.object)
monocle2.fate2.object <- estimateDispersions(monocle2.fate2.object)

diff_test_res <- differentialGeneTest(monocle2.fate2.object,fullModelFormulaStr = "~sm.ns(Pseudotime)")
sig_genes_all <- subset(diff_test_res, qval < 0.1)
sig_genes_all <- sig_genes_all[order(sig_genes_all$qval,decreasing = F),]
sig_gene_names <- row.names(sig_genes_all[1:20,])
# sig_gene_names <- row.names(subset(diff_test_res, qval < 0.1))
enrichment_genes_fate2 <- row.names(sig_genes_all[1:50,])

write.table(enrichment_genes_fate2,"./result/enrichment_genes_fate2.txt",row.names = F,quote = F)
saveRDS(monocle2.fate2.object,"./result/monocle2.fate2.object.rds")

library(pheatmap)
library(colorRamps)

x_monocle <- monocle2.fate2.object
newdata<-data.frame(Pseudotime = seq(min(pData(x_monocle)$Pseudotime), 
                                     max(pData(x_monocle)$Pseudotime), length.out = 100))
m <- genSmoothCurves(x_monocle[sig_gene_names,],
                                     new_data = newdata,
                                     cores = 10,
                                     trend_formula = "~sm.ns(Pseudotime, df=3)")
m = m[!apply(m, 1, sum) == 0, ]
m = log10(m + 1)
m = m[!apply(m, 1, sd) == 0, ]
m = Matrix::t(scale(Matrix::t(m), center = TRUE))
m = m[is.na(row.names(m)) == FALSE, ]
m[is.nan(m)] = 0
m[m > 3] = 3
m[m < -3] = -3
heatmap_matrix <- m
row_dist <- as.dist((1 - cor(Matrix::t(heatmap_matrix)))/2)
row_dist[is.na(row_dist)] <- 1
bks <- seq(-3.1, 3.1, by = 0.1)
hmcols <- blue2green2red(length(bks) - 1)
ph <- pheatmap(heatmap_matrix, useRaster = T, cluster_cols = FALSE, 
        cluster_rows = TRUE, show_rownames = F, show_colnames = F, 
        clustering_distance_rows = row_dist, clustering_method = "ward.D2", 
        cutree_rows = 6, silent = TRUE, filename = NA, 
        breaks = bks, border_color = NA, color = hmcols)
annotation_row <- data.frame(Cluster = factor(cutree(ph$tree_row, 
            6)))
annotation_col <- NA
feature_label <- as.character(fData(x_monocle[sig_gene_names,])[row.names(heatmap_matrix), 
       "gene_short_name"])
feature_label[is.na(feature_label)] <- row.names(heatmap_matrix)
row_ann_labels <- as.character(fData(x_monocle[sig_gene_names,])[row.names(annotation_row), 
       "gene_short_name"])
row_ann_labels[is.na(row_ann_labels)] <- row.names(annotation_row)
row_ann_labels <- row.names(annotation_row)
row.names(heatmap_matrix) <- feature_label
colnames(heatmap_matrix) <- c(1:ncol(heatmap_matrix))
ph_res <- pheatmap(heatmap_matrix[, ], useRaster = T, cluster_cols = FALSE, 
        cluster_rows = TRUE, show_rownames = TRUE, 
        show_colnames = F, clustering_distance_rows = row_dist, 
        clustering_method = "ward.D2", cutree_rows = 6, 
        annotation_row = annotation_row, annotation_col = annotation_col, 
        treeheight_row = 20, breaks = bks, fontsize = 6, color = hmcols, 
        border_color = NA, silent = TRUE, filename = NA)
# grid::grid.rect(gp = grid::gpar("fill", col = NA))
# grid::grid.draw(ph_res$gtable)

ggsave("./result/Plots/fate2-pseudotime-heatmap.pdf",ph_res)

### GO enrichment ###
go_enrichment <- function(genelist){
       library(clusterProfiler)
       library(org.Hs.eg.db)

       genelist = genelist$symbol
       ego = enrichGO(gene          = genelist,
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                keyType = "SYMBOL",
                readable      = TRUE)
       return(ego)
}
library(ggplot2)

enrichment_genes_fate1 <- read.table("./result/enrichment_genes_fate1.txt",header = T)
ego1 <- go_enrichment(enrichment_genes_fate1)
p1 <- clusterProfiler::dotplot(ego1,showCategory = 5) + ggtitle("fate1")

enrichment_genes_fate2 <- read.table("./result/enrichment_genes_fate2.txt",header = T)
ego2 <- go_enrichment(enrichment_genes_fate2)
p2 <- clusterProfiler::dotplot(ego2,showCategory = 5) + ggtitle("fate2")

library(xlsx)
write.xlsx(ego1,"./result/fate1_fate2_GO_enrichment_all_result.xlsx",sheetName = "fate1")
write.xlsx(ego2,"./result/fate1_fate2_GO_enrichment_all_result.xlsx",sheetName = "fate2",append = TRUE)

ggsave("./result/Plots/fate1_fate2_GO_enrichment-top5.pdf",p1+p2,width=10,height=6)

### ###

``
