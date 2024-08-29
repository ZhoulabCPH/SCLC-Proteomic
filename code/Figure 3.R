#####   protein-protein analysis   #####
data <- read.csv("dreamai_134SCLC_7715protein_normalized_and_log2.csv",
                 stringsAsFactors = F,check.names = F,
                 row.names = 1)
data.new <- data[,c(-1,-2)]
data.new <- t(data.new)
protein_cor <- cor(data.new, method = c("spearman"))
#write.csv(protein_cor,"protein_cor_matrix.csv",quote = F)

###  3A and 3B  ###
library(Seurat)
library(tidyverse)
library(ggplot2)
data <- read.csv("protein_cor_matrix.csv",stringsAsFactors = F,
                 row.names = 1,check.names = F)
seurat_object  <- CreateSeuratObject(counts = data)
seurat_object  <- FindVariableFeatures(seurat_object, 
                                       selection.method = "vst")
top10 <- head(VariableFeatures(seurat_object), 10)
plot1 <- VariableFeaturePlot(seurat_object)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
seurat_object@assays$RNA$scale.data <- as.matrix(data[VariableFeatures(seurat_object),])
seurat_object  <- RunPCA(seurat_object , 
                         features = VariableFeatures(seurat_object))
plot1 <- DimPlot(seurat_object , reduction = "pca",group.by="orig.ident")
plot2 <- ElbowPlot(seurat_object , ndims=30, reduction="pca")
plotc <- plot1+plot2
seurat_object  <- FindNeighbors(seurat_object , dims = 1:10)
seurat_object  <- FindClusters(seurat_object , resolution = 0.5)
seurat_object  <- RunUMAP(seurat_object , dims = 1:10)
seurat_object  <- RunTSNE(seurat_object , dims = 1:10)
DimPlot(seurat_object , reduction = "umap",label = TRUE,raster=FALSE)
DimPlot(seurat_object , reduction = "tsne",label = TRUE,raster=FALSE,
        cols = c("#FFD384","#F5C0C0","#966C3B","#A4EBF3",
                 "#A1CAE2","#70AF85","#CD5D7D"))
metadata <- seurat_object@meta.data
cell_names <- rownames(metadata)
clusters <- metadata$seurat_clusters 
combined_data <- data.frame(CellName = cell_names, Cluster = clusters)
data <- read.csv("dreamai_134SCLC_7715protein_raw_intensity.csv",
                 stringsAsFactors = F,check.names = F,
                 row.names = 1)
data.new <- data[,c(-1,-2)]
cluster <- read.csv("top10_newcluster.csv",
                    stringsAsFactors = F,row.names = 1,check.names = F)
data.new <- data.new[,rownames(cluster)]
c1 <- rownames(cluster)[which(cluster$consensus_cluster == "Cluster1")]
c2 <- rownames(cluster)[which(cluster$consensus_cluster == "Cluster2")]
data.new_C1 <- data.new[,c1]
data.new_C2 <- data.new[,c2]
##  calculated FC
FC1_median <- apply(data.new_C1,1,median)/apply(data.new_C2,1,median)
FC1_median[FC1_median > 1.5] <- 1.5
FC2_median <- apply(data.new_C2,1,median)/apply(data.new_C1,1,median)
FC2_median[FC2_median > 1.5] <- 1.5
##  add FC in data
seurat_object <- AddMetaData(seurat_object, metadata = FC1_median, col.name = "FC1_median")
seurat_object <- AddMetaData(seurat_object, metadata = FC2_median, col.name = "FC2_median")
FeaturePlot(seurat_object, reduction = "tsne",features = "FC1_median")+
  scale_colour_gradient2(low = "#0A0A0A", mid = "#f7f7f7", high = "#FAAC7E",
                         midpoint = 1)
FeaturePlot(seurat_object, reduction = "tsne",features = "FC2_median")+
  scale_colour_gradient2(low = "#0A0A0A", mid = "#f7f7f7", high = "#FAAC7E",
                         midpoint = 1)


###   3C  ###
library(org.Hs.eg.db)
library(clusterProfiler)
library(pathview)
library(enrichplot)
library(dplyr)
library(msigdbr)
data <- read.csv("dreamai_134SCLC_7715protein_normalized_and_log2.csv",
                 stringsAsFactors = F,check.names = F,
                 row.names = 1)
protein_lable <- read.csv("module_cluster_singlecell.csv",
                          stringsAsFactors = F,row.names = 1,check.names = F)
protein_lable$symbol <- data[rownames(protein_lable),2]
module1 <- protein_lable[which(protein_lable$Cluster == 1),]
module2 <- protein_lable[which(protein_lable$Cluster == 2),]
module3 <- protein_lable[which(protein_lable$Cluster == 3),]
module4 <- protein_lable[which(protein_lable$Cluster == 4),]
module5 <- protein_lable[which(protein_lable$Cluster == 5),]
module6 <- protein_lable[which(protein_lable$Cluster == 6),]
module7 <- protein_lable[which(protein_lable$Cluster == 7),]
m_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, entrez_gene)
module1 <- bitr(module1$symbol,fromType="SYMBOL",toType=c("ENTREZID"),OrgDb="org.Hs.eg.db")
module2 <- bitr(module2$symbol,fromType="SYMBOL",toType=c("ENTREZID"),OrgDb="org.Hs.eg.db")
module3 <- bitr(module3$symbol,fromType="SYMBOL",toType=c("ENTREZID"),OrgDb="org.Hs.eg.db")
module4 <- bitr(module4$symbol,fromType="SYMBOL",toType=c("ENTREZID"),OrgDb="org.Hs.eg.db")
module5 <- bitr(module5$symbol,fromType="SYMBOL",toType=c("ENTREZID"),OrgDb="org.Hs.eg.db")
module6 <- bitr(module6$symbol,fromType="SYMBOL",toType=c("ENTREZID"),OrgDb="org.Hs.eg.db")
module7 <- bitr(module7$symbol,fromType="SYMBOL",toType=c("ENTREZID"),OrgDb="org.Hs.eg.db")
module1_hall <- enricher(module1$ENTREZID, TERM2GENE=m_t2g,
                         pvalueCutoff =1,pAdjustMethod = "fdr",
                         universe = NULL,qvalueCutoff = 1,
                         minGSSize = 10,maxGSSize = 500)
module2_hall <- enricher(module2$ENTREZID, TERM2GENE=m_t2g,
                         pvalueCutoff =1,pAdjustMethod = "fdr",
                         universe = NULL,qvalueCutoff = 1,
                         minGSSize = 10,maxGSSize = 500)
module3_hall <- enricher(module3$ENTREZID, TERM2GENE=m_t2g,
                         pvalueCutoff =1,pAdjustMethod = "fdr",
                         universe = NULL,qvalueCutoff = 1,
                         minGSSize = 10,maxGSSize = 500)
module4_hall <- enricher(module4$ENTREZID, TERM2GENE=m_t2g,
                         pvalueCutoff =1,pAdjustMethod = "fdr",
                         universe = NULL,qvalueCutoff = 1,
                         minGSSize = 10,maxGSSize = 500)
module5_hall <- enricher(module5$ENTREZID, TERM2GENE=m_t2g,
                         pvalueCutoff =1,pAdjustMethod = "fdr",
                         universe = NULL,qvalueCutoff = 1,
                         minGSSize = 10,maxGSSize = 500)
module6_hall <- enricher(module6$ENTREZID, TERM2GENE=m_t2g,
                         pvalueCutoff =1,pAdjustMethod = "fdr",
                         universe = NULL,qvalueCutoff = 1,
                         minGSSize = 10,maxGSSize = 500)
module7_hall <- enricher(module7$ENTREZID, TERM2GENE=m_t2g,
                         pvalueCutoff =1,pAdjustMethod = "fdr",
                         universe = NULL,qvalueCutoff = 1,
                         minGSSize = 10,maxGSSize = 500)
module1_hall <- as.data.frame(module1_hall)
module2_hall <- as.data.frame(module2_hall)
module3_hall <- as.data.frame(module3_hall)
module4_hall <- as.data.frame(module4_hall)
module5_hall <- as.data.frame(module5_hall)
module6_hall <- as.data.frame(module6_hall)
module7_hall <- as.data.frame(module7_hall)
hall_type <- read.csv("hallmarker_6type.csv",stringsAsFactors = F,
                      row.names = 1,check.names = F)
module_hall_GeneRatio <- matrix(,50,7)
rownames(module_hall_GeneRatio) <- rownames(hall_type)
module_hall_P <- matrix(,50,7)
rownames(module_hall_P) <- rownames(hall_type)
for (i in 1:50){
  weizhi1 <- which(rownames(module_hall)[i] == module1_hall$ID)
  weizhi2 <- which(rownames(module_hall)[i] == module2_hall$ID)
  weizhi3 <- which(rownames(module_hall)[i] == module3_hall$ID)
  weizhi4 <- which(rownames(module_hall)[i] == module4_hall$ID)
  weizhi5 <- which(rownames(module_hall)[i] == module5_hall$ID)
  weizhi6 <- which(rownames(module_hall)[i] == module6_hall$ID)
  weizhi7 <- which(rownames(module_hall)[i] == module7_hall$ID)
  if (length(weizhi1) > 0){
    module_hall_GeneRatio[i,1] <- module1_hall[weizhi1,"Count"]
    module_hall_P[i,1] <- module1_hall[weizhi1,"pvalue"]
  }
  else {
    module_hall_GeneRatio[i,1] <- NA
    module_hall_P[i,1] <- NA
  }
  if (length(weizhi2) > 0){
    module_hall_GeneRatio[i,2] <- module2_hall[weizhi2,"Count"]
    module_hall_P[i,2] <- module2_hall[weizhi2,"pvalue"]
  }
  else {
    module_hall_GeneRatio[i,2] <- NA
    module_hall_P[i,2] <- NA
  }
  if (length(weizhi3) > 0){
    module_hall_GeneRatio[i,3] <- module3_hall[weizhi3,"Count"]
    module_hall_P[i,3] <- module3_hall[weizhi3,"pvalue"]
  }  
  else {
    module_hall_GeneRatio[i,3] <- NA
    module_hall_P[i,3] <- NA
  }
  if (length(weizhi4) > 0){
    module_hall_GeneRatio[i,4] <- module4_hall[weizhi4,"Count"]
    module_hall_P[i,4] <- module4_hall[weizhi4,"pvalue"]
  }
  else {
    module_hall_GeneRatio[i,4] <- NA
    module_hall_P[i,4] <- NA
  }
  if (length(weizhi5) > 0){
    module_hall_GeneRatio[i,5] <- module5_hall[weizhi5,"Count"]
    module_hall_P[i,5] <- module5_hall[weizhi5,"pvalue"]
  }
  else {
    module_hall_GeneRatio[i,5] <- NA
    module_hall_P[i,5] <- NA
  }
  if (length(weizhi6) > 0){
    module_hall_GeneRatio[i,6] <- module6_hall[weizhi6,"Count"]
    module_hall_P[i,6] <- module6_hall[weizhi6,"pvalue"]
  }
  else {
    module_hall_GeneRatio[i,6] <- NA
    module_hall_P[i,6] <- NA
  }
  if (length(weizhi7) > 0){
    module_hall_GeneRatio[i,7] <- module7_hall[weizhi7,"Count"]
    module_hall_P[i,7] <- module7_hall[weizhi7,"pvalue"]
  }
  else {
    module_hall_GeneRatio[i,7] <- NA
    module_hall_P[i,7] <- NA
  }
}
module_hall_P[,1] <- p.adjust(module_hall_P[,1],method = "fdr")
module_hall_P[,2] <- p.adjust(module_hall_P[,2],method = "fdr")
module_hall_P[,3] <- p.adjust(module_hall_P[,3],method = "fdr")
module_hall_P[,4] <- p.adjust(module_hall_P[,4],method = "fdr")
module_hall_P[,5] <- p.adjust(module_hall_P[,5],method = "fdr")
module_hall_P[,6] <- p.adjust(module_hall_P[,6],method = "fdr")
module_hall_P[,7] <- p.adjust(module_hall_P[,7],method = "fdr")
colnames(module_hall_P) <- c("Module1","Module2","Module3","Module4","Module5","Module6","Module7")
colnames(module_hall_GeneRatio) <- c("Module1","Module2","Module3","Module4","Module5","Module6","Module7")
module_hall_P_level <- module_hall_P
module_hall_P_level[module_hall_P > 0.05] <- "Un-enriched"
module_hall_P_level[module_hall_P < 0.05] <- "<0.05"
module_hall_P_level[module_hall_P < 0.01] <- "<0.01"
module_hall_P_level[module_hall_P < 0.001] <- "<0.001"
module_hall_P <- t(module_hall_P)
module_hall_GeneRatio <- t(module_hall_GeneRatio)
module_hall_P_level <- as.data.frame(module_hall_P_level)
annotation_col <- data.frame(Module1 = module_hall_P_level$Module1,
                             Module2 = module_hall_P_level$Module2,
                             Module3 = module_hall_P_level$Module3,
                             Module4 = module_hall_P_level$Module4,
                             Module5 = module_hall_P_level$Module5,
                             Module6 = module_hall_P_level$Module6,
                             Module7 = module_hall_P_level$Module7)
rownames(annotation_col) <- rownames(module_hall_P_level)
ann_colors = list(Module1 = c("<0.001" = "#FAAB78","<0.01" = "#FFDCA9",
                              "<0.05" = "#FCF9BE","Un-enriched" = "#FFFFFF"),
                  Module2 = c("<0.001" = "#FAAB78","<0.01" = "#FFDCA9",
                              "<0.05" = "#FCF9BE","Un-enriched" = "#FFFFFF"),
                  Module3 = c("<0.001" = "#FAAB78","<0.01" = "#FFDCA9",
                              "<0.05" = "#FCF9BE","Un-enriched" = "#FFFFFF"),
                  Module4 = c("<0.001" = "#FAAB78","<0.01" = "#FFDCA9",
                              "<0.05" = "#FCF9BE","Un-enriched" = "#FFFFFF"),
                  Module5 = c("<0.001" = "#FAAB78","<0.01" = "#FFDCA9",
                              "<0.05" = "#FCF9BE","Un-enriched" = "#FFFFFF"),
                  Module6 = c("<0.001" = "#FAAB78","<0.01" = "#FFDCA9",
                              "<0.05" = "#FCF9BE","Un-enriched" = "#FFFFFF"),
                  Module7 = c("<0.001" = "#FAAB78","<0.01" = "#FFDCA9",
                              "<0.05" = "#FCF9BE","Un-enriched" = "#FFFFFF"))
module_hall_P <- module_hall_P[,rownames(hall_type)]
module_hall_GeneRatio <- module_hall_GeneRatio[,rownames(hall_type)]
pheatmap(module_hall_GeneRatio,fontsize=6,cellheight = 10,
         color  = colorRampPalette(c("#FFFFFF","#F8C4B4","#F25F5C"))(100),
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         clustering_method = "ward.D2",
         border_color = "grey60",
         cluster_cols = F, cluster_rows = F,
         show_rownames = T, show_colnames = T)



###  3D  ###
cluster <- read.csv("top10_newcluster.csv",stringsAsFactors = F,
                    row.names = 1)
sub1_name <- rownames(cluster)[which(cluster$consensus_cluster == "Cluster1")]
sub2_name <- rownames(cluster)[which(cluster$consensus_cluster == "Cluster2")]
data_raw <- read.csv("dreamai_134SCLC_7715protein_raw_intensity.csv",
                     stringsAsFactors = F,check.names = F,row.names = 1)
data1 <- read.csv("dreamai_134SCLC_7715protein_normalized_and_log2.csv",
                 stringsAsFactors = F,check.names = F,row.names = 1)
data_raw <- data_raw[,-c(1,2)]
data <- data1[,-c(1,2)]
module <- read.csv("module_cluster_singlecell.csv")
module1 <- module[which(module$Cluster == 1),]
module3 <- module[which(module$Cluster == 3),]
module4 <- module[which(module$Cluster == 4),]
module6 <- module[which(module$Cluster == 6),]
data_raw_module1 <- data_raw[module1$CellName,]
data_raw_module3 <- data_raw[module3$CellName,]
data_raw_module4 <- data_raw[module4$CellName,]
data_raw_module6 <- data_raw[module6$CellName,]
data_module1 <- data[module1$CellName,]
data_module3 <- data[module3$CellName,]
data_module4 <- data[module4$CellName,]
data_module6 <- data[module6$CellName,]
data_raw_module1_subtype1 <- data_raw_module1[,sub1_name]
data_raw_module3_subtype1 <- data_raw_module3[,sub1_name]
data_raw_module4_subtype1 <- data_raw_module4[,sub1_name]
data_raw_module6_subtype1 <- data_raw_module6[,sub1_name]
data_raw_module1_subtype2 <- data_raw_module1[,sub2_name]
data_raw_module3_subtype2 <- data_raw_module3[,sub2_name]
data_raw_module4_subtype2 <- data_raw_module4[,sub2_name]
data_raw_module6_subtype2 <- data_raw_module6[,sub2_name]
data_module1_subtype1 <- data_module1[,sub1_name]
data_module3_subtype1 <- data_module3[,sub1_name]
data_module4_subtype1 <- data_module4[,sub1_name]
data_module6_subtype1 <- data_module6[,sub1_name]
data_module1_subtype2 <- data_module1[,sub2_name]
data_module3_subtype2 <- data_module3[,sub2_name]
data_module4_subtype2 <- data_module4[,sub2_name]
data_module6_subtype2 <- data_module6[,sub2_name]
module1_fc <- c()
module1_p <- c()
module3_fc <- c()
module3_p <- c()
module4_fc <- c()
module4_p <- c()
module6_fc <- c()
module6_p <- c()
for (i in 1:dim(data_module1_subtype1)[1]){
  fc <- rowMeans(data_raw_module1_subtype1[i,])/rowMeans(data_raw_module1_subtype2[i,])
  p <- wilcox.test(as.numeric(data_module1_subtype1[i,]),
                   as.numeric(data_module1_subtype2[i,]))
  module1_p <- c(module1_p,p$p.value)
  module1_fc <- c(module1_fc,fc)
}
for (i in 1:dim(data_module3_subtype1)[1]){
  fc <- rowMeans(data_raw_module3_subtype1[i,])/rowMeans(data_raw_module3_subtype2[i,])
  p <- wilcox.test(as.numeric(data_module3_subtype1[i,]),
                   as.numeric(data_module3_subtype2[i,]))
  module3_p <- c(module3_p,p$p.value)
  module3_fc <- c(module3_fc,fc)
}
for (i in 1:dim(data_module4_subtype1)[1]){
  fc <- rowMeans(data_raw_module4_subtype1[i,])/rowMeans(data_raw_module4_subtype2[i,])
  p <- wilcox.test(as.numeric(data_module4_subtype1[i,]),
                   as.numeric(data_module4_subtype2[i,]))
  module4_p <- c(module4_p,p$p.value)
  module4_fc <- c(module4_fc,fc)
}
for (i in 1:dim(data_module6_subtype1)[1]){
  fc <- rowMeans(data_raw_module6_subtype1[i,])/rowMeans(data_raw_module6_subtype2[i,])
  p <- wilcox.test(as.numeric(data_module6_subtype1[i,]),
                   as.numeric(data_module6_subtype2[i,]))
  module6_p <- c(module6_p,p$p.value)
  module6_fc <- c(module6_fc,fc)
}
module1_result <- data.frame(Protein = rownames(data_module1_subtype1),
                             Symbol = data1[rownames(data_module1_subtype1),2],
                             FC = module1_fc,
                             P = module1_p,
                             adjustP = p.adjust(module1_p,method = "fdr"),
                             sig = module1_p)
module1_subtype1_up <- intersect(which(module1_result$FC > 1.5),
                                 which(module1_result$adjustP < 0.05))
module1_subtype2_up <- intersect(which(module1_result$FC < 2/3),
                                 which(module1_result$adjustP < 0.05))
module1_result$sig[module1_subtype1_up] <- "Subtype1 UP"
module1_result$sig[module1_subtype2_up] <- "Subtype2 UP"
module1_result$sig[-c(module1_subtype1_up,module1_subtype2_up)] <- "no-sig"
module3_result <- data.frame(Protein = rownames(data_module3_subtype1),
                             Symbol = data1[rownames(data_module3_subtype1),2],
                             FC = module3_fc,
                             P = module3_p,
                             adjustP = p.adjust(module3_p,method = "fdr"),
                             sig = module3_p)
module3_subtype1_up <- intersect(which(module3_result$FC > 1.5),
                                 which(module3_result$adjustP < 0.05))
module3_subtype2_up <- intersect(which(module3_result$FC < 2/3),
                                 which(module3_result$adjustP < 0.05))
module3_result$sig[module3_subtype1_up] <- "Subtype1 UP"
module3_result$sig[module3_subtype2_up] <- "Subtype2 UP"
module3_result$sig[-c(module3_subtype1_up,module3_subtype2_up)] <- "no-sig"
module4_result <- data.frame(Protein = rownames(data_module4_subtype1),
                             Symbol = data1[rownames(data_module4_subtype1),2],
                             FC = module4_fc,
                             P = module4_p,
                             adjustP = p.adjust(module4_p,method = "fdr"),
                             sig = module4_p)
module4_subtype1_up <- intersect(which(module4_result$FC > 1.5),
                                 which(module4_result$adjustP < 0.05))
module4_subtype2_up <- intersect(which(module4_result$FC < 2/3),
                                 which(module4_result$adjustP < 0.05))
module4_result$sig[module4_subtype1_up] <- "Subtype1 UP"
module4_result$sig[module4_subtype2_up] <- "Subtype2 UP"
module4_result$sig[-c(module4_subtype1_up,module4_subtype2_up)] <- "no-sig"
module6_result <- data.frame(Protein = rownames(data_module6_subtype1),
                             Symbol = data1[rownames(data_module6_subtype1),2],
                             FC = module6_fc,
                             P = module6_p,
                             adjustP = p.adjust(module6_p,method = "fdr"),
                             sig = module6_p)
module6_subtype1_up <- intersect(which(module6_result$FC > 1.5),
                                 which(module6_result$adjustP < 0.05))
module6_subtype2_up <- intersect(which(module6_result$FC < 2/3),
                                 which(module6_result$adjustP < 0.05))
module6_result$sig[module6_subtype1_up] <- "Subtype1 UP"
module6_result$sig[module6_subtype2_up] <- "Subtype2 UP"
module6_result$sig[-c(module6_subtype1_up,module6_subtype2_up)] <- "no-sig"
library(ggnewscale)
library(tidyverse)
library(ggplot2)
library(ggrepel)
setwd("D:\\北京小细胞肺癌_蛋白质组\\结果3")
module1 <- read.csv("module1_different_protein.csv",stringsAsFactors = F,
                    row.names = 1,check.names = F)
module3 <- read.csv("module3_different_protein.csv",stringsAsFactors = F,
                    row.names = 1,check.names = F)
module4 <- read.csv("module4_different_protein.csv",stringsAsFactors = F,
                    row.names = 1,check.names = F)
module6 <- read.csv("module6_different_protein.csv",stringsAsFactors = F,
                    row.names = 1,check.names = F)
module1$log2FC <- log2(module1$FC)
module3$log2FC <- log2(module3$FC)
module4$log2FC <- log2(module4$FC)
module6$log2FC <- log2(module6$FC)
### top 10 sig
module1.new <- module1[c(which(module1$sig == "Subtype2 UP"),
                         which(module1$sig == "Subtype1 UP")),]
module1.new <- module1.new[order(module1.new$log2FC),]
top10.module1.new <- module1.new[c(1:5,876:880),]
top10.module1.new$label <- rep("Module1",10)
module3.new <- module3[c(which(module3$sig == "Subtype2 UP"),
                         which(module3$sig == "Subtype1 UP")),]
module3.new <- module3.new[order(module3.new$log2FC),]
top10.module3.new <- module3.new[c(1:5,188:192),]
top10.module3.new$label <- rep("Module3",10)

module4.new <- module4[c(which(module4$sig == "Subtype2 UP"),
                         which(module4$sig == "Subtype1 UP")),]
module4.new <- module4.new[order(module4.new$log2FC),]
top10.module4.new <- module4.new[c(1:5,140:144),]
top10.module4.new$label <- rep("Module4",10)
module6.new <- module6[c(which(module6$sig == "Subtype2 UP"),
                         which(module6$sig == "Subtype1 UP")),]
module6.new <- module6.new[order(module6.new$log2FC),]
top10.module6.new <- module4.new[c(1:5,8:12),]
top10.module6.new$label <- rep("Module6",10)
module1.plot <- module1[which(module1$sig != "no-sig"),] 
module1.plot$label <- rep("Module1",dim(module1.plot)[1])
module1.plot$ptname <- rownames(module1.plot)
module3.plot <- module3[which(module3$sig != "no-sig"),] 
module3.plot$label <- rep("Module3",dim(module3.plot)[1])
module3.plot$ptname <- rownames(module3.plot)
module4.plot <- module4[which(module4$sig != "no-sig"),] 
module4.plot$label <- rep("Module4",dim(module4.plot)[1])
module4.plot$ptname <- rownames(module4.plot)
module6.plot <- module6[which(module6$sig != "no-sig"),] 
module6.plot$label <- rep("Module6",dim(module6.plot)[1])
module6.plot$ptname <- rownames(module6.plot)
merge.plot <- rbind(module1.plot,module3.plot,
                    module4.plot,module6.plot)
merge.plot.top10 <- rbind(top10.module1.new,top10.module3.new,
                          top10.module4.new,top10.module6.new)
merge.plot.top10$gene <- rownames(merge.plot.top10)
dbar <- 
  merge.plot.top10 %>% 
  group_by(label) %>% 
  summarise_all(list(min = min, max = max)) %>% 
  select(label, log2FC_min, log2FC_max, label) %>% 
  rename(label = label)
ggplot()+
  geom_col(data = dbar,  # 
           mapping = aes(x = label,y = log2FC_min),
           fill = "#dcdcdc",alpha = 0.6, width = 0.7) +
  geom_col(data = dbar, # 
           mapping = aes(x = label,y = log2FC_max),
           fill = "#dcdcdc",alpha = 0.6, width = 0.7) +
  geom_jitter(data=merge.plot,
              aes(x=label,y=log2FC,color=sig),width =0.35) +
  geom_tile(data = merge.plot.top10, # 
            aes(x = label,y = 0,fill = label),
            height=0.5,color = "black",alpha = 0.6,
            show.legend = F) +
  scale_fill_manual(values = c('#F5C0C0','#A4EBF3','#A1CAE2','#CD5D7D'))+
  ggsci::scale_color_npg() + # 自定义颜色
  scale_color_manual(values = c('#E76F51','#AEB6FF')) +
  ylab("log2 (Fold-change)")



###  3E  ###
library(org.Hs.eg.db)
library(clusterProfiler)
library(pathview)
library(enrichplot)
library(dplyr)
library(msigdbr)
module1 <- read.csv("module1_different_protein.csv",stringsAsFactors = F,
                    row.names = 1,check.names = F)
module3 <- read.csv("module3_different_protein.csv",stringsAsFactors = F,
                    row.names = 1,check.names = F)
module4 <- read.csv("module4_different_protein.csv",stringsAsFactors = F,
                    row.names = 1,check.names = F)
module6 <- read.csv("module6_different_protein.csv",stringsAsFactors = F,
                    row.names = 1,check.names = F)
subtype1_up <- c(module4$Symbol[which(module4$sig == "Subtype1 UP")],
                 module6$Symbol[which(module6$sig == "Subtype1 UP")])
subtype2_up <- c(module1$Symbol[which(module1$sig == "Subtype2 UP")],
                 module3$Symbol[which(module3$sig == "Subtype2 UP")])
###  SCLC subtype1
gene1 <- bitr(subtype1_up, 
              fromType="SYMBOL", 
              toType=c("ENTREZID"),
              OrgDb="org.Hs.eg.db")
kegg_subtype1 <- enrichKEGG(gene = gene1$ENTREZID,
                            keyType = "kegg",
                            organism = 'hsa',
                            minGSSize = 1,
                            maxGSSize = 500,
                            pvalueCutoff = 1,
                            qvalueCutoff = 1,
                            pAdjustMethod = "fdr")
bp_subtype1 <- enrichGO(gene1$ENTREZID, 
                        OrgDb = org.Hs.eg.db, 
                        ont='BP', 
                        keyType = 'ENTREZID',
                        minGSSize = 1,
                        maxGSSize = 500,
                        pvalueCutoff = 1,
                        qvalueCutoff = 1,
                        pAdjustMethod = "fdr")
###  SCLC subtype2
gene2 <- bitr(subtype2_up, 
              fromType="SYMBOL", 
              toType=c("ENTREZID"),
              OrgDb="org.Hs.eg.db")
kegg_subtype2 <- enrichKEGG(gene = gene2$ENTREZID,
                            keyType = "kegg",
                            organism = 'hsa',
                            minGSSize = 1,
                            maxGSSize = 500,
                            pvalueCutoff = 1,
                            qvalueCutoff = 1,
                            pAdjustMethod = "fdr")
bp_subtype2 <- enrichGO(gene2$ENTREZID, 
                        OrgDb = org.Hs.eg.db, 
                        ont='BP', 
                        keyType = 'ENTREZID',
                        minGSSize = 1,
                        maxGSSize = 500,
                        pvalueCutoff = 1,
                        qvalueCutoff = 1,
                        pAdjustMethod = "fdr")
kegg_subtype1 <- as.data.frame(kegg_subtype1)
bp_subtype1 <- as.data.frame(bp_subtype1)
kegg_subtype2 <- as.data.frame(kegg_subtype2)
bp_subtype2 <- as.data.frame(bp_subtype2)
kegg_subtype1 <- kegg_subtype1[which(kegg_subtype1$p.adjust < 0.05),]
bp_subtype1 <- bp_subtype1[which(bp_subtype1$p.adjust < 0.05),]
kegg_subtype2 <- kegg_subtype2[which(kegg_subtype2$p.adjust < 0.05),]
bp_subtype2 <- bp_subtype2[which(bp_subtype2$p.adjust < 0.05),]
kegg_subtype1 <- kegg_subtype1[rev(order(kegg_subtype1$Count)),]
bp_subtype1 <- bp_subtype1[rev(order(bp_subtype1$Count)),]
kegg_subtype2 <- kegg_subtype2[rev(order(kegg_subtype2$Count)),]
bp_subtype2 <- bp_subtype2[rev(order(bp_subtype2$Count)),]
library(ggplot2)
bp <- read.csv("BP_total_result.csv",stringsAsFactors = F,
               row.names = 1,check.names = F)
kegg <- read.csv("KEGG_total_result.csv",stringsAsFactors = F,
                 row.names = 1,check.names = F)
bp$logp <- (-log10(bp$p.adjust))
kegg$logp <- (-log10(kegg$p.adjust))
bp$logp[which(bp$subtype == "subtype1")] <- (-bp$logp[which(bp$subtype == "subtype1")])
kegg$logp[which(kegg$subtype == "subtype1")] <- (-kegg$logp[which(kegg$subtype == "subtype1")])
bp <- bp[order(bp$logp),]
kegg <- kegg[order(kegg$logp),]
bp$name <- paste(bp$Description," (",bp$ID,")",sep = "")
kegg$name <- paste(kegg$Description," (",kegg$ID,")",sep = "")
bp$name <- factor(bp$name,levels = bp$name)
kegg$name <- factor(kegg$name,levels = kegg$name)
##  BP
ggplot(bp, aes(x=name, y=logp, fill=subtype)) + 
  geom_bar(stat="identity",  position=position_dodge()) + 
  scale_fill_manual(values = c("#E76F51","#AEB6FF"))+
  theme_bw() +
  theme(axis.text.x=element_text(hjust = 1,colour="black",size=10,angle = 90), 
        axis.text.y=element_text(size=10,colour="black"), 
        axis.title.y=element_text(size = 10,colour="black"), 
        legend.text=element_text(colour="black", 
                                 size=10),
        legend.title=element_text(colour="black", 
                                  size=10),
        panel.grid.major = element_blank(),   
        panel.grid.minor = element_blank())+
  ylab("-log10(FDR)")
##  KEGG
ggplot(kegg, aes(x=name, y=logp, fill=subtype)) + 
  geom_bar(stat="identity",  position=position_dodge()) + 
  scale_fill_manual(values = c("#E76F51","#AEB6FF"))+
  theme_bw() +
  theme(axis.text.x=element_text(hjust = 1,colour="black",size=10,angle = 90), 
        axis.text.y=element_text(size=10,colour="black"), 
        axis.title.y=element_text(size = 10,colour="black"), 
        legend.text=element_text(colour="black",  
                                 size=10),
        legend.title=element_text(colour="black",
                                  size=10),
        panel.grid.major = element_blank(),   
        panel.grid.minor = element_blank())+
  ylab("-log10(FDR)")


###  3F   ###
library(survival)
library(survminer)
library(ggplot2)
library(pheatmap)
cluster <- read.csv("top10_newcluster.csv",stringsAsFactors = F,
                    row.names = 1,check.names = F)
cluster <- cluster[order(cluster$consensus_cluster),]
pro <- c('ADH1B', 'ALDH2', 'ACSL5', 'ACADL', 'ACSL1', 'ACAA1', 
         'HADHA', 'ACADS', 'HADHB', 'ACSL4', 'ACADSB')
data_pro <- read.csv("dreamai_134SCLC_7715protein_normalized_and_log2.csv",
                     stringsAsFactors = F,check.names = F,
                     row.names = 1)
clinical <- read.csv("134SCLC_clinical.csv",
                     stringsAsFactors = F,row.names = 1)
data_fad <- data_pro[match(pro,data_pro$`Gene name`),]
rownames(data_fad) <- data_fad$`Gene name`
data_fad <- data_fad[,-c(1,2)]
clinical <- clinical[rownames(cluster),]
data_fad <- data_fad[,rownames(cluster)]
clinical <- cbind(clinical,t(data_fad))
linshi <- apply(data_fad,1,scale)
linshi <- t(linshi)
colnames(linshi) <- colnames(data_fad)
linshi[linshi > 2] <- 2
linshi[linshi < (-2)] <- (-2)
annotation_col <- data.frame(label = cluster$consensus_cluster)
rownames(annotation_col) <- rownames(cluster)
ann_colors = list(label = c("subtype1" = "#E76F51","subtype2" = "#AEB6FF"))
pheatmap(linshi,border_color = "NA",
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         color = colorRampPalette(colors = c("#9C4DCC","#FFFFFF","#FFD700"))(1000),
         clustering_method = "ward.D2",
         cluster_cols = F, cluster_rows = F,
         show_rownames = T, show_colnames = F)
# HR
###  OS 
clinical_os <- clinical[which(clinical$OS <= 60),]
cox <- matrix(,,7)
colnames(cox) <- c('lower .95','upper .95','coef', 'exp(coef)','se(coef)','z','Pr(>|z|)' )
for (i in 38:48){
  a <- summary(coxph(Surv(as.numeric(OS),OS.State)~clinical_os[,i],data = clinical_os))
  cox <- rbind(cox,c(a$conf.int[,3:4],a$coefficients[,1:5]))
}
cox <- cbind(cox,paste(round(cox[,4],3)," (",round(cox[,1],3)," - ",round(cox[,2],3),")",sep = ""))
cox[,7] <- round(as.numeric(cox[,7]),3)
cox <- cox[-1,]
rownames(cox) <- colnames(clinical_os)[38:48]
cox_os <- cox
###  DFS 
clinical_dfs <- clinical[which(clinical$DFS <= 57),]
cox <- matrix(,,7)
colnames(cox) <- c('lower .95','upper .95','coef', 'exp(coef)','se(coef)','z','Pr(>|z|)' )
for (i in 38:48){
  a <- summary(coxph(Surv(as.numeric(DFS),DFS.State)~clinical_dfs[,i],data = clinical_dfs))
  cox <- rbind(cox,c(a$conf.int[,3:4],a$coefficients[,1:5]))
}
cox <- cbind(cox,paste(round(cox[,4],3)," (",round(cox[,1],3)," - ",round(cox[,2],3),")",sep = ""))
cox[,7] <- round(as.numeric(cox[,7]),3)
cox <- cox[-1,]
rownames(cox) <- colnames(clinical_dfs)[38:48]
cox_dfs <- cox
