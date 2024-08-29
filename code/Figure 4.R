
###   4A and 4C   ###
library(pheatmap)
estimate <- read.csv("ESTIMATE_result.csv",stringsAsFactors = F,
                     row.names = 1,check.names = F)
xcell <- read.csv("xcell_result.csv",stringsAsFactors = F,
                  row.names = 1,check.names = F)
xcell <- xcell[,-1]
imm <- read.csv("imm_cell_ssgsea.csv",stringsAsFactors = F,
                   row.names = 1,check.names = F)
hall <- read.csv("hallmark_ssgsea.csv",stringsAsFactors = F,
                 row.names = 1,check.names = F)
tcell_related <- read.csv("T_ssgsea.csv",stringsAsFactors = F,
                          row.names = 1,check.names = F)
single <- read.csv("SCLC_single_cell_ssgsea.csv",stringsAsFactors = F,
                          row.names = 1,check.names = F)
fad_score <- read.csv("hsa00071_ssgsea.csv",
                      stringsAsFactors = F,row.names = 1,check.names = F)
fad_score <- as.data.frame(t(fad_score[,1:2]))
cluster <- read.csv("consensusCluster_label.csv",
                    stringsAsFactors = F,row.names = 1,check.names = F)
cluster <- cluster[order(cluster$immune_cluster),]
fad_score <- fad_score[,rownames(cluster)]
estimate <- estimate[,rownames(cluster)]
xcell <- xcell[,rownames(cluster)]
hall <- hall[,rownames(cluster)]
tcell_related <- tcell_related[,rownames(cluster)]
single <- single[,rownames(cluster)]
data <- rbind(estimate,xcell)
# FAD score
linshi <- apply(fad_score,1,scale)
linshi <- t(linshi)
colnames(linshi) <- colnames(fad_score)
annotation_col <- data.frame(imm = cluster$immune_cluster2,
                             sclc = cluster$consensus_cluster,
                             TF = cluster$TF_subtype)
rownames(annotation_col) <- rownames(cluster)
ann_colors = list(imm = c("IEX" = "#99CCFF","IIDS" = "#F580A6","IADS" = "#FF5733"),
                  sclc = c("Cluster1" = "#E76F51","Cluster2" = "#AEB6FF"),
                  TF = c("A" = "#EC8F5E","N" = "#F3B664","P" = "#F1EB90","Y" = "#9FBB73"))
linshi[linshi > 2] <- 2
linshi[linshi < (-2)] <- (-2)
pheatmap(linshi,
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         color  = colorRampPalette(c("#D9D9D9","#FFFFFF","#FF8080"))(100),
         clustering_method = "ward.D",
         border_color = "NA",
         cluster_cols = F, cluster_rows = F,
         show_rownames = T, show_colnames = F)
## immune
linshi <- apply(data,1,scale)
linshi <- t(linshi)
colnames(linshi) <- colnames(data)
annotation_col <- data.frame(imm = cluster$immune_cluster2,
                             sclc = cluster$consensus_cluster,
                             TF = cluster$TF_subtype)
rownames(annotation_col) <- rownames(cluster)
ann_colors = list(imm = c("IEX" = "#99CCFF","IIDS" = "#F580A6","IADS" = "#FF5733"),
                  sclc = c("Cluster1" = "#E76F51","Cluster2" = "#AEB6FF"),
                  TF = c("A" = "#EC8F5E","N" = "#F3B664","P" = "#F1EB90","Y" = "#9FBB73"))
linshi[linshi > 2] <- 2
linshi[linshi < (-2)] <- (-2)
pheatmap(linshi,
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         color  = colorRampPalette(c("#9C4DCC","#FFFFFF","#FFD700"))(100),
         clustering_method = "ward.D",
         border_color = "NA",
         cluster_cols = F, cluster_rows = F,
         show_rownames = T, show_colnames = F)
p <- c()
for (i in 1:dim(data)[1]){
  p_linshi <- kruskal.test(as.numeric(data[i,])~cluster$immune_cluster)
  p <- c(p,p_linshi$p.value)
}
cell_result <- data.frame(p = p,
                          cell = rownames(data))
cell_result$p.adj <- p.adjust(cell_result$p,method = "fdr")
sig_name <- cell_result$cell[which(cell_result$p.adj < 0.001)]
name <- c("ImmuneScore","StromalScore","ESTIMATEScore",
          "CD8+ T-cells","CD8+ naive T-cells",
          "CD4+ T-cells","CD4+ naive T-cells",
          "CD4+ memory T-cells","B-cells",
          "NKT","Macrophages","Macrophages M1",
          "Macrophages M2","DC","aDC","cDC",
          "Tregs","Neurons","Osteoblast","Endothelial cells",
          "Fibroblasts")
pheatmap(linshi[name,],gaps_col = c(63,102),
         #cellheight = 13,cellwidth = 2,
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         color  = colorRampPalette(c("#D9D9D9","#FFFFFF","#FF8080"))(100),
         clustering_method = "ward.D2",
         border_color = "NA",
         cluster_cols = F, cluster_rows = F,
         show_rownames = T, show_colnames = T)
pheatmap(linshi[c(name[4:21],"CD8+ Tcm","CD8+ Tem" ),],gaps_col = c(63,102),
         annotation_col = annotation_col,
         cellheight = 13,cellwidth = 2,
         annotation_colors = ann_colors,
         color  = colorRampPalette(c("#9C4DCC","#FFFFFF","#FFD700"))(100),
         clustering_method = "ward.D2",
         border_color = "NA",
         cluster_cols = F, cluster_rows = F,
         show_rownames = T, show_colnames = F)
cell_result[match(sig_name,cell_result$cell),]
# immune cell
imm <- imm[,colnames(linshi)]
linshi <- apply(imm,1,scale)
linshi <- t(linshi)
colnames(linshi) <- colnames(imm)
annotation_col <- data.frame(imm = cluster$immune_cluster,
                             sclc = cluster$consensus_cluster)
rownames(annotation_col) <- rownames(cluster)
ann_colors = list(imm = c("IS1" = "#99CCFF","IS2" = "#F8B195","IS3" = "#FF5733"),
                  sclc = c("Cluster1" = "#E76F51","Cluster2" = "#AEB6FF"))
linshi[linshi > 2] <- 2
linshi[linshi < (-2)] <- (-2)
pheatmap(linshi,gaps_col = c(63,102),
         annotation_col = annotation_col,
         #cellheight = 13,cellwidth = 2,
         annotation_colors = ann_colors,
         color  = colorRampPalette(c("#9C4DCC","#FFFFFF","#FFD700"))(100),
         clustering_method = "ward.D2",
         border_color = "NA",
         cluster_cols = F, cluster_rows = F,
         show_rownames = T, show_colnames = T)
# single cell
single <- single[,colnames(linshi)]
linshi_single <- apply(single,1,scale)
linshi_single <- t(linshi_single)
colnames(linshi_single) <- colnames(single)
annotation_col <- data.frame(imm = cluster$immune_cluster,
                             sclc = cluster$consensus_cluster)
rownames(annotation_col) <- rownames(cluster)
ann_colors = list(imm = c("IS1" = "#99CCFF","IS2" = "#F8B195","IS3" = "#FF5733"),
                  sclc = c("Cluster1" = "#E76F51","Cluster2" = "#AEB6FF"))
linshi_single[linshi_single > 2] <- 2
linshi_single[linshi_single < (-2)] <- (-2)
pheatmap(linshi_single,gaps_col = c(63,102),
         annotation_col = annotation_col,
         cellheight = 13,cellwidth = 2,
         annotation_colors = ann_colors,
         color  = colorRampPalette(c("#9C4DCC","#FFFFFF","#FFD700"))(100),
         clustering_method = "ward.D2",
         border_color = "NA",
         cluster_cols = F, cluster_rows = T,
         show_rownames = T, show_colnames = F)
p <- c()
for (i in 1:dim(single)[1]){
  p_linshi <- kruskal.test(as.numeric(single[i,])~cluster$immune_cluster)
  p <- c(p,p_linshi$p.value)
}
single_result <- data.frame(p = p,
                          single = rownames(single))
single_result$p.adj <- p.adjust(single_result$p,method = "fdr")
name <- single_result$single[which(single_result$p.adj < 0.05)]
pheatmap(linshi_single[name,],gaps_col = c(63,102),
         annotation_col = annotation_col,
         cellheight = 13,cellwidth = 2,
         annotation_colors = ann_colors,
         color  = colorRampPalette(c("#9C4DCC","#FFFFFF","#FFD700"))(100),
         clustering_method = "ward.D2",
         border_color = "NA",
         cluster_cols = F, cluster_rows = T,
         show_rownames = T, show_colnames = F)
# hallmark
hall <- hall[,colnames(linshi)]
linshi_hall <- apply(hall,1,scale)
linshi_hall <- t(linshi_hall)
colnames(linshi_hall) <- colnames(hall)
annotation_col <- data.frame(imm = cluster$immune_cluster,
                             sclc = cluster$consensus_cluster)
rownames(annotation_col) <- rownames(cluster)
ann_colors = list(imm = c("IS1" = "#99CCFF","IS2" = "#F8B195","IS3" = "#FF5733"),
                  sclc = c("Cluster1" = "#E76F51","Cluster2" = "#AEB6FF"))
linshi_hall[linshi_hall > 2] <- 2
linshi_hall[linshi_hall < (-2)] <- (-2)
pheatmap(linshi_hall[18:19,],gaps_col = c(63,102),
         annotation_col = annotation_col,
         cellheight = 13,cellwidth = 2,
         annotation_colors = ann_colors,
         color  = colorRampPalette(c("#9C4DCC","#FFFFFF","#FFD700"))(100),
         clustering_method = "ward.D2",
         border_color = "NA",
         cluster_cols = F, cluster_rows = F,
         show_rownames = T, show_colnames = F)
p <- c()
for (i in 1:dim(hall)[1]){
  p_linshi <- kruskal.test(as.numeric(hall[i,])~cluster$immune_cluster)
  p <- c(p,p_linshi$p.value)
}
hall_result <- data.frame(p = p,
                          hall = rownames(hall))
hall_result$p.adj <- p.adjust(hall_result$p,method = "fdr")
tcell_related <- tcell_related[,colnames(linshi)]
linshi_tcell_related <- apply(tcell_related,1,scale)
linshi_tcell_related <- t(linshi_tcell_related)
colnames(linshi_tcell_related) <- colnames(tcell_related)
annotation_col <- data.frame(imm = cluster$immune_cluster,
                             sclc = cluster$consensus_cluster)
rownames(annotation_col) <- rownames(cluster)
ann_colors = list(imm = c("IS1" = "#99CCFF","IS2" = "#F8B195","IS3" = "#FF5733"),
                  sclc = c("Cluster1" = "#E76F51","Cluster2" = "#AEB6FF"))
linshi_tcell_related[linshi_tcell_related > 2] <- 2
linshi_tcell_related[linshi_tcell_related < (-2)] <- (-2)
pheatmap(linshi_tcell_related,gaps_col = c(63,102),
         annotation_col = annotation_col,
         cellheight = 13,cellwidth = 2,
         annotation_colors = ann_colors,
         color  = colorRampPalette(c("#9C4DCC","#FFFFFF","#FFD700"))(100),
         clustering_method = "ward.D2",
         border_color = "NA",
         cluster_cols = F, cluster_rows = F,
         show_rownames = T, show_colnames = F)
p <- c()
for (i in 1:dim(tcell_related)[1]){
  p_linshi <- kruskal.test(as.numeric(tcell_related[i,])~cluster$immune_cluster)
  p <- c(p,p_linshi$p.value)
}
tcell_related_result <- data.frame(p = p,
                                   tcell_related = rownames(tcell_related))
tcell_related_result$p.adj <- p.adjust(tcell_related_result$p,
                                       method = "fdr")


###    4B    ###
library(ggplot2)
library(gghalves)
library(ggpubr)
estimate <- read.csv("ESTIMATE_result.csv",stringsAsFactors = F,
                     row.names = 1,check.names = F)
cluster <- read.csv("consensusCluster_label.csv",
                    stringsAsFactors = F,row.names = 1,check.names = F)
cluster <- cluster[order(cluster$immune_cluster),]
estimate <- as.data.frame(t(estimate)[rownames(cluster),])
estimate$imm <- cluster$immune_cluster
# ImmuneScore
a1 <- ggplot() +
  geom_half_boxplot(data = estimate, 
                    aes(x=imm, y=ImmuneScore, fill = imm),
                    width = 0.4, position = position_dodge(0.7), 
                    alpha=1, outlier.shape = NA)+
  geom_half_point(data = estimate, 
                  aes(x=imm, y=ImmuneScore, color=imm), 
                  transformation = position_jitter(width = 0.1,height = 0),
                  shape = 16, alpha=1, size = 0.8)+
  scale_fill_manual(values = c("#99CCFF","#F8B195","#FF5733"))+
  scale_color_manual(values = c("#99CCFF","#F8B195","#FF5733"))+
  stat_compare_means(data = estimate, 
                     aes(x=imm, y=ImmuneScore),
                     method = "kruskal")+   theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 0, vjust = 0.6),
        legend.position="right")+
  labs(y='Immune score')
# StromalScore
a2 <- ggplot() +
  geom_half_boxplot(data = estimate, 
                    aes(x=imm, y=StromalScore, fill = imm),
                    width = 0.4, position = position_dodge(0.7), 
                    alpha=1, outlier.shape = NA)+
  geom_half_point(data = estimate, 
                  aes(x=imm, y=StromalScore, color=imm), 
                  transformation = position_jitter(width = 0.1,height = 0),
                  shape = 16, alpha=1, size = 0.8)+
  scale_fill_manual(values = c("#99CCFF","#F8B195","#FF5733"))+
  scale_color_manual(values = c("#99CCFF","#F8B195","#FF5733"))+
  stat_compare_means(data = estimate, 
                     aes(x=imm, y=StromalScore),
                     method = "kruskal")+   theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 0, vjust = 0.6),
        legend.position="right")+
  labs(y='Stromal score')
# ESTIMATEScore
a3 <- ggplot() +
  geom_half_boxplot(data = estimate, 
                    aes(x=imm, y=ESTIMATEScore, fill = imm),
                    width = 0.4, position = position_dodge(0.7), 
                    alpha=1, outlier.shape = NA)+
  geom_half_point(data = estimate, 
                  aes(x=imm, y=ESTIMATEScore, color=imm), 
                  transformation = position_jitter(width = 0.1,height = 0),
                  shape = 16, alpha=1, size = 0.8)+
  scale_fill_manual(values = c("#99CCFF","#F8B195","#FF5733"))+
  scale_color_manual(values = c("#99CCFF","#F8B195","#FF5733"))+
  stat_compare_means(data = estimate, 
                     aes(x=imm, y=ESTIMATEScore),
                     method = "kruskal")+   theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 0, vjust = 0.6),
        legend.position="right")+
  labs(y='ESTIMATE score')
ggarrange(a1,a2,a3,nrow = 1,ncol = 3)
estimate1 <- estimate[which(estimate$imm != "IS1"),]
wilcox.test(estimate1$ImmuneScore~estimate1$imm)
wilcox.test(estimate1$StromalScore~estimate1$imm)
wilcox.test(estimate1$ESTIMATEScore~estimate1$imm)
estimate2 <- estimate[which(estimate$imm != "IS2"),]
wilcox.test(estimate2$ImmuneScore~estimate2$imm)
wilcox.test(estimate2$StromalScore~estimate2$imm)
wilcox.test(estimate2$ESTIMATEScore~estimate2$imm)
estimate3 <- estimate[which(estimate$imm != "IS3"),]
wilcox.test(estimate3$ImmuneScore~estimate3$imm)
wilcox.test(estimate3$StromalScore~estimate3$imm)
wilcox.test(estimate3$ESTIMATEScore~estimate3$imm)
library(ggplot2)
library(ggpubr)
library(gghalves)
setwd("D:\\北京小细胞肺癌_蛋白质组\\结果4")
xcell <- read.csv("xcell_result.csv",stringsAsFactors = F,
                  row.names = 1,check.names = F)
xcell <- xcell[,-1]
xcell <- as.data.frame(t(xcell))
cluster <- read.csv("免疫结果_consensusCluster_label.csv",
                    stringsAsFactors = F,row.names = 1,check.names = F)
xcell <- xcell[rownames(cluster),]
ada <- c(1:13,16:21)
ina <- c(14,22:34)
cluster$adaptive <- rowMeans(xcell[,ada])
cluster$innate <- rowMeans(xcell[,ina])
cluster$immune_cluster2 <- factor(cluster$immune_cluster2,
                                  levels = c("IEX","IIDS","IADS"))
a1 <- ggplot() +
  geom_half_boxplot(data = cluster, 
                    aes(x=immune_cluster2, y=innate, fill = immune_cluster2),
                    width = 0.4, position = position_dodge(0.7), 
                    alpha=1, outlier.shape = NA)+
  geom_half_point(data = cluster, 
                  aes(x=immune_cluster2, y=innate, color=immune_cluster2), 
                  transformation = position_jitter(width = 0.1,height = 0),
                  shape = 16, alpha=1, size = 0.8)+
  scale_fill_manual(values = c("#99CCFF","#F580A6","#FF5733"))+
  scale_color_manual(values = c("#99CCFF","#F580A6","#FF5733"))+
  stat_compare_means(data = cluster, 
                     aes(x=immune_cluster2, y=innate),
                     transformation = position_jitter(width = 0.1,height = 0),
                     method = "kruskal")+   theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 0, vjust = 0.6),
        legend.position="right")+
  labs(y='innate')
a2 <- ggplot() +
  geom_half_boxplot(data = cluster, 
                    aes(x=immune_cluster2, y=adaptive, fill = immune_cluster2),
                    width = 0.4, position = position_dodge(0.7), 
                    alpha=1, outlier.shape = NA)+
  geom_half_point(data = cluster, 
                  aes(x=immune_cluster2, y=adaptive, color=immune_cluster2), 
                  transformation = position_jitter(width = 0.1,height = 0),
                  shape = 16, alpha=1, size = 0.8)+
  scale_fill_manual(values = c("#99CCFF","#F580A6","#FF5733"))+
  scale_color_manual(values = c("#99CCFF","#F580A6","#FF5733"))+
  stat_compare_means(data = cluster, 
                     aes(x=immune_cluster2, y=adaptive),
                     transformation = position_jitter(width = 0.1,height = 0),
                     method = "kruskal")+   theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 0, vjust = 0.6),
        legend.position="right")+
  labs(y='adaptive')
ggarrange(a1,a2,nrow = 1,ncol = 2)
my_comparisons <- list(c("IEX", "IADS"), c("IEX", "IIDS"),c("IADS", "IIDS"))
ggplot(cluster, aes(x=immune_cluster2, y=innate,color=immune_cluster2)) +
  geom_boxplot(aes(fill=immune_cluster2),alpha=0.1)+ geom_jitter()+  
  stat_compare_means(comparisons=my_comparisons,
                     label.y = c(3000, 3200, 3400),
                     method="wilcox.test")
ggplot(cluster, aes(x=immune_cluster2, y=adaptive,color=immune_cluster2)) +
  geom_boxplot(aes(fill=immune_cluster2),alpha=0.1)+ geom_jitter()+  
  stat_compare_means(comparisons=my_comparisons,
                     label.y = c(3000, 3200, 3400),
                     method="wilcox.test")

###   4D   ###
library(reshape2)
library(ggplot2)
library(grDevices) 
library(RColorBrewer)
library(ggisoband) 
library(ggdensity)
xcell <- read.csv("xcell_result.csv",stringsAsFactors = F,
                  row.names = 1,check.names = F)
xcell <- xcell[,-1]
xcell <- as.data.frame(t(xcell))
cluster <- read.csv("consensusCluster_label.csv",
                    stringsAsFactors = F,row.names = 1,check.names = F)
xcell <- xcell[rownames(cluster),]
ada <- c(1:13,16:21)
ina <- c(14,22:34)
cluster$adaptive <- rowMeans(xcell[,ada])
cluster$innate <- rowMeans(xcell[,ina])
cluster$immune_cluster2 <- factor(cluster$immune_cluster2,
                                  levels = c("IEX","IIDS","IADS"))
colormap <- colorRampPalette(rev(brewer.pal(11,'Spectral')))(32)
ggplot(data = cluster, aes(adaptive , innate, fill = immune_cluster2)) +
  geom_hdr(xlim = c(0, 5000), ylim = c(0, 5000)) +
  theme_bw() + 
  scale_fill_manual(values = c("#99CCFF","#F580A6","#FF5733"))


###   4E   ###
library(survival)
library(survminer)
library(ggplot2)
cluster <- read.csv("consensusCluster_label.csv",
                    stringsAsFactors = F,row.names = 1)
clinical <- read.csv("134SCLC_clinical.csv",
                     stringsAsFactors = F,row.names = 1)
clinical <- clinical[rownames(cluster),]
clinical$cluster2 <- cluster$immune_cluster
cox.zph(coxph(Surv(as.numeric(OS),OS.State)~cluster2,data = clinical))
cox.zph(coxph(Surv(as.numeric(DFS),DFS.State)~cluster2,data = clinical))
# OS
summary(coxph(Surv(as.numeric(OS),OS.State)~cluster2,data = clinical))
surv <- survfit(Surv(as.numeric(OS),OS.State)~cluster2,data = clinical)
surv
survdiff(Surv(as.numeric(OS),OS.State)~cluster2,data = clinical)
summary(surv)
summary(coxph(Surv(as.numeric(OS),OS.State)~cluster2,data = clinical))
ggsurvplot(surv,pval = TRUE,risk.table = TRUE, risk.table.col = "strata", 
           palette = c("#99CCFF","#F8B195","#FF5733"),
           legend = c(2,0.5), legend.title = "",conf.int = F,
           xlab = "Time (month)",ylab = "Overall Survival")
# DFS
summary(coxph(Surv(as.numeric(DFS),DFS.State)~cluster2,data = clinical))
surv <- survfit(Surv(as.numeric(DFS),DFS.State)~cluster2,data = clinical)
surv
survdiff(Surv(as.numeric(DFS),DFS.State)~cluster2,data = clinical)
summary(surv)
summary(coxph(Surv(as.numeric(DFS),DFS.State)~cluster2,data = clinical))
ggsurvplot(surv,pval = TRUE,risk.table = TRUE,risk.table.col = "strata", 
           palette = c("#99CCFF","#F8B195","#FF5733"),
           legend = c(2,0.5), legend.title = "", conf.int = F,
           xlab = "Time (months)",ylab = "DFS")
clinical1 <- clinical[which(clinical$cluster2 != "IS2"),]
summary(coxph(Surv(as.numeric(OS),OS.State)~cluster2,data = clinical1))
summary(coxph(Surv(as.numeric(DFS),DFS.State)~cluster2,data = clinical1))
clinical2 <- clinical[which(clinical$cluster2 != "IS1"),]
summary(coxph(Surv(as.numeric(OS),OS.State)~cluster2,data = clinical2))
summary(coxph(Surv(as.numeric(DFS),DFS.State)~cluster2,data = clinical2))


###    4F  
library(pheatmap)
data <- read.csv("dreamai_134SCLC_7715protein_normalized_and_log2.csv",
                 stringsAsFactors = F,check.names = F,
                 row.names = 1)
MHCI <- c("HLA-A","HLA-B","HLA-C","TAP1","TAP2","B2M")
MHCII <- c("HLA-DPA1","HLA-DPB1","HLA-DPB2","HLA-DQA1","HLA-DQA2",
           "HLA-DQB1","HLA-DQB2","HLA-DRB1","HLA-DRB5","HLA-DRB6")
drug <- c("PDCD1","CD274","PDCD1LG2","CTLA4","LAG3","TIM3",
          "TIGIT","BTLA","IDO1","VSIR","TNFRSF18",
          "CD27","CD28","ICOS","TNFRSF4","TNFRSF9")
exhaustion <- c("TCF7","TBX21","TOX")
cyto <- c("PRF1","GZMA","GZMB","GZMH","GZMK","IFNG","NKG7")
data_1 <- data[match(c(cyto,MHCI,MHCII,drug,exhaustion),data$`Gene name`),]
data_1 <- data_1[which(is.na(data_1$`Gene name`) == F),]
rownames(data_1) <- data_1$`Gene name`
data_1 <- data_1[,-c(1,2)]
data_new <- apply(data_1, 1, scale)
data_new <- t(data_new)
colnames(data_new) <- colnames(data_1)
data_new[data_new > 2] <- 2
data_new[data_new < (-2)] <- (-2)
data_new <- data_new[,colnames(linshi)]
red <- "#FF8080";
blue <- "#8AC6D1";
white <- rgb(255,255,255,maxColorValue = 255)
annotation_col <- data.frame(imm = cluster$immune_cluster,
                             sclc = cluster$consensus_cluster,
                             TF = cluster$TF_subtype)
rownames(annotation_col) <- rownames(cluster)
ann_colors = list(imm = c("IS1" = "#99CCFF","IS2" = "#F580A6","IS3" = "#FF5733"),
                  sclc = c("Cluster1" = "#E76F51","Cluster2" = "#AEB6FF"),
                  TF = c("A" = "#EC8F5E","N" = "#F3B664","P" = "#F1EB90","Y" = "#9FBB73"))
pheatmap(data_new,gaps_col = c(63,102),
         annotation_col = annotation_col,
         cellheight = 13,cellwidth = 2,
         annotation_colors = ann_colors,
         color  = colorRampPalette(c(blue,white,red))(100),
         clustering_method = "ward.D2",
         border_color = "NA",
         cluster_cols = F, cluster_rows = F,
         show_rownames = T, show_colnames = F)
data_1 <- data_1[,colnames(linshi)]
p <- c()
for (i in 1:dim(data_1)[1]){
  p_linshi <- kruskal.test(as.numeric(data_1[i,])~cluster$immune_cluster)
  p <- c(p,p_linshi$p.value)
}
pro_result <- data.frame(p = p,
                         gene = rownames(data_1))
pro_result$p.adj <- p.adjust(pro_result$p,
                             method = "fdr")


###  4G
library(ggplot2)
library(ggpubr)
library(plyr)
library(gghalves)
setwd("D:\\北京小细胞肺癌_蛋白质组\\结果4")
cluster <- read.csv("免疫结果_consensusCluster_label.csv",
                    stringsAsFactors = F,row.names = 1,check.names = F)
cluster <- cluster[order(cluster$immune_cluster),]
setwd("D:\\北京小细胞肺癌_蛋白质组")
data <- read.csv("dreamai_134SCLC_7715protein_normalized_and_log2.csv",
                 stringsAsFactors = F,check.names = F,
                 row.names = 1)
data <- data[match(c("B2M","IDO1","VSIR"),data$`Gene name`),]
rownames(data) <- data$`Gene name`
data <- data[,rownames(cluster)]
data <- t(data)[rownames(cluster),]
cluster <- cbind(cluster,data)
ggplot() +
  geom_half_boxplot(data = cluster, 
                    aes(x=immune_cluster, y=B2M, fill = immune_cluster),
                    width = 0.4, position = position_dodge(0.7), 
                    alpha=1, outlier.shape = NA)+
  geom_half_point(data = cluster, 
                  aes(x=immune_cluster, y=B2M, color=immune_cluster), 
                  transformation = position_jitter(width = 0.1,height = 0),
                  shape = 16, alpha=1, size = 0.8)+
  scale_fill_manual(values = c("#99CCFF", "#F580A6","#FF5733"))+
  scale_color_manual(values = c("#99CCFF", "#F580A6","#FF5733"))+
  stat_compare_means(data = cluster, 
                     aes(x=immune_cluster, y=B2M),
                     method = "kruskal")+   theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 0, vjust = 0.6),
        legend.position="right")+
  labs(y='B2M')
ggplot() +
  geom_half_boxplot(data = cluster, 
                    aes(x=immune_cluster, y=IDO1, fill = immune_cluster),
                    width = 0.4, position = position_dodge(0.7), 
                    alpha=1, outlier.shape = NA)+
  geom_half_point(data = cluster, 
                  aes(x=immune_cluster, y=IDO1, color=immune_cluster), 
                  transformation = position_jitter(width = 0.1,height = 0),
                  shape = 16, alpha=1, size = 0.8)+
  scale_fill_manual(values = c("#99CCFF", "#F580A6","#FF5733"))+
  scale_color_manual(values = c("#99CCFF", "#F580A6","#FF5733"))+
  stat_compare_means(data = cluster, 
                     aes(x=immune_cluster, y=IDO1),
                     method = "kruskal")+   theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 0, vjust = 0.6),
        legend.position="right")+
  labs(y='IDO1')
ggplot() +
  geom_half_boxplot(data = cluster, 
                    aes(x=immune_cluster, y=VSIR, fill = immune_cluster),
                    width = 0.4, position = position_dodge(0.7), 
                    alpha=1, outlier.shape = NA)+
  geom_half_point(data = cluster, 
                  aes(x=immune_cluster, y=VSIR, color=immune_cluster), 
                  transformation = position_jitter(width = 0.1,height = 0),
                  shape = 16, alpha=1, size = 0.8)+
  scale_fill_manual(values = c("#99CCFF", "#F580A6","#FF5733"))+
  scale_color_manual(values = c("#99CCFF", "#F580A6","#FF5733"))+
  stat_compare_means(data = cluster, 
                     aes(x=immune_cluster, y=VSIR),
                     method = "kruskal")+   theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 0, vjust = 0.6),
        legend.position="right")+
  labs(y='VSIR')
