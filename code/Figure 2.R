###  de novo clustering  ###
library(dplyr)
library(ggplot2)
data <- read.csv("dreamai_134SCLC_7715protein_normalized_and_log2.csv",
                 stringsAsFactors = F,check.names = F,
                 row.names = 1)
data <- data[,3:136]
# calculating CV top 10%
cal_cv=function(x){  # 自定义函数 标准差/平均值
  y=na.omit(x)
  return(sd(y)/mean(y))
}
cv_result <- apply(data, 1, cal_cv) # 在每一行上应用自定义cal_cv函数
cv_result <- cv_result[cv_result %>% order %>% rev]
cv_result <- data.frame(x = cv_result)
ggplot(cv_result, aes(x = x)) + 
  geom_density() +
  theme_bw() +
  geom_vline(xintercept = cv_result[772,])+
  labs(title = "Density Plot", x = "Value", y = "Density")
pro_top10 <- rownames(cv_result)[1:772]  ## 前 10% 高变蛋白
###  NMF  ###
library(NMF)
data <- read.csv("dreamai_134SCLC_7715protein_normalized_and_log2.csv",
                 stringsAsFactors = F,check.names = F,
                 row.names = 1)
data1 <- data[pro_top10,3:136]
ranks <- 2:10
estim.coad <- nmf(data1,ranks, nrun=30)
plot(estim.coad)
# NMF,rank=2
nmf.rank4 <- nmf(data1,rank = 2, 
                 nrun=30,seed = 3, 
                 method = "brunet")
group <- predict(nmf.rank4) 
table(group)
### Consensus  ###
library(ConsensusClusterPlus)
data_new <- apply(data1, 1, scale)
data_new <- t(data_new)
colnames(data_new) <- colnames(data1)
results <- ConsensusClusterPlus(data_new, maxK = 10,
                                reps = 1000, pItem = 0.8,
                                pFeature = 1,  
                                clusterAlg = "km", 
                                distance="euclidean",
                                plot = "pdf")  
result_data <- results[[2]]$ml
rownames(result_data) <- names(results[[2]]$consensusClass)
## best k using PAC method
rivzi2016_PAC=c()  
for(i in 2:10){
  cM=results[[i]]$consensusMatrix
  Fn=ecdf(cM[lower.tri(cM)])  
  rivzi2016_PAC=c(rivzi2016_PAC, Fn(0.9) - Fn(0.1))  
}#end for i# The optimal K
rivzi2016_optK=c(2:10)[which.min(rivzi2016_PAC)]
label <- results[[2]][['consensusClass']]
PAC_data <- data.frame(x = 2:10,
                       pac = rivzi2016_PAC)
ggplot(PAC_data, aes(x = x, y = pac,color = "#4F6F52")) + 
  geom_point() + geom_line() + theme_bw() + 
  scale_color_identity() + labs(title = "Consensus Cluster")
### Hierarchical  ###
library(pheatmap)
data_new <- apply(data1, 1, scale)
data_new <- t(data_new)
colnames(data_new) <- colnames(data1)
red <- "#D94E48";
blue <- "#5175A4";
white <- rgb(255,255,255,maxColorValue = 255)
heat <- pheatmap(data_new,
                 color  = colorRampPalette(c(blue,white,red))(100),
                 clustering_method = "ward.D",
                 border_color = "grey60",
                 cluster_cols = T, cluster_rows = T,
                 show_rownames = F, show_colnames = T)
out.dist=dist(t(data_new)) 
out.hclust=hclust(out.dist,method="ward.D")
out.id=cutree(out.hclust,k=2)  


####  2A  ####
library(pheatmap)
cluster <- read.csv("top10_newcluster.csv",stringsAsFactors = F,
                    row.names = 1)
clinical <- read.csv("134SCLC_clinical.csv",
                     stringsAsFactors = F,row.names = 1,check.names = F)
clinical <- clinical[rownames(cluster),]
clinical$cluster1 <- cluster$consensus_cluster
clinical$cluster2 <- cluster$NMF_cluster
clinical$cluster3 <- cluster$hcluster
clinical$age_group <- clinical$Age
clinical$age_group[clinical$Age > 60] <- ">60"
clinical$age_group[clinical$Age <= 60] <- "<=60"
data <- cbind(clinical$`ASCL1 H-score`,clinical$`NeuroD1 H-score`,
              clinical$`POU2F3 H-score`,clinical$`YAP1 H-score`)
rownames(data) <- rownames(clinical)
colnames(data) <- c('ASCL1','NEUROD1','POU2F3','YAP1')
data <- t(data)
scale_data <- apply(data, 2, scale)
rownames(scale_data) <- rownames(data)
scale_data[scale_data > 2] <- 2
scale_data[scale_data < (-2)] <- (-2)
clinical <- clinical[colnames(scale_data),]

data_ki67 <- cbind(clinical$Ki67,clinical$Ki67)
rownames(data_ki67) <- rownames(clinical)
colnames(data_ki67) <- c('Ki67','Ki67')
scale_data_ki67 <- apply(data_ki67, 2, scale)
rownames(scale_data_ki67) <- rownames(data_ki67)
scale_data_ki67[scale_data_ki67 > 1.5] <- 1.5
scale_data_ki67[scale_data_ki67 < (-1.5)] <- (-1.5)
annotation_col <- data.frame(CD56 = clinical$CD56,
                             CgA = clinical$CgA,
                             Syn = clinical$Syn,
                             AJCCstage = clinical$`AJCC stage`,
                             Mcategory = clinical$`M stage`,
                             Ncategory = clinical$`N stage`,
                             Tcategory = clinical$`T stage`,
                             Tumorthrombosis = clinical$`Tumor  thrombosis`,
                             STAS = clinical$`spread through air spaces_STAS`,
                             Nerveinvasion = clinical$`Nerve invasion`,
                             Vascularinvasion = clinical$`Vascular invasion`,
                             Bronchusinvasion = clinical$`Bronchus invasion`,
                             Lymphaticmetastasis = clinical$`Lymphatic metastasis`,
                             Tumorlocation = clinical$`Tumor location`,
                             Smoking = clinical$Smoking,
                             Sex = clinical$Sex,
                             Age = clinical$age_group,
                             hclus = clinical$cluster3,
                             nmf = clinical$cluster2,
                             consensus = clinical$cluster1)
rownames(annotation_col) <- rownames(clinical)
ann_colors = list(CD56 = c("0" = "#F2F3F4","1" = "#DADADA","2" = "#A7A9AC",
                           "3" = "#3A3A3A","Unknown" = "#FFFFFF"),
                  CgA = c("0" = "#F2F3F4","1" = "#DADADA","2" = "#A7A9AC",
                          "3" = "#3A3A3A","Unknown" = "#FFFFFF"),
                  Syn = c("0" = "#F2F3F4","1" = "#DADADA","2" = "#A7A9AC",
                          "3" = "#3A3A3A","Unknown" = "#FFFFFF"),
                  AJCCstage = c("I" = "#F3CCFF","II" = "#D09CFA",
                                "III" = "#A555EC","Unknown" = "#F1F0E8"),
                  Mcategory = c("M0" = "#CDBBA7","Unknown" = "#F1F0E8"),
                  Ncategory = c("N0" = "#D3E4CD","N1" = "#ADC2A9",
                                "N2" = "#99A799","Unknown" = "#F1F0E8"),
                  Tcategory = c("T1" = "#F1F0C0","T2" = "#B7E5DD",
                                "T3" = "#B1BCE6","T4" = "#9A86A4"),
                  Tumorthrombosis = c("Yes" = "#FFCCCC","No" = "#D3E2F2"),
                  STAS = c("Yes" = "#FFCCCC","No" = "#D3E2F2"),
                  Nerveinvasion = c("Yes" = "#FFCCCC","No" = "#D3E2F2"),
                  Vascularinvasion = c("Yes" = "#FFCCCC","No" = "#D3E2F2"),
                  Bronchusinvasion = c("Yes" = "#FFCCCC","No" = "#D3E2F2","Unknown" = "#F1F0E8"),
                  Lymphaticmetastasis = c("Yes" = "#FFCCCC","No" = "#D3E2F2","Unknown" = "#F1F0E8"),
                  Tumorlocation = c("Left" = "#D7C0AE","Right" = "#9BABB8"),
                  Smoking = c("Yes" = "#FFCCCC","No" = "#D3E2F2"),
                  Sex = c("Female" = "#7FB3D5","Male" = "#DCC6E0"),
                  Age = c("<=60" = "#B4BDFF",">60" = "#FFD28F"),
                  nmf = c("Cluster1" = "#F7E0A3","Cluster2" = "#4C8492"),
                  hclus = c("Cluster1" = "#FFB677","Cluster2" = "#98C6E1"),
                  consensus = c("subtype1" = "#E76F51","subtype2" = "#7A6DB0"))
red <- "#FF8080";
blue <- "#8AC6D1";
white <- rgb(255,255,255,maxColorValue = 255)
clusterfun <- function(x){
  which(x == max(x))
}
scale_data1_c1 <- rownames(clinical[which(clinical$cluster1 == "subtype1"),])
scale_data1_c2 <- rownames(clinical[which(clinical$cluster1 == "subtype2"),])
scale_data1 <- scale_data[,scale_data1_c1]
scale_data2 <- scale_data[,scale_data1_c2]
rank1 <- apply(scale_data1, 2, clusterfun)
rank1 <- rank1[order(rank1)]
scale_data1 <- scale_data1[,names(rank1)]
rank2 <- apply(scale_data2, 2, clusterfun)
rank2 <- rank2[order(rank2)]
scale_data2 <- scale_data2[,names(rank2)]
scale_data_new <- cbind(scale_data1,scale_data2)
pheatmap(scale_data_new,fontsize=6,cellheight = 10,
         color  = colorRampPalette(c(blue,white,red))(100),
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         clustering_method = "ward.D2",
         border_color = "grey60",
         cluster_cols = F, cluster_rows = F,
         show_rownames = T, show_colnames = T)
scale_data_ki67 <- t(scale_data_ki67)
scale_data_ki67 <- scale_data_ki67[,colnames(scale_data_new)]
pheatmap(scale_data_ki67,fontsize=6,cellheight = 10,
         color  = colorRampPalette(c("#DADADA",white,"#F25F5C"))(100),
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         clustering_method = "ward.D2",
         border_color = "grey60",
         cluster_cols = F, cluster_rows = F,
         show_rownames = T, show_colnames = T)

# 全蛋白
data <- read.csv("dreamai_134SCLC_7715protein_normalized_and_log2.csv",
                 stringsAsFactors = F,check.names = F,
                 row.names = 1)
data.new <- data[,c(-1,-2)]
data.new <- data.new[,colnames(scale_data_new)]
data.linshi <- apply(data.new,1,scale)
data.linshi <- t(data.linshi)
colnames(data.linshi) <- colnames(data.new)
data.linshi[data.linshi > 2] <- 2
data.linshi[data.linshi < (-2)] <- (-2)
gene <- read.csv("D:\\北京小细胞肺癌_蛋白质组\\结果2\\7715gene变异系数.csv",
                 stringsAsFactors = F,row.names = 1,check.names = F)
data.linshi_high <- data.linshi[rownames(gene)[1:772],]
data.linshi_other <- data.linshi[rownames(gene)[773:7715],]

pheatmap(data.linshi_high,fontsize=6,
         color  = colorRampPalette(c("#9C4DCC",white,"#FFD700"))(100),
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         clustering_method = "ward.D2",
         border_color = "grey60",
         cluster_cols = F, cluster_rows = T,
         show_rownames = F, show_colnames = T)
pheatmap(data.linshi_other,fontsize=6,
         color  = colorRampPalette(c("#9C4DCC",white,"#FFD700"))(100),
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         clustering_method = "ward.D2",
         border_color = "grey60",
         cluster_cols = F, cluster_rows = F,
         show_rownames = F, show_colnames = T)


###  SCLC PS Sankey plot  2B
library(tidyverse)
library(viridis)
library(patchwork)
library(networkD3)
MisLinks <- read.csv("TF-link_matrix.csv",stringsAsFactors = F)
MisNodes <- read.csv("TF-node_matrix.csv",stringsAsFactors = F)
Node2index = list()
Node2index[MisNodes$name] = 0:length(MisNodes$name)
MisLinks = MisLinks %>%
  mutate(source2 = unlist(Node2index[source])) %>%
  mutate(target2 = unlist(Node2index[target]))
color2project = paste(unique(MisNodes$group_color),collapse = '","')
my_color <- paste0('d3.scaleOrdinal().domain(["',color2project,'"]).range(["',color2project,'"])')
sankeyNetwork(Links = MisLinks, 
              Nodes = MisNodes,
              Source = "source2", 
              Target = "target2",
              Value ="value",
              NodeID = "name",
              NodeGroup = "group_color", 
              colourScale = JS(my_color),
              fontSize = 20)


####  2C
library(ggplot2)
library(ggpubr)
library(plyr)
library(circlize)
library(ggforce)
library(highcharter)
setwd("D:\\北京小细胞肺癌_蛋白质组")
cluster <- read.csv("top10_newcluster.csv",stringsAsFactors = F,
                    row.names = 1)
clinical <- read.csv("D:\\北京小细胞肺癌_蛋白质组\\134SCLC_clinical.csv",
                     stringsAsFactors = F,row.names = 1)
clinical <- clinical[rownames(cluster),]
clinical$age_group <- clinical$Age
clinical$age_group[clinical$Age > 60] <- ">60"
clinical$age_group[clinical$Age <= 60] <- "<=60"
clinical$cluster <- cluster$consensus_cluster
clinical_C1 <- clinical[which(clinical$cluster == "Cluster1"),]
clinical_C2 <- clinical[which(clinical$cluster == "Cluster2"),]
data.type1 <- data.frame(ratio = as.numeric(table(clinical_C1$TF_subtype)/dim(clinical_C1)[1]),
                         type = c("A","N","P","Y"))
data.type1$type <- factor(data.type1$type,
                          levels = c("A","N","P","Y"))
data.type2 <- data.frame(ratio = as.numeric(table(clinical_C2$TF_subtype)/dim(clinical_C2)[1]),
                         type = c("A","N","P","Y"))
data.type2$type <- factor(data.type2$type,
                          levels = c("A","N","P","Y"))
## Syn
table(clinical_C1$Syn)
table(clinical_C2$Syn)
data_duiji <- data.frame(type = c(rep("C1_Syn0",1),rep("C1_Syn1",36),
                                  rep("C1_Syn2",10),rep("C1_Syn3",3),rep("C1_NA",2),
                                  rep("C2_Syn0",6),rep("C2_Syn1",40),
                                  rep("C2_Syn2",12),rep("C2_Syn3",15),rep("C2_NA",9)),
                         cluster = c(rep("C1",1),rep("C1",36),
                                     rep("C1",10),rep("C1",3),rep("C1",2),
                                     rep("C2",6),rep("C2",40),
                                     rep("C2",12),rep("C2",15),rep("C2",9)))

data_duiji$number <- 1
data_duiji1 <- ddply(data_duiji,'type',transform,percent = 1/sum(number)*100)
ggplot(data_duiji1)+
  scale_fill_manual(values = c("#FFFFFF","#F2F3F4","#DADADA","#A7A9AC","#3A3A3A",
                               "#FFFFFF","#F2F3F4","#DADADA","#A7A9AC","#3A3A3A"))+ 
  geom_bar(aes(x=cluster,fill=type),position = "fill")+
  theme_bw()+
  theme(axis.ticks.length=unit(0.5,'cm'))+
  guides(fill=guide_legend(title=NULL)) +
  theme(axis.text.x=element_text(hjust = 1,colour="black",size=10,angle = 45), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
        axis.text.y=element_text(size=10,colour="black"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
        axis.title.y=element_text(size = 10,colour="black"), #设置y轴标题的字体属性
        legend.text=element_text(colour="black",  #设置图例的子标题的字体属性
                                 size=10),
        legend.title=element_text(colour="black", #设置图例的总标题的字体属性
                                  size=10),
        panel.grid.major = element_blank(),   #不显示网格线
        panel.grid.minor = element_blank())
## CgA
table(clinical_C1$CgA)
table(clinical_C2$CgA)
data_duiji <- data.frame(type = c(rep("C1_CgA0",1),rep("C1_CgA1",38),
                                  rep("C1_CgA2",7),rep("C1_CgA3",4),rep("C1_NA",2),
                                  rep("C2_CgA0",14),rep("C2_CgA1",39),
                                  rep("C2_CgA2",13),rep("C2_CgA3",8),rep("C2_NA",8)),
                         cluster = c(rep("C1",1),rep("C1",38),
                                     rep("C1",7),rep("C1",4),rep("C1",2),
                                     rep("C2",14),rep("C2",39),
                                     rep("C2",13),rep("C2",8),rep("C2",8)))

data_duiji$number <- 1
data_duiji1 <- ddply(data_duiji,'type',transform,percent = 1/sum(number)*100)
ggplot(data_duiji1)+
  scale_fill_manual(values = c("#F2F3F4","#DADADA","#A7A9AC","#3A3A3A","#FFFFFF",
                               "#F2F3F4","#DADADA","#A7A9AC","#3A3A3A","#FFFFFF"))+ 
  geom_bar(aes(x=cluster,fill=type),position = "fill")+
  theme_bw()+
  theme(axis.ticks.length=unit(0.5,'cm'))+
  guides(fill=guide_legend(title=NULL)) +
  theme(axis.text.x=element_text(hjust = 1,colour="black",size=10,angle = 45), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
        axis.text.y=element_text(size=10,colour="black"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
        axis.title.y=element_text(size = 10,colour="black"), #设置y轴标题的字体属性
        legend.text=element_text(colour="black",  #设置图例的子标题的字体属性
                                 size=10),
        legend.title=element_text(colour="black", #设置图例的总标题的字体属性
                                  size=10),
        panel.grid.major = element_blank(),   #不显示网格线
        panel.grid.minor = element_blank())
## CD56
table(clinical_C1$CD56)
table(clinical_C2$CD56)
data_duiji <- data.frame(type = c(rep("C1_CD560",1),rep("C1_CD561",31),
                                  rep("C1_CD562",2),rep("C1_CD563",14),rep("C1_NA",4),
                                  rep("C2_CD560",5),rep("C2_CD561",35),
                                  rep("C2_CD562",8),rep("C2_CD563",27),rep("C2_NA",7)),
                         cluster = c(rep("C1",1),rep("C1",31),
                                     rep("C1",2),rep("C1",14),rep("C1",4),
                                     rep("C2",5),rep("C2",35),
                                     rep("C2",8),rep("C2",27),rep("C2",7)))

data_duiji$number <- 1
data_duiji1 <- ddply(data_duiji,'type',transform,percent = 1/sum(number)*100)
ggplot(data_duiji1)+
  scale_fill_manual(values = c("#F2F3F4","#DADADA","#A7A9AC","#3A3A3A","#FFFFFF",
                               "#F2F3F4","#DADADA","#A7A9AC","#3A3A3A","#FFFFFF"))+ 
  geom_bar(aes(x=cluster,fill=type),position = "fill")+
  theme_bw()+
  theme(axis.ticks.length=unit(0.5,'cm'))+
  guides(fill=guide_legend(title=NULL)) +
  theme(axis.text.x=element_text(hjust = 1,colour="black",size=10,angle = 45), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
        axis.text.y=element_text(size=10,colour="black"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
        axis.title.y=element_text(size = 10,colour="black"), #设置y轴标题的字体属性
        legend.text=element_text(colour="black",  #设置图例的子标题的字体属性
                                 size=10),
        legend.title=element_text(colour="black", #设置图例的总标题的字体属性
                                  size=10),
        panel.grid.major = element_blank(),   #不显示网格线
        panel.grid.minor = element_blank())


###   2C  Ki67 boxplot  ###
library(ggplot2)
library(ggpubr)
library(plyr)
library(gghalves)
cluster <- read.csv("top10_newcluster.csv",stringsAsFactors = F,
                    row.names = 1)
clinical <- read.csv("134SCLC_clinical.csv",
                     stringsAsFactors = F,row.names = 1)
clinical <- clinical[rownames(cluster),]
clinical$age_group <- clinical$Age
clinical$age_group[clinical$Age > 60] <- ">60"
clinical$age_group[clinical$Age <= 60] <- "<=60"
clinical$cluster <- cluster$consensus_cluster
ggplot() +
  geom_half_boxplot(data = clinical, 
                    aes(x=cluster, y=Ki67, fill = cluster),
                    width = 0.4, position = position_dodge(0.7), 
                    alpha=1, outlier.shape = NA)+
  geom_half_point(data = clinical, 
                  aes(x=cluster, y=Ki67, color=cluster), 
                  transformation = position_jitter(width = 0.1,height = 0),
                  shape = 16, alpha=1, size = 0.8)+
  scale_fill_manual(values = c("#E76F51","#AEB6FF"))+
  scale_color_manual(values = c("#E76F51","#AEB6FF"))+
  stat_compare_means(data = clinical, 
                     aes(x=cluster, y=Ki67),
                     transformation = position_jitter(width = 0.1,height = 0),
                     method = "wilcox.test")+   theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 0, vjust = 0.6),
        legend.position="right")+
  labs(y='Ki67 H-score')



###  2D
library(survival)
library(survminer)
library(ggplot2)
cluster <- read.csv("top10_newcluster.csv",stringsAsFactors = F,
                    row.names = 1)
clinical <- read.csv("134SCLC_clinical.csv",
                     stringsAsFactors = F,row.names = 1)
clinical$age_group <- clinical$Age
clinical$age_group[clinical$Age > 60] <- ">60"
clinical$age_group[clinical$Age <= 60] <- "<=60"
clinical <- clinical[rownames(cluster),]
clinical$cluster1 <- cluster$NMF_cluster
clinical$cluster2 <- cluster$consensus_cluster
clinical$cluster3 <- cluster$hcluster
##  OS KM curv
clinical_os <- clinical[which(clinical$OS <= 60),]
cox.zph(coxph(Surv(as.numeric(OS),OS.State)~cluster2,data = clinical_os))
clinical_os$AJCC.stage[which(clinical_os$AJCC.stage == "Unknown")] <- NA
clinical_os$Lymphatic.metastasis[which(clinical_os$Lymphatic.metastasis == "Unknown")] <- NA
summary(coxph(Surv(as.numeric(OS),OS.State)~cluster2,data = clinical_os))
summary(coxph(Surv(as.numeric(OS),OS.State)~cluster2+Sex+Smoking+
                Lymphatic.metastasis+age_group+AJCC.stage,data = clinical_os))
surv <- survfit(Surv(as.numeric(OS),OS.State)~cluster2,data = clinical_os)
surv
survdiff(Surv(as.numeric(OS),OS.State)~cluster2,data = clinical_os)
summary(surv)
summary(coxph(Surv(as.numeric(OS),OS.State)~cluster2,data = clinical_os))
ggsurvplot(surv,pval = TRUE,risk.table = TRUE, 
           risk.table.col = "strata", palette = c("#FFB677","#98C6E1"),
           legend = c(2,0.5), legend.title = "", conf.int = F,
           xlab = "Time (month)",ylab = "Overall Survival")
# DFS
clinical_dfs <- clinical[which(clinical$DFS <= 57),]
cox.zph(coxph(Surv(as.numeric(DFS),DFS.State)~cluster2,data = clinical_dfs))
clinical_dfs$AJCC.stage[which(clinical_dfs$AJCC.stage == "Unknown")] <- NA
clinical_dfs$Lymphatic.metastasis[which(clinical_dfs$Lymphatic.metastasis == "Unknown")] <- NA
summary(coxph(Surv(as.numeric(DFS),DFS.State)~cluster2,data = clinical_dfs))
summary(coxph(Surv(as.numeric(DFS),DFS.State)~cluster2+Sex+Smoking+
                Lymphatic.metastasis+age_group+AJCC.stage,data = clinical_dfs))
surv <- survfit(Surv(as.numeric(DFS),DFS.State)~cluster2,data = clinical_dfs)
surv
survdiff(Surv(as.numeric(DFS),DFS.State)~cluster2,data = clinical_dfs)
summary(surv)
summary(coxph(Surv(as.numeric(DFS),DFS.State)~cluster2,data = clinical_dfs))
ggsurvplot(surv,pval = TRUE,risk.table = TRUE, risk.table.col = "strata", 
           palette = c("#FFB677","#98C6E1"),legend = c(2,0.5), legend.title = "",
           conf.int = F, xlab = "Time (month)",ylab = "DFS")


###   2E    ###
library(ggplot2)
library(ggpubr)
my_colors <- c("Yes" = "#FFCCCC", "No" = "#D3E2F2")
all_recurrence <- data.frame(percentage = c(26/38,12/38,33/72,39/72),
                                     subtype = factor(c("S1", "S1", "S2", "S2"), 
                                                      levels = c("S1", "S2")),
                                     re =  factor(c("Yes", "No", "Yes", "No"), 
                                                  levels = c("Yes", "No")))
ggplot(all_recurrence, aes(x = subtype, y = percentage*100, fill = re)) +
  geom_bar(stat = "identity") + theme_bw() +
  labs(y = "Percentage", x = "Subtype", fill = "Recurrence",
       title = "All recurrence") +
  scale_fill_manual(values = my_colors) 
