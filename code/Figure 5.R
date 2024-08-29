####   5A   ####
library(ConsensusClusterPlus)
chcams <- read.csv("chcamsxcell_result.csv",
                   stringsAsFactors = F,check.names = F,row.names = 1)
george <- read.csv("georgexcell_result.csv",
                   stringsAsFactors = F,check.names = F,row.names = 1)
jiang <- read.csv("jiangxcell_result.csv",
                  stringsAsFactors = F,check.names = F,row.names = 1)
cai <- read.csv("caixcell_result.csv",
                stringsAsFactors = F,check.names = F,row.names = 1)
Roper <- read.csv("Roperxcell_result.csv",
                  stringsAsFactors = F,check.names = F,row.names = 1)
merge <- cbind(chcams,george,jiang,cai,Roper)
data <- merge
data_new <- apply(data, 1, scale)
data_new <- t(data_new)
colnames(data_new) <- colnames(data)
results <- ConsensusClusterPlus(data_new[c(1:34),], maxK = 10,
                                reps = 1000, pItem = 0.9,
                                pFeature = 1,  
                                clusterAlg = "km", 
                                distance="euclidean",
                                plot = "pdf")  
label <- results[[3]][['consensusClass']]



###   5B  ###
library(pheatmap)
xcell <- read.csv("merge_xcell.csv",stringsAsFactors = F,
                  row.names = 1,check.names = F)
cluster <- read.csv("consensusCluster_label.csv",
                    stringsAsFactors = F,row.names = 1,check.names = F)
xcell <- xcell[,rownames(cluster)]
linshi <- apply(xcell,1,scale)
linshi <- t(linshi)
colnames(linshi) <- colnames(xcell)
annotation_col <- data.frame(imm = cluster$cluster1,
                             cohort = cluster$cohort)
rownames(annotation_col) <- rownames(cluster)
ann_colors = list(imm = c("IEX" = "#ADD8E6","IIDS" = "#F580A6","IADS" = "#FF5733"),
                  cohort = c("Cai" = "#74B28E","CHCAMS" = "#CEB549","George" = "#8184AE",
                             "Jiang" = "#F1DCDA","Roper" = "#7DBBC1"))
linshi[linshi > 2] <- 2
linshi[linshi < (-2)] <- (-2)
weizhi <- c(21,25,12,11,13,14,8,7,6,4,50,32,33,34,20,1,15,64,48,51)
pheatmap(linshi[weizhi,],gaps_col = c(126,169),
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         color  = colorRampPalette(c("#9C4DCC","#FFFFFF","#FFD700"))(100),
         clustering_method = "ward.D",
         border_color = "NA",
         cluster_cols = F, cluster_rows = F,
         show_rownames = T, show_colnames = F)


####       5C and 5D      ####
library(survival)
library(survminer)
cluster <- read.csv("consensusCluster_label.csv",
                    stringsAsFactors = F,row.names = 1,check.names = F)
# George 
george_clinical <- read.csv("clinical.csv",row.names = 1,
                            stringsAsFactors = F,check.names = F)
# Jiang 
GSE60052_clinical <- read.csv("clinical.csv",row.names = 1,
                              stringsAsFactors = F,check.names = F)
# Roper
Roper_clinical <- read.csv("clinical.csv",row.names = 1,
                           stringsAsFactors = F,check.names = F)
# CHCAMS
chc_clinical <- read.csv("clinical.csv",row.names = 1,
                         stringsAsFactors = F,check.names = F)
all_clinical <- data.frame(os_time = c(george_clinical$overall_survival..months.,
                                       GSE60052_clinical$os_time,
                                       Roper_clinical$OS_month,
                                       chc_clinical$OS),
                           os_state = c(george_clinical$os,
                                        GSE60052_clinical$os,
                                        Roper_clinical$os,
                                        chc_clinical$OS.State),
                           source = c(rep("george",dim(george_clinical)[1]),
                                      rep("GSE60052",dim(GSE60052_clinical)[1]),
                                      rep("Roper",dim(Roper_clinical)[1]),
                                      rep("chc",dim(chc_clinical)[1])))
rownames(all_clinical) <- c(rownames(george_clinical),
                            rownames(GSE60052_clinical),
                            rownames(Roper_clinical),
                            rownames(chc_clinical))
cluster1 <- cluster[rownames(all_clinical),]
all_clinical$cluster <- factor(cluster1$cluster1,
                               levels = c("IEX","IIDS","IADS"))
all_clinical_nonICI <- all_clinical[which(all_clinical$source != "Roper"),]
surv <- survfit(Surv(os_time,os_state)~cluster,data = all_clinical_nonICI)
surv
survdiff(Surv(os_time,os_state)~cluster,data = all_clinical_nonICI)
summary(surv)
summary(coxph(Surv(os_time,os_state)~cluster,data = all_clinical_nonICI))
ggsurvplot(surv,pval = TRUE,risk.table = TRUE, risk.table.col = "strata", 
           palette = c("#99CCFF","#F580A6","#FF5733"),
           legend = c(2,0.5),legend.title = "", 
           conf.int = F,xlab = "Time (month)",ylab = "Overall survival")
cli <- all_clinical_nonICI[which(all_clinical_nonICI$cluster != "IIDS"),]
surv1 <- survfit(Surv(os_time,os_state)~cluster,data = cli)
ggsurvplot(surv1,pval = TRUE,risk.table = TRUE)
# icb
clinical_Roper <- all_clinical[which(all_clinical$source == "Roper"),]
surv_Roper <- survfit(Surv(os_time,os_state)~cluster,data = clinical_Roper)
summary(surv_Roper)
ggsurvplot(surv_Roper,pval = TRUE,risk.table = TRUE, risk.table.col = "strata", 
           palette = c("#99CCFF","#F580A6","#FF5733"),
           legend = c(2,0.5), legend.title = "", conf.int = F,
           xlab = "Time (month)",ylab = "Overall survival")
cli1 <- clinical_Roper[which(clinical_Roper$cluster != "IADS"),]
surv1 <- survfit(Surv(os_time,os_state)~cluster,data = cli1)
ggsurvplot(surv1,pval = TRUE,risk.table = TRUE)



###  5D  ###
library(ggplot2)
library(ggpubr)
library(plyr)
Roper_clinical <- read.csv("clinical.csv",row.names = 1,
                           stringsAsFactors = F,check.names = F)
table(Roper_clinical$response,Roper_clinical$immune)
data_duiji <- data.frame(type = c(rep("CB_IADS",0),rep("CB_IEX",1),rep("CB_IIDS",3),
                                  rep("NCB_IADS",1),rep("NCB_IEX",10),rep("NCB_IIDS",1)),
                         cluster = c(rep("IEX",1),rep("IIDS",3),rep("IADS",0),
                                     rep("IEX",10),rep("IIDS",1),rep("IADS",1)))

data_duiji$number <- 1
data_duiji1 <- ddply(data_duiji,'type',transform,percent = 1/sum(number)*100)
ggplot(data_duiji1)+
  scale_fill_manual(values = c("#FF713D","#FF713D",
                               "#CED4D9","#CED4D9","#CED4D9"))+
  geom_bar(aes(x=cluster,fill=type),position = "fill")+
  theme_bw()+
  theme(axis.ticks.length=unit(0.5,'cm'))+
  guides(fill=guide_legend(title=NULL)) +
  theme(axis.text.x=element_text(hjust = 1,colour="black",size=10,angle = 45), 
        axis.text.y=element_text(size=10,colour="black"), 
        axis.title.y=element_text(size = 10,colour="black"), 
        legend.text=element_text(colour="black", 
                                 size=10),
        legend.title=element_text(colour="black", 
                                  size=10),
        panel.grid.major = element_blank(),   
        panel.grid.minor = element_blank())
fisher.test(table(Roper_clinical$response,Roper_clinical$immune))



###   Validation in Roper protein cohort  all protein heatmap ###
library(pheatmap)
Roper_pro1 <- read.csv("data_dsp.csv",
                      stringsAsFactors = F,check.names = F,row.names = 1)
response <- Roper_pro1$response
Roper_pro <- as.data.frame(t(Roper_pro1[,2:31]))
cluster <- read.csv("consensusCluster_label.csv",
                    stringsAsFactors = F,row.names = 1,check.names = F)
Roper_cluster <- cluster[colnames(Roper_pro),]
linshi <- apply(Roper_pro,1,scale)
linshi <- t(linshi)
colnames(linshi) <- colnames(Roper_pro)
annotation_col <- data.frame(imm = Roper_cluster$cluster1)
rownames(annotation_col) <- rownames(Roper_cluster)
ann_colors = list(imm = c("IEX" = "#99CCFF","IIDS" = "#F580A6",
                          "IADS" = "#FF5733"))
linshi[linshi > 1.5] <- 1.5
linshi[linshi < (-1.5)] <- (-1.5)
pheatmap(linshi,
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         color  = colorRampPalette(c("#4682B4","#FFFFFF","#FF8080"))(100),
         clustering_method = "ward.D2",
         border_color = "NA",
         cluster_cols = T, cluster_rows = T,
         show_rownames = T, show_colnames = T)
Roper_pro1 <- Roper_pro1[rownames(Roper_cluster),]
Roper_cluster <- cbind(Roper_cluster,Roper_pro1)
FC <- c()
p <- c()
for (i in 5:34){
  a <- mean(Roper_cluster[5:6,i])/mean(Roper_cluster[1:2,i])
  b <- wilcox.test(Roper_cluster[,i]~Roper_cluster$cluster1)
  FC <- c(FC,a)
  p <- c(p,b$p.value)
}
library(ggplot2)
result <- data.frame(protein = colnames(Roper_cluster)[5:34],
                     type = rep("Roper cohort",length(FC)),
                     FC = FC,
                     p = p,
                     fdr = p.adjust(p,method = "fdr"))
result1 <- rbind(result,result)
result1$type <- c(rep("Roper cohort",length(FC)),
                  rep("No cohort",length(FC)))

ggplot(data = result1, aes(x = type, y = protein)) + 
  geom_point(aes(size = FC,color = FC), alpha = 0.7)+
  scale_color_gradient2(low = "#003366",
                        mid = "#F5F5F5",
                        high = "#C70039",
                        midpoint = 1)


###   Validation of Roper TCR cohort
library(ggplot2)
library(ggpubr)
pre <- data.frame(patient = c("CL0106","CL0107","CL0108","CL0110",
                              "CL0111","CL0116","CL0124","NCI0422"),
                  imm_lab = c("IEX","IEX","IEX","IEX",
                              "IIDS","IEX","IEX","IIDS"),
                  response = c("NCB","NCB","NCB","NCB",
                              "CB","NCB","NCB","CB"),
                  TCR_CDR3 = c(14394,5087,11431,16189,
                               22795,6462,786,23883),
                  TCR_VJ = c(2268,1685,2151,2283,
                             2514,1857,536,2525))
pair <- data.frame(patient = c("CL0116-pre","CL0116-post",
                               "CL0124-pre","CL0124-post",
                               "NCI0422-pre","NCI0422-relapse"),
                   imm_lab = c("IEX","IEX","IEX","IEX",
                               "IIDS","IIDS"),
                   response = c("NCB","NCB","NCB","NCB",
                                "CB","CB"),
                   TCR_CDR3 = c(6462,14020,786,2214,
                                23883,5127),
                   TCR_VJ = c(1857,2329,536,1079,
                              2525,1595))
##  Pre
pre$patient <- factor(pre$patient,
                      levels = c("CL0106","CL0107","CL0108","CL0110",
                                 "CL0116","CL0124","CL0111","NCI0422"))
ggplot(pre, aes(x=patient, y=TCR_CDR3, fill=imm_lab )) + 
  geom_bar(stat="identity",  position=position_dodge()) + 
  scale_fill_manual(values = c("#99CCFF","#F580A6"))+
  theme_bw() + ylab("TCR_CDR3")
ggplot(pre, aes(x=patient, y=TCR_VJ, fill=imm_lab )) + 
  geom_bar(stat="identity",  position=position_dodge()) + 
  scale_fill_manual(values = c("#99CCFF","#F580A6"))+
  theme_bw() + ylab("TCR_VJ")
###  Paired  
pair$patient <- factor(pair$patient,
                       levels = pair$patient)
ggplot(pair, aes(x=patient, y=TCR_CDR3, fill=imm_lab )) + 
  geom_bar(stat="identity",  position=position_dodge()) + 
  scale_fill_manual(values = c("#99CCFF","#F580A6"))+
  theme_bw() + ylab("TCR_CDR3")
ggplot(pair, aes(x=patient, y=TCR_VJ, fill=imm_lab )) + 
  geom_bar(stat="identity",  position=position_dodge()) + 
  scale_fill_manual(values = c("#99CCFF","#F580A6"))+
  theme_bw() + ylab("TCR_VJ")





























