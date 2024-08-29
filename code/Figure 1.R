###  median centering method (intensity/median)  and  log2  ###
data <- read.csv("134SCLC_10065protein.csv",
                  stringsAsFactors = F,
                  check.names = F)
data1 <- data[,4:137]
# Subtract the medians from each protein
median_normalized_row <- function (x) {
  narrays <- NCOL(x)
  if (narrays == 1) 
    return(x)
  cmed <- log(apply(x, 1, median, na.rm = TRUE)) # row
  cmed <- exp(cmed - mean(cmed))
  x/cmed
}
centered_data <- median_normalized_row(data1)
centered_data_log2 <- log2(centered_data)
rownames(centered_data_log2) <- rownames(data1)
centered_data_log2 <- cbind(data2[,1:3],centered_data_log2)


###  >=25% NA remove  ###
data2 <- read.csv("134SCLC_10065protein_normalized_and_log2.csv",
                  stringsAsFactors = F,
                  check.names = F)
weizhi_dele2 <- c()
for (i in 1:dim(data2)[1]){
  if (length(which(is.na(data2[i,4:137]) == T)) >= 34){
    weizhi_dele2 <- c(weizhi_dele2,i)
  }
}
data2.new <- data2[-weizhi_dele2,]


###  1C  ###
library(ggplot2)
library(ggpubr)
mat <- data.frame(value = c(56831,10065,7715),
                  name = factor(c('Identified peptides',
                                'Identified proteins',
                                'Filtered proteins'),
                                levels = c('Identified peptides',
                                           'Identified proteins',
                                           'Filtered proteins')))
ggplot(mat, aes(x=name, y=value, fill=name)) + 
  geom_bar(stat="identity", color="black", show.legend = FALSE) +
  scale_fill_manual(values=c('#ECF2FF','#DAF5FF','#B0DAFF')) +
  theme_bw() + scale_y_log10()

library(ggplot2)
library(ggpubr)
data <- read.csv("134SCLC_7715protein_normalized_and_log2.csv",
                 stringsAsFactors = F,check.names = F,
                 row.names = 1)
pro_num_each_pt <- matrix(,,2)
colnames(pro_num_each_pt) <- c('FALSE','TRUE')
for (i in 3:136){
  tab <- table(is.na(data[,i]))
  pro_num_each_pt <- rbind(pro_num_each_pt,tab)
}
pro_num_each_pt <- pro_num_each_pt[-1,] 
rownames(pro_num_each_pt) <- colnames(data)[3:136]
pro_num_each_pt <- as.data.frame(pro_num_each_pt)
pro_num_each_pt$patient <- rownames(pro_num_each_pt)
pro_num_each_pt <- pro_num_each_pt[order(pro_num_each_pt$`FALSE`),]
pro_num_each_pt$patient <- factor(pro_num_each_pt$patient,
                                  levels = pro_num_each_pt$patient)
ggplot(pro_num_each_pt, aes(x=patient, y=`FALSE`,fill=`FALSE`)) +
  geom_bar(stat="identity") +
  theme_bw() +
  geom_hline(yintercept = median(pro_num_each_pt$`FALSE`), 
             linetype="dashed", color="black") +
  scale_fill_gradient(low="#FFF78A", high="#FFC47E") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x="Patient", y="Count")


###  protein count between clinical information  ###
library(ggplot2)
library(ggpubr)
library(gghalves)
data <- read.csv("134SCLC_10065protein.csv",
                 stringsAsFactors = F,check.names = F,
                 row.names = 1)
pro_num_each_pt <- matrix(,,2)
colnames(pro_num_each_pt) <- c('FALSE','TRUE')
for (i in 3:136){
  tab <- table(is.na(data[,i]))
  pro_num_each_pt <- rbind(pro_num_each_pt,tab)
}
pro_num_each_pt <- pro_num_each_pt[-1,] 
rownames(pro_num_each_pt) <- colnames(data)[3:136]
pro_num_each_pt <- as.data.frame(pro_num_each_pt)
clinical1 <- read.csv("134SCLC_clinical.csv",stringsAsFactors = F,
                      row.names = 1,check.names = F,fileEncoding = "UTF-8")
clinical1 <- clinical1[rownames(pro_num_each_pt),]
pro_num_each_pt$age <- clinical1$Age
pro_num_each_pt$sex <- clinical1$Sex
pro_num_each_pt$smoking <- clinical1$Smoking
pro_num_each_pt$loca <- clinical1$`Tumor location`
pro_num_each_pt$age[clinical1$Age > 60] <- ">60"
pro_num_each_pt$age[clinical1$Age <= 60] <- "<=60"
# age
ggplot() +
  geom_half_boxplot(data = pro_num_each_pt, 
                    aes(x=age, y=`FALSE`, fill = age),
                    width = 0.4, position = position_dodge(0.7), 
                    alpha=1, outlier.shape = NA)+
  geom_half_point(data = pro_num_each_pt, 
                  aes(x=age, y=`FALSE`, color=age), 
                  transformation = position_jitter(width = 0.1,height = 0),
                  shape = 16, alpha=1, size = 0.8)+
  stat_compare_means(data = pro_num_each_pt, 
                     aes(x=age, y=`FALSE`),method = "wilcox.test")+
  scale_fill_manual(values = c("#B4BDFF","#FFD28F"))+
  scale_color_manual(values = c("#B4BDFF","#FFD28F"))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 0, vjust = 0.6),
        legend.position="none")+
  labs(y='Protein counts')
# sex
ggplot() +
  geom_half_boxplot(data = pro_num_each_pt, 
                    aes(x=sex, y=`FALSE`, fill = sex),
                    width = 0.4, position = position_dodge(0.7), 
                    alpha=1, outlier.shape = NA)+
  geom_half_point(data = pro_num_each_pt, 
                  aes(x=sex, y=`FALSE`, color=sex), 
                  transformation = position_jitter(width = 0.1,height = 0),
                  shape = 16, alpha=1, size = 0.8)+
  stat_compare_means(data = pro_num_each_pt, 
                     aes(x=sex, y=`FALSE`),method = "wilcox.test")+
  scale_fill_manual(values = c("#B4BDFF","#FFD28F"))+
  scale_color_manual(values = c("#B4BDFF","#FFD28F"))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 0, vjust = 0.6),
        legend.position="none")+
  labs(y='Protein counts')
# smoking
ggplot() +
  geom_half_boxplot(data = pro_num_each_pt, 
                    aes(x=smoking, y=`FALSE`, fill = smoking),
                    width = 0.4, position = position_dodge(0.7), 
                    alpha=1, outlier.shape = NA)+
  geom_half_point(data = pro_num_each_pt, 
                  aes(x=smoking, y=`FALSE`, color=smoking), 
                  transformation = position_jitter(width = 0.1,height = 0),
                  shape = 16, alpha=1, size = 0.8)+
  stat_compare_means(data = pro_num_each_pt, 
                     aes(x=smoking, y=`FALSE`),method = "wilcox.test")+
  scale_fill_manual(values = c("#D3E2F2","#FFCCCC"))+
  scale_color_manual(values = c("#D3E2F2","#FFCCCC"))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 0, vjust = 0.6),
        legend.position="none")+
  labs(y='Protein counts')
# Tumor location
ggplot() +
  geom_half_boxplot(data = pro_num_each_pt, 
                    aes(x=loca, y=`FALSE`, fill = loca),
                    width = 0.4, position = position_dodge(0.7), 
                    alpha=1, outlier.shape = NA)+
  geom_half_point(data = pro_num_each_pt, 
                  aes(x=loca, y=`FALSE`, color=loca), 
                  transformation = position_jitter(width = 0.1,height = 0),
                  shape = 16, alpha=1, size = 0.8)+
  stat_compare_means(data = pro_num_each_pt, 
                     aes(x=loca, y=`FALSE`),method = "wilcox.test")+
  scale_fill_manual(values = c("#B4BDFF","#FFD28F"))+
  scale_color_manual(values = c("#B4BDFF","#FFD28F"))+
  theme_bw()+ 
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 0, vjust = 0.6),
        legend.position="none")+
  labs(y='Protein counts')


###  DreamAI  ###
library(DreamAI)
library(Rcpp)
require("cluster")
require("survival")
require("randomForest")
require("missForest")
require("glmnet")
require("Rcpp")
require("foreach")
require("itertools")
require("iterators")
require("Matrix")
require("devtools")
require("impute")
require("remotes")
data <- read.csv("134SCLC_7715protein_normalized_and_log2.csv",
                 stringsAsFactors = F,check.names = F,
                 row.names = 1)
data1 <- data[,3:136]
data1 <- as.matrix(data1)
rownames(data1) <- data$`Gene name`

impute <- DreamAI(data1,k= 10,maxiter_MF = 10,ntree = 100,maxnodes = NULL,
                  maxiter_ADMIN=30,tol=10^(-2),gamma_ADMIN=NA,
                  gamma=50,CV=FALSE,fillmethod="row_mean",
                  maxiter_RegImpute=10,conv_nrmse = 1e-6,iter_SpectroFM=40, 
                  method = c("KNN"),out="Ensemble")
data_dreamai <- impute$Ensemble
rownames(data_dreamai) <- rownames(data)
data_dreamai <- cbind(data$`Protein description`,data$`Gene name`,data_dreamai)
colnames(data_dreamai)[1:2] <- c('Protein description','Gene name')


###  spearman analysis 1D  ###
library(ggplot2)
library(pheatmap)
data <- read.csv("dreamai_134SCLC_7715protein_raw_intensity.csv",
                 stringsAsFactors = F,check.names = F,
                 row.names = 1)
data <- data[,c(-1,-2)]
spea.mat <- cor(data, method = "spearman")
spea.mat <- spea.mat[upper.tri(spea.mat)]
data_spear <- data.frame(cor = spea.mat,
                         cor1 = spea.mat)
ggplot(data_spear, aes(x=cor)) + 
  geom_histogram(binwidth = 0.001, aes(y=..density..), 
                 fill="#EEC1EA", color=NA, alpha=0.6) +
  theme_bw() + labs(x="Spearman's Correlation", y="Density")
mean(data_spear$cor)
median(data_spear$cor)
min(data_spear$cor)
max(data_spear$cor)
spea.mat <- cor(data, method = "spearman")
red <- "#D94E48";
blue <- "#5175A4";
white <- rgb(255,255,255,maxColorValue = 255)
color_1 <- rev(brewer.pal(10,"RdYlGn"))
heat <- pheatmap(spea.mat,
                 color  = color_1,
                 clustering_method = "ward.D",
                 border_color = "grey60",
                 cluster_cols = F, cluster_rows = F,
                 show_rownames = T, show_colnames = T)



###   1E   ###
library(ggridges)
library(ggplot2)
library(reshape2)
data <- read.csv("dreamai_134SCLC_7715protein_normalized_and_log2.csv",
                 stringsAsFactors = F,check.names = F,
                 row.names = 1)
data <- data[,c(-1,-2)]
data <- as.matrix(t(data))
clinical1 <- read.csv("134SCLC_clinical.csv",stringsAsFactors = F,
                      row.names = 1,check.names = F,fileEncoding = "UTF-8")
clinical1 <- clinical1[rownames(data),]
clinical1$age_group <- clinical1$Age
clinical1$age_group[clinical1$Age > 60] <- ">60"
clinical1$age_group[clinical1$Age <= 60] <- "<=60"
test_wide <- melt(data, varnames = c("Patient","type"),value.name="exp")
test_wide <- as.data.frame(test_wide)
test_wide$age_group <- rep(clinical1$age_group,7715)
test_wide$sex <- rep(clinical1$Sex,7715)
test_wide$smoking <- rep(clinical1$Smoking,7715)
test_wide$loca <- rep(clinical1$`Tumor location`,7715)
# age
a1 <- ggplot(test_wide, aes(x = exp, y = age_group,fill = age_group)) +
  geom_density_ridges(alpha=0.3) + 
  theme_ridges() + theme_bw() + xlab("log2 Protein indensity")+ 
  theme(legend.position="none")
# sex
a2 <- ggplot(test_wide, aes(x = exp, y = sex,fill = sex)) +
  geom_density_ridges(alpha=0.3) + 
  theme_ridges() + theme_bw() + xlab("log2 Protein indensity")+ 
  theme(legend.position="none")
# smoking
a3 <- ggplot(test_wide, aes(x = exp, y = smoking,fill = smoking)) +
  geom_density_ridges(alpha=0.3) + 
  theme_ridges() + theme_bw() + xlab("log2 Protein indensity")+ 
  theme(legend.position="none")
# location
a4 <- ggplot(test_wide, aes(x = exp, y = loca,fill = loca)) +
  geom_density_ridges(alpha=0.3) + 
  theme_ridges() + theme_bw() + xlab("log2 Protein indensity")+ 
  theme(legend.position="none")
ggarrange(a1,a2,a3,a4,nrow = 4,ncol = 1)



###   1F   ###
library(ggplot2)
library(Rtsne)
data <- read.csv("dreamai_134SCLC_7715protein_normalized_and_log2.csv",
                 stringsAsFactors = F,check.names = F,
                 row.names = 1)
data.new <- data[,c(-1,-2)]
data.linshi <- apply(data.new,1,scale)
data.linshi <- t(data.linshi)
colnames(data.linshi) <- colnames(data.new)
tsne_out = Rtsne(t(data.linshi),dims = 2,pca = T,
                 max_iter = 1000,theta = 0.5,
                 perplexity = 10,verbose = F) # 进行t-SNE降维分析
tsne_result = as.data.frame(tsne_out$Y)
colnames(tsne_result) = c("tSNE1","tSNE2")
tsne_result$patient <- substr(colnames(data.linshi),1,1)
ggplot(tsne_result,aes(tSNE1,tSNE2,color=patient)) +
  geom_point(size = 3)+ theme_bw() + 
  scale_color_manual(values = c("#556B2F"))
