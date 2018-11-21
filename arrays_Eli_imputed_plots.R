source("http://bioconductor.org/biocLite.R")
library(Biobase)
library(reshape2)
library(ggplot2)
library(gplots)
library(ggfortify)
library(VennDiagram)
library(dplyr)
library(tidyr)
library(gridExtra)

load("/home/lucas/ISGlobal/Arrays/Eli_Arrays/imputedData_geneLevel.RData")
load("/home/lucas/ISGlobal/Arrays/Eli_Arrays/normalizedData_geneLevel.RData")

updown_num <- 30

############################### Imputed FC ###############

## 10E
eset_10E <- xgene_impute[,xgene_impute@phenoData@data$subclone == "s10E"]

df_10E_CTL <- exprs(eset_10E)[,eset_10E@phenoData@data$status == "CTL"]
df_10E_HS <- exprs(eset_10E)[,eset_10E@phenoData@data$status == "HS"]

df_DIF_10E <- df_10E_HS-df_10E_CTL
df_DIF_10E <- data.frame(df_DIF_10E)

df_DIF_10E$Max_Dif <- as.numeric(apply(df_DIF_10E, 1, function(x) x[which.max(abs(x))]))

df_DIF_10E$Gene <- rownames(df_DIF_10E)
df_DIF_10E$Annot <- xgene_impute@featureData@data$Annot

imputed_10E_difs <- arrange(df_DIF_10E, -abs(Max_Dif))
updown_10E <- arrange(rbind(arrange(df_DIF_10E, -Max_Dif)[1:updown_num,], arrange(df_DIF_10E, Max_Dif)[1:updown_num,]), -Max_Dif)

## 10G
eset_10G <- xgene_impute[,xgene_impute@phenoData@data$subclone == "s10G"]

df_10G_CTL <- exprs(eset_10G)[,eset_10G@phenoData@data$status == "CTL"]
df_10G_HS <- exprs(eset_10G)[,eset_10G@phenoData@data$status == "HS"]

df_DIF_10G <- df_10G_HS-df_10G_CTL
df_DIF_10G <- data.frame(df_DIF_10G)

df_DIF_10G$Max_Dif <- as.numeric(apply(df_DIF_10G, 1, function(x) x[which.max(abs(x))]))

df_DIF_10G$Gene <- rownames(df_DIF_10G)
df_DIF_10G$Annot <- xgene_impute@featureData@data$Annot

imputed_10G_difs <- arrange(df_DIF_10G, -abs(Max_Dif))
updown_10G <- arrange(rbind(arrange(df_DIF_10G, -Max_Dif)[1:updown_num,], arrange(df_DIF_10G, Max_Dif)[1:updown_num,]), -Max_Dif)

## EK0
eset_EK0 <- xgene_impute[,xgene_impute@phenoData@data$subclone == "EK0"]

df_EK0_CTL <- exprs(eset_EK0)[,eset_EK0@phenoData@data$status == "CTL"]
df_EK0_HS <- exprs(eset_EK0)[,eset_EK0@phenoData@data$status == "HS"]

df_DIF_EK0 <- df_EK0_HS-df_EK0_CTL
df_DIF_EK0 <- data.frame(df_DIF_EK0)

df_DIF_EK0$Max_Dif <- as.numeric(apply(df_DIF_EK0, 1, function(x) x[which.max(abs(x))]))

df_DIF_EK0$Gene <- rownames(df_DIF_EK0)
df_DIF_EK0$Annot <- xgene_impute@featureData@data$Annot

imputed_EK0_difs <- arrange(df_DIF_EK0, -abs(Max_Dif))
updown_EK0 <- arrange(rbind(arrange(df_DIF_EK0, -Max_Dif)[1:updown_num,], arrange(df_DIF_EK0, Max_Dif)[1:updown_num,]), -Max_Dif)

############################### Venn Diagrams ###############################
## Top 50

imputed_10E_difs[1:50,]$Gene
imputed_10G_difs[1:50,]$Gene
imputed_EK0_difs[1:50,]$Gene

s10E_s10G <- intersect(imputed_10E_difs[1:50,]$Gene, imputed_10G_difs[1:50,]$Gene)
s10G_EK0 <- intersect(imputed_10G_difs[1:50,]$Gene, imputed_EK0_difs[1:50,]$Gene)
s10E_EK0 <- intersect(imputed_10E_difs[1:50,]$Gene, imputed_EK0_difs[1:50,]$Gene)
s10E_s10G_EK0 <- intersect(intersect(imputed_10E_difs[1:50,]$Gene, imputed_10G_difs[1:50,]$Gene), imputed_EK0_difs[1:50,]$Gene)

grid.newpage()
draw.triple.venn(area1 = 50, area2 = 50, area3 = 50, n12 = length(s10E_s10G), n23 = length(s10G_EK0), n13 = length(s10E_EK0), 
                 n123 = length(s10E_s10G_EK0), category = c("10E", "10G", "EK0"), lty = "blank", 
                 fill = c("#009688", "#FFC107", "#CDDC39"))

top10E <- imputed_10E_difs[1:50,]

## Only 10E genes (AP2-GHS targets?)
top10E[!top10E$Gene %in% s10E_s10G & !top10E$Gene %in% s10E_EK0,]

## >2 FC genes
top_10E <- imputed_10E_difs[imputed_10E_difs$Max_Dif > 1,]
top_10G <- imputed_10G_difs[imputed_10G_difs$Max_Dif > 1,]
top_EK0 <- imputed_EK0_difs[imputed_EK0_difs$Max_Dif > 1,]


s10E_s10G <- intersect(top_10E$Gene, top_10G$Gene)
s10G_EK0 <- intersect(top_10G$Gene, top_EK0$Gene)
s10E_EK0 <- intersect(top_10E$Gene,top_EK0$Gene)
s10E_s10G_EK0 <- intersect(intersect(top_10E$Gene, top_10G$Gene), top_EK0$Gene)

grid.newpage()
draw.triple.venn(area1 = dim(top_10E)[1], area2 = dim(top_10G)[1], area3 = dim(top_EK0)[1], n12 = length(s10E_s10G), n23 = length(s10G_EK0), n13 = length(s10E_EK0), 
                 n123 = length(s10E_s10G_EK0), category = c("10E", "10G", "EK0"), lty = "blank", 
                 fill = c("#009688", "#FFC107", "#CDDC39"))


############################### Plots Custom #############################

##Plots
#for (i in 1:nrow(exprs(xgene))){
for (i in 1:1){
  print(i)
  graf <- melt(exprs(xgene)[i,])
  graf["Time"] <- xgene_impute@phenoData@data$time
  graf["Type"] <- xgene_impute@phenoData@data$status
  graf["Soque"] <- xgene_impute@phenoData@data$subclone
  graf["group"] <- xgene_impute@phenoData@data$soca
  p <- ggplot(graf, aes(x = Time, y = value, col = Soque, linetype = Type, group = group))
  p <- p + geom_point(aes(color = Soque, shape = Soque)) + geom_line() + scale_linetype_manual(values=c("dashed", "solid"))
  print(p)
}


############################### Heatmaps #################################
## Melt 
heatmap_df_10E <- melt(imputed_10E_difs[1:50,], id.vars = "Gene", measure.vars = c("s10E_T1_HS", "s10E_T2_HS", "s10E_T3_HS"))
## Order
heatmap_df_10E$Gene <- factor(heatmap_df_10E$Gene, levels = rev(imputed_10E_difs[1:50,]$Gene), ordered = TRUE)

## Heatmap
plot_10E <- ggplot(heatmap_df_10E, aes(x = variable, y = Gene, fill = value)) + 
  geom_tile() + 
  scale_fill_gradient2(low = "red", high = "green", mid = "black") +
  scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0,0))
plot_10E

## Updown Heatmaps
## Melt 
heatmap_df_10E <- melt(updown_10E, id.vars = "Gene", measure.vars = c("s10E_T1_HS", "s10E_T2_HS", "s10E_T3_HS"))
heatmap_df_10G <- melt(updown_10G, id.vars = "Gene", measure.vars = c("s10G_T1_HS", "s10G_T2_HS", "s10G_T3_HS"))
heatmap_df_EK0 <- melt(updown_EK0, id.vars = "Gene", measure.vars = c("EK0_T1_HS", "EK0_T2_HS", "EK0_T3_HS"))
## Order
heatmap_df_10E$Gene <- factor(updown_10E$Gene, levels = rev(updown_10E$Gene), ordered = TRUE)
heatmap_df_10G$Gene <- factor(updown_10G$Gene, levels = rev(updown_10G$Gene), ordered = TRUE)
heatmap_df_EK0$Gene <- factor(updown_EK0$Gene, levels = rev(updown_EK0$Gene), ordered = TRUE)

## Heatmap
plot_10E <- ggplot(heatmap_df_10E, aes(x = variable, y = Gene, fill = value)) + 
  geom_tile() + 
  scale_fill_gradient2(low = "red", high = "green", mid = "black") +
  scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0,0)) +
  theme(legend.position="top")
#plot_10E

plot_10G <- ggplot(heatmap_df_10G, aes(x = variable, y = Gene, fill = value)) + 
  geom_tile() + 
  scale_fill_gradient2(low = "red", high = "green", mid = "black") +
  scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0,0)) +
  theme(legend.position="top")
#plot_10G

plot_EK0 <- ggplot(heatmap_df_EK0, aes(x = variable, y = Gene, fill = value)) + 
  geom_tile() + 
  scale_fill_gradient2(low = "red", high = "green", mid = "black") +
  scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0,0)) +
  theme(legend.position="top")
#plot_EK0

grid.arrange(plot_10E, plot_10G, plot_EK0, ncol=3)

############################### All Heatmap #################################
##Create df
updown10E_all <- left_join(updown_10E[,-1], df_DIF_10G[,c(2:4, 6)], by = "Gene")
updown10E_all <- left_join(updown10E_all, df_DIF_EK0[,c(2:4, 6)], by = "Gene")
updown10E_all <- updown10E_all[,c(1,2,3,7,8,9,10,11,12,4,5,6)]

## Melt 
heatmap_df <- melt(updown10E_all, id.vars = "Gene", measure.vars = c("s10E_T1_HS", "s10E_T2_HS", "s10E_T3_HS", 
                                                                     "s10G_T1_HS", "s10G_T2_HS", "s10G_T3_HS",
                                                                     "EK0_T1_HS", "EK0_T2_HS", "EK0_T3_HS"))
## Order
heatmap_df$Gene <- factor(updown10E_all$Gene, levels = rev(updown10E_all$Gene), ordered = TRUE)
## Add strain
heatmap_df$Strain <- c(rep("10E", updown_num*6), rep("10G", updown_num*6), rep("EK0", updown_num*6))

## Heatmap
plot_all <- ggplot(heatmap_df, aes(x = variable, y = Gene, fill = value)) + 
  geom_tile() + 
  scale_fill_gradient2(low = "red", high = "green", mid = "black") +
  scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0,0)) +
  facet_wrap(~Strain, scales="free_x")
plot_all

############################### Clustering ##########################
##Create df -> maximal differences in 10E
updown10E_all <- left_join(updown_10E[,-1], df_DIF_10G[,c(2:4, 6)], by = "Gene")
updown10E_all <- left_join(updown10E_all, df_DIF_EK0[,c(2:4, 6)], by = "Gene")
updown10E_all <- updown10E_all[,c(1,2,3,7,8,9,10,11,12,4,5,6)]

# Heatmap
exp_mtx <- as.matrix(scale(updown10E_all[complete.cases(updown10E_all[,1:9]),1:9]))
rownames(exp_mtx) <- updown10E_all[complete.cases(updown10E_all[,1:9]),]$Gene
heatmap(exp_mtx, Colv=NA, scale='none')

##Create df -> maximal differences anywhere
all_DIF <- cbind(df_DIF_10E[,2:4], df_DIF_10G[,2:4], df_DIF_EK0[,c(2:4, 6:7)])
all_DIF$Max_Dif <- as.numeric(apply(all_DIF[,1:9], 1, function(x) x[which.max(abs(x))]))
top_difs <- arrange(all_DIF, -abs(Max_Dif))[1:50,]

# Heatmap
exp_mtx <- as.matrix(top_difs[complete.cases(top_difs[,1:9]),1:9])
rownames(exp_mtx) <- top_difs[complete.cases(top_difs[,1:9]),]$Annot
heatmap(exp_mtx, Colv=NA, scale='none')
heatmap.2(exp_mtx,dendrogram='none', Rowv=TRUE, Colv=NA,trace='none', key = FALSE, col=redgreen(75),
          lhei=c(0.5,25), 
          lwid=c(0.5,25),
          margins = c(8, 16))

############################### PCA #################################

noNA <- xgene[complete.cases(exprs(xgene))]
df <- t(exprs(noNA))
df<- as.data.frame(df)

df$Strain <- noNA@phenoData@data$subclone
df$Status <- noNA@phenoData@data$status
df$Time <- noNA@phenoData@data$time_teor

autoplot(prcomp(df[,1:5315]), data = df , colour = 'Strain')
autoplot(prcomp(df[,1:5315]), data = df , colour = 'Status')
autoplot(prcomp(df[,1:5315]), data = df , colour = 'Time')

autoplot(prcomp(df[,1:5315]), data = df , colour = 'Strain', size = "Time", shape = "Status")

