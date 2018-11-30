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
library(pca3d)
library(venneuler)
library(ggdendro)

load("/home/lucas/ISGlobal/Arrays/Eli_Arrays/imputedData_geneLevel.RData")
load("/home/lucas/ISGlobal/Arrays/Eli_Arrays/normalizedData_geneLevel.RData")

updown_num <- 50

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

#### Max. Fold-Change (any TP) ####
## >3 FC genes
top_10E_up <- imputed_10E_difs[imputed_10E_difs$Max_Dif > log2(3),]
top_10G_up <- imputed_10G_difs[imputed_10G_difs$Max_Dif > log2(3),]
top_EK0_up <- imputed_EK0_difs[imputed_EK0_difs$Max_Dif > log2(3),]


s10E_s10G <- intersect(top_10E_up$Gene, top_10G_up$Gene)
s10G_EK0 <- intersect(top_10G_up$Gene, top_EK0_up$Gene)
s10E_EK0 <- intersect(top_10E_up$Gene,top_EK0_up$Gene)
s10E_s10G_EK0 <- intersect(intersect(top_10E_up$Gene, top_10G_up$Gene), top_EK0_up$Gene)

grid.newpage()
draw.triple.venn(area1 = dim(top_10E_up)[1], area2 = dim(top_10G_up)[1], area3 = dim(top_EK0_up)[1], n12 = length(s10E_s10G), n23 = length(s10G_EK0), n13 = length(s10E_EK0), 
                 n123 = length(s10E_s10G_EK0), category = c("10E", "10G", "EK0"), lty = "blank", 
                 fill = c("#009688", "#FFC107", "#CDDC39"), alpha = 0.5, euler.d=TRUE, scaled=TRUE, overrideTriple ="any")

v <- venneuler(c("10E" = dim(top_10E_up)[1], "10G" = dim(top_10G_up)[1], "EK0" = dim(top_EK0_up)[1], 
                 "10E&10G" = length(s10E_s10G), 
                 "10G&EK0" = length(s10G_EK0), 
                 "10E&EK0" = length(s10E_EK0), 
                 "10E&10G&EK0" = length(s10E_s10G_EK0)))
plot(v)

## <1/3 FC genes
top_10E_down <- imputed_10E_difs[imputed_10E_difs$Max_Dif < log2(1/3),]
top_10G_down <- imputed_10G_difs[imputed_10G_difs$Max_Dif < log2(1/3),]
top_EK0_down <- imputed_EK0_difs[imputed_EK0_difs$Max_Dif < log2(1/3),]


s10E_s10G <- intersect(top_10E_down$Gene, top_10G_down$Gene)
s10G_EK0 <- intersect(top_10G_down$Gene, top_EK0_down$Gene)
s10E_EK0 <- intersect(top_10E_down$Gene,top_EK0_down$Gene)
s10E_s10G_EK0 <- intersect(intersect(top_10E_down$Gene, top_10G_down$Gene), top_EK0_down$Gene)

grid.newpage()
draw.triple.venn(area1 = dim(top_10E_down)[1], area2 = dim(top_10G_down)[1], area3 = dim(top_EK0_down)[1], n12 = length(s10E_s10G), n23 = length(s10G_EK0), n13 = length(s10E_EK0), 
                 n123 = length(s10E_s10G_EK0), category = c("10E", "10G", "EK0"), lty = "blank", 
                 fill = c("#009688", "#FFC107", "#CDDC39"), alpha = 0.5)

v <- venneuler(c("10E" = dim(top_10E_down)[1], "10G" = dim(top_10G_down)[1], "EK0" = dim(top_EK0_down)[1], 
                 "10E&10G" = length(s10E_s10G), 
                 "10G&EK0" = length(s10G_EK0), 
                 "10E&EK0" = length(s10E_EK0), 
                 "10E&10G&EK0" = length(s10E_s10G_EK0)))
plot(v)

#### Fold-Change by TP ####
## >3 FC in TP 1 genes
top_10E_up <- imputed_10E_difs[imputed_10E_difs$s10E_T1_HS > log2(3),]
top_10G_up <- imputed_10G_difs[imputed_10G_difs$s10G_T1_HS > log2(3),]
top_EK0_up <- imputed_EK0_difs[imputed_EK0_difs$EK0_T1_HS > log2(3),]


s10E_s10G <- intersect(top_10E_up$Gene, top_10G_up$Gene)
s10G_EK0 <- intersect(top_10G_up$Gene, top_EK0_up$Gene)
s10E_EK0 <- intersect(top_10E_up$Gene,top_EK0_up$Gene)
s10E_s10G_EK0 <- intersect(intersect(top_10E_up$Gene, top_10G_up$Gene), top_EK0_up$Gene)

grid.newpage()
draw.triple.venn(area1 = dim(top_10E_up)[1], area2 = dim(top_10G_up)[1], area3 = dim(top_EK0_up)[1], n12 = length(s10E_s10G), n23 = length(s10G_EK0), n13 = length(s10E_EK0), 
                 n123 = length(s10E_s10G_EK0), category = c("10E", "10G", "EK0"), lty = "blank", 
                 fill = c("#009688", "#FFC107", "#CDDC39"), alpha = 0.5, euler.d=TRUE)

v <- venneuler(c("10E" = dim(top_10E_up)[1], "10G" = dim(top_10G_up)[1], "EK0" = dim(top_EK0_up)[1], 
                 "10E&10G" = length(s10E_s10G), 
                 "10G&EK0" = length(s10G_EK0), 
                 "10E&EK0" = length(s10E_EK0), 
                 "10E&10G&EK0" = length(s10E_s10G_EK0)))
plot(v)

## <1/3 FC genes
top_10E_down <- imputed_10E_difs[imputed_10E_difs$s10E_T1_HS < log2(1/3),]
top_10G_down <- imputed_10G_difs[imputed_10G_difs$s10G_T1_HS < log2(1/3),]
top_EK0_down <- imputed_EK0_difs[imputed_EK0_difs$EK0_T1_HS < log2(1/3),]


s10E_s10G <- intersect(top_10E_down$Gene, top_10G_down$Gene)
s10G_EK0 <- intersect(top_10G_down$Gene, top_EK0_down$Gene)
s10E_EK0 <- intersect(top_10E_down$Gene,top_EK0_down$Gene)
s10E_s10G_EK0 <- intersect(intersect(top_10E_down$Gene, top_10G_down$Gene), top_EK0_down$Gene)

grid.newpage()
draw.triple.venn(area1 = dim(top_10E_down)[1], area2 = dim(top_10G_down)[1], area3 = dim(top_EK0_down)[1], n12 = length(s10E_s10G), n23 = length(s10G_EK0), n13 = length(s10E_EK0), 
                 n123 = length(s10E_s10G_EK0), category = c("10E", "10G", "EK0"), lty = "blank", 
                 fill = c("#009688", "#FFC107", "#CDDC39"), alpha = 0.5, euler.d=TRUE)

v <- venneuler(c("10E" = dim(top_10E_down)[1], "10G" = dim(top_10G_down)[1], "EK0" = dim(top_EK0_down)[1], 
                 "10E&10G" = length(s10E_s10G), 
                 "10G&EK0" = length(s10G_EK0), 
                 "10E&EK0" = length(s10E_EK0), 
                 "10E&10G&EK0" = length(s10E_s10G_EK0)))
plot(v)

## >3 FC in TP 2 genes
top_10E_up <- imputed_10E_difs[imputed_10E_difs$s10E_T2_HS > log2(3),]
top_10G_up <- imputed_10G_difs[imputed_10G_difs$s10G_T2_HS > log2(3),]
top_EK0_up <- imputed_EK0_difs[imputed_EK0_difs$EK0_T2_HS > log2(3),]


s10E_s10G <- intersect(top_10E_up$Gene, top_10G_up$Gene)
s10G_EK0 <- intersect(top_10G_up$Gene, top_EK0_up$Gene)
s10E_EK0 <- intersect(top_10E_up$Gene,top_EK0_up$Gene)
s10E_s10G_EK0 <- intersect(intersect(top_10E_up$Gene, top_10G_up$Gene), top_EK0_up$Gene)

grid.newpage()
draw.triple.venn(area1 = dim(top_10E_up)[1], area2 = dim(top_10G_up)[1], area3 = dim(top_EK0_up)[1], n12 = length(s10E_s10G), n23 = length(s10G_EK0), n13 = length(s10E_EK0), 
                 n123 = length(s10E_s10G_EK0), category = c("10E", "10G", "EK0"), lty = "blank", 
                 fill = c("#009688", "#FFC107", "#CDDC39"), alpha = 0.5, euler.d=TRUE)

v <- venneuler(c("10E" = dim(top_10E_up)[1], "10G" = dim(top_10G_up)[1], "EK0" = dim(top_EK0_up)[1], 
                 "10E&10G" = length(s10E_s10G), 
                 "10G&EK0" = length(s10G_EK0), 
                 "10E&EK0" = length(s10E_EK0), 
                 "10E&10G&EK0" = length(s10E_s10G_EK0)))
plot(v)

## <1/3 FC genes
top_10E_down <- imputed_10E_difs[imputed_10E_difs$s10E_T2_HS < log2(1/3),]
top_10G_down <- imputed_10G_difs[imputed_10G_difs$s10G_T2_HS < log2(1/3),]
top_EK0_down <- imputed_EK0_difs[imputed_EK0_difs$EK0_T2_HS < log2(1/3),]


s10E_s10G <- intersect(top_10E_down$Gene, top_10G_down$Gene)
s10G_EK0 <- intersect(top_10G_down$Gene, top_EK0_down$Gene)
s10E_EK0 <- intersect(top_10E_down$Gene,top_EK0_down$Gene)
s10E_s10G_EK0 <- intersect(intersect(top_10E_down$Gene, top_10G_down$Gene), top_EK0_down$Gene)

grid.newpage()
draw.triple.venn(area1 = dim(top_10E_down)[1], area2 = dim(top_10G_down)[1], area3 = dim(top_EK0_down)[1], n12 = length(s10E_s10G), n23 = length(s10G_EK0), n13 = length(s10E_EK0), 
                 n123 = length(s10E_s10G_EK0), category = c("10E", "10G", "EK0"), lty = "blank", 
                 fill = c("#009688", "#FFC107", "#CDDC39"), alpha = 0.5, euler.d=TRUE)

v <- venneuler(c("10E" = dim(top_10E_down)[1], "10G" = dim(top_10G_down)[1], "EK0" = dim(top_EK0_down)[1], 
                 "10E&10G" = length(s10E_s10G), 
                 "10G&EK0" = length(s10G_EK0), 
                 "10E&EK0" = length(s10E_EK0), 
                 "10E&10G&EK0" = length(s10E_s10G_EK0)))
plot(v)

## >3 FC in TP 3 genes
top_10E_up <- imputed_10E_difs[imputed_10E_difs$s10E_T3_HS > log2(3),]
top_10G_up <- imputed_10G_difs[imputed_10G_difs$s10G_T3_HS > log2(3),]
top_EK0_up <- imputed_EK0_difs[imputed_EK0_difs$EK0_T3_HS > log2(3),]


s10E_s10G <- intersect(top_10E_up$Gene, top_10G_up$Gene)
s10G_EK0 <- intersect(top_10G_up$Gene, top_EK0_up$Gene)
s10E_EK0 <- intersect(top_10E_up$Gene,top_EK0_up$Gene)
s10E_s10G_EK0 <- intersect(intersect(top_10E_up$Gene, top_10G_up$Gene), top_EK0_up$Gene)

grid.newpage()
draw.triple.venn(area1 = dim(top_10E_up)[1], area2 = dim(top_10G_up)[1], area3 = dim(top_EK0_up)[1], n12 = length(s10E_s10G), n23 = length(s10G_EK0), n13 = length(s10E_EK0), 
                 n123 = length(s10E_s10G_EK0), category = c("10E", "10G", "EK0"), lty = "blank", 
                 fill = c("#009688", "#FFC107", "#CDDC39"), alpha = 0.5, euler.d=TRUE)

v <- venneuler(c("10E" = dim(top_10E_up)[1], "10G" = dim(top_10G_up)[1], "EK0" = dim(top_EK0_up)[1], 
                 "10E&10G" = length(s10E_s10G), 
                 "10G&EK0" = length(s10G_EK0), 
                 "10E&EK0" = length(s10E_EK0), 
                 "10E&10G&EK0" = length(s10E_s10G_EK0)))
plot(v)

## <1/3 FC genes
top_10E_down <- imputed_10E_difs[imputed_10E_difs$s10E_T3_HS < log2(1/3),]
top_10G_down <- imputed_10G_difs[imputed_10G_difs$s10G_T3_HS < log2(1/3),]
top_EK0_down <- imputed_EK0_difs[imputed_EK0_difs$EK0_T3_HS < log2(1/3),]


s10E_s10G <- intersect(top_10E_down$Gene, top_10G_down$Gene)
s10G_EK0 <- intersect(top_10G_down$Gene, top_EK0_down$Gene)
s10E_EK0 <- intersect(top_10E_down$Gene,top_EK0_down$Gene)
s10E_s10G_EK0 <- intersect(intersect(top_10E_down$Gene, top_10G_down$Gene), top_EK0_down$Gene)

grid.newpage()
draw.triple.venn(area1 = dim(top_10E_down)[1], area2 = dim(top_10G_down)[1], area3 = dim(top_EK0_down)[1], n12 = length(s10E_s10G), n23 = length(s10G_EK0), n13 = length(s10E_EK0), 
                 n123 = length(s10E_s10G_EK0), category = c("10E", "10G", "EK0"), lty = "blank", 
                 fill = c("#009688", "#FFC107", "#CDDC39"), alpha = 0.5, euler.d=TRUE)

v <- venneuler(c("10E" = dim(top_10E_down)[1], "10G" = dim(top_10G_down)[1], "EK0" = dim(top_EK0_down)[1], 
                 "10E&10G" = length(s10E_s10G), 
                 "10G&EK0" = length(s10G_EK0), 
                 "10E&EK0" = length(s10E_EK0), 
                 "10E&10G&EK0" = length(s10E_s10G_EK0)))
plot(v)

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
  p <- p + geom_point(aes(color = Soque, shape = Soque)) + geom_path() + scale_linetype_manual(values=c("dashed", "solid"))
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

############################### Heatmap + Dendogram ##########################
##Create df -> maximal differences anywhere
all_DIF <- cbind(df_DIF_10E[,2:4], df_DIF_10G[,2:4], df_DIF_EK0[,c(2:4, 6:7)])
all_DIF$Max_Dif <- as.numeric(apply(all_DIF[,1:9], 1, function(x) x[which.max(abs(x))]))
top_difs <- arrange(all_DIF, -abs(Max_Dif))[1:100,]
top_difs <- top_difs[complete.cases(top_difs),]

exp_mtx <- as.matrix(top_difs[complete.cases(top_difs[,1:9]),1:9])
rownames(exp_mtx) <- top_difs[complete.cases(top_difs[,1:9]),]$Gene

heatmap(exp_mtx, Colv=NA, scale='none')
heatmap.2(exp_mtx,dendrogram='row', Rowv=TRUE, Colv=NA,trace='none', key = TRUE, col=redgreen(75))

############################### Heatmap_10E hclust + ggplot2 ##########################

##Create df -> maximal differences in 10E
updown10E_all <- left_join(updown_10E[,-1], df_DIF_10G[,c(2:4, 6)], by = "Gene")
updown10E_all <- left_join(updown10E_all, df_DIF_EK0[,c(2:4, 6)], by = "Gene")
updown10E_all <- updown10E_all[,c(1,2,3,7,8,9,10,11,12,4,5,6)]
updown10E_all <- updown10E_all[complete.cases(updown10E_all),]

# Heatmap
data <- scale(updown10E_all[,1:9])
ord <- hclust( dist(data, method = "euclidean"), method = "sinle" )$order
ord

heatmap_df <- melt(updown10E_all, id.vars = "Gene", measure.vars = c("s10E_T1_HS", "s10E_T2_HS", "s10E_T3_HS", 
                                                                "s10G_T1_HS", "s10G_T2_HS", "s10G_T3_HS",
                                                                "EK0_T1_HS", "EK0_T2_HS", "EK0_T3_HS"))
## Order
heatmap_df$Gene <- factor(updown10E_all$Gene, levels = updown10E_all$Gene[ord], ordered = TRUE)
## Add strain
heatmap_df$Strain <- c(rep("10E", table(heatmap_df$variable)[1]+table(heatmap_df$variable)[2]+table(heatmap_df$variable)[3]),
                       rep("10G", table(heatmap_df$variable)[4]+table(heatmap_df$variable)[5]+table(heatmap_df$variable)[6]), 
                       rep("EK0", table(heatmap_df$variable)[7]+table(heatmap_df$variable)[8]+table(heatmap_df$variable)[9]))

## Heatmap
plot_all <- ggplot(heatmap_df, aes(x = variable, y = Gene, fill = value)) + 
  geom_tile() + 
  scale_fill_gradient2(low = "red", high = "green", mid = "black") +
  scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0,0)) +
  facet_wrap(~Strain, scales="free_x") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
plot_all

############################### Heatmap_all hclust + ggplot2 ##########################

##Create df -> maximal differences anywhere
all_DIF <- cbind(df_DIF_10E[,2:4], df_DIF_10G[,2:4], df_DIF_EK0[,c(2:4, 6:7)])
all_DIF$Max_Dif <- as.numeric(apply(all_DIF[,1:9], 1, function(x) x[which.max(abs(x))]))
top_difs <- arrange(all_DIF, -Max_Dif)
top_difs <- top_difs[abs(top_difs$Max_Dif) > log2(4),]
top_difs <- top_difs[complete.cases(top_difs),]

## GGPLOT2 hclust
data <- scale(top_difs[,1:9])
ord <- hclust( dist(data, method = "euclidean"), method = "ward.D" )$order
ord

heatmap_df <- melt(top_difs, id.vars = "Gene", measure.vars = c("s10E_T1_HS", "s10E_T2_HS", "s10E_T3_HS", 
                                                                     "s10G_T1_HS", "s10G_T2_HS", "s10G_T3_HS",
                                                                     "EK0_T1_HS", "EK0_T2_HS", "EK0_T3_HS"))
## Order
heatmap_df$Gene <- factor(top_difs$Gene, levels = top_difs$Gene[ord], ordered = TRUE)
## Add strain
heatmap_df$Strain <- c(rep("10E", table(heatmap_df$variable)[1]+table(heatmap_df$variable)[2]+table(heatmap_df$variable)[3]),
                       rep("10G", table(heatmap_df$variable)[4]+table(heatmap_df$variable)[5]+table(heatmap_df$variable)[6]), 
                       rep("EK0", table(heatmap_df$variable)[7]+table(heatmap_df$variable)[8]+table(heatmap_df$variable)[9]))

## Heatmap
plot_all <- ggplot(heatmap_df, aes(x = variable, y = Gene, fill = value)) + 
  geom_tile() + 
  scale_fill_gradient2(low = "red", high = "green", mid = "black") +
  scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0,0)) +
  facet_wrap(~Strain, scales="free_x") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
plot_all

plot(hclust( dist(data, method = "euclidean"), method = "ward.D" ), cex = 0.6)

############################### PCA #################################

noNA <- xgene[complete.cases(exprs(xgene))]
df <- t(exprs(noNA))
df<- as.data.frame(df)

df$Strain <- noNA@phenoData@data$subclone
df$Status <- noNA@phenoData@data$status
df$Time <- noNA@phenoData@data$time_teor
df$Group <- noNA@phenoData@data$soca

# autoplot(prcomp(df[,1:5315]), data = df , colour = 'Strain')
# autoplot(prcomp(df[,1:5315]), data = df , colour = 'Status')
# autoplot(prcomp(df[,1:5315]), data = df , colour = 'Time')

autoplot(prcomp(df[,1:5315]), data = df , colour = 'Strain', size = "Time", shape = "Status")

pca <- prcomp(df[,1:5315])
df_pca <- as.data.frame(pca$x)
df_pca$Strain <- noNA@phenoData@data$subclone
df_pca$Status <- noNA@phenoData@data$status
df_pca$Time <- noNA@phenoData@data$time_teor
df_pca$Group <- noNA@phenoData@data$soca


p <- ggplot(df_pca, aes(x=PC1,y=PC2, col = Strain, linetype = Status, group = Group))
p <- p + geom_point(aes(color = Strain, shape = Strain), size=2) 
p <- p + geom_path() 
p <- p + scale_linetype_manual(values=c("dashed", "solid")) 
print(p)


pca3d(prcomp(df[,1:5315]), 
      col= c("palegreen1", "palegreen1","seagreen3", "seagreen3", "green4", "green4", "darkgreen", "darkgreen",
             "lightblue", "lightblue","steelblue2", "steelblue2", "dodgerblue3", "dodgerblue3", "blue4", "blue4",
             "yellow", "yellow","tan2", "tan2", "brown2", "brown2", "firebrick4", "firebrick4"),
      shape = rep(c("sphere", "cube"), 12),
      radius = 2.5)

