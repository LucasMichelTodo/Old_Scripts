source("http://bioconductor.org/biocLite.R")
library(Biobase)
library(reshape2)
library(ggplot2)
library(ggfortify)
library(VennDiagram)
library(dplyr)
library(tidyr)

load("/home/lucas/ISGlobal/Arrays/Eli_Arrays/imputedData_geneLevel.RData")
load("/home/lucas/ISGlobal/Arrays/Eli_Arrays/normalizedData_geneLevel.RData")

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
  scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0,0)) +
  theme(axis.text.y = element_text(colour=ifelse(var, 'red', 'black')))
plot_t1

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

