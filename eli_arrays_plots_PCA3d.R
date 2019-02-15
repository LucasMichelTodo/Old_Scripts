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

