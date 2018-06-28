source("http://bioconductor.org/biocLite.R")
#biocLite("Biobase")
library("limma")
#install.packages("Rcpp")
#install.packages("dplyr")
library(dplyr)
library(tidyr)

# eset <- read.delim("/home/lucas/ISGlobal/Arrays/20180417 Read/US10283823_258456110003_S01_GE2_1105_Oct12_2_4.txt",check.names=FALSE,stringsAsFactors=FALSE)
# x <- read.maimages("/home/lucas/ISGlobal/Arrays/20180417 Read/US10283823_258456110003_S01_GE2_1105_Oct12_2_4.txt" ,source="agilent")
# eset <- as(x, "ExpressionSet")
# 
# > as(object, "ExpressionSet")
# 
# y <- backgroundCorrect(x,method="normexp")
# str(x)
# MA <- normalizeWithinArrays(y)
# fit <- lmFit(MA)
# fit <- eBayes(fit)
# topTable(fit)

## Read files
rosetta <- read.table("/home/lucas/ISGlobal/gens_array_rosetta_annotated.txt", sep = "\t")

ind_1_1 <- read.table("/home/lucas/ISGlobal/Arrays/Oriol_heatmaps/US10283823_258456110002/US10283823_258456110002_S01_GE2_1105_Oct12_ACC_2_1_1.txt",sep="\t",stringsAsFactors=FALSE, skip = 9, header = TRUE)
ind_1_2 <- read.table("/home/lucas/ISGlobal/Arrays/Oriol_heatmaps/US10283823_258456110002/US10283823_258456110002_S01_GE2_1105_Oct12_ACC_2_1_2.txt",sep="\t",stringsAsFactors=FALSE, skip = 9, header = TRUE)
ind_1_3<- read.table("/home/lucas/ISGlobal/Arrays/Oriol_heatmaps/US10283823_258456110002/US10283823_258456110002_S01_GE2_1105_Oct12_ACC_2_1_3.txt",sep="\t",stringsAsFactors=FALSE, skip = 9, header = TRUE)

ctl_2_1 <- read.table("/home/lucas/ISGlobal/Arrays/Oriol_heatmaps/US10283823_258456110002/US10283823_258456110002_S01_GE2_1105_Oct12_ACC_2_2_1.txt",sep="\t",stringsAsFactors=FALSE, skip = 9, header = TRUE)
ctl_2_2 <- read.table("/home/lucas/ISGlobal/Arrays/Oriol_heatmaps/US10283823_258456110002/US10283823_258456110002_S01_GE2_1105_Oct12_ACC_2_2_2.txt",sep="\t",stringsAsFactors=FALSE, skip = 9, header = TRUE)
ctl_2_3 <- read.table("/home/lucas/ISGlobal/Arrays/Oriol_heatmaps/US10283823_258456110002/US10283823_258456110002_S01_GE2_1105_Oct12_ACC_2_2_3.txt",sep="\t",stringsAsFactors=FALSE, skip = 9, header = TRUE)

## Create df
df <- cbind(ind_1_1[,c(7,8,11)], ctl_2_1[,11], ind_1_2[,11], ctl_2_2[,11], ind_1_3[,11], ctl_2_3[,11])
df["Gene_name"] <- sapply(df$SystematicName, function(x) roseta[roseta$V1 == x,2])
colnames(df) <- c("ProbeName", "SystematicName", "Ind_1", "Ctl_1", "Ind_2", "Ctl_2", "Ind_3", "Ctl_3", "Gene_name")

## Create eSet
exprsx <- as.matrix(df[,c(-1,-2,-9)])
fdata <- new("AnnotatedDataFrame",df[,c(1:2,9)])
time <- factor(c(1,1,2,2,3,3))
soca <- factor(rep(c("Ind", "Ctl"),3))
pdata <- data.frame(soca=soca,time=time); rownames(pdata) <- colnames(df[,c(-1,-2,-9)])
pdata <- new("AnnotatedDataFrame",pdata)
x <- new("ExpressionSet",exprs=exprsx,featureData=fdata,phenoData=pdata)

## Collapse proves into genes
geneid <- fdata@data$Gene_name
geneid <- as.character(geneid)
genesToRename.sd <- NA

myRma <- function(x) {
  if (class(x)=='numeric') {
    ans <- x
  } else {
    ans <- medpolish(x,trace.iter=FALSE,na.rm=TRUE)
    ans <- ans$overall + ans$col
  }
  return(ans)
}

renameGenesAndSummarize <- function(genesToRename.sd,exprsx,geneid,summaryMethod=myRma) {
  # if (!is.na(genesToRename.sd)) {
  #   convTable <- read.csv(genesToRename.sd,header=TRUE)
  #   myBadOligos <- featureNames(x)[featureNames(x) %in% convTable$OligoID]
  #   for (i in 1:length(myBadOligos)) myBadOligos[i] <- as.character(convTable[convTable$OligoID==myBadOligos[i],'NewGeneID.'])
  #   geneid[featureNames(x) %in% convTable$OligoID] <- myBadOligos
  # }
  
  xgene <- by(exprsx[,c(-1,-2,-9)],geneid,myRma) 
  xgene <- do.call('rbind',xgene)
  
  mysd <- function(x) { ans <- ifelse(sum(!is.na(x))==1,0,sd(x,na.rm=TRUE)); return(ans) }
  sdgene <- aggregate(exprsx[,c(-1,-2,-9)],by=list(geneid),FUN=mysd)
  
  names(sdgene)[1] <- 'geneid'
  xgene <- data.frame(geneid=rownames(xgene),xgene); rownames(xgene) <- NULL
  
  fdata <- data.frame(xgene$geneid, "")
  fdata <- new("AnnotatedDataFrame",data.frame(fdata))
  rownames(fdata) <- as.character(xgene$geneid)
  
  exprsxgene <- as.matrix(xgene[,-1])
  rownames(exprsxgene) <- as.character(xgene$geneid);
  eset <- new("ExpressionSet",exprs=exprsxgene,featureData=fdata,phenoData=pdata)
  return(list(eset=eset,sdgene=sdgene,fdata=fdata,geneid=geneid))
}


tmp <- renameGenesAndSummarize(genesToRename.sd=genesToRename.sd,exprsx=df,geneid=geneid,summaryMethod=myRma)
xgene <- tmp[['eset']]; sdgene <- tmp[['sdgene']]; fdata <- tmp[['fdata']]; geneid <- tmp[['geneid']]
rownames(exprs(xgene))

# Estimate times

bozdechPath <- '/home/lucas/ISGlobal/Arrays/AnalisisR/Files2executeProgram/bozdech_Hb3_clean2.csv'
LemieuxFunctionsPath <- '/home/lucas/ISGlobal/Arrays/AnalisisR/Files2executeProgram/lemieux_et_al_pipeline_functions.r'
figuresPath <- '/home/lucas/ISGlobal/Arrays/Oriol_heatmaps'

getTimeEstimation <- function(x,dataPath,functionsPath,figuresPath,B=100) {
  #  x: the expressionSet for which we want to estimate times (our data).
  #  dataPath: path to data that will be used to estimate timepoints (from Bozdech et al)
  #  functionsPath: path to the script containing the functions from Lemieux's paper.
  #  figuresPath: where we want to save the output plots.
  source(functionsPath)
  #  z <- read.csv(dataPath, as.is = T,sep='\t')
  z <- read.csv(dataPath, as.is = T)
  #  colnames(z)[1] <- 'Name'
  #  oldTime <- as.numeric(as.character(pData(x)$time))
  oldTime <- as.numeric(substr(sampleNames(xgene),nchar(sampleNames(xgene)),nchar(sampleNames(xgene))))
  x <- exprs(x)
  x <- data.frame(Name=as.character(rownames(x)),x,stringsAsFactors=FALSE); rownames(x) <- NULL
  data <- sync_data(x, z)
  x <- data[[1]]
  z <- data[[2]]
  x <- ordinal(x, use.name = T)
  z <- ordinal(z, use.name = T)
  #  z.na <- cbind(z[,1:22], rep(NA, nrow(z)), z[,23:27], rep(NA, nrow(z)), z[,28:56])
  z.na <- cbind(z[,1:22], rep(NA, nrow(z)), z[,23:27], rep(NA, nrow(z)), z[,28:46])  
  z <- t(apply(z.na, 1, smooth.missing))
  sigma.epsilon <- 789.9056
  z.smooth <- smooth.ref(z, method = "spline", spar = 0.5)
  z.smooth.hourly <- z.smooth[,ll.par$hourly]
  #  sigma.eta <- mean(sd(z[,11:ncol(z)] - z.smooth.hourly, na.rm = T), na.rm=T)
  sigma.eta <- mean(sd(z - z.smooth.hourly, na.rm = T), na.rm=T)  
  new.sigma <- sqrt(sigma.eta^2 + sigma.epsilon^2)
  ll <- compute.ll(x = x, z = z.smooth, sigma = new.sigma, bootstrap = T, B = B, sample.rate = 0.50)
  myTimes <- mle(ll)
  png(file.path(figuresPath,'defaultPlots1.png'))
  plot.ll(ll)
  dev.off()
  png(file.path(figuresPath,'defaultPlots2.png'))
  plot.mle(ll)
  dev.off()
  png(file.path(figuresPath,'ownPlots1.png'))
  plot(density(myTimes),main='Estimated times density')
  dev.off()
  png(file.path(figuresPath,'ownPlots2.png'))
  plot(oldTime, as.numeric(myTimes),xlab='Old times',ylab='Estimated times',xlim=c(-5,50),ylim=c(-5,50))
  abline(0,1,col=2,lwd=2)
  abline(v=oldTime,lwd=0.5,lty=3)
  dev.off()
  return(myTimes)
}

estimatedTimes <- getTimeEstimation(xgene,bozdechPath,LemieuxFunctionsPath,file.path(figuresPath),B=100)
estimatedTimes <- estimatedTimes + c(0,0,48,48,48,48)
pData(xgene)$etime <- estimatedTimes

### HEATMAPS

write.table(heatmap_df, "/home/lucas/Documents/oriol_heatmap.csv", sep = "\t")
heatmap_df <- data.frame(c(exprs(xgene)[,1] - exprs(xgene)[,2]), c(exprs(xgene)[,3] - exprs(xgene)[,4]), c(exprs(xgene)[,5] - exprs(xgene)[,6]))
colnames(heatmap_df) <- c("T1", "T2", "T3")
heatmap_df["gene"] <- rownames(heatmap_df)

sel_all <- filter(heatmap_df, abs(T1) > 0.39794 | abs(T2) > 0.39794 | abs(T3) > 0.39794)

sel_top1 <- arrange(heatmap_df, -T1)[c(1:20),]
sel_bottom1 <- arrange(heatmap_df, T1)[c(1:20),]
sel1 <- rbind(sel_top1, sel_bottom1)
sel1 <- arrange(sel1, T1)

# sel1 <- arrange(heatmap_df, -abs(T1))[1:30,]
# sel1 <- arrange(sel1, T1)

sel_top2 <- arrange(heatmap_df, -T2)[c(1:20),]
sel_bottom2 <- arrange(heatmap_df, T2)[c(1:20),]
sel2 <- rbind(sel_top2, sel_bottom2)
sel2 <- arrange(sel2, T2)

sel_top3 <- arrange(heatmap_df, -T3)[c(1:20),]
sel_bottom3 <- arrange(heatmap_df, T3)[c(1:20),]
sel3 <- rbind(sel_top3, sel_bottom3)
sel3 <- arrange(sel3, T3)

heatmap_df_t1 <- melt(sel1)
heatmap_df_t1$gene <- factor(heatmap_df_t1$gene, levels=unique(as.character(heatmap_df_t1$gene)))

heatmap_df_t2 <- melt(sel2)
heatmap_df_t2$gene <- factor(heatmap_df_t2$gene, levels=unique(as.character(heatmap_df_t2$gene)))

heatmap_df_t3 <- melt(sel3)
heatmap_df_t3$gene <- factor(heatmap_df_t3$gene, levels=unique(as.character(heatmap_df_t3$gene)))

plot_t1 <- ggplot(heatmap_df_t1, aes(x = variable, y = gene, fill = value)) + 
  geom_tile() + 
  scale_fill_gradient2(low = "red", high = "green", mid = "black") +
  scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0,0))

plot_t1

plot_t2 <- ggplot(heatmap_df_t2, aes(x = variable, y = gene, fill = value)) + 
  geom_tile() + 
  scale_fill_gradient2(low = "red", high = "green", mid = "black") +
  scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0,0))

plot_t2

plot_t3 <- ggplot(heatmap_df_t3, aes(x = variable, y = gene, fill = value)) + 
  geom_tile() + 
  scale_fill_gradient2(low = "red", high = "green", mid = "black") +
  scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0,0))

plot_t3



## Hierarchical Clustering

hm_all <- melt(sel_all)
plot_all <- ggplot(hm_all, aes(x = variable, y = gene, fill = value)) + 
  geom_tile() + 
  scale_fill_gradient2(low = "red", high = "green", mid = "black") +
  scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0,0))

plot_all

rownames(sel_all) <- sel_all$gene
heatmap(
  as.matrix(sel_all[,1:3]), Rowv=as.dendrogram(hclust(dist(as.matrix(sel_all[,1:3]))),
  Colv=NA, cexRow = 2)
)


plot_grid(plot_t1, plot_t2, plot_t3, align = "v")



