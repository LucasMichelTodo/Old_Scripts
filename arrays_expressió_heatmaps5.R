source("http://bioconductor.org/biocLite.R")
library(limma)
library(dplyr)
library(tidyr)
library(Biobase)
library(reshape2)
library(ggplot2)

dir <- "/home/lucas/ISGlobal/Arrays/Oriol_heatmaps"
figuresPath <- "/home/lucas/ISGlobal/Arrays/Oriol_heatmaps/Plots/"

############################### Read files and create eSet ##################

## Read files
roseta <- read.csv("/home/lucas/ISGlobal/gens_array_rosetta_annotated.txt", sep = "\t", as.is = TRUE, header = FALSE)
dim(roseta)

noBG <- read.table("/home/lucas/ISGlobal/Arrays/20180417 Read/US10283823_258456110002_S01_GE2_1105_Oct12_1_1.txt",sep="\t",stringsAsFactors=FALSE, skip = 9, header = TRUE)
summary(noBG$gProcessedSignal)
summary(ind_1_1$gProcessedSignal)

head(noBG$gProcessedSignal)
head(ind_1_1$gProcessedSignal)

ind_1_1 <- read.table("/home/lucas/ISGlobal/Arrays/Oriol_heatmaps/US10283823_258456110002/US10283823_258456110002_S01_GE2_1105_Oct12_ACC_2_1_1.txt",sep="\t",stringsAsFactors=FALSE, skip = 9, header = TRUE)
ind_1_2 <- read.table("/home/lucas/ISGlobal/Arrays/Oriol_heatmaps/US10283823_258456110002/US10283823_258456110002_S01_GE2_1105_Oct12_ACC_2_1_2.txt",sep="\t",stringsAsFactors=FALSE, skip = 9, header = TRUE)
ind_1_3<- read.table("/home/lucas/ISGlobal/Arrays/Oriol_heatmaps/US10283823_258456110002/US10283823_258456110002_S01_GE2_1105_Oct12_ACC_2_1_3.txt",sep="\t",stringsAsFactors=FALSE, skip = 9, header = TRUE)

a <- ind_1_1[!ind_1_1$SystematicName %in% no,]
b <- a[a$gMedianSignal < 120 & !startsWith(a$SystematicName, "RC"),]
b["name"] <- unlist(lapply(b$SystematicName, function(x) roseta[roseta$V1 == x,]$V3[1])) 
table(b[!is.na(b$name),]$name)

summary(ind_1_1$gProcessedSignal)

no <- c("NegativeControl", "DarkCorner")


table(ind_1_1[ind_1_1$gIsPosAndSignif == 0,]$LogRatio)

ind_1_1[ind_1_1$SystematicName == "PF11_0112",]
table(ind_1_1$gIsPosAndSignif)
table(ind_1_1[ind_1_1$gIsPosAndSignif == 0,]$SystematicName)
table(ind_1_1$rIsPosAndSignif)

table(ind_1_1$gIsSaturated)
table(ind_1_1$rIsSaturated)
ind_1_1[ind_1_1$rIsSaturated == 1,]

table(ind_1_1$gIsWellAboveBG)
ind_1_1[ind_1_1$gIsWellAboveBG == 0,]
table(ind_1_1$rIsWellAboveBG)

ctl_2_1 <- read.table("/home/lucas/ISGlobal/Arrays/Oriol_heatmaps/US10283823_258456110002/US10283823_258456110002_S01_GE2_1105_Oct12_ACC_2_2_1.txt",sep="\t",stringsAsFactors=FALSE, skip = 9, header = TRUE)
ctl_2_2 <- read.table("/home/lucas/ISGlobal/Arrays/Oriol_heatmaps/US10283823_258456110002/US10283823_258456110002_S01_GE2_1105_Oct12_ACC_2_2_2.txt",sep="\t",stringsAsFactors=FALSE, skip = 9, header = TRUE)
ctl_2_3 <- read.table("/home/lucas/ISGlobal/Arrays/Oriol_heatmaps/US10283823_258456110002/US10283823_258456110002_S01_GE2_1105_Oct12_ACC_2_2_3.txt",sep="\t",stringsAsFactors=FALSE, skip = 9, header = TRUE)

## Create df
df <- cbind(ind_1_1[,c(7,8,11)], ctl_2_1[,11], ind_1_2[,11], ctl_2_2[,11], ind_1_3[,11], ctl_2_3[,11])
df["Gene_id"] <- roseta$V2
df["Gene_name"] <- roseta$V4
df["Annot"] <- roseta$V3

df["Annot"] <- gsub("Plasmodium", "Pl.", df$Annot)
df["Annot"] <- gsub("protein", "prot.", df$Annot)
df["Annot"] <- gsub("membrane", "memb.", df$Annot)
df["Annot"] <- gsub("conserved", "cvd.", df$Annot)
df["Annot"] <- gsub("function", "func.", df$Annot)
df["Annot"] <- gsub("unknown", "ukwn.", df$Annot)
df["Annot"] <- gsub("exported", "xptd.", df$Annot)
df["Annot"] <- gsub("pseudogene", "pseudo", df$Annot)
df["Annot"] <- gsub("putative", "put.", df$Annot)
df["Annot"] <- gsub("%2C", "", df$Annot)

colnames(df) <- c("ProbeName", "SystematicName", "Ind_1", "Ctl_1", "Ind_2", "Ctl_2", "Ind_3", "Ctl_3", "Gene_id", "Gene_name", "Annot")
df[df$ProbeName == "HSDHFR_v7.1_P1/1",]$Gene_id <- "HSDHFR"
df[df$ProbeName == "HSDHFR_v7.1_P1/1",]$Gene_name <- "HSDHFR"
df[df$ProbeName == "HSDHFR_v7.1_P1/1",]$Annot <- "hsdhfr"

# Treure gens no únics
no_unics <- read.table("/home/lucas/ISGlobal/Arrays/gens_no_unics.txt", sep = "\t", as.is = TRUE, header = TRUE)
df <- df[!df$SystematicName %in% no_unics$GeneID,]

# Afegir informació sobre gens variants
variant <- read.csv("/home/lucas/ISGlobal/Gen_Referencies/Gens_variants_extended.txt", header = TRUE, sep = "\t")
df["Variant"] <- df$Gene_id %in% variant$ID

# Afegir families i GO terms
go <- read.csv("/home/lucas/ISGlobal/Gen_Referencies/PlasmoDB-37_Pfalciparum3D7_GO.gaf", header = FALSE, sep = "\t", skip = 1)
go["gene_id"] <- gsub("EuPathDB:", "", go$V17)
go$gene_id <- gsub(".1", "", go$gene_id)
table(df$Gene_id[!is.na(df$Gene_id)] %in% go$gene_id)  

## Create eSet
exprsx <- as.matrix(df[,c(-1,-2,-9,-10,-11,-12)])
fdata <- new("AnnotatedDataFrame",df[,c(1:2,9:12)])
time <- factor(c(1,1,2,2,3,3))
soca <- factor(rep(c("Ind", "Ctl"),3))
pdata <- data.frame(soca=soca,time=time); rownames(pdata) <- colnames(df[,c(-1,-2,-9,-10,-11,-12)])
pdata <- new("AnnotatedDataFrame",pdata)
x <- new("ExpressionSet",exprs=exprsx,featureData=fdata,phenoData=pdata)
save(x,file=file.path(dir,'normalizedData_probeLevel.RData'))

############################### Rename and Summarize ##################################

geneid <- fdata@data$Gene_id
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
  
  xgene <- by(exprsx[,c(-1,-2,-9,-10,-11,-12)],geneid,myRma) 
  xgene <- do.call('rbind',xgene)
  
  mysd <- function(x) { ans <- ifelse(sum(!is.na(x))==1,0,sd(x,na.rm=TRUE)); return(ans) }
  sdgene <- aggregate(exprsx[,c(-1,-2,-9,-10,-11,-12)],by=list(geneid),FUN=mysd)
  
  names(sdgene)[1] <- 'geneid'
  xgene <- data.frame(geneid=rownames(xgene),xgene); rownames(xgene) <- NULL
  
  fdata <- by(df[,9:12],geneid,unique)
  genenames <- names(fdata)
  fdata <- do.call('rbind',fdata)
  
  fdata <- new("AnnotatedDataFrame",data.frame(fdata))
  rownames(fdata) <- as.character(xgene$geneid)
  
  exprsxgene <- as.matrix(xgene[,-1])
  rownames(exprsxgene) <- as.character(xgene$geneid);
  eset <- new("ExpressionSet",exprs=exprsxgene,featureData=fdata,phenoData=pdata)
  return(list(eset=eset,sdgene=sdgene,fdata=fdata,geneid=geneid))
}


tmp <- renameGenesAndSummarize(genesToRename.sd=genesToRename.sd,exprsx=df,geneid=geneid,summaryMethod=myRma)
xgene <- tmp[['eset']]; sdgene <- tmp[['sdgene']]; fdata <- tmp[['fdata']]; geneid <- tmp[['geneid']]
write.csv(data.frame(probe=df[,1],gene=geneid),file.path(figuresPath,'geneProbes.csv'),row.names=FALSE) #probes of each gene
rownames(exprs(xgene))

############################### Estimate times ###############

bozdechPath <- '/home/lucas/ISGlobal/Arrays/AnalisisR/Files2executeProgram/bozdech_Hb3_clean2.csv'
LemieuxFunctionsPath <- '/home/lucas/ISGlobal/Arrays/AnalisisR/Files2executeProgram/lemieux_et_al_pipeline_functions.r'

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
save(estimatedTimes,file=file.path(dir,'estimatedTimes.RData'))
load(file.path(dir,'estimatedTimes.RData'))
write.csv(estimatedTimes,file.path(figuresPath,'Estimated_Times.csv'))
pData(xgene)$time <- estimatedTimes

#save ExpressionSet at gene level
save(xgene,file=file.path(dir,'normalizedData_geneLevel.RData'))
save(xgene.noRatio,file=file.path(dir,'normalizedData_geneLevel_noRatio.RData'))

#boxplot after summarization
pdf(file.path(figuresPath,'boxplot_afterSummarization.pdf'))
boxplot(exprs(xgene),main='summarization method: median poslish')
dev.off()

############################### HEATMAPS #########################################################
load(file.path(dir,'xgene_estimated.RData'))

heatmap_df <- data.frame(c(exprs(xgene)[,1] - exprs(xgene)[,2]), c(exprs(xgene)[,3] - exprs(xgene)[,4]), c(exprs(xgene)[,5] - exprs(xgene)[,6]))
colnames(heatmap_df) <- c("T1", "T2", "T3")

newlist = c()
for (x in 1:length(xgene@featureData@data$Gene_name)){
  if (is.na(xgene@featureData@data$Gene_name[x])){
    newlist = c(newlist,xgene@featureData@data$Annot[x])
  } else {
    newlist = c(newlist,xgene@featureData@data$Gene_name[x])
  }
}

heatmap_df["gene"] <- paste0(newlist, ": ", rownames(heatmap_df))

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


ids <- unlist(lapply(heatmap_df_t1$gene, function(x) strsplit(as.character(x), split = ": ", fixed= TRUE)[[1]][2]))
var <- ids %in% variant$ID

plot_t1 <- ggplot(heatmap_df_t1, aes(x = variable, y = gene, fill = value)) + 
  geom_tile() + 
  scale_fill_gradient2(low = "red", high = "green", mid = "black") +
  scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0,0)) +
  theme(axis.text.y = element_text(colour=ifelse(var, 'red', 'black')))

plot_t1

ids <- unlist(lapply(heatmap_df_t2$gene, function(x) strsplit(as.character(x), split = ": ", fixed= TRUE)[[1]][2]))
var <- ids %in% variant$ID

plot_t2 <- ggplot(heatmap_df_t2, aes(x = variable, y = gene, fill = value)) + 
  geom_tile() + 
  scale_fill_gradient2(low = "red", high = "green", mid = "black") +
  scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0,0)) +
  theme(axis.text.y = element_text(colour=ifelse(var, 'red', 'black')))

plot_t2

ids <- unlist(lapply(heatmap_df_t3$gene, function(x) strsplit(as.character(x), split = ": ", fixed= TRUE)[[1]][2]))
var <- ids %in% variant$ID

plot_t3 <- ggplot(heatmap_df_t3, aes(x = variable, y = gene, fill = value)) + 
  geom_tile() + 
  scale_fill_gradient2(low = "red", high = "green", mid = "black") +
  scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0,0)) +
  theme(axis.text.y = element_text(colour=ifelse(var, 'red', 'black')))

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

############################### SD ##########################
#sd
exprsx <- as.matrix(sdgene[,-1])
rownames(exprsx) <- as.character(sdgene$geneid)
colnames(exprsx) <- sampleNames(xgene)
sdgene <- new("ExpressionSet",exprs=exprsx,featureData=fdata,phenoData=pdata)
save(sdgene,file=file.path(dir,'normalizedData_SDgeneLevel.RData'))
xout <- data.frame(fData(xgene),exprs(xgene))
write.csv(xout,file.path(dir,'normalizedData_geneLevel.csv'),row.names=FALSE)
xout <- data.frame(fData(sdgene),exprs(sdgene))
write.csv(xout,file.path(dir,'normalizedData_SDgeneLevel.csv'),row.names=FALSE)

############################### PCAs ################

#PCA (probe level)
load(file.path(dir,'normalizedData_probeLevel.RData'))
xnm <- x[rowSums(is.na(exprs(x)))==0,]
pcdat <- prcomp(t(exprs(xnm)),na.action=na.pass)
n <- sub('\\.gpr','',rownames(pData(x)))
g <- as.numeric(factor(substring(n,nchar(n)-1,nchar(n))))

pdf(file.path(figuresPath,'pca_probeLevel.pdf'))
xlim <- 1.5*range(pcdat$x[,1]); ylim <- 1.25*range(pcdat$x[,2])
plot(pcdat$x[,1],pcdat$x[,2],col=g,xlab='PC1',ylab='PC2',xlim=xlim,ylim=ylim)
text(pcdat$x[,1],pcdat$x[,2],n,pos=1,cex=.8)
dev.off()

#PCA (gene level)
load(file.path(dir,'normalizedData_geneLevel.RData'))
x <- xgene
xnm <- xgene[rowSums(is.na(exprs(xgene)))==0,]
pcdat <- prcomp(t(exprs(xnm)))
n <- sub('\\.gpr','',rownames(pData(x)))
g <- as.numeric(factor(substring(n,nchar(n)-1,nchar(n))))

pdf(file.path(figuresPath,'pca_geneLevel.pdf'))
xlim <- 1.5*range(pcdat$x[,1]); ylim <- 1.25*range(pcdat$x[,2])
plot(pcdat$x[,1],pcdat$x[,2],col=g,xlab='PC1',ylab='PC2',xlim=xlim,ylim=ylim)
text(pcdat$x[,1],pcdat$x[,2],n,pos=1,cex=.8)
dev.off()

############################### Plots at probe level ##########################################

load(file.path(dir,'normalizedData_probeLevel.RData'))
x <- x[,order(pData(x)$soca,pData(x)$time)]

time <- as.numeric(levels(pData(x)$time))
soques <- levels(pData(x)$soca)
ylim <- range(exprs(x)[!is.na(exprs(x))])

myRange <- range(time)
myPoints.y <- round(seq(range(ylim)[1],range(ylim)[2],1))

# Ara només en fa 10
for (i in 1:10){
  fileName <- paste("/home/lucas/ISGlobal/Arrays/Oriol_heatmaps/",gsub('/','',gsub(':','#',featureNames(x)[i])),'.png',sep='') #colon is not supported under windows
  png(fileName)
  cat(paste(i,'\n'))
  plot.new()
  par(las=2)
  plot.window(xlim=myRange,ylim=ylim)
  axis(2,at=myPoints.y)
  axis(1,at=time,lab=as.character(time))
  title(featureNames(x)[i],xlab='Time',ylab='Log ratio')
  for (j in 1:length(soques)) {
    plotMe <- exprs(x[i,which(pData(x[i,])$soca==soques[j])])
    lines(time[!is.na(plotMe)], plotMe[!is.na(plotMe)],col=j,pch=15)
    points(time[!is.na(plotMe)], plotMe[!is.na(plotMe)],pch=15,col=j)
  }
  legend(ifelse(sum(is.na(exprs(x[i,which(pData(x[i,])$time==43)])))==length(soques) | mean(exprs(x[i,which(pData(x[i,])$time==43)]),na.rm=TRUE)<mean(ylim),"topright","bottomright"),soques,pch=15,col=1:length(soques),lty=1)
  dev.off()
}


############################### Point estimations ##########################################

load(file.path(dir,'normalizedData_geneLevel.RData'))

exprsx <- exprs(xgene)
soca <- as.factor(pData(xgene)$soca)
oldTime <- as.factor(substr(sampleNames(xgene),nchar(sampleNames(xgene)),nchar(sampleNames(xgene))))

#for (i in 1:nrow(exprsx)) {
for (i in 1:10){
  cat(paste(i,'\n'))
  y <- as.numeric(exprsx[i,])
  if (sum(is.na(y))>0) {
    soca2correct <- tapply(y,soca,function(x) sum(!is.na(x))>4)
    oldTime2correct <- tapply(y,oldTime,function(x) sum(is.na(x))!=length(x))
    if (sum(soca2correct)>1) { #stop if we do not have soques with at least 5 not NA
      soca2correct <- names(soca2correct[soca2correct])
      oldTime2correct <- names(oldTime2correct[oldTime2correct])
      sel <- soca %in% soca2correct & oldTime %in% oldTime2correct
      if (sum(sel)>0) { #do not select times that are NA on all soques
        #        cat('in','\n')
        myDat <- data.frame(y=y,soca=soca,oldTime=oldTime)
        lm1 <- lm(y ~ soca + oldTime, data=myDat)
        soca2pred <- as.character(soca[sel][is.na(y[sel])])
        oldTime2pred <- as.character(oldTime[sel][is.na(y[sel])])
        y[sel][is.na(y[sel])] <- predict(lm1,newdata=data.frame(soca=soca2pred,oldTime=oldTime2pred))
        exprsx[i,] <- y
      }
    }
  }
}

estim <- is.na(exprs(xgene)) & !is.na(exprsx)
exprs(xgene) <- exprsx
save(xgene,estim,file=file.path(dir,'xgene_estimated.RData'))

############################### Plots at gene level ######

load(file.path(dir,'xgene_estimated.RData'))

myOrder <- order(pData(xgene)$soca,pData(xgene)$time)
xgene <- xgene[,myOrder]
estim <- estim[,myOrder]

time <- pData(xgene)$time
soques <- levels(pData(xgene)$soca)
ylim <- range(exprs(xgene)[!is.na(exprs(xgene))])
if (!file.exists(file.path(figuresPath,'soques_genelevel_estimated/'))) dir.create(file.path(figuresPath,'soques_genelevel_estimated/'))
outDir <- file.path(figuresPath,'soques_genelevel_estimated/')

myRange <- range(time)
myPoints.y <- round(seq(range(ylim)[1],range(ylim)[2],1))
myPoints.x <- round(seq(myRange[1],ceiling(myRange[2]),2))

#for (i in 1:nrow(xgene)) {
for (i in 1:10) {
  cat(paste(i,'\n'))
  fileName <- paste(outDir,gsub('/','',gsub(':','#',featureNames(xgene)[i])),'.png',sep='') #colon is not supported under windows
  png(fileName)
  plot.new()
  par(las=2)
  plot.window(xlim=myRange,ylim=ylim)
  axis(2,at=myPoints.y)
  axis(1,at=myPoints.x)
  #title(paste(featureNames(xgene)[i],ifelse(any(estim[i,]),'(has point estimation(s))','')),ylab='Log ratio',xlab='Time')
  title(newlist[i])
  for (j in 1:length(soques)) {
    plotMe <- exprs(xgene[i,which(pData(xgene[i,])$soca==soques[j])])
    lines(time[which(pData(xgene[i,])$soca==soques[j])][!is.na(plotMe)], plotMe[!is.na(plotMe)],col=j,pch=15)
    myPch <- ifelse(estim[i,which(pData(xgene[i,])$soca==soques[j])],1,15)
    points(time[which(pData(xgene[i,])$soca==soques[j])][!is.na(plotMe)], plotMe[!is.na(plotMe)],pch=myPch,col=j)
  }
  oldTime <- unlist(lapply(strsplit(sampleNames(xgene),','),function(x) x[length(x)]))
  legend(ifelse(sum(is.na(exprs(xgene[i,which(oldTime==43)])))==length(soques) | mean(exprs(xgene[i,which(oldTime==43)]),na.rm=TRUE)<mean(ylim),"topright","bottomright"),soques,pch=15,col=1:length(soques),lty=1)
  dev.off()
}

# Només "diferencials"

if (!file.exists(file.path(figuresPath,'soques_genelevel_nomes_diferencials/'))) dir.create(file.path(figuresPath,'soques_genelevel_nomes_diferencials/'))
outDir <- file.path(figuresPath,'soques_genelevel_nomes_diferencials/')
difs <- rownames(exprs(xgene))[abs(exprs(xgene)[,4] - exprs(xgene)[,1]) > 0.39794 | 
                                 abs(exprs(xgene)[,5] - exprs(xgene)[,2])  > 0.39794 | 
                                 abs(exprs(xgene)[,6] - exprs(xgene)[,3])  > 0.39794]

# for (i in 1:nrow(xgene)) {
#   if (rownames(xgene)[i] %in% difs){
#     cat(paste(i,'\n'))
#     fileName <- paste(outDir,gsub('/','',gsub(':','#',featureNames(xgene)[i])),'.png',sep='') #colon is not supported under windows
#     png(fileName)
#     plot.new()
#     par(las=2)
#     plot.window(xlim=myRange,ylim=ylim)
#     axis(2,at=myPoints.y)
#     axis(1,at=myPoints.x)
#     #title(paste(featureNames(xgene)[i],ifelse(any(estim[i,]),'(has point estimation(s))','')),ylab='Log ratio',xlab='Time')
#     title(newlist[i])
#     for (j in 1:length(soques)) {
#       plotMe <- exprs(xgene[i,which(pData(xgene[i,])$soca==soques[j])])
#       lines(time[which(pData(xgene[i,])$soca==soques[j])][!is.na(plotMe)], plotMe[!is.na(plotMe)],col=j,pch=15)
#       myPch <- ifelse(estim[i,which(pData(xgene[i,])$soca==soques[j])],1,15)
#       points(time[which(pData(xgene[i,])$soca==soques[j])][!is.na(plotMe)], plotMe[!is.na(plotMe)],pch=myPch,col=j)
#     }
#     oldTime <- unlist(lapply(strsplit(sampleNames(xgene),','),function(x) x[length(x)]))
#     legend(ifelse(sum(is.na(exprs(xgene[i,which(oldTime==43)])))==length(soques) | mean(exprs(xgene[i,which(oldTime==43)]),na.rm=TRUE)<mean(ylim),"topright","bottomright"),soques,pch=15,col=1:length(soques),lty=1)
#     dev.off()
#   }
# }
############################### Arees #####################

timeCorrectedEpxrs <- function(exprsx,maxMin.time.est,minMax.time.est,time,soque) {
  selClose <- function(y,x,time.est,pos='before') {
    myTime <- time.est[soque==y]
    if (any(myTime==x)) {
      ans <- rep(FALSE,length(myTime))
    } else {
      if (pos=='previous') ans <- myTime==max(myTime[myTime<x]) else ans <- myTime==min(myTime[myTime>x])
    }
    return(ans)
  }
  getNewVals <- function(x,sel1,sel2) {
    y1 <- exprsx[,sel1]
    y2 <- exprsx[,sel2]
    x1 <- time.est[sel1]
    x2 <- time.est[sel2]  
    a <- matrix((x-x1)/(x2-x1))
    b <- y2-y1
    c <- t(apply(b,1,function(x) x * a))
    ans <- c + y1
    return(ans)
  }
  sel1 <- as.logical(sapply(unique(soque),function(y) selClose(y,maxMin.time.est,time.est,'previous')))
  sel2 <- as.logical(sapply(unique(soque),function(y) selClose(y,maxMin.time.est,time.est,'next')))
  newValsBefore <- getNewVals(maxMin.time.est,sel1,sel2)  
  exprsx[,sel1] <- newValsBefore
  sel1 <- as.logical(sapply(unique(soque),function(y) selClose(y,minMax.time.est,time.est,'previous')))
  sel2 <- as.logical(sapply(unique(soque),function(y) selClose(y,minMax.time.est,time.est,'next')))
  newValsAfter <- getNewVals(minMax.time.est,sel1,sel2)  
  exprsx[,sel2] <- newValsAfter
  return(exprsx)
}

getMaxDif <- function(x) {
  if (sum(!is.na(x))>1) abs(max(x[!is.na(x)])-min(x[!is.na(x)])) else NA
}

compArea <- function(soqueName,soque,exprsx,time.est,from,to) {
  if (missing(from) & missing(to)) {
    sel <- soque == soqueName
    exprsx <- exprsx[,sel]
    x <- time.est[sel]
  } else if (!missing(from) & missing(to)) {
    uniqueTimesInSoque <- as.logical(sapply(unique(soque),function(x) !duplicated(time.est[soque==x]))) #if we have repeated times use only first
    sel <- soque == soqueName & time.est>=from & uniqueTimesInSoque
    prevTime <- c(max(time.est[soque == soqueName & time.est<from]),min(time.est[soque == soqueName & time.est>from]))
    prevSel <- soque == soqueName & time.est %in% prevTime & uniqueTimesInSoque
    firstexprs <- apply(exprsx[,prevSel],1,function(x) ifelse(any(is.na(x)),NA,approxfun(x=prevTime,y=x)(from)))
    exprsx <- cbind(firstexprs,exprsx[,sel])
    x <- c(from,time.est[sel])
  } else if (!missing(from) & !missing(to)) {
    uniqueTimesInSoque <- as.logical(sapply(unique(soque),function(x) !duplicated(time.est[soque==x]))) #if we have repeated times use only first
    sel <- soque == soqueName & (time.est>=from & time.est<=to) & uniqueTimesInSoque
    prevTime <- c(max(time.est[soque == soqueName & time.est<from]),min(time.est[soque == soqueName & time.est>from]))
    prevSel <- soque == soqueName & time.est %in% prevTime & uniqueTimesInSoque
    nextTime <- c(max(time.est[soque == soqueName & time.est<to]),min(time.est[soque == soqueName & time.est>to]))
    nextSel <- soque == soqueName & time.est %in% nextTime & uniqueTimesInSoque
    firstexprs <- apply(exprsx[,prevSel],1,function(x) ifelse(any(is.na(x)),NA,approxfun(x=prevTime,y=x)(from)))
    lastexprs <- apply(exprsx[,nextSel],1,function(x) ifelse(any(is.na(x)),NA,approxfun(x=nextTime,y=x)(to)))
    exprsx <- cbind(firstexprs,exprsx[,sel],lastexprs)
    x <- c(from,time.est[sel],to)
  } else if (missing(from) & !missing(to)) {
    uniqueTimesInSoque <- as.logical(sapply(unique(soque),function(x) !duplicated(time.est[soque==x]))) #if we have repeated times use only first
    sel <- soque == soqueName & time.est<=to & uniqueTimesInSoque
    nextTime <- c(max(time.est[soque == soqueName & time.est<to]),min(time.est[soque == soqueName & time.est>to]))
    nextSel <- soque == soqueName & time.est %in% nextTime & uniqueTimesInSoque
    lastexprs <- apply(exprsx[,nextSel],1,function(x) ifelse(any(is.na(x)),NA,approxfun(x=nextTime,y=x)(to)))        
    exprsx <- cbind(exprsx[,sel],lastexprs)
    x <- c(time.est[sel],to)
  }
  ans <- matrix(NA,nrow(exprsx))
  notNaGenes <- !apply(exprsx,1,function(x) any(is.na(x)))
  y <- exprsx[notNaGenes,]
  idx <- 1:(ncol(exprsx)-1)
  base <- x[idx+1] - x[idx]
  high.square <- sapply(idx,function(x) rowMin(y[,c(x,x+1)]))
  high.triangle <- sapply(idx,function(x) rowMax(y[,c(x,x+1)]))
  square <- rowSums(t(apply(high.square,1,function(x) x * base)))
  triangle <- high.triangle - high.square
  triangle <- rowSums(t(apply(triangle,1,function(x) x * base / 2)))
  area <- square + triangle
  ans[notNaGenes] <- area
  return(ans)
}

getArea <- function(mybreaks,soques,exprsx,time.est,soque) {
  area1 <- sapply(soques,function(x) compArea(x,soque,exprsx,time.est,to=mybreaks[3]))
  area2 <- sapply(soques,function(x) compArea(x,soque,exprsx,time.est,from=mybreaks[2],to=mybreaks[4]))
  area3 <- sapply(soques,function(x) compArea(x,soque,exprsx,time.est,from=mybreaks[3]))
  area4 <- area1 + area3 - area2
  colnames(area1) <- paste('left',colnames(area1),sep='.')
  colnames(area2) <- paste('mid',colnames(area2),sep='.')
  colnames(area3) <- paste('right',colnames(area3),sep='.')
  colnames(area4) <- paste('sides',colnames(area4),sep='.')
  area <- cbind(area1,area2,area3,area4)
  area1.maxDif <- apply(area1,1,getMaxDif)
  area2.maxDif <- apply(area2,1,getMaxDif)
  area3.maxDif <- apply(area3,1,getMaxDif)
  area4.maxDif <- apply(area4,1,getMaxDif)
  area.maxDif <- cbind(area1.maxDif,area2.maxDif,area3.maxDif,area4.maxDif)
  colnames(area.maxDif) <- cbind('area.left','area.mid','area.right','area.sides')
  ans <- list(area=area,area.maxDif=area.maxDif)
  return(ans)
}

getSampledSoque <- function(x,soque,time) {
  ans <- vector('character',length=length(soque))
  for (i in 1:length(unique(time))) {
    if (i==1) {
      ans[time==unique(time)[i]] <- as.character(sample(soque[time==unique(time[i])])) #sample soques (on each time)
    } else {
      tmpsoque <- soque[time==unique(time[i])]
      numSoques <- length(time[time==unique(time)[i]])
      prvTimeLast <- ans[time==unique(time[i-1])][numSoques]
      ans[time==unique(time)[i]][sample(numSoques-1)[1]] <- prvTimeLast
      tmpsoque <- tmpsoque[tmpsoque!=prvTimeLast]
      for (j in 1:numSoques) {
        if (ans[time==unique(time)[i]][j]=='') {
          tmpsoqueIn <- as.character(sample(tmpsoque[tmpsoque!=ans[time==unique(time)[i-1]][j]]))[1]        
          ans[time==unique(time)[i]][j] <- tmpsoqueIn
          tmpsoque <- tmpsoque[tmpsoque!=tmpsoqueIn]
        }
      }
    }
  }
  return(ans)
}

getPval <- function(soqueNames,sampledSoque,exprsx,time,time.est,maxDif,b,mybreaks) {
  myFun <- function(soqueNames,soque,exprsx,time,time.est,mybreaks) {
    myOrder <- order(soque,time.est)
    area.maxDif <- getArea(mybreaks,soqueNames,exprsx[,myOrder],time.est[myOrder],soque[myOrder])[['area.maxDif']]
    maxDif <- apply(area.maxDif,1,function(x) ifelse(any(!is.na(x)),max(x,na.rm=TRUE),NA))
    return(maxDif)
  }
  if (mc.cores==1) {
    maxDif.perm <- lapply(sampledSoque,function(x) myFun(soqueNames,x,exprsx,time,time.est,mybreaks))
  } else {
    maxDif.perm <- mclapply(sampledSoque,function(x) myFun(soqueNames,x,exprsx,time,time.est,mybreaks),mc.preschedule=TRUE,mc.cores=mc.cores)
  }
  maxDif.perm <- do.call(cbind,maxDif.perm)
  maxDif.perm <<- maxDif.perm
  ans <- sapply(1:nrow(maxDif.perm),function(i) sum(maxDif.perm[i,][!is.na(maxDif.perm[i,])] > maxDif[i]) / sum(!is.na(maxDif.perm[i,])) )
  return(ans)
}

#params
set.seed('20101028')
b <- 10000 #number of permutations

#load and order
load(file.path(dir,'xgene_estimated.RData'))
myOrder <- order(pData(xgene)$soca,pData(xgene)$time)
xgene <- xgene[,myOrder]
estim <- estim[,myOrder]
geneDesc <- fData(xgene)$Gene_name

#preprocess
soques <- levels(pData(xgene)$soca)
if (!file.exists(file.path(figuresPath,'soques_genelevel_estimated/'))) dir.create(file.path(figuresPath,'soques_genelevel_estimated/'))
outDir <- file.path(figuresPath,'soques_genelevel_estimated/')
time.est <- pData(xgene)$time
time <- as.numeric(substr(sampleNames(xgene),nchar(sampleNames(xgene)),nchar(sampleNames(xgene))))
soque <- pData(xgene)$soca
exprsx <- exprs(xgene)

#set same min and max time
maxMin.time.est <- max(sapply(unique(soque),function(x) min(time.est[soque==x])))
minMax.time.est <- min(sapply(unique(soque),function(x) max(time.est[soque==x])))
#exprsx <- timeCorrectedEpxrs(exprsx,maxMin.time.est,minMax.time.est,time,soque)
sel.previous2first <- as.logical(sapply(unique(soque),function(x) time.est[soque==x]==max(time.est[soque==x][time.est[soque==x]<=maxMin.time.est])))
sel.next2last <- as.logical(sapply(unique(soque),function(x) time.est[soque==x]==min(time.est[soque==x][time.est[soque==x]>=minMax.time.est])))
pData(xgene)$time[sel.previous2first] <- maxMin.time.est
pData(xgene)$time[sel.next2last] <- minMax.time.est
time.est <- pData(xgene)$time
#solve issue where two points are before first or after last
myFun <- function(x,gb) {
  mytime <- time.est[soque==x]
  afterLast <- mytime[mytime>=minMax.time.est]
  if (length(afterLast)>1) {
    if (gb=='good') ans <- mytime==minMax.time.est else ans <- mytime>minMax.time.est
  } else {
    ans <- rep(FALSE,length(mytime))
  }
  return(ans)
}
# good <- as.logical(sapply(unique(soque),function(x) myFun(x,'good')))
# wrong <- as.logical(sapply(unique(soque),function(x) myFun(x,'wrong')))
# exprsx[,wrong] <- exprsx[,good]
# time.est[wrong] <- time.est[good]
# 
# #set the minimum value of each gene equal to 0
# notAllNa <- apply(exprsx,1,function(x) any(!is.na(x)))
# exprsx[notAllNa,] <- t(apply(exprsx[notAllNa,],1,function(x) x - min(x,na.rm=TRUE))) #

#compute areas
mybreaks <- seq(maxMin.time.est,minMax.time.est,length.out=5)
tmp <- getArea(mybreaks,soques,exprsx,time.est,soque)
area <- tmp[['area']]; area.maxDif <- tmp[['area.maxDif']]
maxDif <- apply(area.maxDif,1,function(x) ifelse(any(!is.na(x)),max(x,na.rm=TRUE),NA))

#get pval
dummy <- as.list(1:b)
if (mc.cores==1) {
  sampledSoque <- lapply(dummy,function(x) getSampledSoque(x,soque,time))
} else {
  sampledSoque <- mclapply(dummy,function(x) getSampledSoque(x,soque,time),mc.set.seed=TRUE,mc.cores=mc.cores)
}
pvalue <- getPval(soques,sampledSoque,exprsx,time,time.est,maxDif,b,mybreaks) #takes long!!
pvalue.adj <- p.adjust(pvalue,'BH')

#export
xout <- data.frame(Name=geneDesc,area.maxDif,maxDif,pvalue,pvalue.adj,estimatedPoints=rowSums(estim))
xout <- xout[order(xout[,'maxDif'],decreasing=TRUE),]
write.csv(xout,file.path(outDir,'areas_estimated.csv'),na='')
sampledSoquesOut <- cbind(time,do.call(cbind,sampledSoque)); colnames(sampledSoquesOut)[-1] <- paste('random',1:(ncol(sampledSoquesOut)-1))
write.csv(sampledSoquesOut,file.path(outDir,'sampledSoques.csv'),row.names=FALSE)
rownames(maxDif.perm) <- featureNames(xgene)
maxDif.perm.out <- cbind(maxDif,maxDif.perm); colnames(maxDif.perm.out)[-1] <- paste('random',1:(ncol(maxDif.perm.out)-1))
write.csv(maxDif.perm.out,file.path(outDir,'maxDifPerm.csv'))
if (b>100) write.csv(maxDif.perm.out[,1:101],file.path(outDir,'maxDifPerm_100.csv'))
xout <- data.frame(Name=geneDesc,area)
rownames(xout) <- featureNames(xgene)
write.csv(xout,file.path(outDir,'areas_subclons.csv'),na='')

