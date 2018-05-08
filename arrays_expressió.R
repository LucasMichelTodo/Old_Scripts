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

roseta <- read.table("/home/lucas/ISGlobal/gens_array_rosetta.txt")

ind_1_1 <- read.table("/home/lucas/ISGlobal/Arrays/20180417 Read/US10283823_258456110002_S01_GE2_1105_Oct12_1_1.txt",sep="\t",stringsAsFactors=FALSE, skip = 9, header = TRUE)
ind_1_2 <- read.table("/home/lucas/ISGlobal/Arrays/20180417 Read/US10283823_258456110002_S01_GE2_1105_Oct12_1_2.txt",sep="\t",stringsAsFactors=FALSE, skip = 9, header = TRUE)
ind_1_3<- read.table("/home/lucas/ISGlobal/Arrays/20180417 Read/US10283823_258456110002_S01_GE2_1105_Oct12_1_3.txt",sep="\t",stringsAsFactors=FALSE, skip = 9, header = TRUE)

ctl_2_1 <- read.table("/home/lucas/ISGlobal/Arrays/20180417 Read/US10283823_258456110002_S01_GE2_1105_Oct12_2_1.txt",sep="\t",stringsAsFactors=FALSE, skip = 9, header = TRUE)
ctl_2_2 <- read.table("/home/lucas/ISGlobal/Arrays/20180417 Read/US10283823_258456110002_S01_GE2_1105_Oct12_2_2.txt",sep="\t",stringsAsFactors=FALSE, skip = 9, header = TRUE)
ctl_2_3 <- read.table("/home/lucas/ISGlobal/Arrays/20180417 Read/US10283823_258456110002_S01_GE2_1105_Oct12_2_3.txt",sep="\t",stringsAsFactors=FALSE, skip = 9, header = TRUE)

df <- cbind(ind_1_1[,c(7,8,11)], ctl_2_1[,11], ind_1_2[,11], ctl_2_2[,11], ind_1_3[,11], ctl_2_3[,11])
df["Gene_name"] <- sapply(df$SystematicName, function(x) roseta[roseta$V1 == x,2])
colnames(df) <- c("ProbeName", "SystematicName", "Ind_1", "Ctl_1", "Ind_2", "Ctl_2", "Ind_3", "Ctl_3")
# ctl <- cbind(ctl_2_1[,c(7,8,11)], ctl_2_2[,11], ctl_2_3[,11])
# ctl["Gene_name"] <- sapply(ind$SystematicName, function(x) roseta[roseta$V1 == x,2])
 
emtx <- as.matrix(df[c(-1),c(-1,-2,-9)])
minimalSet <- ExpressionSet(assayData=emtx)
pData <- read.table("/home/lucas/ISGlobal/Arrays/20180417 Read/pData_prova1.txt", header = TRUE, row.names = 1)
rownames(pData)==colnames(emtx)
summary(pData)
metadata <- data_frame(labelDescription=c("Induced/Control, status", "Time of collection: 4-4.5h/5-10h, 18-20h"), row.names=c("Induced", "Time"))
phenoData <- new("AnnotatedDataFrame", data=pData, varMetadata=metadata)
eSet <- ExpressionSet(assayData=emtx, phenoData=phenoData)
x <- eSet

#boxplot after normalization
boxplot(exprs(eSet))

#PCA (probe level)

xnm <- x[rowSums(is.na(exprs(x)))==0,]
pcdat <- prcomp(t(exprs(xnm)),na.action=na.pass)
n <- sub('\\.gpr','',rownames(pData(x)))
g <- as.numeric(factor(substring(n,nchar(n)-1,nchar(n))))

pdf(file.path(figuresPath,'pca_probeLevel.pdf'))
xlim <- 1.5*range(pcdat$x[,1]); ylim <- 1.25*range(pcdat$x[,2])
plot(pcdat$x[,1],pcdat$x[,2],col=g,xlab='PC1',ylab='PC2',xlim=xlim,ylim=ylim)
text(pcdat$x[,1],pcdat$x[,2],n,pos=1,cex=.8)



