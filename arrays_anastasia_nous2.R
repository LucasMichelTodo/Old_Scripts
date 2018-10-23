source("http://bioconductor.org/biocLite.R")
#biocLite("BiocUpgrade")
library(limma)
library(dplyr)
library(tidyr)
library(Biobase)
library(reshape2)
library(ggplot2)
library(xlsx)
library(ggfortify)

dir <- "/home/lucas/ISGlobal/Arrays/Anastasia_Arrays_Nous"
figuresPath <- "/home/lucas/ISGlobal/Arrays/Anastasia_Arrays_Nous/Plots/"

############################### Read files and create ProbeSet ##################

load(file.path(dir,'gene_list.RData'))

removeLowProbes <- function(df){
  gmedian <- median(sort(df[df$SystematicName %in% genes_list,]$gProcessedSignal)[1:100])
  rmedian <- median(sort(df[df$SystematicName %in% genes_list,]$rProcessedSignal)[1:100])
  df[df$gProcessedSignal < 3*gmedian & df$rProcessedSignal < 3*rmedian, "LogRatio"] <- NA
  return(df)
}

## Read files
#roseta <- read.csv("/home/lucas/ISGlobal/Arrays/Anastasia_Arrays_Nous/probe_names_rosetta_annotated.txt", sep = "\t", as.is = TRUE, header = FALSE)
#dim(roseta)

array_description <- read.csv("/home/lucas/ISGlobal/Arrays/Anastasia_Arrays_Nous/anastasia_array_description.csv", sep = "\t", as.is = TRUE, header = FALSE)
gene_names <- read.csv("/home/lucas/ISGlobal/Arrays/Anastasia_Arrays_Nous/anastasia_genelist_rosetta_annotated.txt", sep = "\t", as.is = TRUE, header = FALSE)

ctl_3 <- read.table("/home/lucas/ISGlobal/Arrays/Anastasia_Arrays_Nous/20180726b Read/US10283823_253249110045_S01_GE2_1105_Oct12_1_1.txt",sep="\t",stringsAsFactors=FALSE, skip = 9, header = TRUE)
ctl_3 <- removeLowProbes(ctl_3)
llcm_3 <- read.table("/home/lucas/ISGlobal/Arrays/Anastasia_Arrays_Nous/20180726b Read/US10283823_253249110045_S01_GE2_1105_Oct12_1_2.txt",sep="\t",stringsAsFactors=FALSE, skip = 9, header = TRUE)
llcm_3 <- removeLowProbes(llcm_3)
ll0_3 <- read.table("/home/lucas/ISGlobal/Arrays/Anastasia_Arrays_Nous/20180726b Read/US10283823_253249110045_S01_GE2_1105_Oct12_1_3.txt",sep="\t",stringsAsFactors=FALSE, skip = 9, header = TRUE)
ll0_3 <- removeLowProbes(ll0_3)

ctl_4 <- read.table("/home/lucas/ISGlobal/Arrays/Anastasia_Arrays_Nous/20180726b Read/US10283823_253249110045_S01_GE2_1105_Oct12_1_4.txt",sep="\t",stringsAsFactors=FALSE, skip = 9, header = TRUE)
ctl_4 <- removeLowProbes(ctl_4)
llcm_4 <- read.table("/home/lucas/ISGlobal/Arrays/Anastasia_Arrays_Nous/20180726b Read/US10283823_253249110045_S01_GE2_1105_Oct12_2_1.txt",sep="\t",stringsAsFactors=FALSE, skip = 9, header = TRUE)
llcm_4 <- removeLowProbes(llcm_4)
ll0_4 <- read.table("/home/lucas/ISGlobal/Arrays/Anastasia_Arrays_Nous/20180726b Read/US10283823_253249110045_S01_GE2_1105_Oct12_2_2.txt",sep="\t",stringsAsFactors=FALSE, skip = 9, header = TRUE)
ll0_4 <- removeLowProbes(ll0_4)

ctl_5 <- read.table("/home/lucas/ISGlobal/Arrays/Anastasia_Arrays_Nous/20180726b Read/US10283823_253249110045_S01_GE2_1105_Oct12_2_3.txt",sep="\t",stringsAsFactors=FALSE, skip = 9, header = TRUE)
ctl_5 <- removeLowProbes(ctl_5)
ll0_5 <- read.table("/home/lucas/ISGlobal/Arrays/Anastasia_Arrays_Nous/20180726b Read/US10283823_253249110045_S01_GE2_1105_Oct12_2_4.txt",sep="\t",stringsAsFactors=FALSE, skip = 9, header = TRUE)
ll0_5 <- removeLowProbes(ll0_5)


## Create probe_df
probe_df <- cbind(ctl_3[,c(7,8,11,14,15)], llcm_3[,c(11,14,15)], ll0_3[,c(11,14,15)],
            ctl_4[,c(11,14,15)], llcm_4[,c(11,14,15)], ll0_4[,c(11,14,15)],
            ctl_5[,c(11,14,15)], ll0_5[,c(11,14,15)])

probe_df["Gene_id"] <- array_description$V4
probe_df["Gene_name"] <- gene_names$V3
probe_df["Annot"] <- array_description$V5

#probe_df["Gene_id"] <- roseta$V2
#probe_df["Gene_name"] <- roseta$V4
#probe_df["Annot"] <- roseta$V3

#probe_df <- probe_df[apply(probe_df, 1, function(x) !all(is.na(x[3:11]))),]

probe_df["Annot"] <- gsub("Plasmodium", "Pl.", probe_df$Annot)
probe_df["Annot"] <- gsub("protein", "prot.", probe_df$Annot)
probe_df["Annot"] <- gsub("membrane", "memb.", probe_df$Annot)
probe_df["Annot"] <- gsub("conserved", "cvd.", probe_df$Annot)
probe_df["Annot"] <- gsub("function", "func.", probe_df$Annot)
probe_df["Annot"] <- gsub("unknown", "ukwn.", probe_df$Annot)
probe_df["Annot"] <- gsub("exported", "xptd.", probe_df$Annot)
probe_df["Annot"] <- gsub("pseudogene", "pseudo", probe_df$Annot)
probe_df["Annot"] <- gsub("putative", "put.", probe_df$Annot)
probe_df["Annot"] <- gsub("%2C", "", probe_df$Annot)

probe_df[is.na(probe_df$Gene_id),]$Gene_id <- probe_df[is.na(probe_df$Gene_id),]$SystematicName

colnames(probe_df) <- c("ProbeName", "SystematicName", 
                        "Ctl_3", "Ctl_3_Gr", "Ctl_3_Rd",
                        "LLCM_3", "LLCM_3_Gr", "LLCM_3_Rd",
                        "LL0_3", "LL0_3_Gr", "LL0_3_Rd",
                        "Ctl_4", "Ctl_4_Gr", "Ctl_4_Rd",
                        "LLCM_4", "LLCM_4_Gr", "LLCM_4_Rd",
                        "LL0_4", "LL0_4_Gr", "LL0_4_Rd",
                        "Ctl_5", "Ctl_5_Gr", "Ctl_5_Rd",
                        "LL0_5", "LL0_5_Gr", "LL0_5_Rd",
                        "Gene_id", "Gene_name", "Annot")

probe_df[probe_df$ProbeName == "HSDHFR_v7.1_P1/1",]$Gene_id <- "HSDHFR"
probe_df[probe_df$ProbeName == "HSDHFR_v7.1_P1/1",]$Gene_name <- "HSDHFR"
probe_df[probe_df$ProbeName == "HSDHFR_v7.1_P1/1",]$Annot <- "hsdhfr"

# Treure gens no únics
#no_unics <- read.table("/home/lucas/ISGlobal/Arrays/gens_no_unics.txt", sep = "\t", as.is = TRUE, header = TRUE)
#probe_df <- probe_df[!probe_df$SystematicName %in% no_unics$GeneID,]

# Treure gens no únics
no_unics <- array_description[!is.na(array_description$Status) & array_description$Status == "drop",]$ProbeName
probe_df <- probe_df[!probe_df$ProbeName %in% no_unics,]

# Afegir sondes noves
plasmid_probes <- read.csv("/home/lucas/ISGlobal/Arrays/sondes_plasmids.csv", header = TRUE, sep = "\t")

for (i in 1:dim(plasmid_probes)[1]){
  probe_df[probe_df$ProbeName == plasmid_probes[i,"NAME"],"Gene_name"] <- as.character(plasmid_probes[i,"ACCESSION_STRING"])
  probe_df[probe_df$ProbeName == plasmid_probes[i,"NAME"],"Gene_id"] <- as.character(plasmid_probes[i,"ACCESSION_STRING"])
  probe_df[probe_df$ProbeName == plasmid_probes[i,"NAME"],"Annot"] <- as.character(plasmid_probes[i,"DESCRIPTION"])
}

#probe_df[probe_df$ProbeName %in% plasmid_probes$NAME,]

# Afegir informació sobre gens variants
variant <- read.csv("/home/lucas/ISGlobal/Gen_Referencies/Gens_variants_extended.txt", header = TRUE, sep = "\t")
probe_df["Variant"] <- probe_df$Gene_id %in% variant$ID

# Afegir families i GO terms
go <- read.csv("/home/lucas/ISGlobal/Gen_Referencies/PlasmoDB-37_Pfalciparum3D7_GO.gaf", header = FALSE, sep = "\t", skip = 1)
go["gene_id"] <- gsub("EuPathDB:", "", go$V17)
go$gene_id <- gsub(".1", "", go$gene_id)
table(probe_df$Gene_id[!is.na(probe_df$Gene_id)] %in% go$gene_id)

control_probes <- c("DarkCorner", "NegativeControl", "GE_BrightCorner", "EQC", "RC", "E1A", "ETG")
no_id <- unique(probe_df[is.na(probe_df$Gene_id) & startsWith(probe_df$SystematicName, control_probes),]$SystematicName)

## Change to Log2 scale
# Log Ratio cols (originally log10)
probe_df[,seq(3,24,3)] <- log2(10**probe_df[,seq(3,24,3)])

# Processed raw signal Cols (originally no log)
probe_df[,c(4,5,7,8,10,11,13,14,16,17,19,20,22,23,25,26)] <- log2(probe_df[,c(4,5,7,8,10,11,13,14,16,17,19,20,22,23,25,26)])

write.csv(probe_df, "/home/lucas/ISGlobal/Arrays/Anastasia_Arrays_Nous/probe_probe_df_noNAs.csv", row.names = FALSE)

############################### Read files and create eSet ##################

df <- probe_df[,c(1:2,seq(3,24,3),27:29)]

## Create eSet
exprsx <- as.matrix(df[,c(-1,-2,-11,-12,-13)])
fdata <- new("AnnotatedDataFrame",df[,c(1:2,11:13)])
time <- factor(c(3,3,3,4,4,4,5,5))
soca <- factor(c("ctl", "llcm", "ll0", "ctl", "llcm", "ll0", "ctl", "ll0"))
pdata <- data.frame(soca=soca,time=time); rownames(pdata) <- colnames(df[,c(-1,-2,-11,-12,-13)])
pdata <- new("AnnotatedDataFrame",pdata)
x <- new("ExpressionSet",exprs=exprsx,featureData=fdata,phenoData=pdata)
save(x,file=file.path(dir,'normalizedData_probeLevel.RData'))


############################### Rename and Summarize ##################################
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
  
  xgene <- by(exprsx[,c(-1,-2,-11,-12,-13)],geneid,myRma) 
  xgene <- do.call('rbind',xgene)
  
  mysd <- function(x) { ans <- ifelse(sum(!is.na(x))==1,0,sd(x,na.rm=TRUE)); return(ans) }
  sdgene <- aggregate(exprsx[,c(-1,-2,-11,-12,-13)],by=list(geneid),FUN=mysd)
  
  names(sdgene)[1] <- 'geneid'
  xgene <- data.frame(geneid=rownames(xgene),xgene); rownames(xgene) <- NULL
  
  fdata <- by(df[,11:13],geneid,unique)
  genenames <- names(fdata)
  fdata <- do.call('rbind',fdata)
  
  fdata <- new("AnnotatedDataFrame",data.frame(fdata))
  rownames(fdata) <- as.character(xgene$geneid)
  
  exprsxgene <- as.matrix(xgene[,-1])
  rownames(exprsxgene) <- as.character(xgene$geneid);
  eset <- new("ExpressionSet",exprs=exprsxgene,featureData=fdata,phenoData=pdata)
  return(list(eset=eset,sdgene=sdgene,fdata=fdata,geneid=geneid))
}

geneid <- fdata@data$Gene_id
geneid <- as.character(geneid)
genesToRename.sd <- NA

tmp <- renameGenesAndSummarize(genesToRename.sd=genesToRename.sd,exprsx=df,geneid=geneid,summaryMethod=myRma)
xgene <- tmp[['eset']]; sdgene <- tmp[['sdgene']]; fdata <- tmp[['fdata']]; geneid <- tmp[['geneid']]

write.csv(data.frame(probe=df[,1],gene=geneid),file.path(figuresPath,'geneProbes.csv'),row.names=FALSE) #probes of each gene

geneid2 <- unique(probe_df$Gene_id)
geneid2 <- as.character(geneid)
genesToRename.sd2 <- NA

probe_df <- probe_df[apply(probe_df, 1, function(x) !all(is.na(x[3:11]))),]
df2 <- probe_df[,c(1:2, seq(5,26,3), 27:29)]

renameGenesAndSummarize2 <- function(genesToRename.sd2,exprsx2,geneid2,summaryMethod=myRma) {

  red_gene <- by(exprsx2[,c(3:10)],df2["Gene_id"],myRma) 
  red_gene <- do.call('rbind',red_gene)
  
  mysd <- function(x) { ans <- ifelse(sum(!is.na(x))==1,0,sd(x,na.rm=TRUE)); return(ans) }
  sdgene <- aggregate(exprsx2[,c(3:10)],by=exprsx2["Gene_id"],FUN=mysd)
  
  names(sdgene)[1] <- 'geneid2'
  red_gene <- data.frame(geneid2=rownames(red_gene),red_gene); rownames(red_gene) <- NULL
  
  fdata <- by(exprsx2[,11:13],exprsx2["Gene_id"],unique)
  genenames <- names(fdata)
  fdata <- do.call('rbind',fdata)
  
  fdata <- new("AnnotatedDataFrame",data.frame(fdata))
  rownames(fdata) <- as.character(red_gene$geneid2)
  
  exprsred_gene <- as.matrix(red_gene[,-1])
  rownames(exprsred_gene) <- red_gene$geneid2;
  colnames(exprsred_gene) <- rownames(pdata@data)
  
  eset <- new("ExpressionSet",exprs=exprsred_gene,featureData=fdata,phenoData=pdata)
  return(list(eset=eset,sdgene=sdgene,fdata=fdata,geneid2=geneid2))
}
tmp2 <- renameGenesAndSummarize2(genesToRename.sd=genesToRename.sd,exprsx=df2,geneid2=geneid2,summaryMethod=myRma)
red_gene <- tmp2[['eset']]; sdgene <- tmp2[['sdgene']]; fdata <- tmp2[['fdata']]; geneid <- tmp2[['geneid']]


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
  oldTime <- as.numeric(time)
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
save(estimatedTimes,file=file.path(dir,'estimatedTimes.RData'))
load(file.path(dir,'estimatedTimes.RData'))
write.csv(estimatedTimes,file.path(figuresPath,'Estimated_Times.csv'))
pData(xgene)$time <- estimatedTimes

#save ExpressionSet at gene level
save(xgene,file=file.path(dir,'normalizedData_geneLevel.RData'))
#save(xgene.noRatio,file=file.path(dir,'normalizedData_geneLevel_noRatio.RData'))

#boxplot after summarization
pdf(file.path(figuresPath,'boxplot_afterSummarization.pdf'))
boxplot(exprs(xgene),main='summarization method: median poslish')
dev.off()

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

############################### PCAs Nous ###########################

pca_df <- t(exprs(x)[complete.cases(exprs(x)),])
autoplot(prcomp(pca_df), data = x@phenoData@data, colour = "time")
autoplot(prcomp(pca_df[3:5,]), data = x@phenoData@data[3:5,], colour = "soca") + ggtitle("Timepoint 2")
autoplot(prcomp(pca_df[6:8,]), data = x@phenoData@data[6:8,], colour = "soca") + ggtitle("Timepoint 3")
autoplot(prcomp(pca_df[9:11,]), data = x@phenoData@data[9:11,], colour = "soca") + ggtitle("Timepoint 4")
autoplot(prcomp(pca_df[12:14,]), data = x@phenoData@data[12:14,], colour = "soca") + ggtitle("Timepoint 5")

############################### Array Plots ######################

## Red
rd_df <- as.data.frame(cbind(probe_df$ProbeName, probe_df$SystematicName, probe_df$Gene_id,
                             probe_df$Ctl_3_Rd, probe_df$LLCM_3_Rd, probe_df$LL0_3_Rd,
                             probe_df$Ctl_4_Rd, probe_df$LLCM_4_Rd, probe_df$LL0_4_Rd,
                             probe_df$Ctl_5_Rd, probe_df$LL0_5_Rd))

colnames(rd_df) <- c(c("ProbeName", "SystematicName", "Gene_id",
                       "Ctl_3", "LLcm_3", "LL0_3",
                       "Ctl_4", "LLcm_4", "LL0_4",
                       "Ctl_5", "LL0_5"))

rd_df[,c(-1,-2,-3)] <- apply(rd_df[,c(-1,-2,-3)], 2, function(x) as.numeric(x))
boxplot(rd_df[,c(-1,-2,-3)], main="Red")

rd_df["X_row"] <- ctl_3[ctl_3$ProbeName %in% probe_df$ProbeName,]$Row
rd_df["Y_col"] <- ctl_3[ctl_3$ProbeName %in% probe_df$ProbeName,]$Col

rd_df_m <- melt(rd_df, id.vars = c(1:3, 12:13))
rd_df_m["Time"] <- unlist(lapply(rd_df_m$variable, function(x) substr(as.character(x), nchar(as.character(x)), nchar(as.character(x)))))
rd_df_m$variable <- as.factor(rd_df_m$variable)
p <- ggplot(rd_df_m, aes(x=variable, y=value, fill=Time)) +
  geom_boxplot() 
p

library(RColorBrewer)
cols <- rev(brewer.pal(11, 'Spectral'))

arrayPlot <- function(df) {
  # Get df name as string
  df_name <- deparse(substitute(df))
  
  p1 <- qplot(Col, Row, data=df, color=log2(rMedianSignal)<7) + scale_color_manual(values=c("aliceblue", "black")) + ggtitle(df_name)
  ggsave(p1, filename = paste0(figuresPath, "Array_Plots/", df_name, "_boolean.jpeg"), device = "jpeg")
  print(p1)
  
  p2 <- qplot(Col, Row, data=df, color=log2(rMedianSignal)) + scale_colour_gradientn(colours = cols) + ggtitle(df_name)
  ggsave(p2, filename = paste0(figuresPath, "Array_Plots/", df_name, ".jpeg"), device = "jpeg")
  print(p2)
  
  p3 <- qplot(Col, Row, data=df, color=is.na(LogRatio)) + scale_color_manual(values=c("aliceblue", "red")) + ggtitle(df_name)
  ggsave(p3, filename = paste0(figuresPath, "Array_Plots/", df_name, "_NAs.jpeg"), device = "jpeg")
  print(p3)
}

# arrayPlot(ctl_3)
# arrayPlot(ctl_4)
# arrayPlot(ctl_5)
#
# arrayPlot(ll0_3)
# arrayPlot(ll0_4)
# arrayPlot(ll0_5)
# 
# arrayPlot(llcm_3)
# arrayPlot(llcm_4)


##Green
gr_df <- as.data.frame(cbind(probe_df$ProbeName, probe_df$SystematicName, probe_df$Gene_id,
                             probe_df$Ctl_3_Gr, probe_df$LLCM_3_Gr, probe_df$LL0_3_Gr,
                             probe_df$Ctl_4_Gr, probe_df$LLCM_4_Gr, probe_df$LL0_4_Gr,
                             probe_df$Ctl_5_Gr, probe_df$LL0_5_Gr))

colnames(gr_df) <- c(c("ProbeName", "SystematicName", "Gene_id", 
                       "Ctl_3", "LLcm_3", "LL0_3",
                       "Ctl_4", "LLcm_4", "LL0_4",
                       "Ctl_5", "LL0_5"))

gr_df[,c(-1,-2,-3)] <- apply(gr_df[,c(-1,-2,-3)], 2, function(x) as.numeric(x))
boxplot(gr_df[,c(-1,-2,-3)], main="Green")

gr_df_m <- melt(gr_df, id.vars = c(1:3))
gr_df_m$variable <- as.factor(gr_df_m$variable)
p <- ggplot(gr_df_m, aes(x=variable, y=value)) +
  geom_boxplot() 
p

############################### Info df ###############################

info_df <- rd_df[,c(1:3,12:13)]
info_df["Gene_name"] <- probe_df$Gene_name
info_df["Annot"] <- probe_df$Annot
info_df["VAriant"] <- probe_df$Variant
#write.xlsx(info_df, file = "/home/lucas/ISGlobal/Arrays/Anastasia_Arrays/info_df.xlsx", sheetName = "Sheet1", row.names = FALSE)

############################### Red signal by Probe #########################

# for (i in 1:nrow(probe_df)){
# #for (i in 1:10){
#   print(i)
#   graf <- melt(probe_df[i,c(1:2, seq(5,26,3), 27:29)])
#   graf["Time"] <- estimatedTimes
#   graf["Soca"] <- c(rep(c("Ctl", "LLCM", "LL0"), 2), "Ctl", "LL0")
#   p <- ggplot(graf, aes(x = Time, y = value, group = Soca, color = Soca))
#   p <- p + geom_point() + geom_line() + coord_cartesian(ylim = c(1.5, 20)) + ggtitle(probe_df$ProbeName[i])
#   ggsave(p, file=paste0(figuresPath,"RedProcessed/Probes_RedProcessed/", as.character(i), "_", gsub("/", ":", probe_df$Gene_id[i]), ".jpeg"), device = "jpeg", width = 14, height = 10, units = "cm")
#   #print(p)
# }

#could not open file '/home/lucas/ISGlobal/Arrays/Anastasia_Arrays/Plots/Probes_red/S2-type_5.8S_p_MAL11/13_5.8S:rRNA.jpeg'
#could not open file '/home/lucas/ISGlobal/Arrays/Anastasia_Arrays/Plots/Probes_red/A-type_5.8S_p_MAL5/7_5.8S:rRNA.jpeg'

############################### Probes by Gene ####################################

# ### Ctl probes by gene
# counter <- 1
# for (i in unique(probe_df$Gene_id)){
#   print(counter)
#   graf <- melt(probe_df[probe_df$Gene_id == i,c(1:2,3,12,21,27:30)])
#   p <- ggplot(graf, aes(x = variable, y = value, group = ProbeName, color = ProbeName))
#   p <- p + geom_point() + geom_line() + coord_cartesian(ylim = c(-10, 6)) 
#   ggsave(p, file=paste0(figuresPath,"Ratio/Probes_Ctl/", gsub("/", ":", i), ".jpeg"), device = "jpeg", width = 14, height = 10, units = "cm")
#   #print(p)
#   counter <- counter+1
# }
# 
# ### LL_0 probes by gene
# counter <- 1
# for (i in unique(probe_df$Gene_id)){
#   print(counter)
#   graf <- melt(probe_df[probe_df$Gene_id == i,c(1:2,9,18,24,27:30)])
#   p <- ggplot(graf, aes(x = variable, y = value, group = ProbeName, color = ProbeName))
#   p <- p + geom_point() + geom_line() + coord_cartesian(ylim = c(-10, 6)) 
#   ggsave(p, file=paste0(figuresPath,"Ratio/Probes_LL0/", gsub("/", ":", i), ".jpeg"), device = "jpeg", width = 14, height = 10, units = "cm")
#   #print(p)
#   counter <- counter+1
# }
# 
# ### LL_CM probes by gene
# counter <- 1
# for (i in unique(probe_df$Gene_id)){
#   print(counter)
#   graf <- melt(probe_df[probe_df$Gene_id == i,c(1:2,6,15,27:30)])
#   p <- ggplot(graf, aes(x = variable, y = value, group = ProbeName, color = ProbeName))
#   p <- p + geom_point() + geom_line() + coord_cartesian(ylim = c(-10, 6)) 
#   ggsave(p, file=paste0(figuresPath,"Ratio/Probes_LLCM/", gsub("/", ":", i), ".jpeg"), device = "jpeg", width = 14, height = 10, units = "cm")
#   #print(p)
#   counter <- counter+1
# }

############################### Probes by Gene: Red processed signal ####################################

# ### Ctl probes by gene
# counter <- 1
# for (i in unique(probe_df$Gene_id)){
#   print(counter)
#   graf <- melt(probe_df[probe_df$Gene_id == i,c(1:2,5,14,23,27:30)])
#   p <- ggplot(graf, aes(x = variable, y = value, group = ProbeName, color = ProbeName))
#   p <- p + geom_point() + geom_line() + coord_cartesian(ylim = c(1.5, 20))
#   ggsave(p, file=paste0(figuresPath,"RedProcessed/Probes_Ctl_RedProcessed/", gsub("/", ":", i), ".jpeg"), device = "jpeg", width = 14, height = 10, units = "cm")
#   #print(p)
#   counter <- counter+1
# }
# 
# ### LL_0 probes by gene
# counter <- 1
# for (i in unique(probe_df$Gene_id)){
#   print(counter)
#   graf <- melt(probe_df[probe_df$Gene_id == i,c(1:2,11,20,26,27:30)])
#   p <- ggplot(graf, aes(x = variable, y = value, group = ProbeName, color = ProbeName))
#   p <- p + geom_point() + geom_line() + coord_cartesian(ylim = c(1.5, 20))
#   ggsave(p, file=paste0(figuresPath,"RedProcessed/Probes_LL0_RedProcessed/", gsub("/", ":", i), ".jpeg"), device = "jpeg", width = 14, height = 10, units = "cm")
#   #print(p)
#   counter <- counter+1
# }
# 
# ### LL_CM probes by gene
# counter <- 1
# for (i in unique(probe_df$Gene_id)){
#   print(counter)
#   graf <- melt(probe_df[probe_df$Gene_id == i,c(1:2,8,17,27:30)])
#   p <- ggplot(graf, aes(x = variable, y = value, group = ProbeName, color = ProbeName))
#   p <- p + geom_point() + geom_line() + coord_cartesian(ylim = c(1.5, 20))
#   ggsave(p, file=paste0(figuresPath,"RedProcessed/Probes_LLCM_RedProcessed/", gsub("/", ":", i), ".jpeg"), device = "jpeg", width = 14, height = 10, units = "cm")
#   #print(p)
#   counter <- counter+1
# }

############################### Plots at probe level ##########################################

# load(file.path(dir,'normalizedData_probeLevel.RData'))
# 
# for (i in 1:nrow(x)){
# #for (i in 1:10){
#   print(i)
#   graf <- melt(exprs(x)[i,])
#   graf["Soca"] <- x@phenoData@data$soca
#   graf["Time"] <- estimatedTimes
#   p <- ggplot(graf, aes(x = Time, y = value, group = Soca, color = Soca))
#   p <- p + geom_point() + geom_line() + coord_cartesian(ylim = c(-10, 6)) + ggtitle(x@featureData@data$ProbeName[i])
#   ggsave(p, file=paste0(figuresPath,"Ratio/Probes/", as.character(i), "_", gsub("/", ":", x@featureData@data$Gene_id[i]), ".jpeg"), device = "jpeg", width = 14, height = 10, units = "cm")
#   #print(p)
# }

############################### Point estimations ##########################################

load(file.path(dir,'normalizedData_geneLevel.RData'))

exprsx <- exprs(xgene)
soca <- as.factor(pData(xgene)$soca)
oldTime <- as.factor(substr(sampleNames(xgene),nchar(sampleNames(xgene)),nchar(sampleNames(xgene))))

for (i in 1:nrow(exprsx)) {
#for (i in 1:10){
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

# #### Ratio Signal
# 
# for (i in 1:nrow(xgene)){
# #for (i in 1:100){
#   print(i)
#   graf <- melt(exprs(xgene)[i,])
#   graf["Soca"] <- xgene@phenoData@data$soca
#   graf["Time"] <- estimatedTimes
#   p <- ggplot(graf, aes(x = Time, y = value, group = Soca, color = Soca))
#   p <- p + geom_point() + geom_line() + coord_cartesian(ylim = c(-10, 6)) + ggtitle(xgene@featureData@data$Gene_id[i])
#   ggsave(p, file=paste0(figuresPath,"Ratio/Gens/", gsub("/", ":", xgene@featureData@data$Gene_id[i]), ".jpeg"), device = "jpeg", width = 14, height = 10, units = "cm")
#   #print(p)
# }
# 
# #### Red Signal
# 
# for (i in 1:nrow(red_gene)){
# #for (i in 1:10){
#   print(i)
#   graf <- melt(exprs(red_gene)[i,])
#   graf["Soca"] <- red_gene@phenoData@data$soca
#   graf["Time"] <- estimatedTimes
#   p <- ggplot(graf, aes(x = Time, y = value, group = Soca, color = Soca))
#   p <- p + geom_point() + geom_line() + coord_cartesian(ylim = c(2, 22)) + ggtitle(red_gene@featureData@data$Gene_id[i])
#   ggsave(p, file=paste0(figuresPath,"RedProcessed/Gens_RedProcessed/", gsub("/", ":", red_gene@featureData@data$Gene_id[i]), ".jpeg"), device = "jpeg", width = 14, height = 10, units = "cm")
#   #print(p)
# }

############################### Fold Changes and stuff ################################

## Probe Level
head(probe_df)
fc_probes <- data.frame(c(probe_df$LLCM_3 - probe_df$Ctl_3),
                        c(probe_df$LL0_3 - probe_df$Ctl_3),
                        c(probe_df$LLCM_3 - probe_df$LL0_3),
                        c(probe_df$LLCM_4 - probe_df$Ctl_4),
                        c(probe_df$LL0_4 - probe_df$Ctl_4),
                        c(probe_df$LLCM_4 - probe_df$LL0_4),
                        c(probe_df$LL0_5 - probe_df$Ctl_5))
            
colnames(fc_probes) <- c("LLCM_3", "LL0_3", "LLcmvsLL03", "LLCM_4", "LL0_4", "LLcmvsLL04", "LL0_5")
fc_probes <- cbind(fc_probes, probe_df[,c(1:2,27:30)])

fcp_llcm_3 <- arrange(fc_probes, -abs(LLCM_3))[1:100,c(1,8:13)]
fcp_llcm_4 <- arrange(fc_probes, -abs(LLCM_4))[1:100,c(4,8:13)]
fcp_ll0_3 <- arrange(fc_probes, -abs(LL0_3))[1:100,c(2,8:13)]
fcp_ll0_4 <- arrange(fc_probes, -abs(LL0_4))[1:100,c(5,8:13)]
fcp_ll0_5 <- arrange(fc_probes, -abs(LL0_5))[1:100,c(7,8:13)]
fcp_llcm_ll0_3 <- arrange(fc_probes, -abs(LLcmvsLL03))[1:100,c(3,8:13)]
fcp_llcm_ll0_4 <- arrange(fc_probes, -abs(LLcmvsLL04))[1:100,c(6,8:13)]

write.xlsx(fcp_llcm_3, file = "/home/lucas/ISGlobal/Arrays/Anastasia_Arrays_Nous/FC/Top_fc_probes_llcm_3.xlsx", sheetName = "Sheet1", row.names = FALSE)
write.xlsx(fcp_llcm_4, file = "/home/lucas/ISGlobal/Arrays/Anastasia_Arrays_Nous/FC/Top_fc_probes_llcm_4.xlsx", sheetName = "Sheet1", row.names = FALSE)
write.xlsx(fcp_ll0_3, file = "/home/lucas/ISGlobal/Arrays/Anastasia_Arrays_Nous/FC/Top_fc_probes_ll0_3.xlsx", sheetName = "Sheet1", row.names = FALSE)
write.xlsx(fcp_ll0_4, file = "/home/lucas/ISGlobal/Arrays/Anastasia_Arrays_Nous/FC/Top_fc_probes_ll0_4.xlsx", sheetName = "Sheet1", row.names = FALSE)
write.xlsx(fcp_ll0_5, file = "/home/lucas/ISGlobal/Arrays/Anastasia_Arrays_Nous/FC/Top_fc_probes_ll0_5.xlsx", sheetName = "Sheet1", row.names = FALSE)
write.xlsx(fcp_llcm_ll0_3, file = "/home/lucas/ISGlobal/Arrays/Anastasia_Arrays_Nous/FC/Top_fc_probes_llcm_ll0_3.xlsx", sheetName = "Sheet1", row.names = FALSE)
write.xlsx(fcp_llcm_ll0_4, file = "/home/lucas/ISGlobal/Arrays/Anastasia_Arrays_Nous/FC/Top_fc_probes_llcm3_ll0_4.xlsx", sheetName = "Sheet1", row.names = FALSE)


## Gene Level
load(file.path(dir,'xgene_estimated.RData'))
head(exprs(xgene))
fc_df <- data.frame(c(exprs(xgene)[,"LLCM_3"]-exprs(xgene)[,"Ctl_3"]),
                      exprs(xgene)[,"LL0_3"]-exprs(xgene)[,"Ctl_3"],
                      exprs(xgene)[,"LLCM_3"]-exprs(xgene)[,"LL0_3"],
                      exprs(xgene)[,"LLCM_4"]-exprs(xgene)[,"Ctl_4"],
                      exprs(xgene)[,"LL0_4"]-exprs(xgene)[,"Ctl_4"],
                      exprs(xgene)[,"LLCM_4"]-exprs(xgene)[,"LL0_4"],
                      exprs(xgene)[,"LL0_5"]-exprs(xgene)[,"Ctl_5"]) 


colnames(fc_df) <- c("LLCM_3", "LL0_3", "LLCMvsLL0_3", "LLCM_4", "LL0_4", "LLCMvsLL0_4", "LL0_5")
fc_df <- cbind(fc_df, xgene@featureData@data[,c(1:3)])
hist(as.matrix(fc_df[,c(1:7)]), breaks = 200)

melted <- melt(fc_df)
top100 <- arrange(melted, -abs(value))[1:300,]
write.xlsx(top100, file = "/home/lucas/ISGlobal/Arrays/Anastasia_Arrays_Nous/FC/top300_gene_alltimes.xlsx", sheetName = "Sheet1", row.names = FALSE)

fc_genes_llcm_3 <- arrange(fc_df, -abs(LLCM_3))[1:100,c(1,8:10)]
fc_genes_llcm_4 <- arrange(fc_df, -abs(LLCM_4))[1:100,c(4,8:10)]
fc_genes_ll0_3 <- arrange(fc_df, -abs(LL0_3))[1:100,c(2,8:10)]
fc_genes_ll0_4 <- arrange(fc_df, -abs(LL0_4))[1:100,c(5,8:10)]
fc_genes_ll0_5 <- arrange(fc_df, -abs(LL0_5))[1:100,c(7,8:10)]
fc_genes_llcm_ll0_3 <- arrange(fc_df, -abs(LLCMvsLL0_3))[1:100,c(3,8:10)]
fc_genes_llcm_ll0_4 <- arrange(fc_df, -abs(LLCMvsLL0_4))[1:100,c(6,8:10)]

write.xlsx(fc_genes_llcm_3, file = "/home/lucas/ISGlobal/Arrays/Anastasia_Arrays_Nous/FC/Top_fc_genes_llcm_3.xlsx", sheetName = "Sheet1", row.names = FALSE)
write.xlsx(fc_genes_llcm_4, file = "/home/lucas/ISGlobal/Arrays/Anastasia_Arrays_Nous/FC/Top_fc_genes_llcm_4.xlsx", sheetName = "Sheet1", row.names = FALSE)
write.xlsx(fc_genes_ll0_3, file = "/home/lucas/ISGlobal/Arrays/Anastasia_Arrays_Nous/FC/Top_fc_genes_ll0_3.xlsx", sheetName = "Sheet1", row.names = FALSE)
write.xlsx(fc_genes_ll0_4, file = "/home/lucas/ISGlobal/Arrays/Anastasia_Arrays_Nous/FC/Top_fc_genes_ll0_4.xlsx", sheetName = "Sheet1", row.names = FALSE)
write.xlsx(fc_genes_ll0_5, file = "/home/lucas/ISGlobal/Arrays/Anastasia_Arrays_Nous/FC/Top_fc_genes_ll0_5.xlsx", sheetName = "Sheet1", row.names = FALSE)
write.xlsx(fc_genes_llcm_ll0_3, file = "/home/lucas/ISGlobal/Arrays/Anastasia_Arrays_Nous/FC/Top_fc_genes_llcm_ll0_3.xlsx", sheetName = "Sheet1", row.names = FALSE)
write.xlsx(fc_genes_llcm_ll0_4, file = "/home/lucas/ISGlobal/Arrays/Anastasia_Arrays_Nous/FC/Top_fc_genes_llcm3_ll0_4.xlsx", sheetName = "Sheet1", row.names = FALSE)

# for (i in 1:nrow(top100)){
#   print(i)
#   graf <- melt(exprs(xgene)[rownames(exprs(xgene)) == top100$Gene_id[i],])
#   graf["Time"] <- estimatedTimes
#   graf["Soque"] <- xgene$soca
#   p <- ggplot(graf, aes(x = Time, y = value, col = Soque))
#   p <- p + geom_point(aes(color = Soque, shape = Soque)) + geom_line() + coord_cartesian(ylim = c(-10, 6)) + ggtitle(top100$Gene_id[i])
#   ggsave(p, file=paste0(dir,"/FC/Plots/",gsub('/','',gsub(':','#',top100$Gene_id[i])),'.jpeg',sep=''), device = "jpeg", width = 14, height = 10, units = "cm")
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
  if (length(idx) == 1){
    square <- apply(high.square,1,function(x) x * base)
  } else {
    square <- rowSums(t(apply(high.square,1,function(x) x * base)))
  }
  triangle <- high.triangle - high.square
  if (length(idx) == 1){
    triangle <- apply(triangle,1,function(x) x * base / 2)
  } else {
    triangle <- rowSums(t(apply(triangle,1,function(x) x * base / 2)))
  }
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
        if (ans[time==unique(time)[i]][j]=="") {
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

xgene <- xgene[,c(-2,-5)]

myOrder <- order(pData(xgene)$soca,pData(xgene)$time)
xgene <- xgene[,myOrder]
estim <- estim[,myOrder]
geneDesc <- fData(xgene)$Gene_name

#preprocess
soques <- levels(droplevels(pData(xgene)$soca))
if (!file.exists(file.path(figuresPath,'soques_genelevel_estimated/'))) dir.create(file.path(figuresPath,'soques_genelevel_estimated/'))
outDir <- file.path(figuresPath,'soques_genelevel_estimated/')
time.est <- as.numeric(pData(xgene)$time)
time <- as.numeric(c(3,4,5,3,4,5)) ## Aquesta linia canvia d'script a script (hem tret LLCM!)
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
good <- as.logical(sapply(unique(soque),function(x) myFun(x,'good')))
wrong <- as.logical(sapply(unique(soque),function(x) myFun(x,'wrong')))
exprsx[,wrong] <- exprsx[,good]
time.est[wrong] <- time.est[good]

#set the minimum value of each gene equal to 0
notAllNa <- apply(exprsx,1,function(x) any(!is.na(x)))
exprsx[notAllNa,] <- t(apply(exprsx[notAllNa,],1,function(x) x - min(x,na.rm=TRUE))) #

#compute areas
mybreaks <- seq(maxMin.time.est,minMax.time.est,length.out=5)
tmp <- getArea(mybreaks,soques,exprsx,time.est,soque)
area <- tmp[['area']]; area.maxDif <- tmp[['area.maxDif']]
maxDif <- apply(area.maxDif,1,function(x) ifelse(any(!is.na(x)),max(x,na.rm=TRUE),NA))

#get pval
# mc.cores <- 4
# dummy <- as.list(1:b)
# if (mc.cores==1) {
#   sampledSoque <- lapply(dummy,function(x) getSampledSoque(x,soque,time))
# } else {
#   sampledSoque <- mclapply(dummy,function(x) getSampledSoque(x,soque,time),mc.set.seed=TRUE,mc.cores=mc.cores)
# }
# 
# start_time <- Sys.time()
# pvalue <- getPval(soques,sampledSoque,exprsx,time,time.est,maxDif,b,mybreaks) #takes long!!
# end_time <- Sys.time()
# end_time - start_time
# 
# pvalue.adj <- p.adjust(pvalue,'BH')

#export
# xout <- data.frame(Name=geneDesc,area.maxDif,maxDif,pvalue,pvalue.adj,estimatedPoints=rowSums(estim))
# xout <- xout[order(xout[,'maxDif'],decreasing=TRUE),]
xout <- data.frame(Name=geneDesc,area.maxDif,maxDif,estimatedPoints=rowSums(estim))
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


############################### Arees without LLCM

############################### Seleccions Arees ###############################


load(file.path(dir,'xgene_estimated.RData'))

area_df <- as.data.frame(area)
area_df["Gene"] <- xgene@featureData$Gene_id

area_FC <- data.frame(c(area_df[,2]-area_df[,1]),
                      c(area_df[,3]-area_df[,1]),
                      c(area_df[,3]-area_df[,2]),
                      c(area_df[,5]-area_df[,4]),
                      c(area_df[,6]-area_df[,4]),
                      c(area_df[,6]-area_df[,5]),
                      c(area_df[,8]-area_df[,7]),
                      c(area_df[,9]-area_df[,7]),
                      c(area_df[,9]-area_df[,8]),
                      c(area_df[,"Gene"]))

colnames(area_FC) <- c("LL_0_left", "LL_CM_left", "CM_0_left", "LL_0_mid", "LL_CM_mid", "CM_0_mid", "LL_0_right", "LL_CM_right", "CM_0_right", "Gene")

marea_left <- melt(area_FC[,c(1,2,3,10)])

top50_areas_left <- arrange(marea_left, -abs(value))[1:50,]
write.xlsx(top50_areas_left, file = paste0(figuresPath, "/Top50areas_left/top_area_left.xlsx"))

for (i in 1:nrow(top50_areas_left)){
#for (i in 1:2){
  print(i)
  graf <- melt(exprs(xgene)[rownames(exprs(xgene)) == top50_areas_left$Gene[i],])
  graf["Time"] <- estimatedTimes
  graf["Soque"] <- xgene$soca
  p <- ggplot(graf, aes(x = Time, y = value, col = Soque))
  p <- p + geom_point(aes(color = Soque, shape = Soque)) + geom_line()
  p <- p + ggtitle(top50_areas_left$Gene[i])
  #print(p)
  ggsave(p, file=paste(figuresPath,"/Top50areas_left/", as.character(i), "_", gsub('/','',gsub(':','#',top50_areas_left$Gene[i])),'.jpeg',sep=''), device = "jpeg", width = 14, height = 10, units = "cm")
}

marea_mid <- melt(area_FC[,c(4,5,6,10)])

top50_areas_mid <- arrange(marea_mid, -abs(value))[1:50,]
write.xlsx(top50_areas_mid, file = paste0(figuresPath, "/Top50areas_mid/top_area_mid.xlsx"))

for (i in 1:nrow(top50_areas_mid)){
  #for (i in 1:2){
  print(i)
  graf <- melt(exprs(xgene)[rownames(exprs(xgene)) == top50_areas_mid$Gene[i],])
  graf["Time"] <- estimatedTimes
  graf["Soque"] <- xgene$soca
  p <- ggplot(graf, aes(x = Time, y = value, col = Soque))
  p <- p + geom_point(aes(color = Soque, shape = Soque)) + geom_line()
  p <- p + ggtitle(top50_areas_mid$Gene[i])
  #print(p)
  ggsave(p, file=paste(figuresPath,"/Top50areas_mid/", as.character(i), "_", gsub('/','',gsub(':','#',top50_areas_mid$Gene[i])),'.jpeg',sep=''), device = "jpeg", width = 14, height = 10, units = "cm")
}

marea_right <- melt(area_FC[,c(7,8,9,10)])

top50_areas_right <- arrange(marea_right, -abs(value))[1:50,]
write.xlsx(top50_areas_right, file = paste0(figuresPath, "/Top50areas_right/top_area_right.xlsx"))

for (i in 1:nrow(top50_areas_right)){
  #for (i in 1:2){
  print(i)
  graf <- melt(exprs(xgene)[rownames(exprs(xgene)) == top50_areas_right$Gene[i],])
  graf["Time"] <- estimatedTimes
  graf["Soque"] <- xgene$soca
  p <- ggplot(graf, aes(x = Time, y = value, col = Soque))
  p <- p + geom_point(aes(color = Soque, shape = Soque)) + geom_line()
  p <- p + ggtitle(top50_areas_right$Gene[i])
  #print(p)
  ggsave(p, file=paste(figuresPath,"/Top50areas_right/", as.character(i), "_", gsub('/','',gsub(':','#',top50_areas_right$Gene[i])),'.jpeg',sep=''), device = "jpeg", width = 14, height = 10, units = "cm")
}


############################### Lipid GO ##############

lipid_go <- read.xlsx(file = paste0(dir,"/Lipid_GO.xlsx"), header = TRUE, sheetIndex = 1)

lipid_go$Gene.ID 
table(top50_areas_left$Gene %in% lipid_go$Gene.ID)
table(top50_areas_mid$Gene %in% lipid_go$Gene.ID)
table(top50_areas_right$Gene %in% lipid_go$Gene.ID)

top50_areas_right[top50_areas_right$Gene %in% lipid_go$Gene.ID,]

table(top50$Gene %in% lipid_go$Gene.ID)

load(file.path(dir,'normalizedData_geneLevel.RData'))
x <- xgene
xnm <- xgene[rowSums(is.na(exprs(xgene)))==0,]
pcdat <- prcomp(exprs(xnm)[rownames(exprs(xnm)) %in% lipid_go$Gene.ID,])
n <- sub('\\.gpr','',rownames(pData(x)))
g <- as.numeric(factor(substring(n,nchar(n)-1,nchar(n))))

pdf(file.path(figuresPath,'pca_lipidGO.pdf'))
xlim <- 1.5*range(pcdat$x[,1]); ylim <- 1.25*range(pcdat$x[,2])
plot(pcdat$x[,1],pcdat$x[,2],col=g,xlab='PC1',ylab='PC2',xlim=xlim,ylim=ylim)
text(pcdat$x[,1],pcdat$x[,2],n,pos=1,cex=.8)
dev.off()



