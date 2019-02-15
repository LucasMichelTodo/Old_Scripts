source("http://bioconductor.org/biocLite.R")
library(limma)
library(dplyr)
library(tidyr)
library(Biobase)
library(reshape2)
library(ggplot2)
library(gplots)

dir <- "/home/lucas/ISGlobal/Arrays/Eli_Arrays"
figuresPath <- "/home/lucas/ISGlobal/Arrays/Eli_Arrays/Plots/"

############################### Read files and create eSet ##################

load(file.path(dir,'gene_list.RData'))

# No-log NA remove -> 3 Median
removeLowProbes <- function(df){
  gmedian <- median(sort(df[df$SystematicName %in% genes_list,]$gProcessedSignal)[1:100])
  rmedian <- median(sort(df[df$SystematicName %in% genes_list,]$rProcessedSignal)[1:100])
  df[df$gProcessedSignal < 3*gmedian & df$rProcessedSignal < 3*rmedian, "LogRatio"] <- NA
  return(df)
}

## Import array data
array_description <- read.csv("/home/lucas/ISGlobal/Arrays/new_array_description_final.csv", sep = "\t", as.is = TRUE, header = TRUE)
gene_names <- read.csv("/home/lucas/ISGlobal/Arrays/Array_Annotation/array_gene_list_rosetta_annotated.txt", sep = "\t", as.is = TRUE, header = FALSE)

# 10E
s10E_T0_CTL <- read.table("/home/lucas/ISGlobal/Arrays/AnalisisR/Per LUCAS_1st analysis_ETF/10E/US10283823_258456110002_S01_GE2_1105_Oct12_1_4.txt",sep="\t",stringsAsFactors=FALSE, skip = 9, header = TRUE)
s10E_T0_CTL <- removeLowProbes(s10E_T0_CTL)
s10E_T0_HS <- read.table("/home/lucas/ISGlobal/Arrays/AnalisisR/Per LUCAS_1st analysis_ETF/10E/US10283823_258456110002_S01_GE2_1105_Oct12_1_4.txt",sep="\t",stringsAsFactors=FALSE, skip = 9, header = TRUE)
s10E_T0_HS <- removeLowProbes(s10E_T0_HS)
s10E_T1_CTL <- read.table("/home/lucas/ISGlobal/Arrays/AnalisisR/Per LUCAS_1st analysis_ETF/10E/US10283823_258456110002_S01_GE2_1105_Oct12_2_4.txt",sep="\t",stringsAsFactors=FALSE, skip = 9, header = TRUE)
s10E_T1_CTL <- removeLowProbes(s10E_T1_CTL)
s10E_T1_HS <- read.table("/home/lucas/ISGlobal/Arrays/AnalisisR/Per LUCAS_1st analysis_ETF/10E/US10283823_258456110003_S01_GE2_1105_Oct12_1_1.txt",sep="\t",stringsAsFactors=FALSE, skip = 9, header = TRUE)
s10E_T1_HS <- removeLowProbes(s10E_T1_HS)
s10G_T2_CTL <- read.table("/home/lucas/ISGlobal/Arrays/AnalisisR/Per LUCAS_1st analysis_ETF/10E/US10283823_258456110003_S01_GE2_1105_Oct12_2_1.txt",sep="\t",stringsAsFactors=FALSE, skip = 9, header = TRUE)
s10G_T2_CTL <- removeLowProbes(s10G_T2_CTL)
s10E_T2_HS <- read.table("/home/lucas/ISGlobal/Arrays/AnalisisR/Per LUCAS_1st analysis_ETF/10E/US10283823_258456110003_S01_GE2_1105_Oct12_1_2.txt",sep="\t",stringsAsFactors=FALSE, skip = 9, header = TRUE)
s10E_T2_HS <- removeLowProbes(s10E_T2_HS)
s10E_T3_CTL <- read.table("/home/lucas/ISGlobal/Arrays/AnalisisR/Per LUCAS_1st analysis_ETF/10E/US10283823_258456110003_S01_GE2_1105_Oct12_2_2.txt",sep="\t",stringsAsFactors=FALSE, skip = 9, header = TRUE)
s10E_T3_CTL <- removeLowProbes(s10E_T3_CTL)
s10E_T3_HS <- read.table("/home/lucas/ISGlobal/Arrays/AnalisisR/Per LUCAS_1st analysis_ETF/10E/US10283823_258456110003_S01_GE2_1105_Oct12_1_3.txt",sep="\t",stringsAsFactors=FALSE, skip = 9, header = TRUE)
s10E_T3_HS <- removeLowProbes(s10E_T3_HS)

#10G
s10G_T0_CTL <- read.table("/home/lucas/ISGlobal/Arrays/AnalisisR/Per LUCAS_1st analysis_ETF/10G/US10283823_258456110003_S01_GE2_1105_Oct12_2_3.txt",sep="\t",stringsAsFactors=FALSE, skip = 9, header = TRUE)
s10G_T0_CTL <- removeLowProbes(s10G_T0_CTL)
s10G_T0_HS <- read.table("/home/lucas/ISGlobal/Arrays/AnalisisR/Per LUCAS_1st analysis_ETF/10G/US10283823_258456110003_S01_GE2_1105_Oct12_2_3.txt",sep="\t",stringsAsFactors=FALSE, skip = 9, header = TRUE)
s10G_T0_HS <- removeLowProbes(s10G_T0_HS)
s10G_T1_CTL <- read.table("/home/lucas/ISGlobal/Arrays/AnalisisR/Per LUCAS_1st analysis_ETF/10G/US10283823_258456110003_S01_GE2_1105_Oct12_1_4.txt",sep="\t",stringsAsFactors=FALSE, skip = 9, header = TRUE)
s10G_T1_CTL <- removeLowProbes(s10G_T1_CTL)
s10G_T1_HS <- read.table("/home/lucas/ISGlobal/Arrays/AnalisisR/Per LUCAS_1st analysis_ETF/10G/US10283823_258456110003_S01_GE2_1105_Oct12_2_4.txt",sep="\t",stringsAsFactors=FALSE, skip = 9, header = TRUE)
s10G_T1_HS <- removeLowProbes(s10G_T1_HS)
s10E_T2_CTL <- read.table("/home/lucas/ISGlobal/Arrays/AnalisisR/Per LUCAS_1st analysis_ETF/10G/US10283823_258456110004_S01_GE2_1105_Oct12_1_1.txt",sep="\t",stringsAsFactors=FALSE, skip = 9, header = TRUE)
s10E_T2_CTL <- removeLowProbes(s10E_T2_CTL)
s10G_T2_HS <- read.table("/home/lucas/ISGlobal/Arrays/AnalisisR/Per LUCAS_1st analysis_ETF/10G/US10283823_258456110004_S01_GE2_1105_Oct12_2_1.txt",sep="\t",stringsAsFactors=FALSE, skip = 9, header = TRUE)
s10G_T2_HS <- removeLowProbes(s10G_T2_HS)
s10G_T3_CTL <- read.table("/home/lucas/ISGlobal/Arrays/AnalisisR/Per LUCAS_1st analysis_ETF/10G/US10283823_258456110004_S01_GE2_1105_Oct12_1_2.txt",sep="\t",stringsAsFactors=FALSE, skip = 9, header = TRUE)
s10G_T3_CTL <- removeLowProbes(s10G_T3_CTL)
s10G_T3_HS <- read.table("/home/lucas/ISGlobal/Arrays/AnalisisR/Per LUCAS_1st analysis_ETF/10G/US10283823_258456110004_S01_GE2_1105_Oct12_2_2.txt",sep="\t",stringsAsFactors=FALSE, skip = 9, header = TRUE)
s10G_T3_HS <- removeLowProbes(s10G_T3_HS)

#EK0
EK0_T0_CTL <- read.table("/home/lucas/ISGlobal/Arrays/AnalisisR/Per LUCAS_1st analysis_ETF/EKO/US10283823_258456110004_S01_GE2_1105_Oct12_1_3.txt",sep="\t",stringsAsFactors=FALSE, skip = 9, header = TRUE)
EK0_T0_CTL <- removeLowProbes(EK0_T0_CTL)
EK0_T0_HS <- read.table("/home/lucas/ISGlobal/Arrays/AnalisisR/Per LUCAS_1st analysis_ETF/EKO/US10283823_258456110004_S01_GE2_1105_Oct12_1_3.txt",sep="\t",stringsAsFactors=FALSE, skip = 9, header = TRUE)
EK0_T0_HS <- removeLowProbes(EK0_T0_HS)
EK0_T1_CTL <- read.table("/home/lucas/ISGlobal/Arrays/AnalisisR/Per LUCAS_1st analysis_ETF/EKO/US10283823_258456110004_S01_GE2_1105_Oct12_2_3.txt",sep="\t",stringsAsFactors=FALSE, skip = 9, header = TRUE)
EK0_T1_CTL <- removeLowProbes(EK0_T1_CTL)
EK0_T1_HS <- read.table("/home/lucas/ISGlobal/Arrays/AnalisisR/Per LUCAS_1st analysis_ETF/EKO/US10283823_258456110004_S01_GE2_1105_Oct12_1_4.txt",sep="\t",stringsAsFactors=FALSE, skip = 9, header = TRUE)
EK0_T1_HS <- removeLowProbes(EK0_T1_HS)
EK0_T2_CTL <- read.table("/home/lucas/ISGlobal/Arrays/AnalisisR/Per LUCAS_1st analysis_ETF/EKO/US10283823_258456110004_S01_GE2_1105_Oct12_2_4.txt",sep="\t",stringsAsFactors=FALSE, skip = 9, header = TRUE)
EK0_T2_CTL <- removeLowProbes(EK0_T2_CTL)
EK0_T2_HS <- read.table("/home/lucas/ISGlobal/Arrays/AnalisisR/Per LUCAS_1st analysis_ETF/EKO/US10283823_258456110005_S01_GE2_1105_Oct12_1_1.txt",sep="\t",stringsAsFactors=FALSE, skip = 9, header = TRUE)
EK0_T2_HS <- removeLowProbes(EK0_T2_HS)
EK0_T3_CTL <- read.table("/home/lucas/ISGlobal/Arrays/AnalisisR/Per LUCAS_1st analysis_ETF/EKO/US10283823_258456110005_S01_GE2_1105_Oct12_2_1.txt",sep="\t",stringsAsFactors=FALSE, skip = 9, header = TRUE)
EK0_T3_CTL <- removeLowProbes(EK0_T3_CTL)
EK0_T3_HS <- read.table("/home/lucas/ISGlobal/Arrays/AnalisisR/Per LUCAS_1st analysis_ETF/EKO/US10283823_258456110005_S01_GE2_1105_Oct12_1_2.txt",sep="\t",stringsAsFactors=FALSE, skip = 9, header = TRUE)
EK0_T3_HS <- removeLowProbes(EK0_T3_HS)

## Merge arrays in 1 df
df <- cbind(s10E_T0_CTL[,c(7,8,11,14,15)], 
            s10E_T0_HS[,c(11,14,15)], 
            s10E_T1_CTL[,c(11,14,15)], 
            s10E_T1_HS[,c(11,14,15)], 
            s10E_T2_CTL[,c(11,14,15)],
            s10E_T2_HS[,c(11,14,15)],
            s10E_T3_CTL[,c(11,14,15)],
            s10E_T3_HS[,c(11,14,15)],
            s10G_T0_CTL[,c(11,14,15)], 
            s10G_T0_HS[,c(11,14,15)], 
            s10G_T1_CTL[,c(11,14,15)], 
            s10G_T1_HS[,c(11,14,15)], 
            s10G_T2_CTL[,c(11,14,15)],
            s10G_T2_HS[,c(11,14,15)],
            s10G_T3_CTL[,c(11,14,15)],
            s10G_T3_HS[,c(11,14,15)],
            EK0_T0_CTL[,c(11,14,15)], 
            EK0_T0_HS[,c(11,14,15)], 
            EK0_T1_CTL[,c(11,14,15)], 
            EK0_T1_HS[,c(11,14,15)], 
            EK0_T2_CTL[,c(11,14,15)],
            EK0_T2_HS[,c(11,14,15)],
            EK0_T3_CTL[,c(11,14,15)],
            EK0_T3_HS[,c(11,14,15)])

colnames(df) <- c("ProbeName", "SystematicName", 
                  "s10E_T0_CTL","s10E_T0_CTL_G","s10E_T0_CTL_R", "s10E_T0_HS", "s10E_T0_HS_G", "s10E_T0_HS_R",
                  "s10E_T1_CTL","s10E_T1_CTL_G","s10E_T1_CTL_R", "s10E_T1_HS", "s10E_T1_HS_G", "s10E_T1_HS_R",
                  "s10E_T2_CTL","s10E_T2_CTL_G","s10E_T2_CTL_R", "s10E_T2_HS", "s10E_T2_HS_G", "s10E_T2_HS_R",
                  "s10E_T3_CTL","s10E_T3_CTL_G","s10E_T3_CTL_R", "s10E_T3_HS", "s10E_T3_HS_G", "s10E_T3_HS_R",
                  "s10G_T0_CTL","s10G_T0_CTL_G","s10G_T0_CTL_R", "s10G_T0_HS", "s10G_T0_HS_G", "s10G_T0_HS_R",
                  "s10G_T1_CTL","s10G_T1_CTL_G","s10G_T1_CTL_R", "s10G_T1_HS", "s10G_T1_HS_G", "s10G_T1_HS_R",
                  "s10G_T2_CTL","s10G_T2_CTL_G","s10G_T2_CTL_R", "s10G_T2_HS", "s10G_T2_HS_G", "s10G_T2_HS_R",
                  "s10G_T3_CTL","s10G_T3_CTL_G","s10G_T3_CTL_R", "s10G_T3_HS", "s10G_T3_HS_G", "s10G_T3_HS_R",
                  "EK0_T0_CTL","EK0_T0_CTL_G","EK0_T0_CTL_R", "EK0_T0_HS", "EK0_T0_HS_G", "EK0_T0_HS_R",
                  "EK0_T1_CTL","EK0_T1_CTL_G","EK0_T1_CTL_R", "EK0_T1_HS", "EK0_T1_HS_G", "EK0_T1_HS_R",
                  "EK0_T2_CTL","EK0_T2_CTL_G","EK0_T2_CTL_R", "EK0_T2_HS", "EK0_T2_HS_G", "EK0_T2_HS_R",
                  "EK0_T3_CTL","EK0_T3_CTL_G","EK0_T3_CTL_R", "EK0_T3_HS", "EK0_T3_HS_G", "EK0_T3_HS_R")


df["Gene_id"] <- array_description$New_Target
df["Gene_name"] <- gene_names$V3
df["Annot"] <- array_description$New_Annot

df[is.na(df$Gene_id),]$Gene_id <- df[is.na(df$Gene_id),]$SystematicName

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

df[df$ProbeName == "HSDHFR_v7.1_P1/1",]$Gene_id <- "HSDHFR"
df[df$ProbeName == "HSDHFR_v7.1_P1/1",]$Gene_name <- "HSDHFR"
df[df$ProbeName == "HSDHFR_v7.1_P1/1",]$Annot <- "hsdhfr"

# Treure gens no únics
no_unics <- array_description[!is.na(array_description$Status) & array_description$Status == "drop",]$ProbeName
df <- df[!df$ProbeName %in% no_unics,]

# Afegir sondes noves
plasmid_probes <- read.csv("/home/lucas/ISGlobal/Arrays/Array_Annotation/sondes_plasmids.csv", header = TRUE, sep = "\t")

for (i in 1:dim(plasmid_probes)[1]){
  df[df$ProbeName == plasmid_probes[i,"NAME"],"Gene_name"] <- as.character(plasmid_probes[i,"ACCESSION_STRING"])
  df[df$ProbeName == plasmid_probes[i,"NAME"],"Gene_id"] <- as.character(plasmid_probes[i,"ACCESSION_STRING"])
  df[df$ProbeName == plasmid_probes[i,"NAME"],"Annot"] <- as.character(plasmid_probes[i,"DESCRIPTION"])
}

#df[df$ProbeName %in% plasmid_probes$NAME,]

# Afegir informació sobre gens variants
variant <- read.csv("/home/lucas/ISGlobal/Gen_Referencies/Gens_variants_extended.txt", header = TRUE, sep = "\t")
df["Variant"] <- df$Gene_id %in% variant$ID

# Afegir families i GO terms
go <- read.csv("/home/lucas/ISGlobal/Gen_Referencies/PlasmoDB-37_Pfalciparum3D7_GO.gaf", header = FALSE, sep = "\t", skip = 1)
go["gene_id"] <- gsub("EuPathDB:", "", go$V17)
go$gene_id <- gsub(".1", "", go$gene_id)
table(df$Gene_id[!is.na(df$Gene_id)] %in% go$gene_id)

## Change to log2 (from log10 or from Processsignal)

# Log Ratio cols (originally log10)
df[,seq(3,72,3)] <- log2(10**df[,seq(3,72,3)])

# Processed raw signal Cols (originally no log)
df[,seq(4,73,3)] <- log2(df[,seq(4,73,3)])
df[,seq(5,74,3)] <- log2(df[,seq(5,74,3)])

## Create eSet
exprsx <- as.matrix(df[,seq(3,72,3)])
fdata <- new("AnnotatedDataFrame",df[,c(1:2,75:78)])
time <- factor(rep(c(0,0,1,1,2,2,3,3),3))
soca <- factor(gsub("(_.*_)", "_",colnames(exprsx), perl = TRUE))
pdata <- data.frame(soca=soca,time=time); rownames(pdata) <- colnames(df[,seq(3,72,3)])
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
  
  xgene <- by(exprsx[,seq(3,72,3)],geneid,myRma) 
  xgene <- do.call('rbind',xgene)
  
  mysd <- function(x) { ans <- ifelse(sum(!is.na(x))==1,0,sd(x,na.rm=TRUE)); return(ans) }
  sdgene <- aggregate(exprsx[,seq(3,72,3)],by=list(geneid),FUN=mysd)
  
  names(sdgene)[1] <- 'geneid'
  xgene <- data.frame(geneid=rownames(xgene),xgene); rownames(xgene) <- NULL
  
  fdata <- by(df[,75:78],geneid,unique)
  genenames <- names(fdata)
  fdata <- do.call('rbind',fdata)
  
  fdata <- new("AnnotatedDataFrame",data.frame(fdata))
  rownames(fdata) <- as.character(xgene$geneid)
  
  exprsxgene <- as.matrix(xgene[,-1])
  rownames(exprsxgene) <- as.character(xgene$geneid);
  eset <- new("ExpressionSet",exprs=exprsxgene,featureData=fdata,phenoData=pdata)
  return(list(eset=eset,sdgene=sdgene,fdata=fdata,geneid=geneid))
}
renameGenesAndSummarize_red <- function(genesToRename.sd,exprsx,geneid,summaryMethod=myRma) {
  # if (!is.na(genesToRename.sd)) {
  #   convTable <- read.csv(genesToRename.sd,header=TRUE)
  #   myBadOligos <- featureNames(x)[featureNames(x) %in% convTable$OligoID]
  #   for (i in 1:length(myBadOligos)) myBadOligos[i] <- as.character(convTable[convTable$OligoID==myBadOligos[i],'NewGeneID.'])
  #   geneid[featureNames(x) %in% convTable$OligoID] <- myBadOligos
  # }
  
  red_gene <- by(exprsx[,seq(5,74,3)],df["Gene_id"],myRma) 
  red_gene <- do.call('rbind',red_gene)
  
  mysd <- function(x) { ans <- ifelse(sum(!is.na(x))==1,0,sd(x,na.rm=TRUE)); return(ans) }
  sdgene <- aggregate(exprsx[,seq(5,74,3)],by=exprsx["Gene_id"],FUN=mysd)
  
  names(sdgene)[1] <- 'geneid'
  red_gene <- data.frame(geneid=rownames(red_gene),red_gene); rownames(red_gene) <- NULL
  
  fdata <- by(df[,75:78],exprsx["Gene_id"],unique)
  genenames <- names(fdata)
  fdata <- do.call('rbind',fdata)
  
  fdata <- new("AnnotatedDataFrame",data.frame(fdata))
  rownames(fdata) <- as.character(red_gene$geneid)
  
  exprsred_gene <- as.matrix(red_gene[,-1])
  rownames(exprsred_gene) <- red_gene$geneid;
  colnames(exprsred_gene) <- rownames(pdata@data)
  eset <- new("ExpressionSet",exprs=exprsred_gene,featureData=fdata,phenoData=pdata)
  return(list(eset=eset,sdgene=sdgene,fdata=fdata,geneid=geneid))
}

tmp <- renameGenesAndSummarize(genesToRename.sd=genesToRename.sd,exprsx=df,geneid=geneid,summaryMethod=myRma)
tmp_red <- renameGenesAndSummarize_red(genesToRename.sd=genesToRename.sd,exprsx=df,geneid=geneid,summaryMethod=myRma)

xgene <- tmp[['eset']]; sdgene <- tmp[['sdgene']]; fdata <- tmp[['fdata']]; geneid <- tmp[['geneid']]
red_gene <- tmp_red[['eset']]; sdgene <- tmp_red[['sdgene']]; fdata <- tmp_red[['fdata']]; geneid <- tmp_red[['geneid']]

write.csv(data.frame(probe=df[,1],gene=geneid),file.path(figuresPath,'geneProbes.csv'),row.names=FALSE) #probes of each gene

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
  myCI <- ci(ll)
  myLL <- loglik(ll)
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
  times <- list(myTimes, myCI)
  return(myTimes)
}

estimatedTimes <- getTimeEstimation(xgene,bozdechPath,LemieuxFunctionsPath,file.path(figuresPath),B=100)

write.table(estimatedTimes, "/home/lucas/ISGlobal/Arrays/Eli_Arrays/estimated_times.csv")

CTLEstimatedTimes <- as.numeric(unlist(lapply(estimatedTimes[rep(c(TRUE, FALSE), 12)], function(x) rep(x, 2))))
names(CTLEstimatedTimes) <- names(estimatedTimes)
estimatedTimes <- CTLEstimatedTimes

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

############################### Plots at probe level ##########################################

load(file.path(dir,'normalizedData_probeLevel.RData'))

## Ratio
#for (i in 1:nrow(x)){
for (i in 1:10){
  print(i)
  graf <- melt(exprs(x)[i,])
  graf["Soca"] <- c(rep("10E", 8), rep("10G", 8), rep("EK0", 8))
  graf["Time"] <- estimatedTimes
  graf["Type"] <- rep(c("CTL", "HS"), 12)
  graf["group"] <- c(rep(c(1,2), 4), rep(c(3,4), 4), rep(c(5,6), 4))
  p <- ggplot(graf, aes(x = Time, y = value, col = Soca, linetype = Type, group = group))
  p <- p + geom_point(aes(color = Soca, shape = Soca)) + geom_line() + scale_linetype_manual(values=c("dashed", "solid")) 
  p <- p + coord_cartesian(ylim = c(-9, 7)) + ggtitle(x@featureData@data$ProbeName[i])
  ggsave(p, file=paste0(figuresPath,"Ratio/Probes/", as.character(i), "_", gsub("/", ":", x@featureData@data$Gene_id[i]), ".jpeg"), device = "jpeg", width = 14, height = 10, units = "cm")
  #print(p)
}

## Red Signal
#for (i in 1:nrow(df)){
for (i in 1:10){
  print(i)
  graf <- melt(df[i,seq(5,74,3)])
  graf["Soca"] <- c(rep("10E", 8), rep("10G", 8), rep("EK0", 8))
  graf["Time"] <- estimatedTimes
  graf["Type"] <- rep(c("CTL", "HS"), 12)
  graf["group"] <- c(rep(c(1,2), 4), rep(c(3,4), 4), rep(c(5,6), 4))
  p <- ggplot(graf, aes(x = Time, y = value, col = Soca, linetype = Type, group = group))
  p <- p + geom_point(aes(color = Soca, shape = Soca)) + geom_line() + scale_linetype_manual(values=c("dashed", "solid")) 
  p <- p + coord_cartesian(ylim = c(-1.5, 18.5)) + ggtitle(df$ProbeName[i])
  ggsave(p, file=paste0(figuresPath,"RedSignal/Probes/", as.character(i), "_", gsub("/", ":", df$Gene_id[i]), ".jpeg"), device = "jpeg", width = 14, height = 10, units = "cm")
  #print(p)
}

############################### Point estimations ##########################################

load(file.path(dir,'normalizedData_geneLevel.RData'))

exprsx <- exprs(xgene)
soca <- as.factor(pData(xgene)$soca)
oldTime <- as.factor(rep(c(0,0,1,1,2,2,3,3),3)) ### OOOOJUUUU ERROOOR

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

## Ratio
#for (i in 1:nrow(exprs(xgene))){
for (i in 1:10){
  print(i)
  graf <- melt(exprs(xgene)[i,])
  graf["Soca"] <- c(rep("10E", 8), rep("10G", 8), rep("EK0", 8))
  graf["Time"] <- estimatedTimes
  graf["Type"] <- rep(c("CTL", "HS"), 12)
  graf["group"] <- c(rep(c(1,2), 4), rep(c(3,4), 4), rep(c(5,6), 4))
  p <- ggplot(graf, aes(x = Time, y = value, col = Soca, linetype = Type, group = group))
  p <- p + geom_point(aes(color = Soca, shape = Soca)) + geom_line() + scale_linetype_manual(values=c("dashed", "solid")) 
  p <- p + coord_cartesian(ylim = c(-9, 7)) + ggtitle(xgene@featureData@data$Gene_id[i])
  #ggsave(p, file=paste0(figuresPath,"Ratio/Genes/", as.character(i), "_", gsub("/", ":", xgene@featureData@data$Gene_id[i]), ".jpeg"), device = "jpeg", width = 14, height = 10, units = "cm")
  print(p)
}

## Red Signal
#for (i in 1:nrow(exprs(red_gene))){
for (i in 1:10){
  print(i)
  graf <- melt(exprs(red_gene)[i,])
  graf["Soca"] <- c(rep("10E", 8), rep("10G", 8), rep("EK0", 8))
  graf["Time"] <- estimatedTimes
  graf["Type"] <- rep(c("CTL", "HS"), 12)
  graf["group"] <- c(rep(c(1,2), 4), rep(c(3,4), 4), rep(c(5,6), 4))
  p <- ggplot(graf, aes(x = Time, y = value, col = Soca, linetype = Type, group = group))
  p <- p + geom_point(aes(color = Soca, shape = Soca)) + geom_line() + scale_linetype_manual(values=c("dashed", "solid")) 
  p <- p + coord_cartesian(ylim = c(-1.5, 18.5)) + ggtitle(red_gene@featureData@data$Gene_id[i])
  ggsave(p, file=paste0(figuresPath,"RedSignal/Genes/", as.character(i), "_", gsub("/", ":", red_gene@featureData@data$Gene_id[i]), ".jpeg"), device = "jpeg", width = 14, height = 10, units = "cm")
  #print(p)
}

############################### Arees #####################

# Per alguna raó 10G el compta desde 25h (primer timepoint) en comptes de 27h (primer timepoint compartit) !!!

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
myOrder <- order(pData(xgene)$soca,pData(xgene)$time)
xgene <- xgene[,myOrder]
estim <- estim[,myOrder]
geneDesc <- fData(xgene)$Gene_name

#preprocess
soques <- levels(pData(xgene)$soca)
if (!file.exists(file.path(figuresPath,'soques_genelevel_estimated/'))) dir.create(file.path(figuresPath,'soques_genelevel_estimated/'))
outDir <- file.path(figuresPath,'soques_genelevel_estimated/')
time.est <- pData(xgene)$time
time <- as.numeric(rep(c(1,2,3,4), 6))
soque <- pData(xgene)$soca
exprsx <- exprs(xgene)

#set same min and max time
maxMin.time.est <- max(sapply(unique(soque),function(x) min(time.est[soque==x])))
minMax.time.est <- min(sapply(unique(soque),function(x) max(time.est[soque==x])))
#exprsx <- timeCorrectedEpxrs(exprsx,maxMin.time.est,minMax.time.est,time,soque)
#sel.previous2first <- as.logical(sapply(unique(soque),function(x) time.est[soque==x]==max(time.est[soque==x][time.est[soque==x]<=maxMin.time.est])))
#sel.next2last <- as.logical(sapply(unique(soque),function(x) time.est[soque==x]==min(time.est[soque==x][time.est[soque==x]>=minMax.time.est])))

sel.previous2first <- time.est <= maxMin.time.est
sel.next2last <- time.est >= minMax.time.est

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
rownames(area) <- rownames(exprsx)

write.csv(area, "/home/lucas/ISGlobal/Arrays/Eli_Arrays/Arres_revisited/arees_corregit.csv")

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

 
############################### Llistes #################################

# Top 50
tops <- c()
strains <- c("s10E", "s10G", "EK0")
for (strain in strains){
  load(file.path(dir,'xgene_estimated.RData'))
  xgene <- xgene[,startsWith(colnames(exprs(xgene)), strain)]
  
  heatmap_df <- data.frame(c(exprs(xgene)[,2] - exprs(xgene)[,1]), 
                           c(exprs(xgene)[,4] - exprs(xgene)[,5]), 
                           c(exprs(xgene)[,6] - exprs(xgene)[,5]), 
                           c(exprs(xgene)[,8] - exprs(xgene)[,7]))
  
  colnames(heatmap_df) <- c("T1", "T2", "T3", "T4")
  
  newlist = c()
  for (x in 1:length(xgene@featureData@data$Gene_name)){
    if (is.na(xgene@featureData@data$Gene_name[x])){
      newlist = c(newlist,xgene@featureData@data$Annot[x])
    } else {
      newlist = c(newlist,xgene@featureData@data$Gene_name[x])
    }
  }
  
  heatmap_df["gene"] <- paste0(newlist, ": ", rownames(heatmap_df))
  #sel_FC25 <- filter(heatmap_df, abs(T1) > 0.39794 | abs(T2) > 0.39794 | abs(T3) > 0.39794)
  max(heatmap_df$T4)
  m <- melt(heatmap_df)
  top50 <- arrange(m, -abs(value))[1:50,]
  sel_FC25 <- m[abs(m$value) > 0.39794,]
  top50 <- arrange(top50, -value)
  sel_FC25 <- arrange(sel_FC25, -value)
  top50["FC"] <- 10**top50$value
  sel_FC25["FC"] <- 10**sel_FC25$value
  write.csv(top50,file.path(figuresPath, paste0(strain,"_top50.csv")))
  write.csv(sel_FC25,file.path(figuresPath, paste0(strain,"_FC25.csv")))
  tops <- c(tops, top50$gene)
}


load(file.path(dir,'xgene_estimated.RData'))
xgene <- xgene[,startsWith(colnames(exprs(xgene)), "s10G")]

heatmap_df <- data.frame(c(exprs(xgene)[,2] - exprs(xgene)[,1]), 
                         c(exprs(xgene)[,4] - exprs(xgene)[,5]), 
                         c(exprs(xgene)[,6] - exprs(xgene)[,5]), 
                         c(exprs(xgene)[,8] - exprs(xgene)[,7]))

colnames(heatmap_df) <- c("T1", "T2", "T3", "T4")

newlist = c()
for (x in 1:length(xgene@featureData@data$Gene_name)){
  if (is.na(xgene@featureData@data$Gene_name[x])){
    newlist = c(newlist,xgene@featureData@data$Annot[x])
  } else {
    newlist = c(newlist,xgene@featureData@data$Gene_name[x])
  }
}

heatmap_df["gene"] <- paste0(newlist, ": ", rownames(heatmap_df))
sel_all <- filter(heatmap_df, abs(T2) > 0.39794 | abs(T3) > 0.39794 | abs(T4) > 0.39794)
sel_all[,c(2:4)] <- 10**sel_all[,c(2:4)]
head(sel_all)

############################### Heatmaps #################################

tops <- unique(tops)

load(file.path(dir,'xgene_estimated.RData'))

newlist = c()
for (x in 1:length(xgene@featureData@data$Gene_name)){
  if (is.na(xgene@featureData@data$Gene_name[x])){
    newlist = c(newlist,xgene@featureData@data$Annot[x])
  } else {
    newlist = c(newlist,xgene@featureData@data$Gene_name[x])
  }
}

short_names <- paste0(newlist, ": ", xgene@featureData@data$Gene_id)
xgene@featureData@data["short_names"] <- short_names
rownames(xgene) <-xgene@featureData@data$short_names

heatmap_df <- data.frame(c(exprs(xgene)[,2] - exprs(xgene)[,1]), 
                         c(exprs(xgene)[,4] - exprs(xgene)[,5]), 
                         c(exprs(xgene)[,6] - exprs(xgene)[,5]), 
                         c(exprs(xgene)[,8] - exprs(xgene)[,7]),
                         c(exprs(xgene)[,10] - exprs(xgene)[,9]), 
                         c(exprs(xgene)[,12] - exprs(xgene)[,11]), 
                         c(exprs(xgene)[,14] - exprs(xgene)[,13]),
                         c(exprs(xgene)[,16] - exprs(xgene)[,15]), 
                         c(exprs(xgene)[,18] - exprs(xgene)[,17]),
                         c(exprs(xgene)[,20] - exprs(xgene)[,19]), 
                         c(exprs(xgene)[,22] - exprs(xgene)[,21]), 
                         c(exprs(xgene)[,24] - exprs(xgene)[,23]))

colnames(heatmap_df) <- gsub("_CTL", "", colnames(exprs(xgene))[rep(c(TRUE, FALSE), 12)])
heatmap_df <- heatmap_df[,c(-1,-5,-9)]
sel_all <- heatmap_df[rownames(heatmap_df) %in% tops, ]
sel_all_T<- sel_all[,c(1,4,7,2,5,8,3,6,9)]

# Hierarchical clustering with cols
heatmap(
  as.matrix(sel_all), Rowv=as.dendrogram(hclust(dist(as.matrix(sel_all))),
                                               Colv=FALSE, cexRow = 2)
)

# Hierarchical clustering without cols

hclustfunc <- function(x) hclust(x, method="complete")
distfunc <- function(x) dist(x,method="euclidean")
# par(mar=c(1,1,1,1))
# png("/home/lucas/ISGlobal/Arrays/Eli_Arrays/Plots/heatmap_tops.png", height = 400, width = 800)
heatmap.2(col=redgreen(75),
          as.matrix(sel_all_T), 
          Colv=FALSE,
          keysize = 1,
          dendrogram="row",
          trace="none", 
          hclust=hclustfunc,
          distfun=distfunc,
          colsep=c(3,6),
          margins = c(6, 16),
          cexRow=0.5,
          lmat = rbind(3:4,2:1), lwid = c(0.3,4), lhei = c(0.5,4))
dev.off()
          

############################### Plots Custom #############################

tops <- unlist(lapply(tops, function(x) strsplit(x, ": ")[[1]][2]))
load(file.path(dir,'normalizedData_geneLevel.RData'))
xgene <- xgene[xgene@featureData@data$Gene_id %in% tops,]
# CTLEstimatedTimes <- as.numeric(unlist(lapply(estimatedTimes[rep(c(TRUE, FALSE), 12)], function(x) rep(x, 2))))
# names(CTLEstimatedTimes) <- names(estimatedTimes)
custom_plots_dir <- paste0(figuresPath, "Custom_plots/")

##Plots
#for (i in 1:nrow(exprs(xgene))){
for (i in 1:10){
  print(i)
  graf <- melt(exprs(xgene)[i,])
  #graf["Time"] <- estimatedTimes
  graf["Time"] <- c("T0", "T0", "T1", "T1", "T2", "T2", "T3", "T3")
  graf["Type"] <- rep(c("CTL", "HS"), 12)
  graf["Soque"] <- c(rep("s10E", 8), rep("s10G", 8), rep("EK0", 8))
  graf["group"] <- c(rep(c(1,2), 4), rep(c(3,4), 4), rep(c(5,6), 4))
  p <- ggplot(graf, aes(x = Time, y = value, col = Soque, linetype = Type, group = group))
  p <- p + geom_point(aes(color = Soque, shape = Soque)) + geom_line() + scale_linetype_manual(values=c("dashed", "solid"))
  ggsave(p, file=paste(custom_plots_dir,gsub('/','',gsub(':','#',featureNames(xgene)[i])),'.jpeg',sep=''), device = "jpeg", width = 14, height = 10, units = "cm")
}

# graf <- melt(exprs(xgene)[1,])
# #graf["Time"] <- CTLEstimatedTimes
# graf["Time"] <- c("T0", "T0", "T1", "T1", "T2", "T2", "T3", "T3")
# graf["Type"] <- rep(c("CTL", "HS"), 12)
# graf["Soque"] <- c(rep("s10E", 8), rep("s10G", 8), rep("EK0", 8))
# graf["group"] <- c(rep(c(1,2), 4), rep(c(3,4), 4), rep(c(5,6), 4))
# p <- ggplot(graf, aes(x = Time, y = value, col = Soque, linetype = Type, group = group))
# p <- p + geom_point(aes(color = Soque, shape = Soque)) + geom_line() + scale_linetype_manual(values=c("dashed", "solid"))
# p






############################### Impute reference Set: Loess ##########################################
ref <- read.csv(file = bozdechPath, as.is = T)

preds <- list()
for (i in 1:nrow(ref)){
  m <- melt(ref[i,])
  m["Time"] <- c(1:22, 24:28, 30:48)
  mloess <- loess(value ~ Time, m, span=0.4)
  pred <- predict(mloess, data.frame(Time=seq(1, 48, 0.1)))
  preds[[i]] <- pred
}
pred_mtx <- do.call("rbind", preds)
pred_df <- as.data.frame(pred_mtx)
colnames(pred_df) <- seq(1, 48, 0.1)
pred_df["Name"] <- ref$Name
write.csv(pred_df, file = "/media/lucas/Disc4T/Projects/Arrays_Eli/bodzek_imputed.csv")

m1 <- melt(ref[1,])
m1["Time"] <- c(1:22, 24:28, 30:48)
m1 <-rbind(m1, data.frame(Name="PF3D7_1318700", variable="TP23", value=NA, Time=23))
m1 <-rbind(m1, data.frame(Name="PF3D7_1318700", variable="TP29", value=NA, Time=29))
m1 <- m1[order(m1$Time),]

m1loess <- loess(value ~ Time, m1, span=0.4)
pred2 <- predict(m1loess, data.frame(Time=seq(1, 48, 0.1)))
pred <- predict(m1loess, data.frame(Time=c(23,29)))

ggplot(m1, aes(x = Time, y = value)) + geom_point() +  geom_smooth(span = 0.4) + geom_smooth(span = 0.1, color = "red")
m1loess <- loess(value ~ Time, m1, span=0.4)

m1[23, "value"] <- pred[1]
m1[29, "value"] <- pred[2]

m1["Status"] <- "real"
m1[c(23,29), "Status"] <- "predicted"

ggplot(m1, aes(x = Time, y = value)) + geom_point(aes(color = Status)) +  geom_smooth(span = 0.4) + coord_cartesian(ylim = c(0, 3))

pred2 <- predict(m1loess, data.frame(Time=seq(1, 48, 0.1)))
pred_df <- data.frame(pred2, "Time" = seq(1, 48, 0.1), "Status" = "pred")
 
ggplot(pred_df, aes(x = Time, y = pred2)) + geom_point(aes(color = Status))+ coord_cartesian(ylim = c(0, 3))

############################### Impute reference Set: Splines ##########################################
ref <- read.csv(file = bozdechPath, as.is = T)
preds <- list()
for (i in 1:nrow(ref)){
  m <- melt(ref[i,])
  m["Time"] <- c(1:22, 24:28, 30:48)
  msplines <- smooth.spline(x=m$Time, y=m$value, spar=0.5)
  pred <- predict(msplines, data.frame(Time=seq(1, 48, 0.1)))
  preds[[i]] <- pred$y
}
pred_mtx <- do.call("cbind", preds)
pred_mtx <- t(pred_mtx)
pred_df <- as.data.frame(pred_mtx)
colnames(pred_df) <- seq(1, 48, 0.1)
pred_df["Name"] <- ref$Name
write.csv(pred_df, file = "/media/lucas/Disc4T/Projects/Arrays_Eli/bodzek_imputed_splines.csv")

m1 <- melt(ref[1,])
m1["Time"] <- c(1:22, 24:28, 30:48)
ggplot(m1, aes(x = Time, y = value)) + geom_point() +  geom_smooth(method = "spline") + coord_cartesian(ylim = c(0, 3))

pdf <- as.data.frame(pred, "Time" = seq(1, 48, 0.1))
ggplot(pdf, aes(x = Time, y = Time.1)) + geom_point() + coord_cartesian(ylim = c(0, 3))



############################## Calculate similarity to Timepoint ##########################################
figPath <- "/media/lucas/Disc4T/Projects/Arrays_Eli/Arrays/Plots/"
ref <- read.csv(file = bozdechPath, as.is = T)
cors <- list()
idx <- 1
for (smpl in rownames(pData(xgene))){
  vals <- c()
  for (i in 1:dim(ref)[1]){
    if (ref[i,]$Name %in% rownames(exprs(xgene))){
      vals <- c(vals, exprs(xgene)[rownames(exprs(xgene)) == ref[i,]$Name, smpl])
    } else {
      vals <- c(vals, NA)
    }
  }
  ref[smpl] <- vals
  noNA_ref <- ref[complete.cases(ref),]
  smpl_cor <- cor(noNA_ref[,2:47], noNA_ref[,smpl])
  cors[[idx]] <- smpl_cor
  max_row <- which.max(smpl_cor)
  maxTime <- rownames(smpl_cor)[max_row]
  png(paste0(figPath, "corplot", smpl))
  plot(noNA_ref[,maxTime], noNA_ref[,smpl])
  dev.off()
  idx <- idx+1
} 

cor_mtx <- do.call("cbind", cors)
cor_df <- as.data.frame(cor_mtx)
colnames(cor_df) <- rownames(pData(xgene))
write.csv(cor_df, file = "/media/lucas/Disc4T/Projects/Arrays_Eli/time_estimation_correlations_logged.csv")

#################

vals <- c()
for (i in 1:dim(ref)[1]){
  if (ref[i,]$Name %in% rownames(exprs(xgene))){
    vals <- c(vals, exprs(xgene)[rownames(exprs(xgene)) == ref[i,]$Name, "s10E_T0_CTL"])
  } else {
    vals <- c(vals, NA)
  }
}
length(vals)

ref["s10E_CTL_27"] <- vals


which.max(x)
rownames(x)[26]

noNA_ref <- ref[complete.cases(ref),]
noNA_ref["nolog_10E"] <- 2**(noNA_ref$s10E_CTL_27)

cor(noNA_ref[,2:47], noNA_ref$s10E_CTL_27)
cor(noNA_ref[,2:48], noNA_ref$nolog_10E)

plot(noNA_ref$TP27, noNA_ref$s10E_CTL_27)
plot(noNA_ref$TP27, noNA_ref$nolog_10E)

summary(noNA_ref)
hist(noNA_ref$nolog_10E, breaks = 200)

noOutliers <- noNA_ref[noNA_ref$nolog_10E < 3,]
dim(noNA_ref)
dim(noOutliers)

cor(noOutliers[,2:48], noOutliers$nolog_10E)
plot(noOutliers$TP27, noOutliers$nolog_10E)

##
imputed_df <- pred_df
vals <- c()
for (i in 1:dim(imputed_df)[1]){
  if (imputed_df[i,]$Name %in% rownames(exprs(xgene))){
    vals <- c(vals, exprs(xgene)[rownames(exprs(xgene)) == imputed_df[i,]$Name, "s10E_T0_CTL"])
  } else {
    vals <- c(vals, NA)
  }
}
length(vals)

imputed_df["s10E_CTL_27"] <- vals

noNA_ref <- imputed_df[complete.cases(imputed_df),]
#noNA_ref["nolog_10E"] <- 2**(noNA_ref$s10E_CTL_27)

dim(noNA_ref)
cor_imputed <- cor(noNA_ref[,c(-472,-473)], noNA_ref$s10E_CTL_27)
#cor(noNA_ref[,2:48], noNA_ref$nolog_10E)

plot(noNA_ref$TP27, noNA_ref$s10E_CTL_27)
plot(noNA_ref$TP27, noNA_ref$nolog_10E)

summary(noNA_ref)
hist(noNA_ref$nolog_10E, breaks = 200)

pData(xgene)

##

imputed_df <- pred_df
vals <- c()
for (i in 1:dim(imputed_df)[1]){
  if (imputed_df[i,]$Name %in% rownames(exprs(xgene))){
    vals <- c(vals, exprs(xgene)[rownames(exprs(xgene)) == imputed_df[i,]$Name, "s10E_T1_CTL"])
  } else {
    vals <- c(vals, NA)
  }
}
length(vals)

imputed_df["s10E_CTL_27"] <- vals

noNA_ref <- imputed_df[complete.cases(imputed_df),]
#noNA_ref["nolog_10E"] <- 2**(noNA_ref$s10E_CTL_27)

dim(noNA_ref)
cor_imputed <- cor(noNA_ref[,c(-472,-473)], noNA_ref$s10E_CTL_27)
which(cor_imputed == max(cor_imputed), arr.ind = T)

pData(xgene)






































