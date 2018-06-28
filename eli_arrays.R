## Import array data
# 10E
s10E_T0_CTL <- read.table("/home/lucas/ISGlobal/Arrays/AnalisisR/Per LUCAS_1st analysis_ETF/10E/US10283823_258456110002_S01_GE2_1105_Oct12_1_4.txt",sep="\t",stringsAsFactors=FALSE, skip = 9, header = TRUE)
s10E_T0_HS <- read.table("/home/lucas/ISGlobal/Arrays/AnalisisR/Per LUCAS_1st analysis_ETF/10E/US10283823_258456110002_S01_GE2_1105_Oct12_1_4.txt",sep="\t",stringsAsFactors=FALSE, skip = 9, header = TRUE)
s10E_T1_CTL <- read.table("/home/lucas/ISGlobal/Arrays/AnalisisR/Per LUCAS_1st analysis_ETF/10E/US10283823_258456110002_S01_GE2_1105_Oct12_2_4.txt",sep="\t",stringsAsFactors=FALSE, skip = 9, header = TRUE)
s10E_T1_HS <- read.table("/home/lucas/ISGlobal/Arrays/AnalisisR/Per LUCAS_1st analysis_ETF/10E/US10283823_258456110003_S01_GE2_1105_Oct12_1_1.txt",sep="\t",stringsAsFactors=FALSE, skip = 9, header = TRUE)
s10E_T2_CTL <- read.table("/home/lucas/ISGlobal/Arrays/AnalisisR/Per LUCAS_1st analysis_ETF/10E/US10283823_258456110003_S01_GE2_1105_Oct12_2_1.txt",sep="\t",stringsAsFactors=FALSE, skip = 9, header = TRUE)
s10E_T2_HS <- read.table("/home/lucas/ISGlobal/Arrays/AnalisisR/Per LUCAS_1st analysis_ETF/10E/US10283823_258456110003_S01_GE2_1105_Oct12_1_2.txt",sep="\t",stringsAsFactors=FALSE, skip = 9, header = TRUE)
s10E_T3_CTL <- read.table("/home/lucas/ISGlobal/Arrays/AnalisisR/Per LUCAS_1st analysis_ETF/10E/US10283823_258456110003_S01_GE2_1105_Oct12_2_2.txt",sep="\t",stringsAsFactors=FALSE, skip = 9, header = TRUE)
s10E_T3_HS <- read.table("/home/lucas/ISGlobal/Arrays/AnalisisR/Per LUCAS_1st analysis_ETF/10E/US10283823_258456110003_S01_GE2_1105_Oct12_1_3.txt",sep="\t",stringsAsFactors=FALSE, skip = 9, header = TRUE)

#10G
s10G_T0_CTL <- read.table("/home/lucas/ISGlobal/Arrays/AnalisisR/Per LUCAS_1st analysis_ETF/10G/US10283823_258456110003_S01_GE2_1105_Oct12_2_3.txt",sep="\t",stringsAsFactors=FALSE, skip = 9, header = TRUE)
s10G_T0_HS <- read.table("/home/lucas/ISGlobal/Arrays/AnalisisR/Per LUCAS_1st analysis_ETF/10G/US10283823_258456110003_S01_GE2_1105_Oct12_2_3.txt",sep="\t",stringsAsFactors=FALSE, skip = 9, header = TRUE)
s10G_T1_CTL <- read.table("/home/lucas/ISGlobal/Arrays/AnalisisR/Per LUCAS_1st analysis_ETF/10G/US10283823_258456110003_S01_GE2_1105_Oct12_1_4.txt",sep="\t",stringsAsFactors=FALSE, skip = 9, header = TRUE)
s10G_T1_HS <- read.table("/home/lucas/ISGlobal/Arrays/AnalisisR/Per LUCAS_1st analysis_ETF/10G/US10283823_258456110003_S01_GE2_1105_Oct12_2_4.txt",sep="\t",stringsAsFactors=FALSE, skip = 9, header = TRUE)
s10G_T2_CTL <- read.table("/home/lucas/ISGlobal/Arrays/AnalisisR/Per LUCAS_1st analysis_ETF/10G/US10283823_258456110004_S01_GE2_1105_Oct12_1_1.txt",sep="\t",stringsAsFactors=FALSE, skip = 9, header = TRUE)
s10G_T2_HS <- read.table("/home/lucas/ISGlobal/Arrays/AnalisisR/Per LUCAS_1st analysis_ETF/10G/US10283823_258456110004_S01_GE2_1105_Oct12_2_1.txt",sep="\t",stringsAsFactors=FALSE, skip = 9, header = TRUE)
s10G_T3_CTL <- read.table("/home/lucas/ISGlobal/Arrays/AnalisisR/Per LUCAS_1st analysis_ETF/10G/US10283823_258456110004_S01_GE2_1105_Oct12_1_2.txt",sep="\t",stringsAsFactors=FALSE, skip = 9, header = TRUE)
s10G_T3_HS <- read.table("/home/lucas/ISGlobal/Arrays/AnalisisR/Per LUCAS_1st analysis_ETF/10G/US10283823_258456110004_S01_GE2_1105_Oct12_2_2.txt",sep="\t",stringsAsFactors=FALSE, skip = 9, header = TRUE)

#EK0
EK0_T0_CTL <- read.table("/home/lucas/ISGlobal/Arrays/AnalisisR/Per LUCAS_1st analysis_ETF/EKO/US10283823_258456110004_S01_GE2_1105_Oct12_1_3.txt",sep="\t",stringsAsFactors=FALSE, skip = 9, header = TRUE)
EK0_T0_HS <- read.table("/home/lucas/ISGlobal/Arrays/AnalisisR/Per LUCAS_1st analysis_ETF/EKO/US10283823_258456110004_S01_GE2_1105_Oct12_1_3.txt",sep="\t",stringsAsFactors=FALSE, skip = 9, header = TRUE)
EK0_T1_CTL <- read.table("/home/lucas/ISGlobal/Arrays/AnalisisR/Per LUCAS_1st analysis_ETF/EKO/US10283823_258456110004_S01_GE2_1105_Oct12_2_3.txt",sep="\t",stringsAsFactors=FALSE, skip = 9, header = TRUE)
EK0_T1_HS <- read.table("/home/lucas/ISGlobal/Arrays/AnalisisR/Per LUCAS_1st analysis_ETF/EKO/US10283823_258456110004_S01_GE2_1105_Oct12_1_4.txt",sep="\t",stringsAsFactors=FALSE, skip = 9, header = TRUE)
EK0_T2_CTL <- read.table("/home/lucas/ISGlobal/Arrays/AnalisisR/Per LUCAS_1st analysis_ETF/EKO/US10283823_258456110004_S01_GE2_1105_Oct12_2_4.txt",sep="\t",stringsAsFactors=FALSE, skip = 9, header = TRUE)
EK0_T2_HS <- read.table("/home/lucas/ISGlobal/Arrays/AnalisisR/Per LUCAS_1st analysis_ETF/EKO/US10283823_258456110005_S01_GE2_1105_Oct12_1_1.txt",sep="\t",stringsAsFactors=FALSE, skip = 9, header = TRUE)
EK0_T3_CTL <- read.table("/home/lucas/ISGlobal/Arrays/AnalisisR/Per LUCAS_1st analysis_ETF/EKO/US10283823_258456110005_S01_GE2_1105_Oct12_2_1.txt",sep="\t",stringsAsFactors=FALSE, skip = 9, header = TRUE)
EK0_T3_HS <- read.table("/home/lucas/ISGlobal/Arrays/AnalisisR/Per LUCAS_1st analysis_ETF/EKO/US10283823_258456110005_S01_GE2_1105_Oct12_1_2.txt",sep="\t",stringsAsFactors=FALSE, skip = 9, header = TRUE)

## Merge arrays in 1 df
eli_all_df <- cbind(s10E_T0_CTL[,c(7,8,11)], 
                    s10E_T0_HS[,11], 
                    s10E_T1_CTL[,11], 
                    s10E_T1_HS[,11], 
                    s10E_T2_CTL[,11],
                    s10E_T2_HS[,11],
                    s10E_T3_CTL[,11],
                    s10E_T3_HS[,11],
                    s10G_T0_CTL[,11], 
                    s10G_T0_HS[,11], 
                    s10G_T1_CTL[,11], 
                    s10G_T1_HS[,11], 
                    s10G_T2_CTL[,11],
                    s10G_T2_HS[,11],
                    s10G_T3_CTL[,11],
                    s10G_T3_HS[,11],
                    EK0_T0_CTL[,11], 
                    EK0_T0_HS[,11], 
                    EK0_T1_CTL[,11], 
                    EK0_T1_HS[,11], 
                    EK0_T2_CTL[,11],
                    EK0_T2_HS[,11],
                    EK0_T3_CTL[,11],
                    EK0_T3_HS[,11])

colnames(eli_all_df)[c(-1:-3)] <- lapply(colnames(eli_all_df)[c(-1:-3)], function(x) substr(x,1,nchar(x)-6))
colnames(eli_all_df)[3] <- "s10E_T0_CTL"

## Import rosseta and translate gene_names:
roseta <- read.table("/home/lucas/ISGlobal/gens_array_rosetta.txt")
eli_all_df["Gene_name"] <- sapply(eli_all_df$SystematicName, function(x) roseta[roseta$V1 == x,2])

## Unifiy rows by gene_name
unified_df <- eli_all_df %>% group_by(Gene_name) %>% summarise_all(mean)

# Create eset
emtx <- as.matrix(unified_df[,c(-1,-2,-3)])
rownames(emtx) <- unified_df$Gene_name
table(is.na(rownames(emtx)))
emtx <- emtx[!(is.na(rownames(emtx))),]

xgene <- ExpressionSet(assayData=emtx)
## Get time estimations

bozdechPath <- '/home/lucas/ISGlobal/Arrays/AnalisisR/Files2executeProgram/bozdech_Hb3_clean2.csv'
LemieuxFunctionsPath <- '/home/lucas/ISGlobal/Arrays/AnalisisR/Files2executeProgram/lemieux_et_al_pipeline_functions.r'

# Time recalibration function (Lemieux)

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
  oldTime <- rep(c("T0","T0","T1","T1","T2","T2","T3","T3"), 3)
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

# Time recalibration

outDir <- file.path("/home/lucas/ISGlobal/Arrays/AnalisisR/Eli_plots/")
figuresPath <- outDir

if (!file.exists(file.path(figuresPath,'timeAdjusting/'))) dir.create(file.path(figuresPath,'timeAdjusting/'))
estimatedTimes <- getTimeEstimation(xgene,bozdechPath,LemieuxFunctionsPath,file.path(figuresPath,'timeAdjusting/'),B=100)
save(estimatedTimes,file=file.path(dir,'estimatedTimes.RData'))
load(file.path(dir,'estimatedTimes.RData'))
pData(xgene)$time <- estimatedTimes
write.csv(estimatedTimes,file.path(figuresPath,'Estimated_Times.csv'))

##Plots

library(reshape2)
library(ggplot2)


for (i in 1:nrow(unified_df)){
  graf <- melt(unified_df[i,c(-1:-3)])
  graf["Time"] <- rep(estimatedTimes, 3)
  graf["Type"] <- rep(c("CTL", "HS"), 12)
  graf["Soque"] <- c(rep("s10E", 8), rep("s10G", 8), rep("EK0", 8))
  graf["group"] <- c(rep(c(1,2), 4), rep(c(3,4), 4), rep(c(5,6), 4))
  p <- ggplot(graf, aes(x = Time, y = value, col = Soque, linetype = Type, group = group))
  p <- p + geom_point(aes(color = Soque, shape = Soque)) + geom_line() + scale_linetype_manual(values=c("dashed", "solid"))
  ggsave(p, file=paste(outDir,unified_df$Gene_name[i],'.png',sep=''), width = 14, height = 10, units = "cm")
}


graf <- melt(unified_df[1,c(-1:-3)])
graf["Time"] <- rep(estimatedTimes)
graf["Type"] <- rep(c("CTL", "HS"), 12)
graf["Soque"] <- c(rep("s10E", 8), rep("s10G", 8), rep("EK0", 8))
graf["group"] <- c(rep(c(1,2), 4), rep(c(3,4), 4), rep(c(5,6), 4))
p <- ggplot(graf, aes(x = Time, y = value, col = Soque, linetype = Type, group = group))
p <- p + geom_point(aes(color = Soque, shape = Soque)) + geom_line() + scale_linetype_manual(values=c("dashed", "solid"))
p
#ggsave(p, file=paste(outDir,unified_df$Gene_name[i],'.png',sep=''), width = 14, height = 10, units = "cm")



unified_df[unified_df$Gene_name == "PF3D7_0103700",]




