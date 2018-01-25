library(ggplot2)
library(ggfortify)

####### Crear Taula ######

met_cov <- read.csv("/home/lucas/ISGlobal/Chip_Seq/DATA/Aligns/q5/Coverage/10G_met_cov.csv", header = TRUE,  sep ="\t", quote = "", row.names = NULL, stringsAsFactors = FALSE)
met_cov <- met_cov[,c(1:5)]
colnames(met_cov) <- c("ID", "met_ORF", "met_5", "met_3", "Chr")
ac_cov <- read.csv("/home/lucas/ISGlobal/Chip_Seq/DATA/Aligns/q5/Coverage/10G_ac_cov.csv", header = TRUE,  sep ="\t", quote = "", row.names = NULL, stringsAsFactors = FALSE)
ac_cov <- ac_cov[,c(1:5)]
colnames(ac_cov) <- c("ID", "ac_ORF", "ac_5", "ac_3", "Chr")

cov_df <- cbind(met_cov, ac_cov)
cov_df <- cov_df[,c(1:4,7:9,5)]

trans_df <- readWorksheetFromFile("/home/lucas/ISGlobal/Chip_Seq/Transcripció_CSV/3D7_Variantome_AllData_withGam.xls", sheet = 1)
trans_df <- trans_df[,c(1,26)]

rosetta <- read.table("/home/lucas/ISGlobal/Gen_Referencies/Gene_references_rosetta.txt", fill = TRUE, fileEncoding = "UTF-8")
for (i in 1:length(trans_df[,1])){ # Traduïr noms.
  if (trans_df[i,1] %in% rosetta[,3]){
    trans_df[i,"ID"] <- rosetta[rosetta[,3] %in% trans_df[i,1],1][1]
  } else if (trans_df[i,1] %in% rosetta[,4]){
    trans_df[i,"ID"] <- rosetta[rosetta[,4] %in% trans_df[i,1],1][1]
  } else if (trans_df[i,1] %in% rosetta[,5]){
    trans_df[i,"ID"] <- rosetta[rosetta[,5] %in% trans_df[i,1],1][1]
  }
}

ggplot(trans_df, aes(log(trans_df$Aver.2Higher10G.))) + geom_histogram(bins = 200)

met_df <- read.table("/home/lucas/ISGlobal/Chip_Seq/DATA/Aligns/q5/Coverage/10G_met_cov.csv", header = TRUE, sep = "\t", fileEncoding = "UTF-8")
met_df <- met_df[,c(1:5)]
ggplot(met_df, aes(log(met_df$Gene.cov))) + geom_histogram(bins = 200)


silenced <- met_df[log(met_df$Gene.cov) > 3,]$Gene
noexprs <- trans_df[log(trans_df$Aver.2Higher10G.) < 4,]$ID
table(silenced %in% noexprs)
table(noexprs %in% silenced)

cov_df["Silenced"] <- FALSE
cov_df[cov_df$ID %in% silenced,"Silenced"] <- TRUE

cov_df["Expressed"] <- TRUE
cov_df[cov_df$ID %in% noexprs,"Expressed"] <- FALSE

variant <- read.csv("/home/lucas/ISGlobal/Gen_Referencies/Gens_variants3.txt", header = FALSE,  sep ="\t", quote = "", row.names = NULL, stringsAsFactors = FALSE)
variant <- variant$V1

cov_df["Variant"] <- FALSE
cov_df[cov_df$ID %in% variant,"Variant"] <- TRUE

Trans2 <- read.table("/home/lucas/ISGlobal/Chip_Seq/Transcripció_CSV/Trans2.csv", sep="\t", header = TRUE)
Trans2["expressedin_10G"] <- Trans2$Trans < 0
difexp_10G <- Trans2[1:20, c(4:6)]
upexp <- difexp_10G[difexp_10G$expressedin_10G,]$ID
downexp <- difexp_10G[!difexp_10G$expressedin_10G,]$ID

cov_df["Dif_exp"] <- FALSE
cov_df[cov_df$ID %in% upexp,"Dif_exp"] <- "UP"
cov_df[cov_df$ID %in% downexp,"Dif_exp"] <- "DOWN"

####### Gràfics ######
cov_df["Class"] <- "Regular"
cov_df[cov_df$Variant == TRUE & cov_df$Silenced == TRUE,]$Class <- "Variant-Silenced"
cov_df[cov_df$Variant == TRUE & cov_df$Silenced == FALSE,]$Class <- "Variant-Active"
table(cov_df$Class)

df <- rbind(sample_n(cov_df[cov_df$Class == "Regular",], 528), cov_df[cov_df$Class != "Regular",]) #Sampling

ggplot(df, aes(x = df$met_5, y = df$ac_5, color = df$Class)) +
  scale_x_log10() + scale_y_log10() +
  geom_point(size=1)

##---------------------------##
df <- cov_df[cov_df$Dif_exp != FALSE,]
df["Class"] <- "Regular"
df[df$Variant == TRUE & df$Dif_exp =="UP",]$Class <- "Variant-Active"
df[df$Variant == TRUE & df$Dif_exp =="DOWN",]$Class <- "Variant-Inactive"
ggplot(df, aes(x = df$met_ORF, y = df$ac_ORF, color = df$Class)) +
  scale_x_log10() + scale_y_log10() +
  geom_point(size=0.2)

pca_df <- log(df[,c(2:7)]+0.01)
autoplot(prcomp(pca_df), data = df, colour = "Class")

df <- log(cov_df[,c(2:7)]+0.01)

autoplot(prcomp(df), data = cov_df, colour = "Class")
