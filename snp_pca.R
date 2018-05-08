snp_df <- read.csv("/home/lucas/ISGlobal/Chip_Seq/DATA/Aligns/q5/Variant_Calling/Run_Eliprotocol_12_02_18/Annotation/All_SNP_table_anotated_difs.txt", header = TRUE, sep = "\t")

head(iris)

library(ggfortify)
df <- as.data.frame(t(snp_df[,c(10:14)]))
samples <- rownames(df)
df["Sample"] <- samples
autoplot(prcomp(df[,-194]), data=df, colour="Sample")

