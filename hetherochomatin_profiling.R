library(xlsx)

gam <- read.table("/home/lucas/ISGlobal/Chip_Seq/Hetherochomation_Voss/Gam_II_III_pf2004.bed", header = FALSE, sep = "\t")
schi <- read.table("/home/lucas/ISGlobal/Chip_Seq/Hetherochomation_Voss/Schizont_pf2004.bed", header = FALSE, sep = "\t")
schi <- read.table("/home/lucas/ISGlobal/Chip_Seq/Hetherochomation_Voss/Schizont_pf2004.bed", header = FALSE, sep = "\t")

df <- read.xlsx(file = "/home/lucas/ISGlobal/Chip_Seq/Hetherochomation_Voss/mmc7.xlsx", sheetIndex = 1)
df["Dif2"] <- (df$gam.II.II.HP1.ChIP_RPKM/df$gam.II.II.Input_RPKM) / (df$schizonts.HP1.ChIP_RPKM/df$schizonts.Input_RPKM)
head(df)
df[,c(1:3,29)]

write.table(df[,c(1:3,29)], file = "/home/lucas/ISGlobal/Chip_Seq/Hetherochomation_Voss/gam_over_schi.bed", row.names = FALSE, quote = FALSE, sep = "\t")
df["Dif3"] <- df$gam.II.III.z.score/df$schizonts.z.score
head(df)
plot(log2(df$Dif2), df$Dif3)
summary(df$Dif2)

log2(df$Dif2)
