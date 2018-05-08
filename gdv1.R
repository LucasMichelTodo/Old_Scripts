library(dplyr)

gdv1_chip <- read.csv("/home/lucas/ISGlobal/gdv1_chip_Voss.csv", header = TRUE)

s10GvsA7 <- read.csv("/home/lucas/ISGlobal/Chip_Seq/DATA/Aligns/q5/Narrow_fe15/XLS_contrasts/Overlapped_and_filetred/calcFE/bed/annotated2/10G_A7K9_peaks_overlappandfilter_calcFE_annotated.csv", sep = "\t", header = TRUE)
s10GvsE5 <- read.csv("/home/lucas/ISGlobal/Chip_Seq/DATA/Aligns/q5/Narrow_fe15/XLS_contrasts/Overlapped_and_filetred/calcFE/bed/annotated2/10G_E5K9_peaks_overlappandfilter_calcFE_annotated.csv", sep = "\t", header = TRUE)
s12BvsA7 <- read.csv("/home/lucas/ISGlobal/Chip_Seq/DATA/Aligns/q5/Narrow_fe15/XLS_contrasts/Overlapped_and_filetred/calcFE/bed/annotated2/1.2B_A7K9_peaks_overlappandfilter_calcFE_annotated.csv", sep = "\t", header = TRUE)
s12BvsE5 <- read.csv("/home/lucas/ISGlobal/Chip_Seq/DATA/Aligns/q5/Narrow_fe15/XLS_contrasts/Overlapped_and_filetred/calcFE/bed/annotated2/1.2B_E5K9_peaks_overlappandfilter_calcFE_annotated.csv", sep = "\t", header = TRUE)

OFFvsON <- unique(c(as.character(s10GvsA7$Gene), as.character(s10GvsE5$Gene), as.character(s12BvsA7$Gene), as.character(s12BvsE5$Gene)))
OFF_10G <- intersect(as.character(s10GvsA7$Gene), as.character(s10GvsE5$Gene))
OFF_12B <- intersect(as.character(s12BvsA7$Gene), as.character(s12BvsE5$Gene))
OFF <- intersect(OFF_10G, OFF_12B)

sA7vs10G <- read.csv("/home/lucas/ISGlobal/Chip_Seq/DATA/Aligns/q5/Narrow_fe15/XLS_contrasts/Overlapped_and_filetred/calcFE/bed/annotated2/A7K9_10G_peaks_overlappandfilter_calcFE_annotated.csv", sep = "\t", header = TRUE)
sA7vs12B <- read.csv("/home/lucas/ISGlobal/Chip_Seq/DATA/Aligns/q5/Narrow_fe15/XLS_contrasts/Overlapped_and_filetred/calcFE/bed/annotated2/A7K9_1.2B_peaks_overlappandfilter_calcFE_annotated.csv", sep = "\t", header = TRUE)
sE5vs10G <- read.csv("/home/lucas/ISGlobal/Chip_Seq/DATA/Aligns/q5/Narrow_fe15/XLS_contrasts/Overlapped_and_filetred/calcFE/bed/annotated2/E5K9_10G_peaks_overlappandfilter_calcFE_annotated.csv", sep = "\t", header = TRUE)
sE5vs12B <- read.csv("/home/lucas/ISGlobal/Chip_Seq/DATA/Aligns/q5/Narrow_fe15/XLS_contrasts/Overlapped_and_filetred/calcFE/bed/annotated2/E5K9_1.2B_peaks_overlappandfilter_calcFE_annotated.csv", sep = "\t", header = TRUE)

ONvsOFF <- unique(c(as.character(sA7vs10G$Gene), as.character(sA7vs12B$Gene), as.character(sE5vs10G$Gene), as.character(sE5vs12B$Gene)))
ON_A7 <- intersect(as.character(sA7vs10G$Gene), as.character(sA7vs12B$Gene))
ON_E5 <- intersect(as.character(sE5vs10G$Gene), as.character(sE5vs12B$Gene))
ON <- intersect(ON_A7, ON_E5)


arrange(gdv1_chip, -abs(TP3.FC.ON.OFF.HP1..log2.))[1:100,]
gdv1_chip[gdv1_chip$gene.ID=="PF3D7_1222600",]

gdv1_chip["FC"] <- 2^gdv1_chip$TP3.FC.ON.OFF.HP1..log2.
gdv1_df <- gdv1_chip[,c(1:5,20)]
arrange(gdv1_df, -FC)

summary(gdv1_df[gdv1_df$gene.ID %in% OFF,]$FC)
summary(gdv1_df[gdv1_df$gene.ID %in% ON,]$FC)
