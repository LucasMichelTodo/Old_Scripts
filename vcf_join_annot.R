library(dplyr)

vcf <- read.csv("/home/lucas/ISGlobal/Chip_Seq/DATA/Aligns/q5/Variant_Calling/a7_e5_all_recal_table.txt", header = TRUE,  sep ="\t", quote = "", row.names = NULL, stringsAsFactors = FALSE, comment.char = "#")
annot <- read.csv("/home/lucas/ISGlobal/Chip_Seq/DATA/Aligns/q5/Variant_Calling/annotated_a7e5_recall.txt", header = FALSE,  sep ="\t", quote = "", row.names = NULL, stringsAsFactors = FALSE, comment.char = "#")
annot["Chrom"] <- unlist(lapply(annot$V2, function(x) unlist(strsplit(x, ":"))[1]))
annot["Pos"] <- unlist(lapply(annot$V2, function(x) unlist(strsplit(x, ":"))[2]))
annot["Pos"] <- unlist(lapply(annot$Pos, function(x) unlist(strsplit(x, "-"))[1]))
annot$Pos <- as.numeric(annot$Pos)
out <- left_join(vcf, annot, by=c("Chrom","Pos"))

write.csv(out, file = "/home/lucas/ISGlobal/Chip_Seq/DATA/Aligns/q5/Variant_Calling/annotated_a7e5_recall_vcf.txt")
