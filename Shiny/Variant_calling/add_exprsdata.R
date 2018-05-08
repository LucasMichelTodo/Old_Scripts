snp_table <- read.csv("/home/lucas/ISGlobal/Chip_Seq/DATA/Aligns/q5/Variant_Calling/Run_Eliprotocol_12_02_18/Annotation/All_SNP_table_anotated_difs_vepannotated_merged.txt", header = TRUE,  sep ="\t", quote = "", row.names = NULL, stringsAsFactors = FALSE)
dif_exp <-  read.csv("/home/lucas/ISGlobal/Gen_Referencies/gens_diferencials_Rovira_Graells_rosetta.txt", header = FALSE,  sep ="\t", quote = "", row.names = NULL, stringsAsFactors = FALSE)
snp_table["Exp_dif"] <- snp_table$Gene %in% dif_exp$V2
snp_table[snp_table$Var_Type == ".",]$Var_Type <- "non-coding"
write.table(snp_table, "/home/lucas/ISGlobal/Chip_Seq/DATA/Aligns/q5/Variant_Calling/Run_Eliprotocol_12_02_18/Annotation/All_SNP_table_anotated_difs_vepannotated_merged_exprs.txt", quote = FALSE,  sep ="\t")
