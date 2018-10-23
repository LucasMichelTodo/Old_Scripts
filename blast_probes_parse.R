library(xlsx)

################### Load files #########################################
blast_probes <- read.csv("/home/lucas/ISGlobal/Arrays/parsed_result_multitranscript.csv", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
df <- blast_probes[,c(1:5)]
colnames(df) <- c("Probe", "Hit_1", "Score_1", "Hit_2", "Score_2")

df["Target"] <- rep(NA, dim(df)[1])
df["Target_ID"] <- rep(NA, dim(df)[1])

df$Probe <- unlist(lapply(df$Probe, function(x) gsub(":", "/", x, fixed = TRUE)))

array_description <- read.csv("/home/lucas/ISGlobal/Arrays/array_description.txt", sep = "\t", as.is = TRUE, header = TRUE)
df[!df$Probe %in% array_description$ID,]

for (i in 1:dim(df)[1]){
  if (df[i,]$Probe %in% array_description$ID){
    df[i,"Target"] <- unique(array_description[array_description$ID == df[i,]$Probe,]$TopHit) 
  }
}

roseta <- read.csv("/home/lucas/ISGlobal/Gen_Referencies/Rosetta.txt", sep = "\t", as.is = TRUE, header = FALSE)

for (i in 1:dim(df)[1]){
  if (df[i,]$Target %in% roseta$V2){
    df[i,]$Target_ID <- capture.output(cat(roseta[roseta$V2 == df[i,]$Target,]$V1, sep = ", "))
  }
  if (df[i,]$Target %in% roseta$V3){
    df[i,]$Target_ID <- capture.output(cat(roseta[roseta$V3 == df[i,]$Target,]$V1, sep = ", "))
  }
  if (df[i,]$Target %in% roseta$V4){
    df[i,]$Target_ID <- capture.output(cat(roseta[roseta$V4 == df[i,]$Target,]$V1, sep = ", "))
  }
  if (df[i,]$Target %in% roseta$V5){
    df[i,]$Target_ID <- capture.output(cat(roseta[roseta$V5 == df[i,]$Target,]$V1, sep = ", "))
  }
}


kasfak <- read.csv(file = "/home/lucas/ISGlobal/Arrays/SupTable1_Kafsack_MalJ2012.csv", header = TRUE, sep = "\t")
kasfak$OligoID <- as.character(kasfak$OligoID)
kasfak$OligoID <- unlist(lapply(kasfak$OligoID, function(x) gsub(":", "/", x, fixed = TRUE)))


df["Kasfak"] <- rep(NA, dim(df)[1])

for (i in 1:dim(df)[1]){
  if (df[i,]$Probe %in% kasfak$OligoID){
    df[i,]$Kasfak <- as.character(kasfak[kasfak$OligoID == df[i,]$Probe,]$Uniqueness)
  }
}

df["Transcripts"] <- sapply(blast_probes$V6, function(x) gsub(",", ";", x))
write.csv(df[!df$Transcripts == "",], file = "/home/lucas/ISGlobal/Arrays/multitranscripts_table.csv", row.names = FALSE, quote = FALSE)

roseta <- read.csv("/home/lucas/ISGlobal/gens_array_rosetta_annotated.txt", sep = "\t", as.is = TRUE, header = FALSE)
anot <- c()
for (i in df$Target_ID){
  if (!is.na(i)){
    anot <- c(anot, roseta[!is.na(roseta$V2) & roseta$V2 == i,]$V3[1])
  } else {
    anot <- c(anot, NA)
  }
}

df["Anot"] <- anot
df$Target_ID <- sapply(df$Target_ID, function(x) gsub(",", ";", x))

################### Seleccions ###############################

# dif_hit <- apply(df, 1, function(x) !grepl(x["Hit_1"], x["Hit_2"]))
# not_unique <- df[dif_hit & df$Score_2 >= 80,]
# 
# #write.csv(not_unique[not_unique$Kasfak == "UNIQUE",], file = "/home/lucas/ISGlobal/Arrays/probe_blast_table_NOUNIQUE.csv", row.names = FALSE, quote = FALSE)
# 
# not_unique <- df[df$Score_2 >= 80,]
# 
# write.csv(df[df$Kasfak == "NOT_UNIQUE" & df$Score_2 < 80 & !is.na(df$Kasfak),], file = "/home/lucas/ISGlobal/Arrays/kasfak_notunique_goodhits", row.names = FALSE, quote = FALSE)
# df[df$Kasfak == "NOT_UNIQUE" & df$Score_2 < 80 & !is.na(df$Kasfak),]
# 
# wrong_target <- df[apply(df, 1, function(x) !grepl(x["Target_ID"], x["Hit_1"])) & !is.na(df$Target_ID) & df$Score_2 < 80,]
# bad_hit1 <- df[df$Score_1 < 80 & !df$Target %in% keep_list,]


################### Final Selection ######################

keep_list <- c("HA_tag", "BSD", "Cas9", "DiCre", "CsA", "DsRed", "eGFP_mut2", "eYFP", "DDdomain", 
               "FKBP", "Fluciferase", "GFP_pL6", "GFP", "HsDHFR", "mCherry", "Mp_rRNA", "ncRNA_thermreg", 
               "neoR", "Track_RNA", "tetR", "ScUPRT", "ScGlucanSynthase", "ScCD", "Rluciferase", "MAL1_18s")

keep_list <- c(keep_list, df[grepl("rRNA", df$Target),]$Target, df[grepl("tRNA", df$Target),]$Target)

a <- df[df$Score_1 < 80 & !df$Target %in% keep_list,]

keep_list <- c(keep_list, a[sapply(a$Target, function(x) grepl("TR", x)),]$Target)
keep_list <- c(keep_list, a[sapply(a$Target, function(x) grepl("Pfa", x)),]$Target)

keep_list2 <- c(a[sapply(a$Probe, function(x) grepl("RNAZ", x)),]$Probe)

a <- df[df$Score_1 < 80 & !df$Target %in% keep_list & !df$Probe %in% keep_list2,]

write.csv(a, file = "/home/lucas/ISGlobal/Arrays/bad_hit1.csv", row.names = FALSE, quote = FALSE)

keep_list3 <-  df[grepl("S ribosomal RNA", df$Anot),]$Anot


bad_probes <- df[(df$Score_1 < 100 | df$Score_2 > 80) & !df$Target %in% keep_list & !df$Probe %in% keep_list2 & !df$Anot %in% keep_list3,]
good_probes <- df[(df$Score_1 >= 100 & df$Score_2 <= 80) & !df$Target %in% keep_list & !df$Probe %in% keep_list2 & !df$Anot %in% keep_list3,]
exceptions_probes <- df[df$Target %in% keep_list | df$Probe %in% keep_list2 | df$Anot %in% keep_list3,]

################### Reanotate and write files ###########################

good_probes["New_Name"] <- sapply(good_probes$Hit_1, function(x) gsub("\\.[0-9]", "", (x)))
exceptions_probes["New_Name"] <- exceptions_probes$Target
exceptions_probes[!is.na(exceptions_probes$Anot),]$New_Name <- exceptions_probes[!is.na(exceptions_probes$Anot),]$Anot

prova <- good_probes[good_probes$Target_ID != good_probes$New_Anot & !is.na(good_probes$Target_ID),]
write.csv(prova, file = "/home/lucas/ISGlobal/Arrays/a_reanotar.csv", row.names = FALSE, quote = FALSE)

write.csv(bad_probes, file = "/home/lucas/ISGlobal/Arrays/bad_probes.csv", row.names = FALSE, quote = FALSE)
write.csv(good_probes, file = "/home/lucas/ISGlobal/Arrays/good_probes.csv", row.names = FALSE, quote = FALSE)
write.csv(exceptions_probes, file = "/home/lucas/ISGlobal/Arrays/excepcions_probes.csv", row.names = FALSE, quote = FALSE)
