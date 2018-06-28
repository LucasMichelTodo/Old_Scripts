blast_probes <- read.csv("/home/lucas/ISGlobal/Arrays/parsed_result_bits.csv", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
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


for (i in 1:dim(df)[1]){
  if (endsWith(df[i,"Hit_1"], ".1")){
    df[i,"Hit_1"] <- gsub('.{2}$', '', df[i,"Hit_1"])
  }
}

for (i in 1:dim(df)[1]){
  if (endsWith(df[i,"Hit_2"], ".1")){
    df[i,"Hit_2"] <- gsub('.{2}$', '', df[i,"Hit_2"])
  }
}

write.csv(df, file = "/home/lucas/ISGlobal/Arrays/probe_blast_table.csv", row.names = FALSE, quote = FALSE)

write.csv(df[is.na(df$Score_1),], file = "/home/lucas/ISGlobal/Arrays/probe_blast_nohit.csv", row.names = FALSE, quote = FALSE)



