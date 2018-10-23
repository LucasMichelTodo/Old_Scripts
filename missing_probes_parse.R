library(xlsx)

################### Load files #########################################

blast_probes <- read.csv("/home/lucas/ISGlobal/Arrays/missing_probes_parsed_result_multitranscript.csv", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
df <- blast_probes[,c(1:5)]
colnames(df) <- c("Probe", "Hit_1", "Score_1", "Hit_2", "Score_2")

df["Target"] <- rep(NA, dim(df)[1])
df["Target_ID"] <- rep(NA, dim(df)[1])

df$Probe <- unlist(lapply(df$Probe, function(x) gsub(":", "/", x, fixed = TRUE)))

array_description <- read.csv("/home/lucas/ISGlobal/Arrays/Array_Annotation/array_description.txt", sep = "\t", as.is = TRUE, header = TRUE)
df[!df$Probe %in% array_description$ID,]

for (i in 1:dim(df)[1]){
  if (sub('/', ':', df[i,]$Probe) %in% array_description$ID){
    df[i,"Target"] <- unique(array_description[array_description$ID == sub('/', ':', df[i,]$Probe),]$TopHit) 
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


kasfak <- read.csv(file = "/home/lucas/ISGlobal/Arrays/Array_Annotation/SupTable1_Kafsack_MalJ2012.csv", header = TRUE, sep = "\t")
kasfak$OligoID <- as.character(kasfak$OligoID)
kasfak$OligoID <- unlist(lapply(kasfak$OligoID, function(x) gsub(":", "/", x, fixed = TRUE)))


df["Kasfak"] <- rep(NA, dim(df)[1])

for (i in 1:dim(df)[1]){
  if (df[i,]$Probe %in% kasfak$OligoID){
    df[i,]$Kasfak <- as.character(kasfak[kasfak$OligoID == df[i,]$Probe,]$Uniqueness)
  }
}

df["Transcripts"] <- sapply(blast_probes$V6, function(x) gsub(",", ";", x))

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

save(df,file="/home/lucas/ISGlobal/Arrays/Array_Annotation/missing_probes_new_array_annotations.RData")

## Repartir en good, bad, i excepcions
load(file.path("/home/lucas/ISGlobal/Arrays/Array_Annotation/missing_probes_new_array_annotations.RData"))

good_probes <- df[(df$Score_1 >= 100 & df$Score_2 <= 80),]
df <- df[!df$Probe %in% good_probes$Probe,] 
exception_probes <- df[grepl("RNAZ", df$Probe) | grepl("PFC10_API.*RNA", df$Probe),]
bad_probes <- df[!df$Probe %in% exception_probes$Probe,]

## Anotar
good_probes["New_Name"] <- sapply(good_probes$Hit_1, function(x) gsub("\\.[0-9]", "", (x)))
exception_probes["New_Name"] <- exception_probes$Target
exception_probes[startsWith(exception_probes$Probe, "PFC10_API"),]$New_Name <- "PF3D7_API"
RNAZs <- exception_probes[startsWith(exception_probes$Probe, "RNAZ"),]$Probe
exception_probes[exception_probes$Probe %in% RNAZs,]$New_Name <- sapply(RNAZs, function(x) gsub("_v7.1_P1/1", "", (x)))

colnames(good_probes) <- c("Probe", "Hit_1", "Score_1", "Hit_2", "Score_2", "Target_Kasfak", "Target_Kasfak_New_ID", "Kasfak", "Transcripts", "Anot", "New_Target")
colnames(exception_probes) <- colnames(good_probes)
colnames(bad_probes) <- c("Probe", "Hit_1", "Score_1", "Hit_2", "Score_2", "Target_Kasfak", "Target_Kasfak_New_ID", "Kasfak", "Transcripts", "Anot")

good_probes$Probe <- sub('/', ':', good_probes$Probe)

## Escriure files
write.csv(bad_probes, file = "/home/lucas/ISGlobal/Arrays/missing_bad_probes_definitiu.csv", row.names = FALSE, quote = FALSE)
write.csv(good_probes, file = "/home/lucas/ISGlobal/Arrays/missing_good_probes_definitiu.csv", row.names = FALSE, quote = FALSE)
write.csv(exception_probes, file = "/home/lucas/ISGlobal/Arrays/missing_excepcions_probes_definitiu.csv", row.names = FALSE, quote = FALSE)





















