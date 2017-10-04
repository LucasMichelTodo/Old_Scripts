library(readxl)
source(file = "~/ISGlobal/Scripts/Rosetta.R")

rosetta <- read.table("/home/lucas/ISGlobal/Gen_Referencies/Gene_references_rosetta.txt", fill = TRUE, fileEncoding = "UTF-8")
trans_df_raw <- read_excel("/home/lucas/ISGlobal/Chip_Seq/Transcripció_CSV/3D7_Variantome_AllData_withGam.xls", sheet = 2, range = "A1:U5268")

trans_df <- trans_df_raw[,c("X__1", "left.1.2b", "left.10g", "mid.1.2b", "mid.10g", "right.1.2b", "right.10g", "sides.1.2b", "sides.10g")]
colnames(trans_df)[1] <- "Gene"
trans_df <- as.data.frame(trans_df)
trans_df[,2:9] <- sapply(trans_df[,2:9], as.numeric)

for(column in trans_df[,2:9]){
  print(column[1])
}

max_dif <- c()
max_loc <- c()
for(i in 1:nrow(trans_df[,2:9])){
  row <- trans_df[,2:9][i,]
  l <- row[1]-row[2]
  m <- row[3]-row[4]
  r <- row[5]-row[6]
  s <- row[7]-row[8]
  max_dif[i] <- max(l,m,r,s)
  if (length(c("l", "m", "r", "s")[which.max(c(l,m,r,s))]) == 1){ #Evita que peti quan max = NA
    max_loc[i] <- c("l", "m", "r", "s")[which.max(c(l,m,r,s))]
  } else {
    max_loc[i] <- "NULL"
  }
  
}

trans_df["max_dif"] <- max_dif
trans_df["max_loc"] <- max_loc

trans_df[order(-max_dif),]
a <- rosetize(trans_df$Gene)

df_10G <- read.table("/home/lucas/ISGlobal/Chip_Seq/DATA/Aligns/q5/Narrow_fe30/Annotation_hits/Tables/10G_1.2B_annotation_table.csv", header = TRUE) 
df_1.2B <- read.table("/home/lucas/ISGlobal/Chip_Seq/DATA/Aligns/q5/Narrow_fe30/Annotation_hits/Tables/1.2B_10G_annotation_table.csv", header = TRUE) 
