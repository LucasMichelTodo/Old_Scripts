
rosetta <- read.table("/home/lucas/ISGlobal/Gen_Referencies/Gene_references_rosetta.txt", fill = TRUE, fileEncoding = "UTF-8")
families <- read.table("/home/lucas/ISGlobal/Chip_Seq/Transcripció_CSV/Taula_families_simplificada.csv", fill = TRUE, fileEncoding = "UTF-8")
families <- families[,1:2]


for (i in 1:length(families[,1])){
  print(i)
  if (families[i,1] %in% rosetta[,3]){
    if (length(rosetta[rosetta[,3] %in% families[i,1],1]) > 1){
      for (a in 1:length(rosetta[rosetta[,3] %in% families[i,1],1])){
        families[i,2+a] <- rosetta[rosetta[,3] %in% families[i,1],1][a]
      } 
    } else {
        families[i,3] <- rosetta[rosetta[,3] %in% families[i,1],1]
    }
  } else if (families[i,1] %in% rosetta[,4]){
    if (length(rosetta[rosetta[,4] %in% families[i,1],1]) > 1){
      for (a in 1:length(rosetta[rosetta[,4] %in% families[i,1],1])){
        families[i,2+a] <- rosetta[rosetta[,4] %in% families[i,1],1][a]
      } 
    } else {
      families[i,3] <- rosetta[rosetta[,4] %in% families[i,1],1]
    }
  } else if (families[i,1] %in% rosetta[,5]){
    if (length(rosetta[rosetta[,5] %in% families[i,1],1]) > 1){
      for (a in 1:length(rosetta[rosetta[,5] %in% families[i,1],1])){
        families[i,2+a] <- rosetta[rosetta[,5] %in% families[i,1],1][a]
      } 
    } else {
      families[i,3] <- rosetta[rosetta[,5] %in% families[i,1],1]
    }
  }
}

write.table(families, file = "/home/lucas/ISGlobal/Gen_Referencies/Families.csv", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)


