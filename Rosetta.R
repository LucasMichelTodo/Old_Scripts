rosetta <- read.table("/home/lucas/ISGlobal/Gen_Referencies/Gene_references_rosetta.txt", fill = TRUE, fileEncoding = "UTF-8")

rosetize <- function(x){
  result <- c()
  for (i in 1:length(x)){
     if (x[i] %in% rosetta[,3]){
       result[i] <- rosetta[rosetta[,3] == x[i],1][1]
     } else if (x[i] %in% rosetta[,4]){
       result[i] <- rosetta[rosetta[,4] == x[i],1][1]
     } else if (x[i] %in% rosetta[,5]){
       result[i] <- rosetta[rosetta[,5] == x[i],1][1]
     }
  }
  return(result)
  }


rosetta[rosetta[,3] == "MAL1_28s",1][1]
