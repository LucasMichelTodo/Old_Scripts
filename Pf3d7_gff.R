
source("/home/lucas/ISGlobal/Scripts/Funcions_1.R")

gff <- read.table("/home/lucas/ISGlobal/Gen_Referencies/PlasmoDB-31_Pfalciparum3D7.gff", sep ="\t", quote = "") 
colnames(gff) <- c("Chrom", "DB", "Type", "Start", "Stop", "Score", "Strand", "Frame", "Annot")
head(gff)
split1 <- sapply(gff$Annot, function(x) strsplit(as.character(x), split = ";", fixed = TRUE))
split2 <- sapply(split1, function(x) strsplit(as.character(x), split = "=", fixed = TRUE))


# x <- c()
# for (element in split2){
#   for (item in element){
#     print(item[1])
#     x <- c(x,item[1])
#   }
# }
# 
# ano_kinds <- unique(x)

annotations <- c("ID","description", "Parent", "Ontology_term", "protein_source_id", "Note")

apply(gff[1:3,], 1, function(x) print(x))

a <- gff$Annot[6]

for (element in strsplit(as.character(a), split = ";", fixed = TRUE)){
  for(item in element){
    print(item)
    if (grepl("ID", item, fixed = TRUE)){
      id <- substring(item, 4)
      print(id)
    } else {
      print("nope")
  }}
}

get_id <- function(df, row){
  id <- df[row,] 
}

str(split2)

ano <- gff$Annot[1:10]

sapply(ano, function(x) strsplit(as.character(x), split = ";", fixed = TRUE))
strsplit(ano[1], split = ";", fixed = TRUE)
