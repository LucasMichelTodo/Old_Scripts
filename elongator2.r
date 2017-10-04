#! /usr/bin/Rscript

# Load gff and chromosome lengths
df <- read.table("/home/lucas/ISGlobal/Gen_Referencies/PlasmoDB-31_Pfalciparum3D7_Sorted_filtered.gff", sep = "\t")
colnames(df) <-  c("chrom", "source", "feature", "start", "stop", "score", "strand", "frame", "attributes")
chrom_lengths <- read.table("/home/lucas/ISGlobal/Gen_Referencies/Pf3D7.sizes")
df["5prima"] <- rep(0, dim(df)[1])
df["3prima"] <- rep(0, dim(df)[1])

# Main for loop. Must handle separately first,last, and rows where chromosome changes.
for (i in 1:length(df$chrom)){
  if (i == 1) { #Handle first row
    if (df[i,"start"] > 1000) {
      df[i,"start"] <- df[i,"start"]-1000
      if (df[i,"strand"] == "+"){
        df[i,"5prima"] <- 1000
      } else {
        df[i,"3prima"] <- 1000
      }
    } else {
      if (df[i,"strand"] == "+"){
        df[i,"5prima"] <- df[i,"start"]
      } else {
        df[i,"3prima"] <- df[i,"start"]
      }
      df[i,"start"] <- 1
    }} 
  if (i < length(df$chrom)){ #Avoid breaking in last row
    if (df[i,"chrom"] == df[i+1,"chrom"]){ #Avoid comparing two different chromosomes when chromosome changes.
      dif <- df[i+1,"start"] - df[i,"stop"] 
      if ( dif < 1000 && dif > 0){
        if (df[i,"strand"] == "+"){
          df[i, "3prima"] <- dif
        } else {
          df[i, "5prima"] <- dif
        }
        if (df[i+1,"strand"] == "+"){
          df[i+1, "5prima"] <- dif
        } else {
          df[i+1, "3prima"] <- dif
        }
        df[i,"stop"] <- df[i,"stop"] + dif
        df[i+1, "start"] <- df[i+1,"start"] - dif
      } else {
        df[i,"stop"] <- df[i,"stop"] + 1000
        df[i+1, "start"] <- df[i+1, "start"] -1000
        if (df[i,"strand"] == "+"){
          df[i, "3prima"] <- 1000
        } else {
          df[i, "5prima"] <- 1000
        }
        if (df[i+1,"strand"] == "+"){
          df[i+1, "5prima"] <- 1000
        } else {
          df[i+1, "3prima"] <- 1000
        }
      }
      chrom <- df[]
    } else { #Handle rows where there is a chromosome change (last row of a chrom)
      if (df[i+1,"start"] > 1000) {
        df[i+1,"start"] <- df[i+1,"start"]-1000
        if (df[i+1,"strand"] == "+"){
          df[i+1,"5prima"] <- 1000
        } else {
          df[i+1,"3prima"] <- 1000
        }
        } else { #Handle start of next chromosome
          if (df[i+1,"strand"] == "+"){
            df[i+1,"5prima"] <- df[i+1,"start"]
          } else {
            df[i+1,"3prima"] <- df[i+1,"start"]
          }
          df[i+1,"start"] <- 1
        } 
      chrend <- chrom_lengths[chrom_lengths$V1 %in% df[i,"chrom"],2]
      if (df[i,"stop"] + 1000 > chrend) {
        if (df[i,"strand"] == "+"){
          df[i,"3prima"] <- chrend-df[i,"stop"]
        } else {
          df[i,"5prima"] <- chrend-df[i,"stop"]
        }
        df[i,"stop"] <- chrend
      } else { #Handle end of current chromosome
        df[i,"stop"] <- df[i,"stop"] + 1000
        if (df[i,"strand"] == "+"){
          df[i,"3prima"] <- 1000
        } else {
          df[i,"5prima"] <- 1000
        }
      } 
    }
  } else { #Handle last row 
    chrend <- chrom_lengths[chrom_lengths$V1 %in% df[i,"chrom"],2]
    if (df[i,"stop"] + 1000 > chrend) {
      if (df[i,"strand"] == "+"){
        df[i,"3prima"] <- chrend-df[i,"stop"]
      } else {
        df[i,"5prima"] <- chrend-df[i,"stop"]
      }
      df[i,"stop"] <- chrend
    } else {
      df[i,"stop"] <- df[i,"stop"] + 1000
      if (df[i,"strand"] == "+"){
        df[i,"3prima"] <- 1000
      } else {
        df[i,"5prima"] <- 1000
      }
    } 
  }              
}

write.table(df, file = "/home/lucas/ISGlobal/Gen_Referencies/Elongated_genes2.gff", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)


