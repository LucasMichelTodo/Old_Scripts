#! /usr/bin/Rscript

# Load gff and chromosome lengths
df <- read.csv("/home/lucas/ISGlobal/Gen_Referencies/Pdb41_sorted_filtered.gff", sep = "\t")
colnames(df) <-  c("chrom", "source", "feature", "start", "stop", "score", "strand", "frame", "attributes")
chrom_lengths <- read.table("/home/lucas/ISGlobal/Gen_Referencies/Pf3D7.sizes")
df["pre"] <- rep(0, dim(df)[1])
df["post"] <- rep(0, dim(df)[1])

# Main for loop. Must handle separately first,last, and rows where chromosome changes.
for (i in 1:length(df$chrom)){ 
  
  if (i == 1) { #Handle first row
    if (df[i,"start"] > 1000) {
      df[i,"start"] <- df[i,"start"]-1000
      df[i,"pre"] <- 1000
    } else {
        df[i,"pre"] <- df[i,"start"] 
        df[i,"start"] <- 1
    }}
  
  if (i < length(df$chrom)){ #Avoid breaking in last row
    if (df[i,"chrom"] == df[i+1,"chrom"]){ #Avoid comparing two different chromosomes when chromosome changes.
      dif <- df[i+1,"start"] - df[i,"stop"] 
      if ( dif < 1000 && dif > 0){
        df[i,"post"] <- dif
        df[i+1,"pre"] <- dif
        df[i,"stop"] <- df[i,"stop"] + dif
        df[i+1, "start"] <- df[i+1,"start"] - dif
      } else {
        df[i,"stop"] <- df[i,"stop"] + 1000
        df[i+1, "start"] <- df[i+1, "start"] -1000
        df[i,"post"] <- 1000
        df[i+1,"pre"] <- 1000
      }
      chrom <- df[]
    } else { #Handle rows where there is a chromosome change (last row of a chrom)
      if (df[i+1,"start"] > 1000) { #Handle start of next chromosome
        df[i+1,"start"] <- df[i+1,"start"]-1000
        df[i+1, "pre"] <- 1000
      } else {
          df[i+1, "pre"] <- df[i+1, "start"]
          df[i+1,"start"] <- 1
      }
      
      chrend <- chrom_lengths[chrom_lengths$V1 %in% df[i,"chrom"],2]
      if (df[i,"stop"] + 1000 > chrend) { #Handle end of current chromosome
        df[i,"post"] <- chrend - df[i,"stop"]
        df[i,"stop"] <- chrend
      } else {
          df[i,"stop"] <- df[i,"stop"] + 1000
          df[i,"post"] <- 1000
      } 
    }
  } else {
    chrend <- chrom_lengths[chrom_lengths$V1 %in% df[i,"chrom"],2]
    if (df[i,"stop"] + 1000 > chrend) { #Handle last row 
      df[i,"post"] <- chrend - df[i,"stop"]
      df[i,"stop"] <- chrend
    } else {
      df[i,"stop"] <- df[i,"stop"] + 1000
      df[i,"post"] <- 1000
    } 
  }              
}

write.table(df, file = "/home/lucas/ISGlobal/Gen_Referencies/Elongated_genes_3.gff", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
head(df, 50)

