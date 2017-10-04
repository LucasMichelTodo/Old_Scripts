
## FUNCIONS PER A PROJECTE ##

library(Rsamtools)

#### Seleccions "n" últims caràcters ####

substrR <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

#### Seleccionar tot menys "n" últims caràcters ####

substrR2 <- function(x, n){
  substr(x, 1, nchar(x)-n)
}

#### Contar ACGT i percentatge en un BAM ####

# - bam ha de ser un string amb el path al BAM file 

BAM_ACGT_count <- function (bam){
  raw_bam <- scanBam(bam)
  reads_seq <- paste0(unlist(raw_bam[[1]]$seq))
  ACGT_table <- as.data.frame(table(strsplit(reads_seq,"*")))
  ACGT_table$Frac <- as.numeric(sapply(ACGT_table$Freq, function(x) { x/sum(ACGT_table$Freq)*100}))
  ACGT_table
}

# selecciona els arxius de "dir" que acabin amb "suf"

select.files.end <- function(dir, suf){
  files <- list.files(dir)
  files <- files[substrR(files, nchar(suf)) == suf]
  return(files)
}


