library(Rsamtools)
source(file = "~/R/Funcions_1.R")

# Carregar arxius
fls <- select.files.end("/media/lucas/Elements/malaria/data/bowtie_align_bam_bai", ".sort.bam")
fls <- "i10-10Gme.sort.bam"
inpath <- "/media/lucas/Elements/malaria/data/bowtie_align_bam_bai/"
outpath <- "~/ISGlobal/Chip_Seq/Aln_Unaln/"

# Crear BAM de aligned i unaligned
comand <- "samtools view -b -f -4 %in% > %out%"
for (x in fls){
  cmd <- comand
  cmd <- gsub("%in%", paste0(inpath, x), cmd)
  cmd <- gsub("%out%", paste0(outpath, gsub("sort.bam", "un.bam", x)), cmd)
  system(cmd)
}

comand <- "samtools view -b -F -4 %in% > %out%"
for (x in fls){
  cmd <- comand
  cmd <- gsub("%in%", paste0(inpath, x), cmd)
  cmd <- gsub("%out%", paste0(outpath, gsub("sort.bam", "al.bam", x)), cmd)
  system(cmd)
}

# Contar ACGT
un_fls <- select.files.end("~/ISGlobal/Chip_Seq/Aln_Unaln/", ".un.bam")
al_fls <- select.files.end("~/ISGlobal/Chip_Seq/Aln_Unaln/", ".al.bam")
BAM_ACGT_count(paste0(outpath, un_fls))
