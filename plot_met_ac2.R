library(ggplot2)

cov_met <- read.csv("/home/lucas/ISGlobal/Chip_Seq/DATA/Aligns/q5/10G_me_sort_q5_fullcoverage.bed", header = FALSE,  sep ="\t", quote = "", row.names = NULL, stringsAsFactors = FALSE)
cov_ac <- read.csv("/home/lucas/ISGlobal/Chip_Seq/DATA/Aligns/q5/10G_ac_sort_q5_fullcoverage.bed", header = FALSE,  sep ="\t", quote = "", row.names = NULL, stringsAsFactors = FALSE)

ref <- read.csv("/home/lucas/ISGlobal/Gen_Referencies/Elongated_genes2.gff", header = FALSE,  sep ="\t", quote = "", row.names = NULL, stringsAsFactors = FALSE)

cov <- cbind(cov_met, cov_ac[,4])
colnames(cov) <- c("Chom", "Start", "Stop", "Met", "Ac")
cov <- cov[order(cov$Chom),]
rownames(cov) <- NULL 
cov["Type"] <- "intergenic"
dim(ref)

for (i in 1:dim(ref)[1]){
#for (i in 3200:3201){
  print(ref[i,])
  chr <- ref[i,1]
  print(chr)
  pre_start <- ref[i,4]
  pre_stop <- pre_start + ref[i,10]
  post_start <- ref[i,5] - ref[i,11]
  post_stop <- ref[i,5]
  
  # el problema es que retorna els idx dins dels subset (i després seleccionem els index al df sencer)
  if (ref[i,7] == "+") {
    cov[cov$Chom == chr,][which.min(abs(cov[cov$Chom == chr,2] - pre_start)): which.min(abs(cov[cov$Chom == chr,2] - pre_stop)),6] <- "5prima"
    cov[cov$Chom == chr,][which.min(abs(cov[cov$Chom == chr,2] - pre_stop)): which.min(abs(cov[cov$Chom == chr,2] - post_start)),6] <- "ORF"
    cov[cov$Chom == chr,][which.min(abs(cov[cov$Chom == chr,2] - post_start)): which.min(abs(cov[cov$Chom == chr,2] - post_stop)),6] <- "3prima"
  
  } else {
    cov[cov$Chom == chr,][which.min(abs(cov[cov$Chom == chr,2] - pre_start)): which.min(abs(cov[cov$Chom == chr,2] - pre_stop)),6] <- "3prima"
    cov[cov$Chom == chr,][which.min(abs(cov[cov$Chom == chr,2] - pre_stop)): which.min(abs(cov[cov$Chom == chr,2] - post_start)),6] <- "ORF"
    cov[cov$Chom == chr,][which.min(abs(cov[cov$Chom == chr,2] - post_start)): which.min(abs(cov[cov$Chom == chr,2] - post_stop)),6] <- "5prima"
    }
}

#which.min(abs(vect - value)) 

cod_cov <- cov[cov$Type != "intergenic",]

### Plots
ggplot(cod_cov, aes(x = cod_cov$Met, y = cod_cov$Ac, color = cod_cov$Type)) +
  geom_point(size=0.1)

ggplot(cov, aes(x = cov$Met, y = cov$Ac, color = cov$Type)) +
  scale_x_continuous(limits = c(0, 50000)) + scale_y_continuous(limits = c(0, 25000)) +
  geom_point(size=0.1)

ggplot(cov, aes(x = cov$Met, y = cov$Ac, color = cov$Type)) +
  scale_x_continuous(limits = c(0, 20000)) + scale_y_continuous(limits = c(0, 10000)) +
  geom_point(size=0.1)


