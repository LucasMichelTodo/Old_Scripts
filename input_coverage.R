input_10G <- read.csv("/home/lucas/ISGlobal/Chip_Seq/DATA/Aligns/q5/10G_in_cov.csv",  sep ="\t", quote = "", row.names = NULL, stringsAsFactors = FALSE)
input_1.2B <- read.csv("/home/lucas/ISGlobal/Chip_Seq/DATA/Aligns/q5/1.2B_in_cov.csv",  sep ="\t", quote = "", row.names = NULL, stringsAsFactors = FALSE)

summary(input_10G$Gene.cov)
summary(input_1.2B$Gene.cov)

hist(input_10G[input_10G$Gene.cov < 10,]$Gene.cov, breaks = 500)
hist(input_1.2B[input_1.2B$Gene.cov < 10,]$Gene.cov, breaks = 500)

del_10G <- input_10G[input_10G$Gene.cov < 1,]$Gene
del_1.2B <- input_1.2B[input_1.2B$Gene.cov < 1,]$Gene

setdiff(del_1.2B, del_10G)
setdiff(del_10G, del_1.2B)


input_diff <- input_10G[,2:4] - input_1.2B[,2:4]
input_diff["gene"] <- input_10G$Gene

summary(input_diff$Gene.cov)
hist(input_diff[input_diff$Gene.cov <= 5,]$Gene.cov, breaks = 500)
input_diff[abs(input_diff$Gene.cov) > 1.5,]


x <- input_cov[input_cov$Gene %in% dif_10G_cov$Gene,1:2]
write.csv(x, file = "/home/lucas/ISGlobal/input_cov.csv", quote = FALSE)
