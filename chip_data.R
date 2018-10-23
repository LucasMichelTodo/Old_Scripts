library(ggfortify)
library(neuralnet)

#### Import data ####
data <- read.csv(file = "/home/lucas/ISGlobal/Chip_Seq/DATA/Aligns/q5/10G_me_cov_intervals.csv", header = TRUE, sep = "\t", row.names = 1)
data_ac <- read.csv(file = "/home/lucas/ISGlobal/Chip_Seq/DATA/Aligns/q5/10G_ac_cov_intervals.csv", header = TRUE, sep = "\t", row.names = 1)

variant <- read.table("/home/lucas/ISGlobal/Gen_Referencies/Gens_variants_gff2.txt", header = FALSE,  sep ="\t", quote = "", row.names = NULL, stringsAsFactors = FALSE)

#### Merge ####
# Create and scale an acetilation mtx and a methylation mtx. 
# Create two composite matrices, one dividing Ac/Met and one leaving one matx beside the other.


# Create training and test set
datatrain = data[ index, ]
datatest = data[ -index, ]


ac.mtx <- data_ac[,c(-9,-10)] + 1*10**-6
scaled_ac <- scale(ac.mtx, center = TRUE, scale = TRUE)

met.mtx <- data[,c(-9,-10)] + 1*10**-6
scaled_met <- scale(met.mtx, center = TRUE, scale = TRUE)

cancor(t(ac.mtx), t(met.mtx))
scaled_mtx <- scaled_ac/scaled_met

side_mtx <- cbind(scaled_ac, scaled_met)

df <- data.frame(scaled_mtx)
df$Chrom <- data$Chrom
df$Annotations <- data$Annotations
head(df)

side_df <- data.frame(side_mtx)
side_df$Chrom <- data$Chrom
side_df$Annotations <- data$Annotations
side_df$Variant <- rownames(side_df) %in% variant$V1
head(side_df)

# PCA
autoplot(prcomp(df[,1:8]))
autoplot(prcomp(ac.mtx))
autoplot(prcomp(met.mtx))
autoplot(prcomp(side_mtx))

#### Sampling data ####

samplesize = 0.60 * nrow(side_df)
set.seed(80)
index = sample( seq_len ( nrow ( side_df ) ), size = samplesize )

# Create training and test set
datatrain = side_df[ index, ]
datatest = side_df[ -index, ]

head(datatrain)

set.seed(2)
NN = neuralnet(Variant ~ Seg1+Seg2+Seg3+Seg4+Seg5+Seg6+Seg7+Seg8+Seg1.1+Seg2.1+Seg3.1+Seg4.1+Seg5.1+ Seg6.1+Seg7.1+Seg8.1, datatrain, hidden = 3 , act.fct = "logistic", linear.output = F )
NNmet = neuralnet(Variant ~ Seg1.1+Seg2.1+Seg3.1+Seg4.1+Seg5.1+ Seg6.1+Seg7.1+Seg8.1, datatrain, hidden = 3 , act.fct = "logistic", linear.output = F )

# plot neural network
plot(NN)

predict_testNN = compute(NN, datatest[,c(1:16)])

hist(predict_testNN$net.result)

datatest$Predicted <- predict_testNN$net.result > 0.5


table(datatest$Variant == datatest$Predicted)
table(datatest[datatest$Variant,]$Variant == datatest[datatest$Variant,]$Predicted)

# Hist
df$Total <- rowSums(df[,1:8], na.rm = TRUE)
hist(df$Total, breaks = 2000, xlim = c(-200,200))

het <- scaled_mtx[abs(df$Total) > 50,]
autoplot(prcomp(het[,1:8]))

# create heatmap and don't reorder columns
heatmap(het[,1:8], Colv = FALSE)
heatmap(as.matrix(side_df[side_df$Variant,1:16]), Colv = FALSE)

# Calculating totals
data$Total <- rowSums(data[,c(-9,-10)], na.rm = TRUE)
hist(data$Total, breaks = 2000, xlim = c(0,400))
hist(data$Total, breaks = 2000, xlim = c(0,50))

nomet <- data[data$Total < 40,]
met <- data[data$Total >= 40,]

autoplot(prcomp(nomet[,c(-9,-10,-11)]))
autoplot(prcomp(met[,c(-9,-10,-11)]))

# Calculating totals
df$Total <- rowSums(df[,c(-9,-10)], na.rm = TRUE)
hist(df$Total, breaks = 2000)
hist(df$Total, breaks = 2000, xlim = c(0,500))

nomet <- data[data$Total < 40,]
met <- data[data$Total >= 40,]

autoplot(prcomp(nomet[,c(-9,-10,-11)]))
autoplot(prcomp(met[,c(-9,-10,-11)]))


# scale data to mean=0, sd=1 and convert to matrix
met_scaled <- as.matrix(scale(met[,c(-9,-10, -11)]))
nomet_scaled <- as.matrix(scale(nomet[,c(-9,-10, -11)]))

# create heatmap and don't reorder columns
heatmap(met_scaled, Colv = F)


# Cluster rows
hc.rows <- hclust(dist(met_scaled))
plot(hc.rows)
