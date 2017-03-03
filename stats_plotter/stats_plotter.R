source(file = "~/R/Funcions_1.R")

mapq <- read.csv2(file = "~/ISGlobal/TestSet/align_tests/params_1/A7K9_14456_TTAGGC_MAPQ.csv", sep = "\t", header = FALSE)
lens <- read.csv2(file = "~/ISGlobal/TestSet/align_tests/params_1/A7K9_14456_TTAGGC_lengths.csv", sep = "\t", header = FALSE)

hist(mapq)
hist(abs(as.numeric(lens)))
max(as.numeric(lens))