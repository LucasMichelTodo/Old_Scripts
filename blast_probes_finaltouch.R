

good_probes <- read.csv("/home/lucas/ISGlobal/Arrays/Array_Annotation/good_probes_definitiu_withSeq_wtihAnnot.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)
bad_probes <- read.csv("/home/lucas/ISGlobal/Arrays/Array_Annotation/bad_probes_definitiu_withSeq.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)
exceptions_probes <- read.csv("/home/lucas/ISGlobal/Arrays/Array_Annotation/excepcions_probes_definitiu_withSeq_curated.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)

missing_probes <- read.csv("/home/lucas/ISGlobal/Arrays/missing_probes_parsed_result_multitranscript.csv", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
missing_list <- gsub(":", "/", missing_probes$V1)

good_probes[good_probes$Probe %in% missing_list,]$Probe <- sub('/', ':', good_probes[good_probes$Probe %in% missing_list,]$Probe)
bad_probes[bad_probes$Probe %in% missing_list,]$Probe <- sub('/', ':', bad_probes[bad_probes$Probe %in% missing_list,]$Probe)
exceptions_probes[exceptions_probes$Probe %in% missing_list,]$Probe <- sub('/', ':', exceptions_probes[exceptions_probes$Probe %in% missing_list,]$Probe)


write.csv(bad_probes, file = "/home/lucas/ISGlobal/Arrays/bad_probes_final.csv", row.names = FALSE, quote = FALSE)
write.csv(good_probes, file = "/home/lucas/ISGlobal/Arrays/good_probes_final.csv", row.names = FALSE, quote = FALSE)
write.csv(exceptions_probes, file = "/home/lucas/ISGlobal/Arrays/excepcions_probes_final.csv", row.names = FALSE, quote = FALSE)

