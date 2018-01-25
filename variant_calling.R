library(tibble)


var_df <- read.csv("/home/lucas/ISGlobal/Chip_Seq/DATA/Aligns/q5/Variant_Calling/filtered_results.table",  sep ="\t", quote = "", row.names = NULL, stringsAsFactors = FALSE)

colnames(var_df) <- c("Chr", "Pos", "Ref", "Alt", 
                      "1.2B_GT", "1.2B_AD", "1.2B_GQ", 
                      "10G_GT", "10G_AD", "10G_GQ", 
                      "A7_GT", "A7_AD", "A7_GQ", 
                      "C2_GT", "C2_AD", "C2_GQ", 
                      "E5_GT", "E5_AD", "E5_GQ")


get_freq <- function(x){
  freq <- c()
  for (i in x){
    if (is.na(i)){
      freq <- c(freq, NA)
    } else {
      vals <- unlist(strsplit(i, ","))
      freq <- c(freq, as.numeric(vals[2]) /  (as.numeric(vals[1]) + as.numeric(vals[2])))
    } 
  }
  freq[freq == Inf] <- 0
  return(freq)
}

var_df <- add_column(var_df, `1.2B_freq` = get_freq(var_df$`1.2B_AD`), .after = "1.2B_AD")
var_df <- add_column(var_df, `10G_freq` = get_freq(var_df$`10G_AD`), .after = "10G_AD")
var_df <- add_column(var_df, A7_freq = get_freq(var_df$A7_AD), .after = "A7_AD")
var_df <- add_column(var_df, C2_freq = get_freq(var_df$C2_AD), .after = "C2_AD")
var_df <- add_column(var_df, E5_freq = get_freq(var_df$E5_AD), .after = "E5_AD")

var_0.3 <- na.omit(var_df[var_df$`1.2B_freq` > 0.3 | 
                          var_df$`10G_freq` > 0.3 | 
                          var_df$A7_freq > 0.3 | 
                          var_df$C2_freq > 0.3 | 
                          var_df$E5_freq > 0.3, ])
dim(var_0.3)

######## Comparacions ########

## 1.2B vs 10G

na.omit(var_df[abs(var_df$`1.2B_freq` - var_df$`10G_freq`) > 0.3 &
         var_df$`1.2B_GQ` > 20 & var_df$`10G_GQ` > 20,])

var_df[is.na(var_df$`1.2B_freq`) & !is.na(var_df$`10G_freq`) & var_df$`10G_GQ` > 20 |
         is.na(var_df$`10G_freq`) & !is.na(var_df$`1.2B_freq`) & var_df$`1.2B_GQ` > 20,]






