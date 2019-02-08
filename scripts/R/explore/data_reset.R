setwd("C:/Users/wbowers/Documents/tcga_replication/data")

#exp.data <- read.csv("TCGA_SARC_data_raw.csv", row.names = 1, stringsAsFactors = FALSE)
exp.info <- read.csv("TCGA_SARC_info_raw.csv", row.names = 1, stringsAsFactors = FALSE)

colnames(exp.info) <- exp.info[1,]

twelvid <- c() # 12-character shortened IDs
for (i in 1:length(colnames(exp.info))){
  twelvid <- c(twelvid, substr(colnames(exp.info)[i], 1, 12))
}

head(exp.info)[,1:5]
exp.info[1,] <- twelvid

write.csv(exp.info, "TCGA_SARC_mrna_info_2.csv", row.names = FALSE)
