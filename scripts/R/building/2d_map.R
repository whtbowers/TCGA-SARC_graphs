setwd("C:/Users/wbowers/Documents/tcga_replication/data")

# Create directory if doesn't exist
date <- Sys.Date()

if (!file.exists(paste("../figs/", date, sep = ''))){
  dir.create(file.path(paste("../figs/", date ,sep = "")))
}

runpath <- paste("../figs/", date, "/", format(Sys.time(), "%H_%M_%S"), sep = "")
dir.create(file.path(runpath)) 

exp.data <- read.csv("TCGA_SARC_data_raw.csv", row.names = 1, stringsAsFactors = FALSE)
info <- read.csv("TCGA_SARC_info_raw.csv", row.names = 1, stringsAsFactors = FALSE)

exp.data.1000 <- exp.data[1:1000,]

# #Get gene names alone, without pipe
# gnames <- c()
# for (i in 1:nrow(exp.data)){
#   gnames <- c(gnames, strsplit(rownames(exp.data)[i], "\\|")[[1]][1]) # Have to treat pipe as escape condition
# }
# Unfortunately have to keep this way due to repatition of gene names, symbols, etc. Maybe use when pared down to sammer number of genes.

# Try with all patients, then exclude later? Unclear where list of exclused patients is.
# Pairwise Spearman correlation
exp.cor <- cor(exp.data, method="spearman")
