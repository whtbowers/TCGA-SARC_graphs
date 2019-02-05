setwd("C:/Users/wbowers/Documents/tcga_replication/data")

# Create directory if doesn't exist
date <- Sys.Date()

if (!file.exists(paste("../figs/", date, sep = ''))){
  dir.create(file.path(paste("../figs/", date ,sep = "")))
}

runpath <- paste("../figs/", date, "/", format(Sys.time(), "%H_%M_%S"), sep = "")
dir.create(file.path(runpath)) 


exp.data <- read.csv("TCGA_SARC_data_raw.csv", row.names = 1)
info <- read.csv("TCGA_SARC_info_raw.csv", row.names = 1)