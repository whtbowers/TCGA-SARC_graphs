setwd("C:/Users/wbowers/Documents/tcga_replication/data")

exp.data <- read.csv("TCGA_SARC_mrna_data_lnorm_medc_Amatch.csv", row.names = 1, stringsAsFactors = FALSE)

# Create time- and datestamped directories if don't already exist
date <- Sys.Date()

if (!file.exists(paste("../figs/", date, sep = ''))){
dir.create(file.path(paste("../figs/", date ,sep = "")))
}

runpath <- paste("../figs/", date, "/", format(Sys.time(), "%H_%M_%S"), sep = "")
dir.create(file.path(runpath))

library(ConsensusClusterPlus)

exp.clust <- ConsensusClusterPlus(d = as.matrix(exp.data),
                     reps = 1000,
                     maxK = 10,
                     pItem = 1,
                     pFeature = 0.8,
                     clusterAlg = "hc",
                     distance = "pearson",
                     seed = 123.456,
                     title = runpath,
                     plot = "png"
                     )

icl <- calcICL(exp.clust, title=runpath, plot = "png")
icl[["itemConsensus"]]
