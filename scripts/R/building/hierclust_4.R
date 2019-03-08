setwd("C:/Users/wbowers/Documents/tcga_replication_2/data")
set.seed(123.456)

# Clear out empty figure directories
for (item in list.dirs("../figs", recursive = FALSE)){
  if (length(list.files(item, recursive = TRUE)) == 0) {
    unlink(item, recursive = TRUE)
  } else {
    for (lower_item in list.dirs(item, recursive = FALSE)){
      if (length(list.files(lower_item, recursive = FALSE)) == 0){
        unlink(paste(lower_item, sep = ""), recursive = TRUE)
      }
    }
  }
}

remove(list = ls())
# Create directory if doesn't exist
date <- Sys.Date()

if (!file.exists(paste("../figs/", date, sep = ''))){
  dir.create(file.path(paste("../figs/", date ,sep = "")))
}

runpath <- paste("../figs/", date, "/", format(Sys.time(), "%H_%M_%S"), sep = "")
dir.create(file.path(runpath))

# Data have already been transformed, centred, and filtered

exp.data <- read.csv("TCGA_SARC_mrna_data.csv", row.names = 1, stringsAsFactors = FALSE)
exp.info <- read.csv("TCGA_SARC_mrna_info.csv", row.names = 1, stringsAsFactors = FALSE)

library(heatmap3)
library(RColorBrewer)

#Sort data by cluster
exp.info <- exp.info[,order(as.vector(as.numeric(exp.info[5,])))]
exp.data <- exp.data[,order(as.vector(as.numeric(exp.info[5,])))]

col.hmap <- colorRampPalette(brewer.pal(10, "RdBu"))(nrow(exp.data))

# Must be numeric factor
col.clust <- brewer.pal(3, "Set1")[as.vector(as.numeric(exp.info[5,]))]

# Factorise histological subtypes          
hist.fac <- factor(as.vector(as.character(exp.info[4,])),
                   levels = c("DDLPS", "UPS", "STLMS", "ULMS", "MPNST", "MFS", "SS"),
                   labels = c(1, 2, 3, 4, 5, 6, 7)
)

col.hist <- brewer.pal(8, "Dark2")[hist.fac]

sidecols <- cbind(Cluster=col.clust, HistType = col.hist)

# Table displaying colour of each phenotype
coltab <- unique(data.frame(as.vector(as.character(exp.info[4,])), col.hist))

# Draw heatmap

#
