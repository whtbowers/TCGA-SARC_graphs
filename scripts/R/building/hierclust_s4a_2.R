setwd("C:/Users/wbowers/Documents/tcga_replication/data")
set.seed(123.456)

library(tidyverse)

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

# Create directory if doesn't exist
date <- Sys.Date()

if (!file.exists(paste("../figs/", date, sep = ''))){
  dir.create(file.path(paste("../figs/", date ,sep = "")))
}

runpath <- paste("../figs/", date, "/", format(Sys.time(), "%H_%M_%S"), sep = "")
dir.create(file.path(runpath))

# Data have already been transformed, centred, and filtered

exp.data <- read.csv("TGCA_SARC_mrna_data.csv", row.names = 1, stringsAsFactors = FALSE)
exp.info <- read.csv("TGCA_SARC_mrna_info.csv", row.names = 1, stringsAsFactors = FALSE)


library(ConsensusClusterPlus)

con.clust <- ConsensusClusterPlus(d = as.matrix(exp.data),
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

clust3 <- con.clust[[3]]
clust4 <- con.clust[[4]]
clust5 <- con.clust[[5]]

library(RColorBrewer)

col.hmap <- colorRampPalette(brewer.pal(10, "RdBu"))(nrow(clust3$consensusMatrix))

# Must be numeric factor
col.clust <- brewer.pal(3, "Set1")[as.vector(as.numeric(clust3$consensusClass))]

as.vector(as.character(exp.info[4,]))

# Factorise histological subtypes          
hist.fac <- factor(as.vector(as.character(exp.info[4,])),
                   levels = c("DDLPS", "UPS", "STLMS", "ULMS", "MPNST", "MFS", "SS"),
                   labels = c(1, 2, 3, 4, 5, 6, 7)
)

col.hist <- brewer.pal(8, "Dark2")[hist.fac]
library(heatmap.plus)
library(heatmap3)

sidecols <- cbind(Cluster=col.clust, HistType = col.hist)

# Table displaying colour of each phenotype
coltab <- unique(data.frame(as.vector(as.character(exp.info[4,])), col.hist))

png(paste(runpath, "/hmp.png" ,sep="")) 
heatmap3(
  clust3$consensusMatrix,
  col = col.hmap,
  ColSideColors = sidecols,
  showColDendro = F,
  showRowDendro = F,
  legendfun = function()showLegend(legend=as.vector(as.character(coltab[,1])),col=as.vector(as.character(coltab[,2])),cex=1.5)
)
text(x=c(0.18, 0.40, 0.65), y = c(0.85, 0.85, 0.85), labels = c("C1", "C2", "C3"))

dev.off()

exp.info.withclust <- rbind(exp.info, clust3$consensusClass)

write.csv(exp.info.withclust, "TGCA_SARC_mrna_info.csv")
