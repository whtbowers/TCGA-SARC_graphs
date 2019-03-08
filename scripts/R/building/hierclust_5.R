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

# Reduce number of gened for testing

# exp.data <- exp.data[1:500, ]

library(ConsensusClusterPlus)

#Cluster by patient

patient.con.clust <- ConsensusClusterPlus(d = as.matrix(exp.data),
                                  reps = 1000,
                                  maxK = 10,
                                  pItem = 1,
                                  pFeature = 0.8,
                                  clusterAlg = "hc",
                                  distance = "pearson",
                                  seed = 123.456, 
                                  title = paste(runpath, "/patient_clust", sep=""),
                                  plot = "png"
                                  )

# 3 clusters appears optimal

gene.con.clust <- ConsensusClusterPlus(d = as.matrix(t(exp.data)),
                                      reps = 1000,
                                      maxK = 10,
                                      pItem = 1,
                                      pFeature = 0.8,
                                      clusterAlg = "hc",
                                      distance = "pearson",
                                      seed = 123.456, 
                                      title = paste(runpath, "/gene_clust", sep=""),
                                      plot = "png"
                                      )

# 3 clusters appears optimal

# library(heatmap3)
library(RColorBrewer)

col.hmap <- colorRampPalette(brewer.pal(10, "RdBu"))(nrow(exp.data))

# Must be numeric factor
col.clust <- brewer.pal(3, "Set1")[as.vector(as.numeric(exp.info[5,]))]

# Factorise histological subtypes          
hist.fac <- factor(as.vector(as.character(exp.info[4,])),
                   levels = c("DDLPS", "UPS", "STLMS", "ULMS", "MPNST", "MFS", "SS"),
                   labels = c(1, 2, 3, 4, 5, 6, 7)
)

#Order by gene cluster 
# Do not bind, just order, but save table of which genes fell into which cluster
gene.clusts <- as.vector(gene.con.clust[[3]]$consensusClass) # 3 clusters
# gene.clusts <- as.vector(gene.con.clust[[4]]$consensusClass) # 4 clusters
gene.clusts.ind.ord <- order(gene.clusts)
gene.clusts.ord <- gene.clusts[order(gene.clusts)]
exp.data <- exp.data[gene.clusts.ind.ord,]

# Order by patient cluster
exp.info <- exp.info[,order(as.vector(as.numeric(exp.info[5,])))]
exp.data <- exp.data[,order(as.vector(as.numeric(exp.info[5,])))]

# Create colour vector

col.hist <- brewer.pal(8, "Dark2")[hist.fac]

col.gene <- brewer.pal(5, "Set2")[gene.clusts.ord]

topcols <- cbind(Cluster=col.clust, HistType = col.hist)

# Table displaying colour of each phenotype
coltab <- unique(data.frame(as.vector(as.character(exp.info[4,])), col.hist))

# ComplexHeatmap
library(ComplexHeatmap)
library(circlize)

# Create data frame for easier annotation
ann_data <- data.frame(
  HistType = as.vector(as.character(exp.info[4,])),
  TumorClust = as.vector(as.character(exp.info[5,]))
)

# Vectors of histype and associated colour
un.hist <- as.vector(as.character(coltab[, 1]))
un.hist.col <- as.vector(as.character(coltab[, 2]))

# Annotations

ha1 <- rowAnnotation(
  df = data.frame(GeneClust = as.character(gene.clusts.ord)),
  col = list(GeneClust = c("1" = "#ff03af", "2" = "#ffea03", "3" = "#5303ff"))
)

# ha1 <- rowAnnotation(
#   df = data.frame(GeneClust = as.character(gene.clusts.ord)),
#   col = list(GeneClust = c("1" = "#ff03af", "2" = "#ffea03", "3" = "#5303ff", "4" = "#02ff03"))
# )


ha2 <- HeatmapAnnotation(
  df = ann_data,
  col = list(
    HistType = c(
      "DDLPS" = "#1B9E77", 
      "UPS" = "#D95F02",
      "ULMS" = "#E7298A",
      "STLMS" = "#7570B3",
      "MPNST" = "#66A61E",
      "MFS" = "#E6AB02",
      "SS" = "#A6761D"),
    TumorClust = c(
      "1" = "#59ff00",
      "2" = "#ff5100",
      "3" = "#a200ff")
  )
) # TODO Must be way to automate
# If no colours specified, will autocolour, but can pick poor schemes

# png(paste(runpath, "/hmp_3clust.png" ,sep=""), width = 1000, height = 500, units = "px")
# png(paste(runpath, "/hmp_4clust.png" ,sep=""), width = 1000, height = 500, units = "px") png(paste(runpath, "/hmp.png" ,sep=""), width = 1000, height = 500, units = "px")

ha1 + Heatmap(
  matrix = exp.data,
  name = "Median centred mRNA expression",
  top_annotation = ha2,
  row_title = paste(nrow(exp.data), "genes"),
  row_title_side = "right",
  column_title = paste(ncol(exp.data), "patients"),
  column_title_side = "bottom",
  show_row_names = FALSE,
  show_column_names = FALSE,
  cluster_columns = FALSE,
  cluster_rows = FALSE
)

dev.off()

