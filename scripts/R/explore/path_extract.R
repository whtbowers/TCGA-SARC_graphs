setwd("C:/Users/wbowers/Documents/tcga_replication_2/data")
set.seed(123.456)

full.clin <- read.csv("SARC.merged_only_clinical_clin_format.csv", row.names = 1)

# row names
full.clin[,1]

# Find for all rows relating to pathology

path.ind <- c(11) # Preload with barcode row

for (i in 1:nrow(full.clin)) {
  if (grepl("pathology", as.character(full.clin[i,1]))) {
    path.ind <- c(path.ind, i)
  }
}

path.data <- full.clin[path.ind,]

# Set fields as rownames so not altered
rownames(path.data) <- as.vector(as.character(path.data[,1]))
path.data <- path.data[,-1]

# Set patient barcode as rowname
colnames(path.data) <- as.vector(sapply(path.data[1,], toupper))
path.data <- path.data[-1 ,]

# Remove rows >90% NA

for (i in rev(1:nrow(path.data))){
  if (sum(is.na(path.data[i,]))/ncol(path.data) > 0.9) {
    path.data <- path.data[-i,]
  }
}

write.csv(path.data, "sarc_clin_path_ext.csv")

