library(tidyverse)

# Import TCGA data as headerless tibble
reads <- read_tsv(file = "tcga-sarc_FPKM-UQ.txt", col_names = FALSE)
colnames(reads)<- c("id", "fpkm")

# Query Ensembl for gene id
library(biomaRt)
ensembl <- useEnsembl(biomart="ensembl",
                      dataset = "hsapiens_gene_ensembl",
                      GRCh = 37)
# Check first 10 available attributes
#head(listAttributes(ensembl), 10)

# Select attributes
esbl_attr <- getBM(attributes = c("ensembl_gene_id_version",
                                  "hgnc_symbol"),
                   values = 1,
                   mart=ensembl)

# Check first 10 rows query output
head(esbl_attr)

#reads.head <-head(reads, 10)

# Pare down both data frames to just intersects - may seem a bit hacky, but ensure order maintained.

intersection <- intersect(reads$id, esbl_attr$ensembl_gene_id_version)

ind.reads <- which(intersection %in% reads$id)
ind.query <- which(intersection %in% esbl_attr$ensembl_gene_id_version)

mut.reads <- reads[ind.reads,]
mut.query <- esbl_attr[ind.query,]




check <- esbl_attr$ensembl_gene_id_version == query.success$id

fpkms <- reads[which(match.reads$id %in% query.success$ensembl_gene_id_version),]

head(data.frame(query.success$ensembl_gene_id_version, fpkms$id))


# for (read_id in head(reads$id, 10)) {
#   if (read_id %in% esbl_attr$ensembl_gene_id_version || read_id %in% esbl_attr$ensembl_gene_id){
#     id.ind <- 
#     # Add to query.success tibble
#     
#     if (query.success == 0){
#       query.success <- as.tibble(c(read_id))
#     }
#   } else {
#     
#     # Add to query.failure tibble
#     
#     print(0)
#   }
# }
