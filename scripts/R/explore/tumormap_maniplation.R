setwd("C:/Users/wbowers/Documents/tcga_replication/data")
set.seed(123.456)

library(tidyverse)

exp.data <- read.csv("TCGA_SARC_mrna_data_lnorm_medc_nosdfilt.csv", row.names = 1, stringsAsFactors = FALSE)
exp.info <- read.csv("TCGA_SARC_mrna_info_2.csv", stringsAsFactors = FALSE)
clin.data <- read.csv("TCGA_SARC_clinical.csv", stringsAsFactors = FALSE)

# Simplify labels in info data


# Account for duplicate IDs in expression data
all.ord <- exp.data[,match(as.character(match.id[1,]), clin.data$submitter_id)]

match.id <- exp.info[,which(exp.info[1,] %in% intersect(as.character(exp.info[1,]), as.character(clin.data$submitter_id)))]

exp.ord <- as.numeric(exp.goi[match(as.character(match.id[1,]), clin.data$submitter_id)])

all.ord <- exp.data[,match(as.character(match.id[1,]), clin.data$submitter_id)]
inf.ord <- exp.info[,match(as.character(match.id[1,]), clin.data$submitter_id)]

int.clin <- c()
for (i in 1:ncol(match.id)){
  int.clin <- c(int.clin, match(match.id[1,i], clin.data$submitter_id))
}

clin.data <- clin.data[int.clin,]

# Recategorise tumor types for easier sorting
# prim.diag <- clin.data$primary_diagnosis
# prim.diag[(prim.diag == "Malignant fibrous histiocytoma") | (prim.diag == "Giant cell sarcoma") | (prim.diag == "Aggressive fibromatosis") | (prim.diag == "Liposarcoma, well differentiated") | (prim.diag == "Abdominal fibromatosis")] <- NA
# prim.diag[(prim.diag == "Synovial sarcoma, NOS") | (prim.diag == "Synovial sarcoma, biphasic") | (prim.diag == "Synovial sarcoma, spindle cell")] <- "SS"
# prim.diag[(prim.diag == "Leiomyosarcoma, NOS") & (clin.data$tissue_or_organ_of_origin == "Uterus, NOS")] <- "ULMS"
# prim.diag[(prim.diag == "Undifferentiated sarcoma") | (prim.diag == "Pleomorphic liposarcoma")] <- "UPS"
# prim.diag[(prim.diag == "Leiomyosarcoma, NOS") | (prim.diag == "Myxoid leiomyosarcoma")] <- "STLMS"
# prim.diag[prim.diag == "Malignant peripheral nerve sheath tumor"] <- "MPNST"
# prim.diag[prim.diag == "Dedifferentiated liposarcoma"] <- "DDLPS"
# prim.diag[prim.diag == "Fibromyxosarcoma"] <- "MFS"

interm <- cbind(clin.data, prim.diag)

# transpose so can get tumormap by gene and by patient

clin.data.t <- as.data.frame(t(clin.data))

all.ord.t <- as.data.frame(t(all.ord))

write.csv(inf.ord, "TCGA_SARC_mrna_info_2_aln.csv")
write.csv(all.ord, "TCGA_SARC_mrna_data_aln.csv")
write.csv(clin.data.t, "TCGA_SARC_clinical_aln.csv")

write.table(all.ord, "TCGA_SARC_mrna_data_aln.tsv", sep = "\t")
write.table(clin.data.t, "TCGA_SARC_clinical_aln.tsv", sep = "\t")
write.table(all.ord.t, "TCGA_SARC_mrna_data_aln_t.tsv", sep = "\t")
write.table(clin.data, "TCGA_SARC_clinical_aln_t.tsv", sep = "\t")
