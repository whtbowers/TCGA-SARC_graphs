setwd("C:/Users/wbowers/Documents/tcga_replication/data")

clin.data <- read.csv("TCGA_SARC_clinical.tsv", sep = "\t", stringsAsFactors = FALSE)
exp.info <- read.csv("TCGA_SARC_info_raw.csv", row.names = 1)

# Check first 12 digits of id in mRNA data against submitter id of clinical data

clin.id <- as.vector(exp.info[1,]) # Clinical ids from mRNA file

twelvid <- c() # 12-character shortened IDs
for (i in 1:length(clin.id)){
  twelvid <- c(twelvid, substr(toString(exp.info[1,i]), 1, 12))
}

clinid <- t(twelvid)
colnames(clinid) <- colnames(exp.info)
rbind(exp.info[1:2,], clinid, exp.info[-(1:2),])
exp.info <- rbind(clinid, exp.info)

# Ids which appear in mRNA and clinical dataset
id.match <- intersect(twelvid, clin.data$submitter_id)

# Remove rows from clinical data which aren't in intersect of IDs and arrange so matches order of id.match

clin.data.match <- clin.data[which(clin.data$submitter_id %in% id.match),]
clin.data.match <- clin.data.match[match(id.match, clin.data.match$submitter_id),]

sum(clin.data.match$submitter_id == id.match, na.rm = TRUE)/nrow(clin.data.match) # Check order correctly matched (Should equal 1)
# write.csv(clin.data.match, "TCGA_SARC_clinical_Amatch.csv")

# Remove rows from mRNA data which aren't in intersect of IDs and arrange so matches order of id.match

exp.data <- read.csv("TCGA_SARC_data_raw.csv", row.names = 1)

ind.exp.info <- c()
# exp.data.match <- exp.data[,(which(twelvid %in% id.match))]

for (i in 1:ncol(exp.info)){
  if (!(toString(exp.info[1,i]) %in% id.match)){
    ind.exp.info <- c(ind.exp.info, i)
  }
}

exp.info.match <- exp.info[,-ind.exp.info]
exp.data.match <- exp.data[,-ind.exp.info]

write.csv(exp.info.match, "TCGA_SARC_mrna_info_Amatch.csv")
write.csv(exp.data.match, "TCGA_SARC_mrna_data_Amatch.csv")

                               