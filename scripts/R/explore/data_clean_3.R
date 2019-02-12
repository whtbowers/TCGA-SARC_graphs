setwd("C:/Users/wbowers/Documents/tcga_replication/data")
set.seed(123.456)

exp.data <- read.csv("TCGA_SARC_data_raw.csv", row.names = 1, stringsAsFactors = FALSE)
exp.info <- read.csv("TCGA_SARC_info_raw.csv", row.names = 1, stringsAsFactors = FALSE)
clin.data <- read.csv("TCGA_SARC_clinical.csv", row.names = 1, stringsAsFactors = FALSE)

# Filter out all genes with < 90% nonzero expression
ind.90filt <- c()
for (i in 1:nrow(exp.data)){
  if ((sum(exp.data[i,] == 0, na.rm = TRUE)/ncol(exp.data))>=0.1){
    ind.90filt <- c(ind.90filt, i)
  } 
}

exp.data.90filt <- exp.data[-ind.90filt,]

# Check distribution of first gene

library(ggplot2)

ggplot() +
  geom_histogram(aes(x=as.numeric(exp.data.90filt[1,])))

#log2 transform add 0.05 to prevent -inf
exp.data.log2 <- log2(exp.data.90filt+0.05)

# median centre for each gene across all tumours
exp.data.c1 <- apply(exp.data.log2,2,function(x){
  x-median(x)
})

# median centre for each tumour across all genes
# Transpose otherwise automatically transposes
exp.data.c2 <- as.data.frame(t(apply(exp.data.c1,1,function(x) {
  x-median(x)
})))

# write.csv(exp.data.c2, "TCGA_SARC_mrna_data_lnorm_medc_nosdfilt_aln.csv")

# Check new distribution
ggplot() +
  geom_histogram(aes(x=as.numeric(exp.data.c2[1,])))

# Remove genes with std < 2
ind.std2filt <- c()
for (i in 1:nrow(exp.data.c2)){
  if(sd(exp.data.c2[i,]) < 2){
    ind.std2filt <- c(ind.std2filt, i)
  }
}

exp.data.sdfilt <- exp.data.c2[-ind.std2filt,]

dim(exp.data.sdfilt)

# Create submitter ID row for expression info

colnames(exp.info) <- exp.info[1,]

twelvid <- c() # 12-character shortened IDs
for (i in 1:length(colnames(exp.info))){
  twelvid <- c(twelvid, substr(colnames(exp.info)[i], 1, 12))
}

exp.info[1,] <- twelvid

# info cleanup
hist.abbrv <- as.character(as.vector(exp.info[2,]))

hist.abbrv[hist.abbrv == "Malignant Peripheral Nerve Sheath Tumors (MPNST)"] <- "MPNST"
hist.abbrv[(hist.abbrv == "Inflammatory \030MFH\031 / Undifferentiated pleomorphic sarcoma with prominent inflamm") | (hist.abbrv == "Undifferentiated Pleomorphic Sarcoma (UPS); NOS") | (hist.abbrv == "Pleomorphic MFH / Undifferentiated pleomorphic sarcoma") | (hist.abbrv == "Undifferentiated Pleomorphic Sarcoma (UPS)") | (hist.abbrv == "Pleomorphic \030MFH\031/ Undifferentiated pleomorphic sarcoma")] <- "UPS"
hist.abbrv[(hist.abbrv == "Synovial Sarcoma; Poorly differentiated") | (hist.abbrv == "Synovial Sarcoma; Monophasic") | (hist.abbrv == "Synovial Sarcoma; Biphasic")] <- "SS"
hist.abbrv[hist.abbrv == "Dedifferentiated liposarcoma"] <- "DDLPS"
hist.abbrv[hist.abbrv == "Leiomyosarcoma (LMS)"] <- "LMS"
hist.abbrv[hist.abbrv == "Myxofibrosarcoma"] <- "MFS"
hist.abbrv[(hist.abbrv == "0") | (hist.abbrv == "Desmoid Sarcoma") | (hist.abbrv == "Normal Tissue")] <- "other"
exp.info <- rbind(exp.info, hist.abbrv)
unique(as.vector(as.character(exp.info[4,])))

# Only change in this and next table should be LMS -> STLMS and ULMS

table(as.character(exp.info[4,]))

## Ensure STLMS and ULMS distinctly labelled ##

# All uterine tumours in clinical data turn out to be LMS
clin.ids <- clin.data$submitter_id[which(clin.data$tissue_or_organ_of_origin == "Uterus, NOS")]

# Change LMS to ULMS if uterine
exp.info[match(clin.ids, as.vector(as.character(exp.info[1,])))] <- "ULMS"

#Relabel other LMS as STLMS
exp.info[4, exp.info[4,] == "LMS"] <- "STLMS"

table(as.character(exp.info[4,]))

## Remove patients who have type 'other' ##
ncol(exp.info)

exp.data <- exp.data.sdfilt[,-which(exp.info[4,] == "other")]
exp.info <- exp.info[,-which(exp.info[4,] == "other")]

write.csv(exp.info, "TGCA_SARC_mrna_info")
write.csv(exp.data, "TGCA_SARC_mrna_data")