setwd("C:/Users/wbowers/Documents/tcga_replication/data")

library(tidyverse)

# Create directory if doesn't exist
# date <- Sys.Date()
# 
# if (!file.exists(paste("../figs/", date, sep = ''))){
#   dir.create(file.path(paste("../figs/", date ,sep = "")))
# }
# 
# runpath <- paste("../figs/", date, "/", format(Sys.time(), "%H_%M_%S"), sep = "")
# dir.create(file.path(runpath))

exp.data <- read.csv("TCGA_SARC_mrna_data_lnorm_medc_nosdfilt.csv", row.names = 1, stringsAsFactors = FALSE)
exp.info <- read.csv("TCGA_SARC_mrna_info_2.csv", stringsAsFactors = FALSE)
clin.data <- read.csv("TCGA_SARC_clinical_Amatch.csv", stringsAsFactors = FALSE)

## Check if number of patients with difficult diagnoses same as number of excluded patients ##


# Specify gene of interest
goi <- "CCT2"

gnames <- c()
for (i in 1:nrow(exp.data)){
  gnames <- c(gnames, strsplit(rownames(exp.data)[i], "\\|")[[1]][1])
}

# Find row relating to gene of interest
exp.goi <- exp.data[match(goi, gnames),]

# Ensure expression levels same order as binarised censored death
exp.ord <- as.numeric(exp.goi[which(clin.data$submitter_id %in% exp.info[1,])]) ## Wrong but this way round for purpose of illustration

all.ord <- exp.data[,which(clin.data$submitter_id %in% exp.info[1,])] ## Wrong but this way round for purpose of illustration

# Binarise days to death for OS
death.bin <- ifelse(clin.data$days_to_death == '--', 0, 1)

# Separate expression values by median
exp.medsplit <- ifelse(exp.ord >= median(exp.ord), 1, 2)

# Separate expression values by tertile
exp.tertiles <- quantile(exp.ord, probs = seq(0, 1, 1/3))

exp.tersplit <- ifelse(
  exp.ord <= exp.tertiles[2], 1,
  ifelse((exp.ord <= exp.tertiles[3]), 2, NA)
)

# Recategorise tumor types for easier sorting
prim.diag <- clin.data$primary_diagnosis
prim.diag[(prim.diag == "Malignant fibrous histiocytoma") | (prim.diag == "Giant cell sarcoma") | (prim.diag == "Aggressive fibromatosis") | (prim.diag == "Liposarcoma, well differentiated") | (prim.diag == "Abdominal fibromatosis")] <- NA
prim.diag[(prim.diag == "Synovial sarcoma, NOS") | (prim.diag == "Synovial sarcoma, biphasic") | (prim.diag == "Synovial sarcoma, spindle cell")] <- "SS"
prim.diag[(prim.diag == "Leiomyosarcoma, NOS") & (clin.data$tissue_or_organ_of_origin == "Uterus, NOS")] <- "ULMS"
prim.diag[(prim.diag == "Undifferentiated sarcoma") | (prim.diag == "Pleomorphic liposarcoma")] <- "UPS"
prim.diag[(prim.diag == "Leiomyosarcoma, NOS") | (prim.diag == "Myxoid leiomyosarcoma")] <- "STLMS"
prim.diag[prim.diag == "Malignant peripheral nerve sheath tumor"] <- "MPNST"
prim.diag[prim.diag == "Dedifferentiated liposarcoma"] <- "DDLPS"
prim.diag[prim.diag == "Fibromyxosarcoma"] <- "MFS"


# Create df for easier graphing
surv.inf <- data.frame(
  os.times = as.numeric(clin.data$days_to_death),
  pfs.times = (as.numeric(clin.data$days_to_death) - as.numeric(clin.data$days_to_last_follow_up)),
  diag.cluster = prim.diag,
  patient.barcode = clin.data$submitter_id,
  patient.vital_status = death.bin,
  CCT2 = exp.ord,
  med.cutoff = exp.medsplit,
  ter.cutoff = exp.tersplit
)

#install.packages("pbdMPI")
#install.packages("kazaam")
# library(ggfortify)
# 
# autoplot(prcomp(all.ord), colour = surv.inf$diag.cluster)

# Dispay all graphically

cutoffs <- c(med.cutoff, ter.cutoff)

surv.cuts <- list(
  list( "CCT2", "CCT2",  "os.times", "os", "Overall Survival"),
  list("CCT2", "CCT2", "pfs.times", "pfs", "Progression-free Survival"),
  list("CCT2", c("CCT2", "diagnosis"),  "os.times", "dss", "Disease-specific Survival")
)

# Test cutoff plots for median and tertile cutoffs
library(survival)
library(survminer)

fit <- survfit(Surv(pfs.times, patient.vital_status) ~ ter.cutoff, data = surv.inf)
ggsurvplot(
  fit,
  pval = TRUE,
  conf.int = TRUE
  )
# for (cutoff in cutoffs){
#   print(as.character(cutoff))
# }
