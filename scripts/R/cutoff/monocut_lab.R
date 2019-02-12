setwd("C:/Users/wbowers/Documents/tcga_replication/data")
set.seed(123.456)

library(tidyverse)

# Create directory if doesn't exist
date <- Sys.Date()

if (!file.exists(paste("../figs/", date, sep = ''))){
  dir.create(file.path(paste("../figs/", date ,sep = "")))
}

runpath <- paste("../figs/", date, "/", format(Sys.time(), "%H_%M_%S"), sep = "")
dir.create(file.path(runpath))

exp.data <- read.csv("TCGA_SARC_mrna_data_lnorm_medc_nosdfilt.csv", row.names = 1, stringsAsFactors = FALSE)
exp.info <- read.csv("TCGA_SARC_mrna_info_2.csv", stringsAsFactors = FALSE)
clin.data <- read.csv("TCGA_SARC_clinical.csv", stringsAsFactors = FALSE)

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


# Account for duplicate IDs in expression data
match.id <- exp.info[,which(exp.info[1,] %in% intersect(as.character(exp.info[1,]), as.character(clin.data$submitter_id)))]

exp.ord <- as.numeric(exp.goi[match(as.character(match.id[1,]), clin.data$submitter_id)])

all.ord <- exp.data[,match(as.character(match.id[1,]), clin.data$submitter_id)]
inf.ord <- exp.info[,match(as.character(match.id[1,]), clin.data$submitter_id)]

int.clin <- c()
for (i in 1:ncol(match.id)){
  int.clin <- c(int.clin, match(match.id[1,i], clin.data$submitter_id))
}

clin.data <- clin.data[int.clin,] # Clinical data with duplicate rows in order of expression data.

# Indices of id matches in clinical data to allow for duplicates

# exp.ord <- as.numeric(exp.goi[which(exp.info[1,] %in% clin.data$submitter_id)])
# 
# all.ord <- exp.data[,which(exp.info[1,] %in% clin.data$submitter_id)]
# inf.ord <- exp.info[,which(exp.info[1,] %in% clin.data$submitter_id)]

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
sum(exp.tersplit == 1, na.rm = TRUE) # Number of patients in lower tertile
sum(exp.tersplit == 2, na.rm = TRUE) # Number of patients in upper tertile

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
  sex = clin.data$gender,
  race = clin.data$race,
  patient.barcode = clin.data$submitter_id,
  patient.vital_status = death.bin,
  CCT2 = exp.ord,
  med.cutoff = exp.medsplit,
  ter.cutoff = exp.tersplit
)

# Dispay all graphically

# cutoffs <- list(
#   list(surv.inf$med.cutoff, "medcut", "median cutoff"),
#   list(surv.inf$ter.cutoff, "tercut", "tertile cutoff")
# )
# 
# endpoints <- list(
#   list(goi,  surv.inf$os.times, "os", "Overall Survival"),
#   list(goi, surv.inf$pfs.times, "pfs", "Progression-free Survival"),
#   list(goi,  surv.inf$os.times, "dss", "Disease-specific Survival")
# )

# Test cutoff plots for median and tertile cutoffs
library(survival)
library(survminer)

# Create survival based on chosen stat time and event

surv_object <- Surv(time = surv.inf$os.times, event = surv.inf$patient.vital_status)

# Fit survival curve (Stratify by chosen category)

fit1 <- survfit(surv_object ~ med.cutoff, data = surv.inf)

# Plot curve visually

ggsurv <- ggsurvplot(fit1, 
           data = surv.inf, 
           pval = TRUE,
           risk.table = TRUE,
           conf.int = TRUE
           ) 
ggsurv$plot <- ggsurv$plot +
  ggplot2::annotate("text",
                    x = 1200, y = 0.9, # x/y coords of text
                    label = paste("Upper: ", sum(exp.tersplit == 1, na.rm = TRUE),
                                  "\nLower:", sum(exp.tersplit == 2, na.rm = TRUE), sep = "")
                    )

print(ggsurv)

# Run multivariate cox

fit.coxph <- coxph(surv_object ~ med.cutoff + diag.cluster, 
                   data = surv.inf)
ggforest(fit.coxph, data = surv.inf)
