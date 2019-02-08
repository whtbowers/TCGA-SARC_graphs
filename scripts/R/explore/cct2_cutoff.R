setwd("C:/Users/wbowers/Documents/tcga_replication/data")

library(tidyverse)
# Create directory if doesn't exist
date <- Sys.Date()

if (!file.exists(paste("../figs/", date, sep = ''))){
  dir.create(file.path(paste("../figs/", date ,sep = "")))
}

runpath <- paste("../figs/", date, "/", format(Sys.time(), "%H_%M_%S"), sep = "")
dir.create(file.path(runpath))

exp.data <- read.csv("TCGA_SARC_mrna_data_lnorm_medc_nosdfilt_Amatch.csv", row.names = 1, stringsAsFactors = FALSE)
# exp.data <- read.csv("TCGA_SARC_mrna_data_Amatch.csv", row.names = 1, stringsAsFactors = FALSE)
exp.info <- read.csv("TCGA_SARC_mrna_info_Amatch.csv", row.names = 1, stringsAsFactors = FALSE)
clin.data <- read.csv("TCGA_SARC_clinical_Amatch.csv", row.names = 1, stringsAsFactors = FALSE)

# Specify gene of interest
goi <- "CCT2"

gnames <- c()
for (i in 1:length(rownames(exp.data))){
  gnames <- c(gnames, strsplit(rownames(exp.data)[i], "\\|")[[1]][1])
}

# Find row relating to gene of interest
exp.goi <- exp.data[match(goi, gnames),]

# Ensure expression levels same order as binarised censored death
exp.ord <- as.numeric(exp.goi[which(clin.data$submitter_id %in% exp.info[1,])])

# Binarise days to death for OS
death.bin <- ifelse(clin.data$days_to_death == '--', 0, 1)

# create df of relevent info
os.inf <- data.frame(
  times = as.numeric(clin.data$days_to_death),
  patient.barcode = clin.data$submitter_id,
  patient.vital_status = death.bin,
  CCT2 = exp.ord
)

# Find cutoff point
library(survminer)

os.cct2.cut <- surv_cutpoint(
  data = os.inf,
  time = "times",
  event = "patient.vital_status",
  variables = "CCT2"
)

summary(os.cct2.cut)

#Plot CCT2 cutpoint
png(paste(runpath, "/os_cct2_cutpoint.png", sep = ""))
plot(os.cct2.cut, "CCT2", pallette ="npg")
dev.off()

# Categorise CCT2 based upon cutpoint
os.cct2.cat <- surv_categorize(os.cct2.cut)

# Create survival object and plot curves
library(survival)

os.fit <- survfit(Surv(times, patient.vital_status) ~ CCT2, data = os.cct2.cat)

os.survplot <- ggsurvplot(
  os.fit, 
  title = "Overall Survival",
  conf.int = TRUE,
  pval = TRUE)
ggsave(paste(runpath, "/os_cct2_survplot.png", sep = ""),
       plot = print(os.survplot))
