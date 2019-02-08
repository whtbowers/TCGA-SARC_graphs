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
surv.inf <- data.frame(
  os.times = as.numeric(clin.data$days_to_death),
  pfs.times = as.numeric(clin.data$days_to_death) - as.numeric(clin.data$days_to_last_follow_up),
  diagnosis = clin.data$primary_diagnosis,
  patient.barcode = clin.data$submitter_id,
  patient.vital_status = death.bin,
  CCT2 = exp.ord
)

# Cutpoint for PFS
library(survminer)

os.cct2.cut <- surv_cutpoint(
  data = surv.inf,
  time = "os.times",
  event = "patient.vital_status",
  variables = "CCT2"
)

pfs.cct2.cut <- surv_cutpoint(
  data = surv.inf,
  time = "pfs.times",
  event = "patient.vital_status",
  variables = "CCT2"
)

dss.cct2.cut <- surv_cutpoint(
  data = surv.inf,
  time = "os.times",
  event = "patient.vital_status",
  variables = c("CCT2", "diagnosis")
)

surv.cuts <- list(
  list(os.cct2.cut, "CCT2", "CCT2",  "os.times", "os", "Overall Survival"),
  list(pfs.cct2.cut, "CCT2", "CCT2", "pfs.times", "pfs", "Progression-free Survival"),
  list(dss.cct2.cut, "CCT2", c("CCT2", "diagnosis"),  "os.times", "dss", "Disease-specific Survival")
)

library(survival)

# Cutoff and KM plots

for (i in 1:length(surv.cuts)){
  cut <- surv.cuts[[i]][[1]]
  split.variable <- surv.cuts[[i]][[2]]
  surv.variables <- surv.cuts[[i]][[3]]
  times <- surv.cuts[[i]][[4]]
  abrv <- surv.cuts[[i]][[5]]
  title <- surv.cuts[[i]][[6]]
  
  #Plot CCT2 cutpoint
  png(paste(runpath, "/", abrv, "_cct2_cutpoint.png", sep = ""))
  print(plot(cut, split.variable, pallette = "npg"))
  dev.off()
  
  cct2.cat <- surv_categorize(cut)
  
  
  # Create survival object and plot curves
  
  fit <- survfit(Surv(cct2.cat[,1], patient.vital_status) ~ sum(surv.variables), data = cct2.cat)
  
  survplot <- ggsurvplot(
    fit, 
    title = title,
    conf.int = TRUE,
    pval = TRUE)
  ggsave(paste(runpath, "/", abrv,"_cct2_survplot.png", sep = ""),
         plot = print(survplot))
  
  
  
  
}

