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

# mRNA expression data (transformed and centred but not filtered for SD)
exp.data <- read.csv("TCGA_SARC_mrna_data_nosdfilt.csv", row.names = 1, stringsAsFactors = FALSE)

# mRNA expression metadata
exp.info <- read.csv("TCGA_SARC_mrna_info.csv", row.names = 1, stringsAsFactors = FALSE)

# survival data
ends.data <- read.csv("liu_endpoints_sarc.csv", row.names = 1, stringsAsFactors = FALSE)

# old clinical data
clin.tcga <- read.csv("TCGA_SARC_clinical.csv", stringsAsFactors = FALSE)

# Improved clinical data
clin.liu <- read.csv("liu_clin_sarc.csv", row.names = 1, stringsAsFactors = FALSE)

# Specify gene of interest
goi <- "CCT2"

# Isolate gene symbol from row label to find data relating to gene of interest only
gnames <- c()
for (i in 1:nrow(exp.data)){
  gnames <- c(gnames, strsplit(rownames(exp.data)[i], "\\|")[[1]][1])
}
exp.goi <- exp.data[match(goi, gnames),]


## Ensure expression, metainfo, survival and clinical data all in same order

exp.info.ord <- exp.info[,which(as.vector(as.character(exp.info[1,])) %in% clin.liu$bcr_patient_barcode)]
exp.goi.ord <- as.numeric(exp.goi[,which(as.vector(as.character(exp.info[1,])) %in% clin.liu$bcr_patient_barcode)])
clin.tcga.ord <- clin.tcga[match(as.vector(as.character(exp.info.ord[1,])), clin.tcga$submitter_id),]
clin.liu.ord <- clin.liu[match(as.vector(as.character(exp.info.ord[1,])), clin.liu$bcr_patient_barcode),]
ends.data.ord <- ends.data[match(as.vector(as.character(exp.info.ord[1,])), ends.data$bcr_patient_barcode),]

#Start w/decile splitpoint space; increase to centile later
cutpoints <- quantile(exp.goi.ord, seq(0.1, 0.9, by=0.1)) 

# Initially median split
exp.medsplit <- ifelse(exp.goi.ord >= median(exp.goi.ord),
                       paste("High", goi, "expression"), 
                       paste("Low", goi, "expression")
                       )

med.age <- median(clin.liu.ord$age_at_initial_pathologic_diagnosis)
ages <- ifelse(
  clin.liu.ord$age_at_initial_pathologic_diagnosis >= median(clin.liu.ord$age_at_initial_pathologic_diagnosis),
  paste("Older than", med.age, "at diagnosis"),
  paste("Younger than", med.age, "at diagnosis")
)





library(survival)
library(survminer)

#Collect chisquared values and associated pvals

chisqs <- c()
pvas <- c()

for (point in cutpoints){
  
}
# Create df for easier graphing
surv.inf<- data.frame(
  pfi.cen = ends.data.ord$PFI.1,
  pfi.times = ends.data.ord$PFI.time.1,
  exp.goi = exp.medsplit,
  stage = clin.liu.ord$ajcc_pathologic_tumor_stage,
  sex = clin.liu.ord$gender,
  age = ages
)

surv.obj <- Surv(time = surv.inf$pfi.times, event = surv.inf$pfi.cen)

surv.fit <- survfit(surv.obj ~ exp.goi, data = surv.inf) 

ggsurvplot(fit = surv.fit,
           data = surv.inf,
           conf.int = TRUE,
           risk.table = TRUE
           )

# Cox's Proportional hazards
coxph.fit <- coxph(surv.obj ~ age + sex + exp.goi,
                   data = surv.inf
                   )

ggforest(coxph.fit, 
         data = surv.inf)
ggcoxfunctional(fit = coxph.fit,
                data = surv.inf)
ggcoxdiagnostics(fit = coxph.fit,
                 data = surv.inf)

