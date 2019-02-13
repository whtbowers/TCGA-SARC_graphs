setwd("C:/Users/wbowers/Documents/tcga_replication/data")
set.seed(123.456)

library(tidyverse)

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

# Create directory if doesn't exist
date <- Sys.Date()

if (!file.exists(paste("../figs/", date, sep = ''))){
  dir.create(file.path(paste("../figs/", date ,sep = "")))
}

runpath <- paste("../figs/", date, "/", format(Sys.time(), "%H_%M_%S"), sep = "")
dir.create(file.path(runpath))

exp.data <- read.csv("TCGA_SARC_mrna_data_nosdfilt.csv", row.names = 1, stringsAsFactors = FALSE) # Run on non-sd filtered data as sd filter removes CCT2, which here is gene of interest
exp.info <- read.csv("TCGA_SARC_mrna_info.csv", row.names = 1, stringsAsFactors = FALSE)
clin.data <- read.csv("TCGA_SARC_clinical.csv", stringsAsFactors = FALSE)

# Specify gene of interest
goi <- "CCT2"

gnames <- c()
for (i in 1:nrow(exp.data)){
  gnames <- c(gnames, strsplit(rownames(exp.data)[i], "\\|")[[1]][1])
}

# Find row relating to gene of interest
exp.goi <- exp.data[match(goi, gnames),]


## Ensure expression levels same order as binarised censored death

exp.info.ord <- exp.info[,which(as.vector(as.character(exp.info[1,])) %in% clin.data$submitter_id)]

exp.goi.ord <- exp.goi[,which(as.vector(as.character(exp.info[1,])) %in% clin.data$submitter_id)]

clin.data.ord <- clin.data[match(as.vector(as.character(exp.info.ord[1,])), clin.data$submitter_id),]

# View(data.frame(as.vector(as.character(exp.info.ord[1,])), testord$submitter_id))

# Stratified Test/validation split

library(caret)

train.ind <- createDataPartition(
  as.vector(as.character(exp.info.ord[4,])),
  p = 0.7,
  list = FALSE
  )

train.data <- exp.goi.ord[,train.ind]
train.info <- exp.info.ord[,train.ind]
train.clin <- clin.data.ord[train.ind,]

test.data <- exp.goi.ord[,-train.ind]
test.info <- exp.info.ord[,-train.ind]
test.clin <- clin.data.ord[-train.ind,]

# table(as.vector(as.character(test.info[4,])))#/ncol(test.info)
# 
# table(as.vector(as.character(train.info[4,])))#/ncol(train.info)
# 
# table(as.vector(as.character(exp.info.ord[4,])))#/ncol(exp.info.ord)

# cutpoint percentiles

quantiles <- quantile(train.data, seq(0.1, 1, by=0.1)) #Start w/perdecile, then increase to percentile

# Binarise days to death for OS
death.bin <- ifelse(train.clin$days_to_death == '--', 0, 1)

# Separate expression values by median
exp.medsplit <- ifelse(train.data >= median(as.vector(as.numeric(train.data))), 1, 2)

# Create df for easier graphing
surv.inf <- data.frame(
  os.times = as.numeric(train.clin$days_to_death),
  pfs.times = as.vector((as.numeric(train.clin$days_to_death)) - as.vector(as.numeric(train.clin$days_to_last_follow_up))),
  diag.cluster = as.vector(as.character(train.info[4,])),
  sex = train.clin$gender,
  race = train.clin$race,
  patient.barcode = train.clin$submitter_id,
  patient.vital_status = death.bin,
  CCT2 = as.vector(as.numeric(train.data)),
  med.cutoff = as.vector(as.numeric(exp.medsplit))
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
           #risk.table = TRUE,
           conf.int = TRUE,
           log.rank.weights = "survdiff"
           )
ggsurv$plot <- ggsurv$plot +
  ggplot2::annotate("text",
                    x = 1200, y = 0.9, # x/y coords of text
                    label = paste("Upper: ", sum(exp.medsplit == 1, na.rm = TRUE),
                                  "\nLower:", sum(exp.medsplit == 2, na.rm = TRUE), sep = "")
                    )

print(ggsurv)

# Log-rank comparison of curves

lrank <- survdiff(surv_object ~ med.cutoff, data = surv.inf )

lrank$chisq


# 
# # Run multivariate cox
# 
# fit.coxph <- coxph(surv_object ~ med.cutoff + diag.cluster, 
#                    data = surv.inf)
# ggforest(fit.coxph, data = surv.inf)
