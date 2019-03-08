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

valsplit <- ifelse(exp.goi.ord >= median(as.vector(as.numeric(exp.goi.ord))), 1, 2)

death.bin <- ifelse(clin.data.ord$days_to_death == '--', 0, 1)

ages <- as.vector(as.numeric(clin.data.ord$year_of_death)) - as.vector(as.numeric(clin.data.ord$year_of_birth))
age.split <- median((ages), na.rm = TRUE)
age.bin <- ifelse(ages < age.split, paste("Older than", age.split), paste(age.split, "or older"))

surv.inf <- data.frame(
  os.times = as.vector(as.numeric(clin.data.ord$days_to_death)),
  pfs.times = as.vector((as.numeric(clin.data.ord$days_to_death)) - as.vector(as.numeric(clin.data.ord$days_to_last_follow_up))),
  hist.cluster = as.vector(as.character(exp.info.ord[4,])),
  sex = as.vector(as.character(clin.data.ord$gender)),
  age = age.bin,
  patient.barcode = clin.data.ord$submitter_id,
  patient.vital_status = death.bin,
  CCT2 = as.vector(as.numeric(exp.goi.ord)),
  cutoff = as.vector(as.numeric(valsplit))
)

surv_object <- Surv(time = surv.inf$os.times, event = surv.inf$patient.vital_status)

fit.coxph <- coxph(surv_object ~ age + sex + hist.cluster,
                   data = surv.inf)

ggfor <- ggforest(fit.coxph, data = surv.inf, main = "Hazard ratios for Overall surival")

ggsave(paste(runpath, "/os_coxforplot.png", sep=""), width = 6, height = 4, units = "in",device = "png")

print("Cox PH for overall survival")
print(summary(fit.coxph))


