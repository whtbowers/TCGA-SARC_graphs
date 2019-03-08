setwd("C:/Users/wbowers/Documents/tcga_replication/data")
set.seed(123.456)

# Prune empty figure directories from incomplete runs 
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

# Create dated and timestamped figure directories for easier sorting
date <- Sys.Date()

if (!file.exists(paste("../figs/", date, sep = ''))){
  dir.create(file.path(paste("../figs/", date ,sep = "")))
}

runpath <- paste("../figs/", date, "/", format(Sys.time(), "%H_%M_%S"), sep = "")
dir.create(file.path(runpath))

# Run on non-sd filtered data as sd filter removes CCT2, which here is gene of interest
exp.data <- read.csv("TCGA_SARC_mrna_data_nosdfilt.csv", row.names = 1, stringsAsFactors = FALSE)
exp.info <- read.csv("TCGA_SARC_mrna_info.csv", row.names = 1, stringsAsFactors = FALSE)
clin.data <- read.csv("TCGA_SARC_clinical.csv", stringsAsFactors = FALSE)

# Specify gene of interest
goi <- "CCT2"

# Isolate gene symbol from row label to find data relating to gene of interest only
gnames <- c()
for (i in 1:nrow(exp.data)){
  gnames <- c(gnames, strsplit(rownames(exp.data)[i], "\\|")[[1]][1])
}
exp.goi <- exp.data[match(goi, gnames),]


## Ensure expression, metainfo, and clinical data all in same order

exp.info.ord <- exp.info[,which(as.vector(as.character(exp.info[1,])) %in% clin.data$submitter_id)]
exp.goi.ord <- exp.goi[,which(as.vector(as.character(exp.info[1,])) %in% clin.data$submitter_id)]
clin.data.ord <- clin.data[match(as.vector(as.character(exp.info.ord[1,])), clin.data$submitter_id),]


# Stratified Test/validation split

library(caret)

train.ind <- createDataPartition(
  as.vector(as.character(exp.info.ord[4,])),
  p = 0.5, # Start with 0.5 split as per original X-Tile
  list = FALSE
)

train.data <- exp.goi.ord[,train.ind]
train.info <- exp.info.ord[,train.ind]
train.clin <- clin.data.ord[train.ind,]

test.data <- exp.goi.ord[,-train.ind]
test.info <- exp.info.ord[,-train.ind]
test.clin <- clin.data.ord[-train.ind,]

#Start w/decile splitpoint space; increase to centile later
sppoints <- quantile(train.data, seq(0.1, 0.9, by=0.1)) 

# Binarise mortality status
death.bin.train <- ifelse(train.clin$days_to_death == '--', 0, 1)
death.bin.test <- ifelse(test.clin$days_to_death == '--', 0, 1)

library(survival)
library(survminer)

# Collect chi squared values
## TODO: convert for loop to package

chisqs.train <- c()
for (point in sppoints) {
  
  #separate by quantile split values
  
  valsplit <- ifelse(train.data >= point, 1, 2)
  
  surv.inf <- data.frame(
    os.times = as.numeric(train.clin$days_to_death),
    pfs.times = as.vector((as.numeric(train.clin$days_to_death)) - as.vector(as.numeric(train.clin$days_to_last_follow_up))),
    diag.cluster = as.vector(as.character(train.info[4,])),
    sex = train.clin$gender,
    race = train.clin$race,
    patient.barcode = train.clin$submitter_id,
    patient.vital_status = death.bin.train,
    CCT2 = as.vector(as.numeric(train.data)),
    cutoff = as.vector(as.numeric(valsplit))
  )
  
  # Create survival based on chosen stat time and event
  
  surv_object <- Surv(time = surv.inf$os.times, event = surv.inf$patient.vital_status)
  
  # Fit survival curve (Stratify by chosen category)
  
  surv_fit <- survfit(surv_object ~ cutoff, data = surv.inf)
  
  # Log-rank to get chisq
  
  lrank <- survdiff(surv_object ~ cutoff, data = surv.inf )
  chisqs.train <- c(chisqs.train, lrank$chisq)
  
  ggsurv <- ggsurvplot(surv_fit,
                       data = surv.inf,
                       pval = TRUE,
                       #risk.table = TRUE,
                       conf.int = TRUE,
                       #log.rank.weights = "survdiff",
                       title = paste("Overall survival (train), split value = ", format(round(point, 2), nsmall = 2), ", \U1D6D8 =", format(round(lrank$chisq, 2), nsmall = 2), sep = "")
  )
  ggsurv$plot <- ggsurv$plot +
    ggplot2::annotate("text",
                      x = 1200, y = 0.9, # x/y coords of text
                      label = paste("Upper: ", sum(valsplit == 1, na.rm = TRUE),
                                    "\nLower: ", sum(valsplit == 2, na.rm = TRUE), sep = "")
                      
    )
  ggsave(paste(runpath, "/train_os_survival_split_", format(round(point, 2), nsmall = 2),"_chisq_", format(round(lrank$chisq, 2), nsmall = 2), ".png", sep=""), width = 6, height = 4, units = "in",device = "png")
  
}

max.chisq.train <- max(abs((chisqs.train)))
opt.split.train <- sppoints[which(chisqs.train == max.chisq.train)]
cat("From Training data:\nmax chisq: ", max.chisq.train, "\noptimal split value: ", format(round(as.numeric(opt.split.train), 2), nsmall = 2), sep = "")

# Use on test set and see if similar cutpoint
# Copypaste for now, create function/package when happy with functionality

chisqs.test <- c()
for (point in sppoints) {
  
  #separate by quantile split values
  
  valsplit <- ifelse(test.data >= point, 1, 2)
  
  surv.inf <- data.frame(
    os.times = as.numeric(test.clin$days_to_death),
    pfs.times = as.vector((as.numeric(test.clin$days_to_death)) - as.vector(as.numeric(test.clin$days_to_last_follow_up))),
    diag.cluster = as.vector(as.character(test.info[4,])),
    sex = test.clin$gender,
    race = test.clin$race,
    patient.barcode = test.clin$submitter_id,
    patient.vital_status = death.bin.test,
    CCT2 = as.vector(as.numeric(test.data)),
    cutoff = as.vector(as.numeric(valsplit))
  )
  
  # Create survival based on chosen stat time and event
  
  surv_object <- Surv(time = surv.inf$os.times, event = surv.inf$patient.vital_status)
  
  # Fit survival curve
  
  surv_fit <- survfit(surv_object ~ cutoff, data = surv.inf)
  
  # Log-rank to get chisq
  
  lrank <- survdiff(surv_object ~ cutoff, data = surv.inf )
  chisqs.test <- c(chisqs.test, lrank$chisq)
  
  ggsurv <- ggsurvplot(surv_fit,
                       data = surv.inf,
                       pval = TRUE,
                       #risk.table = TRUE,
                       conf.int = TRUE,
                       #log.rank.weights = "survdiff",
                       title = paste("Overall survival (test), split value = ", format(round(point, 2), nsmall = 2), ", \U1D6D8 =", format(round(lrank$chisq, 2), nsmall = 2), sep = "")
  )
  ggsurv$plot <- ggsurv$plot +
    ggplot2::annotate("text",
                      x = 1200, y = 0.9, # x/y coords of text
                      label = paste("Upper: ", sum(valsplit == 1, na.rm = TRUE),
                                    "\nLower: ", sum(valsplit == 2, na.rm = TRUE), sep = "")
                      
    )
  ggsave(paste(runpath, "/test_os_survival_split_", format(round(point, 2), nsmall = 2),"_chisq_", format(round(lrank$chisq, 2), nsmall = 2), ".png", sep=""), width = 6, height = 4, units = "in",device = "png")
  
}

max.chisq.test <- max(abs((chisqs.test)))
opt.split.test <- sppoints[which(chisqs.test == max.chisq.test)]
cat("From Test data:\nmax chisq: ", max.chisq.test, "\noptimal split value: ", format(round(as.numeric(opt.split.test), 2), nsmall = 2), sep = "")

if (opt.split.test == opt.split.train) {
  print("Agreement in split value between test and train sets")
} else {
  print("Test/train sets split values do not agree")
}
