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

# clinical data with pathology fields
clin.path <- read.csv("sarc_clin_path_ext.csv", row.names = 1, stringsAsFactors = FALSE)

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
clin.path.ord <- clin.path[,match(as.vector(as.character(exp.info.ord[1,])), gsub("[.]", "-", colnames(clin.path)))]
clin.liu.ord <- clin.liu[match(as.vector(as.character(exp.info.ord[1,])), clin.liu$bcr_patient_barcode),]
ends.data.ord <- ends.data[match(as.vector(as.character(exp.info.ord[1,])), ends.data$bcr_patient_barcode),]


# 3x split for cross validation, stratified by PFI survival and histology

library(caret)

# test w/just stratification by histological subtype first

train.ind <- createDataPartition(
  as.vector(as.character(exp.info.ord[4,])),
  times = 3,
  p = 1/3,
  list = FALSE
)

# Cross-val test cutpoint
for (i in 1:ncol(train.ind)) {
  
  exp.info.ord.test <- exp.info.ord[, train.ind[,i]]
  
}

exp.info.ord[, train.ind[,2]]

#Start w/decile splitpoint space; increase to centile later
cutpoints <- quantile(exp.goi.ord, seq(0.1, 0.9, by=0.1)) 

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
pvals <- c()

for (point in cutpoints){
  
  exp.split <- ifelse(exp.goi.ord >= point,
                      paste("High", goi, "expression"), 
                      paste("Low", goi, "expression")
                      
  )
  
  # Create df for easier graphing
  surv.inf<- data.frame(
    pfi.cen = ends.data.ord$PFI.1,
    pfi.times = ends.data.ord$PFI.time.1,
    exp.goi = exp.split,
    stage = clin.liu.ord$ajcc_pathologic_tumor_stage,
    sex = clin.liu.ord$gender,
    age = ages
  )
  
  surv.obj <- Surv(time = surv.inf$pfi.times, event = surv.inf$pfi.cen)
  surv.fit <- survfit(surv.obj ~ exp.goi, data = surv.inf) 
  
  # Log-rank to get chisq
  lrank <- survdiff(surv.obj ~ exp.goi, data = surv.inf)
  
  chisqs <- c(chisqs, lrank$chisq)
  pvals <- c(pvals, 1-pchisq(lrank$chisq, 1))
  
  ggsurv <- ggsurvplot(surv.fit,
                       data = surv.inf,
                       pval = TRUE,
                       #risk.table = TRUE,
                       conf.int = TRUE,
                       #log.rank.weights = "survdiff",
                       title = paste("Progression-free interval, split value = ", format(round(point, 2), nsmall = 2), ", \U1D6D8 =", format(round(lrank$chisq, 2), nsmall = 2), sep = "")
  )
  
  ggsave(file = paste(runpath,
                      "/test_os_survival_split_",
                      format(round(point, 2), nsmall = 2),
                      "_chisq_",
                      format(round(lrank$chisq, 2), nsmall = 2),
                      ".png",
                      sep=""),
         plot = print(ggsurv),
         width = 6, height = 4, units = "in",device = "png")
  
}

# Create melt table for better plotting (have to create manually in this case)

valtype <- c(rep("chisq", length(chisqs)), rep("pval", length(pvals)))
vals <- c(chisqs, pvals)
cuts <- c(cutpoints, cutpoints)

line.melt <- data.frame(metric = valtype, value = vals, cuts)

chiplot <- ggplot(data = line.melt) +
  geom_line(aes(x = cuts, y = value, color = metric), size = 2) +
  scale_color_manual(values = c("steelblue", "red")) +
  geom_hline(aes(yintercept = 0.05), col = "darkgray", linetype = 2, size = 2) +
  geom_vline(aes(xintercept = cutpoints[which(chisqs == max(chisqs))]), col = "green")

ggsave(file = paste(runpath, "/chi_pval_lines.png", sep=""), 
       plot = print(chiplot),
       width = 6, height = 4, units = "in",device = "png")

# 
# 
# 
# 
# # Create df for easier graphing
# surv.inf<- data.frame(
#   pfi.cen = ends.data.ord$PFI.1,
#   pfi.times = ends.data.ord$PFI.time.1,
#   exp.goi = exp.medsplit,
#   stage = clin.liu.ord$ajcc_pathologic_tumor_stage,
#   sex = clin.liu.ord$gender,
#   age = ages
# )
# 
# surv.obj <- Surv(time = surv.inf$pfi.times, event = surv.inf$pfi.cen)
# 
# surv.fit <- survfit(surv.obj ~ exp.goi, data = surv.inf) 
# 
# ggsurv <- ggsurvplot(fit = surv.fit,
#            data = surv.inf,
#            conf.int = TRUE,
#            risk.table = TRUE
#            )
# 
# # Cox's Proportional hazards
# coxph.fit <- coxph(surv.obj ~ age + sex + exp.goi,
#                    data = surv.inf
#                    )
# 
# ggforest(coxph.fit, 
#          data = surv.inf)
# ggcoxfunctional(fit = coxph.fit,
#                 data = surv.inf)
# ggcoxdiagnostics(fit = coxph.fit,
#                  data = surv.inf)
# 
