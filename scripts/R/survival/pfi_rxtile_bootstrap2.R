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

# Abeshouse clinical data
clin.abes <- read.csv("tcga_sarc_sample_features.csv", row.names = 1, stringsAsFactors = FALSE)

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

exp.info.ord <- exp.info[,which(as.vector(as.character(exp.info[1,])) %in% clin.liu$bcr_patient_barcode)] # Barcodes x
exp.goi.ord <- as.numeric(exp.goi[,which(as.vector(as.character(exp.info[1,])) %in% clin.liu$bcr_patient_barcode)]) # Barcodes x
clin.path.ord <- clin.path[,match(as.vector(as.character(exp.info.ord[1,])), gsub("[.]", "-", colnames(clin.path)))] # Barcodes x
clin.liu.ord <- clin.liu[match(as.vector(as.character(exp.info.ord[1,])), clin.liu$bcr_patient_barcode),] # Barcodes y
ends.data.ord <- ends.data[match(as.vector(as.character(exp.info.ord[1,])), ends.data$bcr_patient_barcode),] # Barcodes y

#Start w/decile splitpoint space; increase to centile later
cutpoints <- quantile(exp.goi.ord, seq(0.1, 0.9, by=0.1)) 

# 3x split for cross validation, stratified by PFI survival and histology

# library(sampling)
library(caret)
library(survival)
library(survminer)

# Sample by histype where patient survives
live.strat.ind <- createDataPartition(
  as.vector(as.character(exp.info.ord[4,]))[which(ends.data.ord$PFI.1 == 0)],
  p = 0.7,
  list = FALSE
)

live.strat <- exp.info.ord[which(ends.data.ord$PFI.1 == 0)][,live.strat.ind]

# Sample by histype where patient dies

death.strat.ind <- createDataPartition(
  as.vector(as.character(exp.info.ord[4,]))[which(ends.data.ord$PFI.1 == 1)],
  p = 0.7,
  list = FALSE
)

death.strat <- exp.info.ord[which(ends.data.ord$PFI.1 == 1)][,death.strat.ind]

# Combine vector of barcodes
strat.barcodes <- c(as.vector(as.character(live.strat[1,])), as.vector(as.character(death.strat[1,])))

# which to find indices 
train.ind <- which(as.vector(as.character(exp.info.ord[1,])) %in% strat.barcodes)

# Train split
exp.info.ord.train <- exp.info.ord[, train.ind]
exp.goi.ord.train <- exp.goi.ord[train.ind]
clin.path.ord.train <- clin.path.ord[, train.ind]
clin.liu.ord.train <- clin.liu.ord[train.ind,]
ends.data.ord.train <- ends.data.ord[train.ind,]

# Vaidation split
exp.info.ord.val <- exp.info.ord[, -train.ind]
exp.goi.ord.val <- exp.goi.ord[-train.ind]
clin.path.ord.val <- clin.path.ord[, -train.ind]
clin.liu.ord.val <- clin.liu.ord[-train.ind,]
ends.data.ord.val <- ends.data.ord[-train.ind,]

# test w/just stratification by histological subtype first
n.resample <- 100
resample.ind <- createResample(
  y = as.vector(as.character(exp.info.ord.train[4,])), 
  times = n.resample, 
  list = FALSE
                               )

# Initiate variable for graphing data frame

chisqsdata.list <- list()

# Collect max chisq for each resample and associated cutpoints

maxchisqs <- c()
cuts.maxchisq <- c()

# Collect number of resample values greater than or equal to the original sample to calculate pval

p.numerator <- 0

# (1+sum(exp.goi.ord >= exp.goi.ord[resample.ind[,2]]))/(1+length(exp.goi.ord))

# Bootstrap test cutpoint
for (i in 1:ncol(resample.ind)) {
  
  print(paste("Assessing optimal splitpoint for resample", i))
  
  # Generate resample sets of all data
  
  exp.info.ord.res <- exp.info.ord.train[, resample.ind[,i]]
  exp.goi.ord.res <- exp.goi.ord.train[resample.ind[,i]]
  clin.path.ord.res <- clin.path.ord.train[, resample.ind[,i]]
  clin.liu.ord.res <- clin.liu.ord.train[resample.ind[,i],]
  ends.data.ord.res <- ends.data.ord.train[resample.ind[,i],]
  
  # Get chisq value and associated 
  chisqs <- c()
  pvals <- c()
  
  for (point in cutpoints){
    
    exp.split <- ifelse(exp.goi.ord.res >= point,
                        paste("High", goi, "expression"), 
                        paste("Low", goi, "expression")
                        
    )
    
    # Bin age
    med.age <- median(clin.liu.ord.res$age_at_initial_pathologic_diagnosis)
    ages <- ifelse(
      clin.liu.ord.res$age_at_initial_pathologic_diagnosis >= median(clin.liu.ord.res$age_at_initial_pathologic_diagnosis),
      paste("Older than", med.age, "at diagnosis"),
      paste("Younger than", med.age, "at diagnosis")
    )
    
    # Create df for easier graphing
    surv.inf<- data.frame(
      pfi.cen = ends.data.ord.res$PFI.1,
      pfi.times = ends.data.ord.res$PFI.time.1,
      exp.goi = exp.split
    )
    
    # Fit survival objects
    surv.obj <- Surv(time = surv.inf$pfi.times, event = surv.inf$pfi.cen)
    surv.fit <- survfit(surv.obj ~ exp.goi, data = surv.inf)
    
    # Log-rank to get chisq
    lrank <- survdiff(surv.obj ~ exp.goi, data = surv.inf)
    chisqs <- c(chisqs, lrank$chisq)
    pvals <- c(pvals, 1-pchisq(lrank$chisq, 1))
  
  }
  
  #Add max chisq and associated cutpoints to arrays
  maxchisq <- max(chisqs)
  maxchisqs <- c(maxchisqs, maxchisq)
  cut.maxchisq <- as.numeric(cutpoints[which(chisqs == maxchisq)])
  cuts.maxchisq <- c(cuts.maxchisq, cut.maxchisq)
  
  # Collect values to build melt table
  chisqsdata.list[[i]] <- data.frame(
    chisqs = chisqs,
    pvals = pvals,
    resample = rep(as.character(i),length(chisqs)),
    cuts = cutpoints,
    maxchisq = rep(maxchisq,length(chisqs)),
    cut.maxchisq = rep(cut.maxchisq,length(chisqs))
    
                             )

}

chisqs.frame <- do.call(rbind, chisqsdata.list)

## Show distribution of chi-squared 95% CI
ind.ord <- (1:length(cuts.maxchisq))[order(cuts.maxchisq)] # Maintain order of indices
cuts.maxchisq.ord <- cuts.maxchisq[order(cuts.maxchisq)]
limcut <- round(quantile(1:length(cuts.maxchisq), c(0.05, 0.95)))
discard <- rep("Within 95% CI", length(cuts.maxchisq))
discard[c(1:(limcut[1]-1), (limcut[2]+1):length(cuts.maxchisq.ord))] <- "Outside 95% CI"

# 95% ci region
cuts.maxchisq.95ci <- cuts.maxchisq.ord[which(discard == "Within 95% CI")]

# Display distribution of cutpoints
dist.cuts <- ggplot()+
  geom_histogram(aes(x = cuts.maxchisq.ord), fill = "steelblue", bins = 30) +
  xlab("Cutoff value") +
  ylab("Count") +
  labs(title = "Maximium chi-squared-associated cutoff value distribution") + 
  scale_fill_discrete(name = "Cutoff")

ggsave(file = paste(runpath, "/pfi_max_cut_dist.png", sep=""), 
       plot = print(dist.cuts),
       width = 6, height = 4, units = "in",device = "png")

chiplot <- ggplot(data = chisqs.frame) +
    geom_line(aes(x = cuts, y = chisqs, color = resample), size = 1.5) +
    geom_vline(aes(xintercept = median(cuts.maxchisq.95ci)), col = "darkblue", linetype = 2, size = 1, alpha = 0.7) +
    geom_vline(aes(xintercept = min(cuts.maxchisq.95ci)), col = "darkred", linetype = 2, size = 1, alpha = 0.7) +
    geom_vline(aes(xintercept = max(cuts.maxchisq.95ci)), col = "darkred", linetype = 2, size = 1, alpha = 0.7) +
  labs(title = paste0(n.resample, "-fold resampled chi-squared values for survival comparison at varying high/low cutoff for PFI")) +
  xlab("Transformed cutoff value") +
  ylab("Chi-squared") +
  theme(text = element_text(size=15))+
  guides(color = FALSE)
  
ggsave(file = paste(runpath, "/pfi_chi_lines.png", sep=""), 
       plot = print(chiplot),
       width = 12, height = 8, units = "in",device = "png")


# Test cutpoint on training set

opt.cut <- median(cuts.maxchisq.95ci)

opt.cut.train.split <- ifelse(exp.goi.ord.train >= opt.cut,
                            paste("High", goi, "expression"), 
                            paste("Low", goi, "expression")
)

surv.inf.train <- data.frame(
  pfi.cen = ends.data.ord.train$PFI.1,
  pfi.times = ends.data.ord.train$PFI.time.1,
  exp.goi = opt.cut.train.split
)
  
# Create survival objects
surv.obj.train <- Surv(time = surv.inf.train$pfi.times, event = surv.inf.train$pfi.cen)

# Fit and plot single and paired PFI survival curves
surv.fit.train <- survfit(surv.obj.train ~ exp.goi, data = surv.inf.train)
surv.fit.train.sgl <- survfit(surv.obj.train ~ 1, data = surv.inf.train)

ggsurv.train.sgl <- ggsurvplot(fit = surv.fit.train.sgl,
                           data = surv.inf.train,
                           conf.int = TRUE,
                           risk.table = FALSE,
                           pval = TRUE,
                           title = "PFI training set event curve"
)

print(ggsurv.train.sgl)

ggsave(file = paste0(runpath,"/pfi_train_survival_sgl.png"),
       plot = print(ggsurv.train.sgl),
       width = 12, height = 8, units = "in",device = "png")

ggsurv.train <- ggsurvplot(fit = surv.fit.train,
                     data = surv.inf.train,
                     conf.int = TRUE,
                     risk.table = FALSE,
                     pval = TRUE
)

print(ggsurv.train)

######################
##### Validation #####
######################

# Create high/low split using selected optimal cutpoint

opt.cut <- median(cuts.maxchisq.95ci)

opt.cut.val.split <- ifelse(exp.goi.ord.val >= opt.cut,
                    paste("High", goi, "expression"), 
                    paste("Low", goi, "expression")
)

sel.exp.cutoff <- ggplot() +
  geom_histogram(aes(x = exp.goi.ord.val, fill = opt.cut.val.split), bins = 50) +
  xlab(paste("Transformed", goi,"expression")) +
  ylab("Patient count") +
  labs(fill = paste(goi, "expression level"), title = paste0("Selected ", goi, " expression cutoff for PFI"))

ggsave(file = paste0(runpath,"/pfi_val_sel_goi_cutoff.png"),
       plot = print(sel.exp.cutoff),
       width = 12, height = 8, units = "in",device = "png")

#Bin ages in vaidation set
med.age <- median(clin.liu.ord.val$age_at_initial_pathologic_diagnosis)
ages <- ifelse(
  clin.liu.ord.val$age_at_initial_pathologic_diagnosis >= median(clin.liu.ord.val$age_at_initial_pathologic_diagnosis),
  paste("Older than", med.age, "at diagnosis"),
  paste("Younger than", med.age, "at diagnosis")
)
 
# Create df for easier graphing
surv.inf.val <- data.frame(
  pfi.cen = ends.data.ord.val$PFI.1,
  pfi.times = ends.data.ord.val$PFI.time.1,
  exp.goi = opt.cut.val.split,
  stage = clin.liu.ord.val$ajcc_pathologic_tumor_stage,
  sex = clin.liu.ord.val$gender,
  age = ages,
  tumor.width = ifelse(as.vector(as.numeric(clin.path.ord.val["patient.primary_pathology.tumor_sizes.tumor_size.pathologic_tumor_width",])) < median(as.vector(as.numeric(clin.path.ord.val["patient.primary_pathology.tumor_sizes.tumor_size.pathologic_tumor_width",])), na.rm = TRUE), "lower than median width", "higher than median width"),
  tumor.length = ifelse(as.vector(as.numeric(clin.path.ord.val["patient.primary_pathology.tumor_sizes.tumor_size.pathologic_tumor_length",])) < median(as.vector(as.numeric(clin.path.ord.val["patient.primary_pathology.tumor_sizes.tumor_size.pathologic_tumor_length",])), na.rm = TRUE), "lower than median length", "higher than median length"),
  tumor.depth = ifelse(as.vector(as.numeric(clin.path.ord.val["patient.primary_pathology.tumor_sizes.tumor_size.pathologic_tumor_depth",])) < median(as.vector(as.numeric(clin.path.ord.val["patient.primary_pathology.tumor_sizes.tumor_size.pathologic_tumor_depth",])), na.rm = TRUE), "lower than median depth", "higher than median depth")
)

# Create survival objects
surv.obj.val <- Surv(time = surv.inf.val$pfi.times, event = surv.inf.val$pfi.cen)

# Fit and plot single and paired PFI survival curves
surv.fit.val <- survfit(surv.obj.val ~ surv.inf.val$exp.goi, data = surv.inf.val)
surv.fit.val.sgl <- survfit(surv.obj.val ~ 1, data = surv.inf.val)

ggsurv.val.sgl <- ggsurvplot(fit = surv.fit.val.sgl,
                         data = surv.inf.val,
                         conf.int = TRUE,
                         risk.table = FALSE,
                         pval = TRUE,
                         title = "PFI validation set event curve"
)

print(ggsurv.val.sgl)
ggsave(file = paste0(runpath,"/pfi_val_survival_sgl.png"),
       plot = print(ggsurv.val.sgl),
       width = 12, height = 8, units = "in",device = "png")

ggsurv.val <- ggsurvplot(fit = surv.fit.val,
           data = surv.inf.val,
           conf.int = TRUE,
           risk.table = FALSE,
           pval = TRUE
           )
 
print(ggsurv.val)

ggsave(file = paste0(runpath,"/pfi_val_survival_split.png"),
                             plot = print(ggsurv.val),
                             width = 12, height = 8, units = "in",device = "png")
