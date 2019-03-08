setwd("C:/Users/wbowers/Documents/tcga_replication_2/data")
set.seed(123.456)

write.csv(data.frame(gnames), "gene_names.csv")

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

death.bin <- ifelse(clin.liu.ord$vital_status == "Dead", 1, 0)

# Merge times to death or last follow-up

os.times <- clin.liu.ord$death_days_to
os.times[which(clin.liu.ord$death_days_to == "#N/A")] <- clin.liu.ord$last_contact_days_to[which(clin.liu.ord$death_days_to == "#N/A")]

#Start w/ median
exp.split <- ifelse(exp.goi.ord >= median(exp.goi.ord),
                    paste("High", goi, "expression"), 
                    paste("Low", goi, "expression")
                    
)

# Create df for easier graphing
surv.inf<- data.frame(
  os.cen <- death.bin,
  os.times = as.numeric(os.times),
  exp.goi = exp.split,
  stage = clin.liu.ord$ajcc_pathologic_tumor_stage,
  sex = clin.liu.ord$gender#,
  #age = ages
)

surv.obj <- Surv(time = surv.inf$os.times, event = surv.inf$os.cen)
surv.fit <- survfit(surv.obj ~ exp.goi, data = surv.inf) 

ggsurvplot(surv.fit,
           data = surv.inf,
           pval = TRUE,
           #risk.table = TRUE,
           conf.int = TRUE
)
      