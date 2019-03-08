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

## Select subtypes and draw DSS curve with median split
# To do: implement X-Tile to find optimal split for each subtype

subtypes <- unique(as.vector(as.character(exp.info.ord[4,])))


for (subtype in subtypes) {
  print(subtype)
  
  exp.info.sub <- exp.info.ord[, which(exp.info.ord[4,] == subtype)]
  exp.goi.sub <- exp.goi.ord[, which(exp.info.ord[4,] == subtype)]
  clin.data.sub <- clin.data.ord[which(exp.info.ord[4,] == subtype), ]
  
  # Binarise vital status for each clinical subset
  
  death.bin <- ifelse(clin.data.sub$days_to_death == '--', 0, 1)
  
  # Split by median
  
  valsplit <- ifelse(exp.goi.sub >= median(as.vector(as.numeric(exp.goi.sub))), 1, 2)
  
  # Get available ages at death (NB many right-censored)
  # TODO: For each group, find number of censored patients
  
  ages <- as.vector(as.numeric(clin.data.sub$year_of_death)) - as.vector(as.numeric(clin.data.sub$year_of_birth))
  age.split <- median((ages), na.rm = TRUE)
  age.bin <- ifelse(ages < age.split, 1, 2)
  
  
  surv.inf <- data.frame(
        os.times = as.numeric(clin.data.sub$days_to_death),
        pfs.times = as.vector((as.numeric(clin.data.sub$days_to_death)) - as.vector(as.numeric(clin.data.sub$days_to_last_follow_up))),
        hist.cluster = as.vector(as.character(exp.info.sub[4,])),
        sex = ifelse(as.vector(as.character(clin.data.sub$gender)) == "female", 1, 2),
        age = age.bin,
        patient.barcode = clin.data.sub$submitter_id,
        patient.vital_status = death.bin,
        CCT2 = as.vector(as.numeric(exp.goi.sub)),
        cutoff = as.vector(as.numeric(valsplit))
      )
  
  # Create survival curves based

  surv_object <- Surv(time = surv.inf$os.times, event = surv.inf$patient.vital_status)
  
  surv_fit <- survfit(surv_object ~ cutoff, data = surv.inf)
  
  ggsurv <- ggsurvplot(surv_fit,
                       data = surv.inf,
                       pval = TRUE,
                       risk.table = TRUE,
                       conf.int = TRUE,
                       #log.rank.weights = "survdiff",
                       title = paste("Disease specific survival for ", subtype, sep = "")
  )
  ggsurv$plot <- ggsurv$plot +
    ggplot2::annotate("text",
                      x = 1200, y = 0.9, # x/y coords of text
                      label = paste("Upper: ", sum(valsplit == 1, na.rm = TRUE),
                                    "\nLower: ", sum(valsplit == 2, na.rm = TRUE), sep = "")

    )
  ggsave(filename = paste(runpath, "/", subtype, "_dss_kmplot.png", sep=""), plot = print(ggsurv), width = 6, height = 4, units = "in",device = "png")
  
  fit.coxph <- coxph(surv_object ~ age + sex,
                     data = surv.inf)
  
  # Cox PH forest plot
  
  ggfor <- ggforest(fit.coxph, data = surv.inf, main = paste("Hazard ratios for DSS of", subtype))
  
  ggsave(paste(runpath, "/", subtype, "_coxforplot.png", sep=""), plot = print(ggfor), width = 6, height = 4, units = "in", device = "png")
  # dev.off()
  
  # Display mv Cox PH table
  
  print(paste("Cox PH for disease-specific survival for", subtype))
  print(summary(fit.coxph))
  
}