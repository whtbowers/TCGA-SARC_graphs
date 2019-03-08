setwd("C:/Users/wbowers/Documents/tcga_replication_2/data")
set.seed(123.456)

library(RSQLite)

# Import data
exp.data <- read.csv(
  "TCGA_SARC_mrna_data.csv",
  row.names = 1,
  stringsAsFactors = FALSE
  )

exp.info <- read.csv(
  "TCGA_SARC_mrna_info.csv",
  row.names = 1,
  stringsAsFactors = FALSE
  )

# # Updated clinical data
# clin.data <- read.csv(
#   "upd_clin_data.csv",
#   row.names = 1
# )
# 
# # Clinical endpoint data
# clin.endp <- read.csv(
#   "upd_clin_data.csv",
#   row.names = 1
# )

# Updated clinical data, sarcoma only
clin.data <- read.csv(
  "liu_clin_sarc.csv",
  row.names = 1
)

# Clinical endpoint data, sarcoma only
clin.endp <- read.csv(
  "liu_endpoints_sarc.csv",
  row.names = 1
)

#Extract just sarcoma data from clinical data

# clin.data.sarc <- clin.data[which(clin.data$type=="SARC"),]
# clin.endp.sarc <- clin.endp[(clin.endp$type == "SARC"),]
# 
# write.csv(clin.data.sarc, "liu_clin_sarc.csv")
# write.csv(clin.endp.sarc, "liu_endpoints_sarc.csv")

# survival.data

# Transpose data to comply with SQL format
exp.data.t <- as.data.frame(t(exp.data))
exp.info.t <- as.data.frame(t(exp.info))

# rename info columns
colnames(exp.info.t) <- c("barcode_id", "hist_verbose", "checksum", "hist_abbrv", "3k_ConClust")

# Rename columns with rownames - All gene names unique after filtering
gnames <- c()
for (i in 1:nrow(exp.data)){
  gnames <- c(gnames, strsplit(rownames(exp.data)[i], "\\|")[[1]][1])
}
colnames(exp.data.t) <- gnames


# exp.data.t <- exp.data.t[-1,]# Remove first row
exp.info.t


# Add column of long-form IDs
longform.ids <- rownames(exp.data.t)
exp.data.t <- cbind(long.ids = longform.ids, exp.data.t)
exp.info.t <- cbind(longform.ids, exp.info.t)

# Replace rownames with numbers
rownames(exp.data.t) <- 1:nrow(exp.data.t)
rownames(exp.info.t) <- 1:nrow(exp.info.t)

#Initiate empty database
sts.db <- dbConnect(RSQLite::SQLite(), "TCGA-SARC.db")

dbWriteTable(sts.db, "expdata", exp.data.t)
dbWriteTable(sts.db, "expinfo", exp.info.t)
dbWriteTable(sts.db, "clinical", clin.data)
dbWriteTable(sts.db, "endpoints", clin.endp)

dbListTables(sts.db)
