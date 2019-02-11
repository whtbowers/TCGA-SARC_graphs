setwd("C:/Users/wbowers/Documents/tcga_replication/data")


exp.data <- read.csv("TCGA_SARC_mrna_data_lnorm_medc_nosdfilt.csv", row.names = 1, stringsAsFactors = FALSE)
exp.info <- read.csv("TCGA_SARC_mrna_info_2.csv", stringsAsFactors = FALSE)
clin.data <- read.csv("TCGA_SARC_clinical.csv", stringsAsFactors = FALSE)
# manifest <- read.csv("gdc_manifest.2019-02-08.txt", sep = "\t", stringsAsFactors = FALSE)
# 
# intersect(as.character(manifest$md5), as.character(clin.data$submitter_id))
# 
# exp.info[1,]
# clin.data$submitter_id
int <- intersect(as.character(exp.info[1,]), as.character(clin.data$submitter_id))

all.inf <- exp.info[,which(exp.info[1,] %in% int)]


int.clin <- c()

for (i in 1:ncol(all.inf)){
  int.clin <- c(int.clin, match(all.inf[1,i], clin.data$submitter_id))
}

clin.match <- clin.data[int.clin,]

if (sum(duplicated(as.character(exp.info[1,])), na.rm = TRUE)>0){
  print("AYE")
}

exp.int <- as.character(exp.info[1,]) %in% int
sum(exp.int, na.rm = TRUE)
as.character(exp.info[1,][which(duplicated(as.character(exp.info[1,])))])

match(as.character(exp.info[1,]), clin.data$submitter_id)

       