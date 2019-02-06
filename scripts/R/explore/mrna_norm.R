setwd("C:/Users/wbowers/Documents/tcga_replication/data")

exp.data <- read.csv("TCGA_SARC_mrna_data_Amatch.csv", row.names = 1, stringsAsFactors = FALSE)
# exp.data <- exp.data[1:1000,]

# Filter out all genes with < 90% nonzero expression
ind.90filt <- c()
for (i in 1:nrow(exp.data)){
 if ((sum(exp.data[i,] == 0, na.rm = TRUE)/ncol(exp.data))>=0.1){
   ind.90filt <- c(ind.90filt, i)
 } 
}

exp.data.90filt <- exp.data[-ind.90filt,]

# Check distribution of first gene

library(ggplot2)

ggplot() +
  geom_histogram(aes(x=as.numeric(exp.data.90filt[1,])))

#log2 transform add 0.05 to prevent -inf
exp.data.log2 <- log2(exp.data.90filt+0.05)

# median centre for each gene across all tumours
exp.data.c1 <- apply(exp.data.log2,2,function(x){
  x-median(x)
})

# median centre for each tumour across all genes
exp.data.c2 <- as.data.frame(t(apply(exp.data.c1,1,function(x){
  x-median(x)
})))

write.csv(exp.data.c2, "TCGA_SARC_mrna_data_lnorm_medc_nosdfilt_Amatch.csv")

# Check new distribution
ggplot() +
  geom_histogram(aes(x=as.numeric(exp.data.c2[1,])))

# Remove genes with std < 2
ind.std2filt <- c()
for (i in 1:nrow(exp.data.c2)){
  if(sd(exp.data.c1[i,]) < 2){
    ind.std2filt <- c(ind.std2filt, i)
  }
}

exp.data.sdfilt <- exp.data.c2[-ind.std2filt,]

write.csv(exp.data.sdfilt, "TCGA_SARC_mrna_data_lnorm_medc_Amatch.csv")
