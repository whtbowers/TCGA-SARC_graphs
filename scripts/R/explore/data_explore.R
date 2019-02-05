setwd("C:/Users/wbowers/Documents/tcga_replication/data")
df <- read.csv("TCGA.SARC.264.matrix.txt", header=FALSE, sep = "\t",stringsAsFactors = FALSE)

# Separate labels from numerical data
hist.lab <- df[1:3,] # Histological labels and patient tags
hist.lab <- hist.lab[,-1] # drop col only containing 'NAME'

data <- df[-c(1:3),]
rownames(data) <- data[,1]
data <- data[-1]
 
#Set patient labels as colnames for easy alignment
colnames(data) <- hist.lab[1,]

# Save info and data as csv files
write.csv(data, "TCGA_SARC_data_raw.csv")
write.csv(hist.lab, "TCGA_SARC_info_raw.csv")


# Check same number of columns 
ncol(hist.lab)==ncol(data)

# Look at distribution of UBC
library(tidyverse)

ggplot()+
  geom_histogram(aes(y = data[1,]))
