setwd("C:/Users/wbowers/Documents/tcga_replication_2/data")

library(samr)

#Two-class example
twoclass <- read.csv(paste(path.package("samr"), "/excel/twoclass.csv", sep = ""), header = FALSE, stringsAsFactors = FALSE)

# Separate data and labels
twoclass.geneid <- as.character(twoclass[-1, 1])
twoclass.genename <- as.character(twoclass[-1, 2])
twoclass.response <- as.vector(as.character(twoclass[1,-c(1, 2)]))
twoclass.data <- as.matrix(twoclass[-1, -c(1, 2)])
class(twoclass.data) <- "numeric"


## Block permutation
# 1block2 = treatment 1, batch 2
# If 2 batches and 2 treatments, 4! = 24 permutations

# SAM automatically imputes via KNN, default k = 10

twoclass.object <- list(x = twoclass.data,
                        y = twoclass.response,
                        geneid = twoclass.geneid,
                        genenames = twoclass.genename,
                        logged2 = TRUE)

twoclass.sam <- samr(data = twoclass.object,
                     resp.type = "Two class unpaired", 
                     nperms = 100, 
                     nresamp = 20
                     )

#Estimate FDR, FNR, Type I error, max/min significant genes
twoclass.samsize <- samr.assess.samplesize(samr.obj = twoclass.sam, data = twoclass.object, dif = 0)

# Plot how FNR and Type I error change with number of genes selected
png(filename = "../figs/samr/sam_twoclass_estsamplesize.png", res=150, width = 2000, height = 2000)
samr.assess.samplesize.plot(twoclass.samsize, logx = TRUE)
dev.off()

# Default min.foldchange = 0, so no rule applied
twoclass.delta <- as.data.frame(samr.compute.delta.table(samr.obj = twoclass.sam))

# Select delta which provides maximum number of genes for minimal fdr
# If multiple deltas give same number of genes, choose moost stringent (minimal) delta
opt.delta <- min(twoclass.delta$delta[which(twoclass.delta$`# called` == max(twoclass.delta$`# called`[which(twoclass.delta$`median FDR` == min(twoclass.delta$`median FDR`, na.rm = TRUE))]))])

twoclass.siggenes <- samr.compute.siggenes.table(samr.obj = twoclass.sam, del = opt.delta, data = twoclass.object, delta.table = twoclass.delta)

# SAM Q-Q plot
png(filename = "../figs/samr/sam_twoclass_qq.png", res=150, width = 2000, height = 2000)
samr.plot(samr.obj = twoclass.sam)
dev.off()
