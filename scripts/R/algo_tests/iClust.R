library(iClusterPlus)
library(GenomicRanges)
library(gplots)
library(lattice)
# Data from level 3 normed/segmented data (gbm.seg), mutations (gbm.mut), and expression (gbm.exp)
data(gbm)
dim(gbm.mut) # Dimensions of mutation data (84*306)

# Collect mean over columns
mut.rate <- apply(gbm.mut, 2, mean) 

# Select columns where mean mutation rate > 0.02
gbm.mut2 <- gbm.mut[, which(mut.rate > 0.02)]
# gbm.mut2[1:10, 1:8]

# By filtering for this variance, get top variable genes to improve analysis

# Dimension of expression data (84, 1740)
dim(gbm.exp)
gbm.exp[1:3, 1:8] #First three samples, 8 genes

# Difficult to incorporate raw/normed copy number data for iCluster analysis as high dim and spatial correlation of data
# Probably better to use DNAcopy
# BiocManager::install("DNAcopy")

# Dimensions of raw/normalised data (16295, 6)
dim(gbm.seg)
gbm.seg[1:3,] # If use DNAcopy, will be gbm.cn

# Reduce gbm copy number to 5K by removing redundant regions using CNregions

data("variation.hg18.v10.nov.2010")
gbm.cn <- CNregions(seg = gbm.seg,
                    epsilon = 0, 
                    adaptive = FALSE, 
                    rmCNV = TRUE,
                    cnv = variation.hg18.v10.nov.2010[,3:5],
                    frac.overlap = 0.5,
                    rmSmallseg = TRUE,
                    nProbes = 5)

gbm.cn[1:3, 1:5]

# Ensure samples all in correct order:
gbm.cn <- gbm.cn[order(rownames(gbm.cn)),]

# Check samples all in same order for all 3 data sets
all(rownames(gbm.cn)==rownames(gbm.exp) || rownames(gbm.cn==rownames(gbm.mut)))

## Integrative cluster analysis ##
# iClusterPlus fits regularised latent model with integrated cluster assignment
# Run with k desired eigenfeatures, giving k+1 clusters
# Normally have to optimise k and lamda by model tuning

fit.single <- iClusterPlus(dt1 = gbm.mut2,
                           dt2 = gbm.cn,
                           dt3 = gbm.exp,
                           type = c("binomial", "gaussian", "gaussian"),
                           lambda = c(0.04, 0.61, 0.90),
                           K = 2,
                           maxiter = 10)

## Model tuning ##
# If cluster number unknown, test range of k from 1 to estimated reasonable k

#Time tuning process
t <- Sys.time()

set.seed(123)
date()
for(k in 1:5){
  cv.fit <- tune.iClusterPlus(cpus = 1,
    # cpus = detectCores()-1, # Run on all but one core
   dt1 = gbm.mut2,
   dt2 = gbm.cn,
   dt3 = gbm.cn,
   type = c("binomial", "gaussian", "gaussian"),
   K = k,
   n.lambda = 185,
   scale.lambda = c(1, 1, 1),
   maxiter = 20)
  
  save(cv.fit, file = paste("cv.fit.k", k, ".Rdata", sep=""))
}
Sys.time()-t
#####################################################################
# As windows (at least on assigned laptop) can't run this           #
# on more than one core, I'll end it  here before I die of old age. #
#####################################################################