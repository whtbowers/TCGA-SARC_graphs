set.seed(120)
library(pamr)
data(khan)

gene.ids <- khan[,1]
gene.names <- khan[,2]
#khan <- khan[,-c(1,2)]
khan.data <- as.matrix(khan[2:nrow(khan), 3:ncol(khan)])

class(khan.data) <- "numeric"

khan.clabel <- as.character(data.frame(lapply(khan[1,], as.character), stringsAsFactors = FALSE))[-c(1,2)]# class label, ensure remove first 2 empty cols to match length

khan.object <- list(x = khan.data, y = khan.clabel, gene.names = gene.names, gene.ids = gene.ids)

#train model
pamtrain <- pamr.train(data = khan.object)

# Adapt thresholds
newthresh <- pamr.adaptthresh(pamtrain)

# retrain with new thresholds
pamtrain2 <- pamr.train(khan.object, threshold.scale=newthresh)

# Cross validate
pamres <- pamr.cv(pamtrain2, khan.object)

# select threshold with minimum error after CV

opt.thresh <- pamres$threshold[which(pamres$error == min(pamres$error))]

#Plot shrunken class centroids
png(filename = "C:/Users/wbowers/Documents/tcga_replication_2/figs/pamr/pam_centroids.png", res=150, width = 2000, height = 2000)
pamr.plotcen(pamtrain2, khan.object, threshold = opt.thresh)
dev.off()

# Plot CV probability
png(filename = "C:/Users/wbowers/Documents/tcga_replication_2/figs/pamr/pam_cvprobcurves.png", res=150, width = 2000, height = 2000)
pamr.plotcvprob(pamres, khan.object, threshold = opt.thresh)
dev.off()

# Plot cv curves
png(filename = "C:/Users/wbowers/Documents/tcga_replication_2/figs/pamr/pam_cvcurves.png", res=150, width = 2000, height = 2000)
pamr.plotcv(pamres)
dev.off()

# Find false discovery rate
pamfdr <- pamr.fdr(pamtrain2, khan.object)

# Plot FDR
png(filename = "C:/Users/wbowers/Documents/tcga_replication_2/figs/pamr/pam_fdr.png", res=150, width = 2000, height = 2000)
pamr.plotfdr(pamfdr)
dev.off()

# List genes selected as significant
pamr.listgenes(fit = pamtrain2, data = khan.object, threshold = opt.thresh, fitcv = pamres)

# Visualise confusion matrix
confusion <- pamr.confusion(pamres, threshold = opt.thresh)

# Show geneplot (fit has to be result of pamr.train)
png(filename = "C:/Users/wbowers/Documents/tcga_replication_2/figs/pamr/pam_geneplot.png", res=150, width = 10000, height = 10000)
pamr.geneplot(pamtrain2, khan.object, opt.thresh)
dev.off()
