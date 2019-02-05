library(ALL)

data(ALL)
d=exprs(ALL) # Retrieve expression datafrom expression set (eSet) object

# d[1:5, 1:5] # Simply to show data layout

# Find 5000 most informative genes for class detection; use median absolute deviation (MAD)
mads = apply(d, 1, mad)
d=d[rev(order(mads))[1:5000],]

# If want to transform or normalise data, possible using other bioconductor methods here

# Using default settings of hierarchical clustering algo with Pearson correlation distance

d = sweep(d, 1, apply(d, 1, median, na.rm=T))

# Data now ready for analysis

# Selected 80% item resampling (pItem), 80% gene resampling (pFeature), and max k of 6, 50 resamplings, agglomerative hclustering of 1-Pearson correlation distance

# In practice, higher reps (eg 1000) and higher k (eg 20) recommended

library(ConsensusClusterPlus)
title = tempdir()
results = ConsensusClusterPlus(d,
                               maxK= 6,
                               reps= 50, 
                               pItem= 0.8, 
                               pFeature = 1,
                               title = title,
                               clusterAlg = "hc",
                               distance = "pearson",
                               seed = 1262118388.71279,
                               plot = "png")
# Display consensus matrix

results[[2]][["consensusMatrix"]][1:5, 1:5]
heatmap(results[[2]][["consensusMatrix"]])

# Display consensus tree - hclust tree object

results[[2]][["consensusTree"]]

# Display sample classifications

results[[2]][["consensusClass"]][1:5]

# (Optional) calculate cluster-consensus and item consensus

icl = calcICL(results, title=title, plot="png")

# Display ICL

icl[["clusterConsensus"]][1:5,] # Return clusters with top 5 hughest consensuses

# Alternatively can provide custom distance matrix as input
dt <- as.dist(1-cor(d, method="pearson"))

ConsensusClusterPlus(dt,
                     maxK= 4,
                     reps= 10,
                     pItem= 0.8,
                     pFeature = 1,
                     title = title,
                     clusterAlg = "hc",
                     distance = "pearson",
                     seed = 1262118388.71279)



 