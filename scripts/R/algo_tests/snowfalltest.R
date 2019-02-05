library(parallel)
library(snowfall)

sfInit(parallel = TRUE, cpus = detectCores()-1)
