rm(list = ls())
set.seed(1603)
library(tidyverse)
# Functions
source('1-3-prediction-function.R')
forests = readRDS(file = '2-1-2-4sp-forests.RData')

start = 1
stop  = 2

for(i in start:stop){
    result.sim = data.frame(matrix(data = NA,
                                   nrow = length(forests$design),
                                   ncol = 12))
    colnames(result.sim) = c("design", "total.litter", "mean.litter",
                             "sd.litter", "mean.sp.rich.litter",
                             "sd.sp.rich.litter",   "mean.decomp.rate.C",
                             "sd.decomp.rate.C",    "mean.decomp.rate.N",
                             "sd.decomp.rate.N",    "total.decomp.C",
                             "total.decomp.N")
    
    for(j in 1:length(forests[['design']])){
      s.p = simul(forest = forests$species[[j]][[i]])
      write.csv(s.p[['final_forest']],
                paste0('simul/4sp/simul-4sp_design-', j, '_rep-',i,'.csv'),
                row.names = F)
      result.sim[j,] = c(j, s.p[['summary']])
    }
    write.csv(result.sim,
              paste0('simul/4sp/output-4sp_rep-',i,'.csv'),
              row.names = F)
}