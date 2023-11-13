simul = function(forest, restrict.dispersal = T) {
  
  # Define forest
  tree.grid = expand.grid(x = 1:18,
                          y = 1:18)
  
  # Define grid
  point.grid = {
    expand.grid(x = seq(1, 18, 0.1),
                y = seq(1, 18, 0.1)) |>
      # add one column per species for distance, biomass and litter
      mutate(bio.S1 = 0) |>
      mutate(bio.S2 = 0) |>
      mutate(bio.S3 = 0) |>
      mutate(bio.S4 = 0) |>
      mutate(bio.S5 = 0) |>
      mutate(bio.S6 = 0) |>
      mutate(bio.S7 = 0) |>
      mutate(bio.S8 = 0) |>
      mutate(bio.S9 = 0) |>
      mutate(bio.S10 = 0) |>
      mutate(bio.S11 = 0) |>
      mutate(bio.S12 = 0) |>
      mutate(dist.S1 = 0) |>
      mutate(dist.S2 = 0) |>
      mutate(dist.S3 = 0) |>
      mutate(dist.S4 = 0) |>
      mutate(dist.S5 = 0) |>
      mutate(dist.S6 = 0) |>
      mutate(dist.S7 = 0) |>
      mutate(dist.S8 = 0) |>
      mutate(dist.S9 = 0) |>
      mutate(dist.S10 = 0) |>
      mutate(dist.S11 = 0) |>
      mutate(dist.S12 = 0) |>
      mutate(biodist.S1 = 0) |>
      mutate(biodist.S2 = 0) |>
      mutate(biodist.S3 = 0) |>
      mutate(biodist.S4 = 0) |>
      mutate(biodist.S5 = 0) |>
      mutate(biodist.S6 = 0) |>
      mutate(biodist.S7 = 0) |>
      mutate(biodist.S8 = 0) |>
      mutate(biodist.S9 = 0) |>
      mutate(biodist.S10 = 0) |>
      mutate(biodist.S11 = 0) |>
      mutate(biodist.S12 = 0) |>
      mutate(litter.S1 = 0) |>
      mutate(litter.S2 = 0) |>
      mutate(litter.S3 = 0) |>
      mutate(litter.S4 = 0) |>
      mutate(litter.S5 = 0) |>
      mutate(litter.S6 = 0) |>
      mutate(litter.S7 = 0) |>
      mutate(litter.S8 = 0) |>
      mutate(litter.S9 = 0) |>
      mutate(litter.S10 = 0) |>
      mutate(litter.S11 = 0) |>
      mutate(litter.S12 = 0)
  }
  
  # Calculate distance from trees
  sp.biomass = data.frame(
    sp = paste0('S',1:12), 
    biomass = c(0.0107,0.0104,0.000589,0.00139,
                0.00291, 0.0268,0.00288,0.00611,
                0.00153,0.00525,0.00389,0.00250)
  )
  
  for (i in 1:nrow(point.grid)) {
    for (j in 1:nrow(tree.grid)) {
      dist = sqrt((point.grid$x[i] - tree.grid$x[j]) ^ 2 +
                    (point.grid$y[i] - tree.grid$y[j]) ^ 2)
      if (restrict.dispersal == T) {
        if (dist < sqrt(1.5 ^ 2 + 1) & dist > 0) {
          sp = forest[j]
          point.grid[i, paste0('bio.', sp)] = point.grid[i, paste0('bio.', sp)] + sp.biomass$biomass[sp.biomass$sp == sp] # here biomass = .1 for each tree
          point.grid[i, paste0('dist.', sp)] = point.grid[i, paste0('dist.', sp)] + 1 /
            dist
          point.grid[i, paste0('biodist.', sp)] = point.grid[i, paste0('biodist.', sp)] +  sp.biomass$biomass[sp.biomass$sp == sp] /
            dist
        }
      }else{
        if (dist > 0) {
          sp = forest[j]
          point.grid[i, paste0('bio.', sp)] = point.grid[i, paste0('bio.', sp)] + sp.biomass$biomass[sp.biomass$sp == sp] # here biomass = .1 for each tree
          point.grid[i, paste0('dist.', sp)] = point.grid[i, paste0('dist.', sp)] + 1 /
            dist
          point.grid[i, paste0('biodist.', sp)] = point.grid[i, paste0('biodist.', sp)] + sp.biomass$biomass[sp.biomass$sp == sp] /
            dist
        }
      }
    }
  }
  
  #### 1. litter distribution ####
  #### 1.1. Extraction model fit ####
  load(file = "fit-litter.RData")
  
  parameters.b =
    as.matrix(fit) |>
    as.data.frame() |>
    select(paste0('bio[', 1:12, ']')) |>
    colMeans()
  
  parameters.d =
    as.matrix(fit) |>
    as.data.frame() |>
    select(paste0('d[', 1:12, ']')) |>
    colMeans()
  
  parameters.db =
    as.matrix(fit) |>
    as.data.frame() |>
    select(paste0('db[', 1:12, ']')) |>
    colMeans()
  
  #### 1.2. Estimate species specific fall ####
  for (i in 1:nrow(point.grid)) {
    for (sp in 1:12) {
      est =
        parameters.b[sp] * point.grid[i, paste0('bio.S', sp)] +
        parameters.d[sp] * point.grid[i, paste0('dist.S', sp)] +
        parameters.db[sp] * point.grid[i, paste0('biodist.S', sp)]
      if (est > 0) {
        point.grid[i, paste0('litter.S', sp)] = est
      }
    }
  }
  
  #### 1.3. Aggregate indices ####
  point.grid$litter.total = rowSums(point.grid[, grep('litter.S', colnames(point.grid))])
  point.grid$litter.sp.rich = rowSums(point.grid[, grep('litter.S', colnames(point.grid))] !=
                                        0)
  
  #### 2. Predicting decomposition ####
  #### 2.1. Extraction model fit ####
  load(file = "fit-decomposition.RData")
  
  parameters.decom.C = fit_carbon |>
    as.matrix() |>
    colMeans()
  
  parameters.decom.N = fit_nitrogen |>
    as.matrix() |>
    colMeans()
  
  #### 4.2. Predict C and N loss ####
  point.grid =
    bind_cols(point.grid,
              point.grid[,
                         grep('litter.S',
                              colnames(point.grid))] /
                point.grid$litter.total)
  
  point.grid$decomp.C = NA
  point.grid$decomp.N = NA
  
  for (i in 1:nrow(point.grid)) {
    new.data = point.grid |>
      as.matrix()
    newdata = new.data[i, c(53:64, 51, 52)]
    par.C = fit_carbon |>
      as.matrix() |>
      colMeans()
    par.N = fit_nitrogen |>
      as.matrix() |>
      colMeans()
    pred.c =
      sum(par.C[1:12] * newdata[1:12]) + # identity effect
      sum(par.C[17:28] * newdata[1:12] * (1 - newdata[1:12])) + # interaction effect
      par.C[33] * newdata[13] + # biomass effect
      par.C[34] * newdata[14]
    
    if (pred.c > 100) {
      point.grid$decomp.C[i] =  100
    } else if (pred.c < 0) {
      point.grid$decomp.C[i] =  0
    } else{
      point.grid$decomp.C[i] = pred.c
    }
    
    point.grid$decomp.N[i] =
      sum(par.N[1:12] * newdata[1:12]) + # identity effect
      sum(par.N[17:28] * newdata[1:12] * (1 - newdata[1:12])) + # interaction effect
      par.N[33] * newdata[13] + # biomass effect
      par.N[34] * newdata[14]
  }
  #### 3. Export ####
  return(list(
    final_forest = point.grid,
    summary = summarise(
      point.grid,
      total.litter = sum(litter.total),
      mean.litter = mean(litter.total),
      sd.litter = sd(litter.total),
      mean.sp.rich.litter = mean(litter.sp.rich),
      sd.sp.rich.litter = sd(litter.sp.rich),
      mean.decomp.rate.C = mean(decomp.C),
      sd.decomp.rate.C = sd(decomp.C),
      mean.decomp.rate.N = mean(decomp.N),
      sd.decomp.rate.N = sd(decomp.N),
      total.decomp.C = sum(decomp.C * litter.total * 0.5  / 100),
      total.decomp.N = sum(decomp.N * litter.total * 0.04 / 100)
    )
  )
  )
}