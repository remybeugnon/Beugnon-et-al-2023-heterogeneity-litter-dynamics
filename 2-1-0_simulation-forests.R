rm(list = ls())
set.seed(1603)
# Functions
library(tidyverse)
library(purrr)

ncol = 18
nrow = 18

grid = rep(LETTERS[1:9], ncol*nrow/9)

nb.same.type = function(grid){
  # return the number of neighboors that have the same type
  sum.neighboors = function(c){
    i = c[2] # row
    j = c[1] # col
    x = c(i-1, i, i, i+1) 
    y = c(j, j+1, j-1, j)
    
    to.keep = x > 0 & x <= nrow & y > 0 & y <= ncol
    x = x[to.keep]
    y = y[to.keep]
    
    # move back to flat coordinates
    coords.flat = (x-1)*ncol + y
    return(sum(grid[(i-1)*ncol + j] == grid[coords.flat]))
  }
  mat.coords = expand.grid(1:ncol, 1:nrow)
  apply(mat.coords, MARGIN = 1, FUN = sum.neighboors)
}

expected.hypergeom = function(nb.n, nb.per.species){
  return(((nb.per.species-1)*nb.n)/(nrow*ncol-1))
}

to.optimise = function(grid){
  nb.neighbours = c(c(2, rep(3, ncol-2), 2), rep(c(3,rep(4, ncol-2), 3), nrow-2), c(2, rep(3, ncol-2), 2)) 
  expectations.hypergeom = vapply(nb.neighbours, expected.hypergeom, FUN.VALUE = numeric(1), nb.per.species = nrow*ncol/length(unique(grid)))
  observed = nb.same.type(grid) 
  return(sum(observed-expectations.hypergeom))
}

permute = function(vect_forest, nb){
  for(i in 1:nb){
    permutation = sample(1:length(vect_forest), 2)
    a = vect_forest[permutation[1]]
    b = vect_forest[permutation[2]]
    vect_forest = 
      vect_forest %>% 
      replace(., permutation[1], b) %>%
      replace(., permutation[2], a)
  }
  return(vect_forest)
}


#### Design forests ####
designs.9 = {
  list(
    design.bloc = c(
      c(rep('A',6),rep('B',6),rep('C',6)),
      c(rep('A',6),rep('B',6),rep('C',6)),
      c(rep('A',6),rep('B',6),rep('C',6)),
      c(rep('A',6),rep('B',6),rep('C',6)),
      c(rep('A',6),rep('B',6),rep('C',6)),
      c(rep('A',6),rep('B',6),rep('C',6)),
      c(rep('D',6),rep('E',6),rep('F',6)),
      c(rep('D',6),rep('E',6),rep('F',6)),
      c(rep('D',6),rep('E',6),rep('F',6)),
      c(rep('D',6),rep('E',6),rep('F',6)),
      c(rep('D',6),rep('E',6),rep('F',6)),
      c(rep('D',6),rep('E',6),rep('F',6)),
      c(rep('G',6),rep('H',6),rep('I',6)),
      c(rep('G',6),rep('H',6),rep('I',6)),
      c(rep('G',6),rep('H',6),rep('I',6)),
      c(rep('G',6),rep('H',6),rep('I',6)),
      c(rep('G',6),rep('H',6),rep('I',6)),
      c(rep('G',6),rep('H',6),rep('I',6))),
    
    design.block.2 =
      c(
        c(rep('A',3),rep('B',3),rep('C',3)),c(rep('D',3),rep('E',3),rep('F',3)),
        c(rep('A',3),rep('B',3),rep('C',3)),c(rep('D',3),rep('E',3),rep('F',3)),
        c(rep('A',3),rep('B',3),rep('C',3)),c(rep('D',3),rep('E',3),rep('F',3)),
        c(rep('G',3),rep('H',3),rep('I',3)),c(rep('A',3),rep('B',3),rep('C',3)),
        c(rep('G',3),rep('H',3),rep('I',3)),c(rep('A',3),rep('B',3),rep('C',3)),
        c(rep('G',3),rep('H',3),rep('I',3)),c(rep('A',3),rep('B',3),rep('C',3)),
        c(rep('D',3),rep('E',3),rep('F',3)),c(rep('G',3),rep('H',3),rep('I',3)),
        c(rep('D',3),rep('E',3),rep('F',3)),c(rep('G',3),rep('H',3),rep('I',3)),
        c(rep('D',3),rep('E',3),rep('F',3)),c(rep('G',3),rep('H',3),rep('I',3)),
        c(rep('A',3),rep('B',3),rep('C',3)),c(rep('D',3),rep('E',3),rep('F',3)),
        c(rep('A',3),rep('B',3),rep('C',3)),c(rep('D',3),rep('E',3),rep('F',3)),
        c(rep('A',3),rep('B',3),rep('C',3)),c(rep('D',3),rep('E',3),rep('F',3)),
        c(rep('G',3),rep('H',3),rep('I',3)),c(rep('A',3),rep('B',3),rep('C',3)),
        c(rep('G',3),rep('H',3),rep('I',3)),c(rep('A',3),rep('B',3),rep('C',3)),
        c(rep('G',3),rep('H',3),rep('I',3)),c(rep('A',3),rep('B',3),rep('C',3)),
        c(rep('D',3),rep('E',3),rep('F',3)),c(rep('G',3),rep('H',3),rep('I',3)),
        c(rep('D',3),rep('E',3),rep('F',3)),c(rep('G',3),rep('H',3),rep('I',3)),
        c(rep('D',3),rep('E',3),rep('F',3)),c(rep('G',3),rep('H',3),rep('I',3))),
    
    design.2lines = c(
      rep("A", 36), rep("B", 36), rep("C", 36),
      rep("D", 36), rep("E", 36), rep("F", 36),
      rep("G", 36), rep("H", 36), rep("I", 36)),
    
    design.lines = c(
      rep("A", 18), rep("B", 18), rep("C", 18),
      rep("D", 18), rep("E", 18), rep("F", 18),
      rep("G", 18), rep("H", 18), rep("I", 18),
      rep("A", 18), rep("B", 18), rep("C", 18),
      rep("D", 18), rep("E", 18), rep("F", 18),
      rep("G", 18), rep("H", 18), rep("I", 18)),
    
    design.demilines = c(
      rep("A", 9), rep("B", 9), rep("C", 9),
      rep("D", 9), rep("E", 9), rep("F", 9),
      rep("G", 9), rep("H", 9), rep("I", 9),
      rep("A", 9), rep("B", 9), rep("C", 9),
      rep("D", 9), rep("E", 9), rep("F", 9),
      rep("G", 9), rep("H", 9), rep("I", 9),
      rep("A", 9), rep("B", 9), rep("C", 9),
      rep("D", 9), rep("E", 9), rep("F", 9),
      rep("G", 9), rep("H", 9), rep("I", 9),
      rep("A", 9), rep("B", 9), rep("C", 9),
      rep("D", 9), rep("E", 9), rep("F", 9),
      rep("G", 9), rep("H", 9), rep("I", 9))
    
  )
}
designs.4 = {
  list(
    design.bloc = c(
      rep(c(rep('A',9),rep('B',9)),9),
      rep(c(rep('C',9),rep('D',9)),9))
    )
}
designs.2 = {
  list(
    design.bloc = c(rep('A',18*9),rep('B',18*9)))
}

designs.9 = c(
  designs.9, 
  map(.x = 1:400, .f = function(x){
    forest = permute(designs.9[['design.bloc']], x)
    forest})
  )
  
designs.4 = c(
  designs.4, 
  map(.x = 1:500, .f = function(x){
    forest = permute(designs.4[['design.bloc']], x)
    forest})
  )

designs.2 = c(
  designs.2, 
  map(.x = 1:500, .f = function(x){
    forest = permute(designs.2[['design.bloc']], x)
    forest})
  )

#### Forests selections ####
forest.designs.9 = designs.9[c(1:5,15,25,44,56,75,105,155,305)]
forest.designs.4 = designs.4[c(1,12,30,47,75,100,165,500)]
forest.designs.2 = designs.2[c(1,12,26,47,75,114,165,400)]

# Heterogeneity calculation
heterogeneity.forest.9 =
  data.frame(
    forests.9 = 1:13,
    description = c(
      'Blocks', 'Small blocks', 
      'Two rows', 'One row', 'Half rows',
      paste(c(15,25,44,56,75,105,155,305), 'permutations')
    ) %>% factor(., levels = c(
      'Blocks', 'Small blocks', 
      'Two rows', 'One row', 'Half rows',
      paste(c(15,25,44,56,75,105,155,305), 'permutations')
    )),
    homogeneity = map_dbl(.x = 1:13, 
        .f = function(x){
          to.optimise(forest.designs.9[[x]])})
)

heterogeneity.forest.4 =
  data.frame(
    forests.4 = 1:8,
    description = c(
      'Blocks', 
      paste(c(12,30,47,75,100,165,500), 'permutations')
    ) %>% factor(., levels = c(
      'Blocks',
      paste(c(12,30,47,75,100,165,500), 'permutations')
    )),
    homogeneity = map_dbl(.x = 1:8, 
                          .f = function(x){
                            to.optimise(forest.designs.4[[x]])})
  )

heterogeneity.forest.2 =
  data.frame(
    forests.2 = 1:8,
    description = c(
      'Blocks', 
      paste(c(12,26,47,75,114,165,400), 'permutations')
    ) %>% factor(., levels = c(
      'Blocks',
      paste(c(12,26,47,75,114,165,400), 'permutations')
    )),
    homogeneity = map_dbl(.x = 1:8, 
                          .f = function(x){
                            to.optimise(forest.designs.2[[x]])})
  )

#### Saving ####
saveRDS(forest.designs.2, file = '2-1-1-2sp-forests.RData')
saveRDS(forest.designs.4, file = '2-1-2_4sp-forests.RData')
saveRDS(forest.designs.9, file = '2-1-3_9sp-forests.RData')
#### END ####