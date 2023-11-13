rm(list = ls())
set.seed(1603)

library(tidyverse)
library(purrr)
# Function
{
  ncol = 18
  nrow = 18
  nb_species = 9
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
  to.optimise = function(grid, expected){
    nb.neighbours = c(c(2, rep(3, ncol-2), 2), rep(c(3,rep(4, ncol-2), 3), nrow-2), c(2, rep(3, ncol-2), 2)) 
    expectations.hypergeom = vapply(nb.neighbours, expected.hypergeom, FUN.VALUE = numeric(1), nb.per.species = nrow*ncol/9) # length(unique(grid)))
    observed = nb.same.type(grid) 
    return(sum(observed-expectations.hypergeom))
  }
}


#### > 1. 9 sp mix ####
forests = readRDS(file = '2-2-3_9sp-predictions.RData')

list.records = 
  list.files('simul/')
list.records = list.records[grep('9sp', list.records)]
x = list.records[1]
# com.dist = read.csv('simul/simul-summary/community-distance.csv')

df = map_df(.x = list.records[-1], 
            function(x){
              data.frame(
                rep = x, 
                read.csv(paste0('simul/',x)))
            }) 
name.short = colnames(df)
colnames(df) = c(
  'Replicate','Design' ,
  'Total litterfall (g)', 'Mean litterfall (g/dm2)', 'Litterfall SD (g/dm2)',
  'Mean litter species richness (#/dm2)','Litter species richness SD (#/dm2)',
  'Mean carbon decomposition rate (%)', 'Carbon decomposition rate SD (%)',
  'Mean nitrogen decomposition rate (%)', 'Nitrogen decomposition rate SD (%)',
  'Total C loss (g)','Total N loss (g)'
)

df$Replicate = 
  df$Replicate |> 
  str_remove_all('output-2sp_rep-') |>
  str_remove_all('output-4sp_rep-') |>
  str_remove_all('output-9sp_rep-') |> 
  str_remove_all('.csv') |> 
  as.numeric()

df.1 = df |>
  pivot_longer(cols = 3:13)

df.1$name = df.1$name |>
  factor(levels = c(
    'Total litterfall (g)', 'Mean litterfall (g/dm2)', 'Litterfall SD (g/dm2)',
    'Mean litter species richness (#/dm2)','Litter species richness SD (#/dm2)',
    'Mean carbon decomposition rate (%)', 'Carbon decomposition rate SD (%)',
    'Mean nitrogen decomposition rate (%)', 'Nitrogen decomposition rate SD (%)',
    'Total C loss (g)','Total N loss (g)'
  ))

df.1$heterogeneity = NA
for(i in 1:length(forests$design)){
  df.1$heterogeneity[df.1$Design == i] = -(to.optimise(forests$design[[i]])+1)
}

df.1 = 
  df.1 |> 
  filter(Design != 5)
df.1$value[df.1$name == 'Total C loss (g)'] = 
  df.1$value[df.1$name == 'Total C loss (g)']/(17*17*100)

df.2 = df.1 |> filter(name %in% c('Mean litterfall (g/dm2)', 
                                  'Litterfall SD (g/dm2)',
                                  'Mean litter species richness (#/dm2)',
                                  'Mean carbon decomposition rate (%)', 
                                  'Carbon decomposition rate SD (%)',
                                  'Total C loss (g)'))
df.2$name = as.character(df.2$name)
df.2$name[df.2$name == 'Total C loss (g)'] = 'Mean C loss (g/dm2)'
df.2$name = factor(df.2$name, levels = c('Mean litterfall (g/dm2)', 
                                         'Litterfall SD (g/dm2)',
                                         'Mean litter species richness (#/dm2)',
                                         'Mean carbon decomposition rate (%)', 
                                         'Carbon decomposition rate SD (%)',
                                         'Mean C loss (g/dm2)'))

df.3 = df.2 |> 
  filter(name %in% c('Mean litterfall (g/dm2)', 
                     'Litterfall SD (g/dm2)',
                     'Mean litter species richness (#/dm2)'))

p.litter = 
  ggplot(data = df.3, 
         aes(x = heterogeneity, 
             y = value, 
             color = factor(Design))) + 
  geom_jitter(data = df.3, 
                     aes(x = heterogeneity, 
                         y = value, 
                         color = factor(Design)),
                     size = .5, 
                     alpha = .5, 
              width = 10) +
  geom_smooth(data = df.3|> 
                filter(!(Design %in% 2:5)), 
              aes(x = heterogeneity, 
                  y = value, 
                  color = factor(Design)),
              color = 'black', 
              method="loess", 
              level=0.90,
              size = .5) + 
  facet_wrap(vars(name), scales = 'free',ncol =  3) + 
  scale_color_manual(
    values = c(
      '1' = 'darkred',
      '2' = 'red', 
      '3' = 'darkblue',
      '4' = 'blue', 
      '5' = 'lightgreen',
      '6' = 'gray','7' = 'gray','8' = 'gray',
      '9' = 'gray','10' = 'gray','11' = 'gray',
      '12' = 'gray','13' = 'gray','14' = 'black')) + 
  scale_y_continuous(labels = scales::label_number_si()) +
  scale_x_continuous(breaks = c(-950, -730, -805, -480, 0), 
                     labels = c("Blocks", 'Mini blocks',
                                "2 lines", 'Lines', 'Random')) + 
  labs(x = "", 
       y = "Value", 
       title = "A. Litterfall") +
  theme_bw() + 
  theme(legend.position = 'none', 
        axis.text.x = element_text(angle = 45, hjust = 1))
p.litter 

df.4 = df.2 |> filter(name %in% c('Mean carbon decomposition rate (%)', 
                                  'Carbon decomposition rate SD (%)',
                                  'Mean C loss (g/dm2)'))

p.decomp = 
  ggplot(data = df.4, 
         aes(x = heterogeneity, 
             y = value, 
             color = factor(Design))) + 
  geom_jitter(data = df.4, 
              aes(x = heterogeneity, 
                  y = value, 
                  color = factor(Design)),
              size = .5, 
              alpha = .5,
              width = 10) +
  geom_smooth(data = df.4 |> 
                filter(!(Design %in% 2:5)), 
              aes(x = heterogeneity, 
                  y = value, 
                  color = factor(Design)),
              color = 'black', 
              method="loess", 
              level=0.90,
              size = .5) + 
  facet_wrap(vars(name), scales = 'free',ncol =  3) + 
  scale_color_manual(
    values = c(
      '1' = 'darkred',
      '2' = 'red', 
      '3' = 'darkblue',
      '4' = 'blue', 
      '5' = 'lightgreen',
      '6' = 'gray','7' = 'gray','8' = 'gray',
      '9' = 'gray','10' = 'gray','11' = 'gray',
      '12' = 'gray','13' = 'gray','14' = 'black')) + 
  scale_y_continuous(labels = scales::label_number_si()) +
  scale_x_continuous(breaks = c(-950, -730, -805, -480, 0), 
                     labels = c("Blocks", 'Mini blocks',
                                "2 lines", 'Lines', 'Random')) + 
  labs(x = "Tree species spatial heterogeneity", 
       y = "Value", 
       title = "B. Decomposition") +
  theme_bw() + 
  theme(legend.position = 'none', 
        axis.text.x = element_text(angle = 45, hjust = 1))
p.decomp  


Fig.2 = 
  ggpubr::ggarrange(p.litter,
            p.decomp, 
            nrow = 2,
            align = 'hv')
Fig.2

ggsave(plot = Fig.2, 
       filename = '3-1_Fig2.png',
       height = 15, width = 20, units = 'cm')

ggsave(plot = Fig.2, 
       filename = '3-1_Fig2.svg',
       height = 18, width = 20, units = 'cm')