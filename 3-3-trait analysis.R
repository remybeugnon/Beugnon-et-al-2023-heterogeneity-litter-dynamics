rm(list = ls())
set.seed(1603)

library(tidyverse)
library(purrr)
library(randomForest)
library(caret)
library(datasets)

colors.reponses = 
  c(`Total litterfall (g)` = '#106022',                
    `Mean litterfall (g/dm2)` = '#50A162',             
    `Litterfall SD (g/dm2)` = '#80C18F',               
    `Mean litter species richness (#/dm2)` = '#737D15',
    `Litter species richness SD (#/dm2)` = '#9CA738',
    `Litter composition mean distance` = '#C6C068',
    `Mean carbon decomposition rate (%)` = '#D4B86A',  
    `Carbon decomposition rate SD (%)` = '#806315',    
    `Mean nitrogen decomposition rate (%)` = '#83877B',
    `Nitrogen decomposition rate SD (%)` = '#6CA299',  
    `Total C loss (g)` = '#553E00',                    
    `Total N loss (g)` = '#0D5146'
  )

colors.labels = 
  c(`Total litterfall (g)` = 'Total litterfall (g)',                
    `Mean litterfall (g/dm2)` = 'Mean litterfall (g/dm2)',             
    `Litterfall SD (g/dm2)` = 'Litterfall SD (g/dm2)',               
    `Mean litter species richness (#/dm2)` = 'Mean litter species richness (#/dm2)',
    `Litter species richness SD (#/dm2)` = 'Litter species richness SD (#/dm2)',
    `Litter composition mean distance` = 'Litter composition mean distance',
    `Mean carbon decomposition rate (%)` = 'Mean carbon decomposition rate (%)',  
    `Carbon decomposition rate SD (%)` = 'Carbon decomposition rate SD (%)',    
    `Mean nitrogen decomposition rate (%)` = 'Mean nitrogen decomposition rate (%)',
    `Nitrogen decomposition rate SD (%)` = 'Nitrogen decomposition rate SD (%)',  
    `Total C loss (g)` = 'Mean C loss (g/dm2)',                    
    `Total N loss (g)` = 'Mean N loss (g/dm2)'
  )

datasets = readRDS('2-2-3_9sp-predictions.RData')
df.1 = datasets[[2]]
df.1$value[df.1$name == 'Litter composition mean distance'] = 
  df.1$value[df.1$name == 'Litter composition mean distance'] |> 
  log()
forests = readRDS(file = '2-1-3-9sp-forests.RData')
species = c('S01', paste0('S',2:11))
leaf_traits = 
  read.csv(file = 'leaves-trait.csv') %>% 
  group_by(Species, Diversity) %>%
  summarise(SLA = mean(SLA..mm2.mg., na.rm = T),
            LDMC = mean(LDMC..mg.g., na.rm = T),
            CN = mean(C.N, na.rm = T), 
            C = mean(C...., na.rm = T), 
            N = mean(N...., na.rm = T), 
            Mg = mean(Mg..mg.g., na.rm = T),
            Ca = mean(Ca..mg.g., na.rm = T),
            K = mean(K..mg.g., na.rm = T),
            P = mean(P..mg.g., na.rm = T))

leaf_traits = leaf_traits %>% 
  filter(Diversity == 8)

leaf_traits$Species = 
  leaf_traits$Species %>%
  str_replace_all( . , pattern = "Castanea henryi", 'S01') %>%
  str_replace_all( . , pattern = "Nyssa sinensis", 'S2') %>%
  str_replace_all( . , pattern = "Quercus fabri", 'S3') %>%
  str_replace_all( . , pattern = "Quercus serrata", 'S4') %>%
  str_replace_all( . , pattern = "Castanopsis sclerophylla", 'S5') %>%
  str_replace_all( . , pattern = "Choerospondias axillaris", 'S6') %>%
  str_replace_all( . , pattern = "Sapium sebiferum", 'S7') %>%
  str_replace_all( . , pattern = "Lithocarpus glaber", 'S8') %>%
  str_replace_all( . , pattern = "Koelreuteria bipinnata", 'S9') %>%
  str_replace_all( . , pattern = "Liquidambar formosana", 'S10') %>%
  str_replace_all( . , pattern = "Sapindus mukorossi", 'S11') %>%
  str_replace_all( . , pattern = "Cyclobalanopsis glauca", 'S12')

leaf_traits = 
  leaf_traits %>%
  filter(Species %in% species) %>%
  arrange(Species)

library(FactoMineR)
leaf_traits$chemistry = PCA(leaf_traits[, 6:11])$ind$coord[,1]
leaf_traits$thoughness = PCA(leaf_traits[, 4:5])$ind$coord[,1]
leaf_traits$trait.1 = PCA(leaf_traits[, 3:11])$ind$coord[,1]
leaf_traits$trait.2 = PCA(leaf_traits[, 3:11])$ind$coord[,2]

lt.DF = leaf_traits %>% 
  select(Species, Diversity, trait.1, trait.2)  %>% 
  data.frame()
rownames(lt.DF) = lt.DF$Species
lt.DF = lt.DF%>% select(-Diversity, - Species)

outputs = unique(df.1$name)[c(2,3,4,6,7,10)]
var.tested = 'Total C loss (g)'
forest.analyses = function(var.tested, dataset, forests, species, lt.DF){
  dd = dataset %>% 
    filter(name == var.tested) |> 
    mutate(heterogeneity = heterogeneity + abs(min(heterogeneity))) |> 
    filter(!(Design %in% 2:5))
  slope = 
    dd %>%
    arrange(Replicate) %>%
    split(.$Replicate) %>%
    map(~ lm(value ~ heterogeneity, data = .x)) %>%
    map_dfr(~ as.data.frame(t(as.matrix(coef(.))))) 

  compo  = 
    forests[[2]][[1]] %>%
    map(~unique(.x)) %>%
    map(~.x[order(.x)])
  
  results.df = 
    data_frame(
      compo = compo[1:949] %>% 
        map(~ paste0(.x, collapse = '-')) %>%
        unlist(.),
      intercept = slope[,1],
      slope = slope[,2]
    )
  
  results.df$compo = 
    results.df$compo %>% 
    str_replace_all(., 'S1-', 'S01-')
  
  df.compo = 
    data.frame(
      matrix(data=0, nrow = 949, ncol = 11)
    )
  
  for(j in 1:949){
    for(i in 1:11){
      if(results.df[j,1] %>% str_detect(species[i]) == T){
        df.compo[j,i] = 1
      }
    }
  }

  rownames(df.compo) = paste0('R-', 1:949)
  colnames(df.compo) = species
  df.compo = df.compo[,species[order(species)]]

  traits.communities = FD::dbFD(a = df.compo, x = lt.DF)

  results.df = 
    bind_cols(results.df, traits.communities$CWM) %>%
    mutate(FRic = traits.communities$FRic) %>%
    mutate(FDis = traits.communities$FDis) %>%
    mutate(FEve = traits.communities$FEve) %>%
    mutate(RaoQ = traits.communities$RaoQ)
  
  ind <- sample(2, nrow(results.df), 
                replace = TRUE, 
                prob = c(.8, .2))
  train <- results.df[ind==1,]
  test <- results.df[ind==2,]
  
  rf.intercept <- randomForest(
    intercept~ trait.1 + trait.2 + FRic + FDis, 
    data=train, 
    proximity=TRUE
  ) 

  rf.slope <- randomForest(
    slope ~ trait.1 + trait.2 + FRic + FDis, 
    data=train, 
    proximity=TRUE
  ) 

  importance.intercept = 
    importance(rf.intercept)/sum(importance(rf.intercept))
  varImpPlot(rf.intercept,
             sort = T)
  
  importance.slope = 
    importance(rf.slope)/sum(importance(rf.slope))
  varImpPlot(rf.slope,
             sort = T)
  
  r2.intercept = 1-(var((results.df$intercept) - predict(rf.intercept, results.df))/var(results.df$intercept))
  r2.slope = 1-(var((results.df$slope) - predict(rf.slope, results.df))/var(results.df$slope))

  data_frame(y = predict(rf.slope, test), x = test$slope) |> 
  ggplot(aes(x,y)) + 
  geom_point() + 
  geom_smooth()

  data_frame(y = predict(rf.intercept, test), x = test$intercept) |> 
  ggplot(aes(x,y)) + 
  geom_point() + 
  geom_smooth()

  r = 
    list(
    model.out = 
      bind_rows(
        r.intercept = 
          data.frame(
            variable = var.tested, 
            resp = 'intercept', 
            exp = rownames(importance.intercept),
            value = importance.intercept |> unname()
          ),
        r.slope = 
          data.frame(
            variable = var.tested, 
            resp = 'slope', 
            exp = rownames(importance.slope),
            value = importance.slope |> unname()
          )
      ), 
    rf = list(var.tested, rf.intercept, rf.slope, train, r2.intercept, r2.slope)
  )
  return(r)
}

forest.analyses(var.tested = 'Total litterfall (g)',
                dataset = df.1,
                forests = forests,
                species = species,
                lt.DF = lt.DF)

forest.analyses(var.tested = 'Total C loss (g)',
                dataset = df.1,
                forests = forests,
                species = species,
                lt.DF = lt.DF)

results.rf = 
  outputs |>
  map(~ forest.analyses(var.tested = .x, 
                           dataset = df.1,
                           forests = forests,
                           species = species,
                           lt.DF = lt.DF))

summary.rf = 
  results.rf |>
  map_df(~.x[[1]]) |>
  group_by(exp, resp) |>
  summarise(sd.value = sd(value, na.rm = T),
            value = mean(value, na.rm = T))

ggplot(summary.rf, aes(x = exp, y = value, color = resp)) + 
  geom_point() + 
  coord_flip()

d.p = summary.rf |> filter(resp == 'intercept')
d.p$exp = d.p$exp |> 
  str_replace_all('trait.1', 'CWM 1') |>
  str_replace_all('trait.2', 'CWM 2')
d.p$exp = factor(d.p$exp, levels = d.p$exp[order(d.p$value)])

p.intercept = 
  ggplot(data = d.p, aes(x = exp, y = value)) + 
  geom_point(data = d.p, aes(x = exp, y = value),size = 3) +
  geom_errorbar(data = d.p, 
                aes(x = exp,
                    ymin = value - 1.96 * sd.value, 
                    ymax = value + 1.96 * sd.value),
                width = 0, size = .5) + 
  geom_jitter(data = results.rf |> 
                map_df(~.x[[1]]) |> 
                filter(resp == 'intercept') |> 
                mutate(exp = exp |> 
                         str_replace_all('trait.1', 'CWM 1') |>
                         str_replace_all('trait.2', 'CWM 2')),
              aes(x = exp, y = value, color = variable),
              width = .1) + 
  labs(title = "Intercept drivers", y = 'Variable importance', 
       x = 'Explanatory variables', color = '') +
  scale_color_manual(values = colors.reponses, labels = colors.labels) + 
  lims(y = c(0.16, 0.30)) + 
  geom_hline(yintercept = 0.25, linetype = 2) + 
  theme_bw() + 
  coord_flip()
p.intercept

d.s = summary.rf |> filter(resp == 'slope')
d.s$exp = d.s$exp |> 
  str_replace_all('trait.1', 'CWM 1') |>
  str_replace_all('trait.2', 'CWM 2')
d.s$exp = factor(d.s$exp, levels = d.s$exp[order(d.s$value)])

p.slope = 
  ggplot(data = d.s, aes(x = exp, y = value)) + 
  geom_point(data = d.s, aes(x = exp, y = value), size = 3) +
  geom_errorbar(data = d.s,
                aes(x = exp,
                    ymin = value - 1.96 * sd.value,
                    ymax = value + 1.96 * sd.value),
                width = 0, size = .5) +
  geom_jitter(data = results.rf |> 
               map_df(~.x[[1]]) |> 
               filter(resp == 'slope') |> 
               mutate(exp = exp |> 
                        str_replace_all('trait.1', 'CWM 1') |>
                        str_replace_all('trait.2', 'CWM 2')),
             aes(x = exp, y = value, color = variable),
             width = .1) + 
  scale_color_manual(values = colors.reponses, labels = colors.labels) + 
  labs(title = "Slope drivers", y = 'Variable importance', 
       x = '', color = ' ') + 
  lims(y = c(0.16, 0.30)) + 
  geom_hline(yintercept = 0.25, linetype = 2) + 
  coord_flip() + 
  theme_bw() 
p.slope

p = 
  ggpubr::ggarrange(
  p.intercept,
  p.slope, 
  common.legend = T
)
p

d.R2.intercepts = 
  1:6 |> 
  map(~results.rf[[.x]][[2]][[5]]) |>
  unlist() %>% 
  data.frame(name = outputs, R2 = .)

# d.R2.intercepts$name |> levels() = outputs
p.r2.intercept =
  ggplot(data = NULL) + 
  annotate(geom = 'point', y = mean(d.R2.intercepts$R2,na.rm =T), x = 1, size = 3) +
  annotate(geom = 'segment', 
           y = mean(d.R2.intercepts$R2,na.rm =T) + 
             2*sd(d.R2.intercepts$R2,na.rm =T), 
           yend = mean(d.R2.intercepts$R2,na.rm =T) - 
             2*sd(d.R2.intercepts$R2,na.rm =T), 
           x = 1, xend = 1,
           size = .5) +
  geom_point(data = d.R2.intercepts,
             aes(x = 1, y = R2, color = name)) + 
  scale_color_manual(values = colors.reponses, labels = colors.labels) + 
  labs(title = bquote("Intercept models R"^2), y = bquote("R"^2), 
       x = '', color = ' ') + 
  lims(y = c(.1,.25)) + 
  coord_flip() + 
  theme_bw() + 
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank())

d.R2.slopes = 
  1:6 |> 
  map(~results.rf[[.x]][[2]][[6]]) |>
  unlist() %>% 
  data.frame(name = outputs, R2 = .)

p.r2.slope = 
  ggplot(data = NULL) + 
  annotate(geom = 'point', y = mean(d.R2.slopes$R2,na.rm =T), x = 1, size = 3) +
  annotate(geom = 'segment', 
           y = mean(d.R2.slopes$R2,na.rm =T) + 
             2*sd(d.R2.slopes$R2,na.rm =T), 
           yend = mean(d.R2.slopes$R2,na.rm =T) - 
             2*sd(d.R2.slopes$R2,na.rm =T), 
           x = 1, xend = 1,
           size = .5) +
  geom_point(data = d.R2.slopes,
             aes(x = 1, y = R2, color = name)) + 
  scale_color_manual(values = colors.reponses, labels = colors.labels) + 
  labs(title = bquote("Slope models R"^2), y = bquote("R"^2), 
       x = '', color = ' ') + 
  lims(y = c(.1,.25)) + 
  coord_flip() + 
  theme_bw()  + 
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank())