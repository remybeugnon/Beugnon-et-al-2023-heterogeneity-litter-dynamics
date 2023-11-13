rm(list = ls())
set.seed(1603)

library(tidyverse)
library(purrr)

df.2.sp = readRDS('2-2-1_2sp-predictions.RData')[[2]] |> 
  mutate(sp = 2)

df.4.sp = readRDS('2-2-2_4sp-predictions.RData')[[2]] |> 
  mutate(sp = 4)

df.9.sp = readRDS('2-2-3_9sp-predictions')[[2]] |> 
  mutate(sp = 9)


df = df.2.sp |> 
  add_row(df.4.sp) |> 
  add_row(df.9.sp)

response = df$name |> unique()
r = response[5]

p = 
  ggplot(data = df, 
         aes(x = heterogeneity, 
             y = value, 
             color = factor(sp))) + 
  geom_smooth(method="loess",
              linewidth = .5, 
              alpha = .1) + 
  labs(x = "Tree spatial heterogeneity", 
       y = 'Value', 
       color = 'Tree Species Richness') +
  facet_wrap(vars(name), scales = 'free') + 
  theme_bw() + 
  scale_color_viridis_d() + 
  theme(legend.position = 'top')

df.3 = df |> 
  filter(name %in% c('Mean litterfall (g/dm2)', 
                     'Litterfall SD (g/dm2)',
                     'Mean litter species richness (#/dm2)'))

p.litter =
  ggplot(data = df.3, 
         aes(x = heterogeneity, 
             y = value)) + 
  geom_smooth(data = df.3 |> 
                 filter(!(sp == 9 & Design %in% 2:5)),
              method="loess", 
              level=0.95,
              size = .5, 
              color = 'darkgreen',
              aes(linetype = factor(sp))) + 
  facet_wrap(vars(name), scales = 'free', ncol =  3) + 
  scale_y_continuous(labels = scales::label_number_si()) +
  scale_x_continuous(breaks = c(-1100, -500, -180,0), 
                     labels = c('Blocks',
                                '2-sp. random',
                                '4-sp. random',
                                '9-sp. random')) + 
  scale_linetype_manual(values = c(1:3)) +
  labs(x = "Tree species spatial heterogeneity", 
       y = "Value", 
       title = "A. Litterfall") +
  theme_bw() + 
  theme(legend.position = 'none', 
        axis.text.x = element_text(angle = 45, hjust = 1))

p.litter  

df.4 = df |> filter(name %in% c('Mean carbon decomposition rate (%)', 
                                  'Carbon decomposition rate SD (%)',
                                  'Total C loss (g)')) 
df.4$value[df.4$name == 'Total C loss (g)'] =
  df.4$value[df.4$name == 'Total C loss (g)'] / (17*17*100)
df.4$name = as.character(df.4$name)
df.4$name[df.4$name == 'Total C loss (g)'] = 'Mean C loss (g/dm2)'
df.4$name = factor(df.4$name,
                   levels = c('Mean carbon decomposition rate (%)', 
                              'Carbon decomposition rate SD (%)',
                              'Mean C loss (g/dm2)'))

p.decomp =
  ggplot(data = df.4, 
         aes(x = heterogeneity, 
             y = value)) + 
  geom_smooth(data = df.4 |> 
                filter(!(sp == 9 & Design %in% 2:5)),
              method="loess", 
              level=0.95,
              size = .5,
              color = 'brown',
              aes(linetype = factor(sp))) + 
  facet_wrap(vars(name), scales = 'free', ncol =  3) + 
  scale_y_continuous(labels = scales::label_number_si()) +
  scale_x_continuous(breaks = c(-1100, -500, -180,0), 
                     labels = c('Blocks',
                                '2-sp. random',
                                '4-sp. random',
                                '9-sp. random')) + 
  scale_linetype_manual(values = c(1:3)) +
  labs(x = "Tree species spatial heterogeneity", 
       y = "Value", 
       title = "B. Decomposition") +
  theme_bw() + 
  theme(legend.position = 'none', 
        axis.text.x = element_text(angle = 45, hjust = 1))

p.decomp

df.5 = df.4 |> 
  filter(name == 'Mean carbon decomposition rate (%)') |> 
  left_join(
    data.frame(
      sp = c(2,4,9),
      hetero.block = c(-1056.3684, -1020.3684, -948.3684211),
      hetero.random = c(-496.3684, -180.3684, -0.3684211), 
      hetero.mean = c(-746.3684,-512.3684,-476.3684211)
    )
    ) |>
  filter(round(heterogeneity) == round(hetero.block) |
           round(heterogeneity) == round(hetero.mean) |
           round(heterogeneity) == round(hetero.random))

df.5$hetero = NA
df.5$hetero[round(df.5$heterogeneity) == round(df.5$hetero.block)] = 'Block'
df.5$hetero[round(df.5$heterogeneity) == round(df.5$hetero.mean)] = 'Average'
df.5$hetero[round(df.5$heterogeneity) == round(df.5$hetero.random)] = 'Random'

df.5$hetero = factor(df.5$hetero, levels = c('Block', 'Average', 'Random'))

d = 
  df.5 |> 
  group_by(hetero, sp) |> 
  summarise(m = mean(value),
            se = sd(value))
p.div = 
  ggplot(data = df.5) + 
  geom_point(data = d, aes(x = sp, y = m, color = hetero)) +
  geom_line(data = d, aes(x = sp, y = m, color = hetero)) +
  geom_jitter(data = df.5,
              aes(x = sp, y = value, group = hetero, color = hetero),
              alpha = .1, size = .2, width = .1) +
  annotate(geom = 'segment', x = 9.3, xend = 9.3, y = 35, yend = 48, 
           arrow = arrow(length = unit(.2, 'cm'), type = 'closed'), 
           size = 1) +
  annotate(geom = 'text', x = 10.2, y = 42, angle = -90, label = 'Heterogeneity') +
  labs(x = "Tree species richness", 
       y = "Decomposition rate (%)",
       color = "Heterogeneity",
       title = "C. BEF relationship") +
  scale_x_continuous(trans = 'log2', breaks = c(2,4,9)) +
  scale_color_manual(values = c('gray20', 'darkblue', 'darkred'))+
  theme_bw()

p.div

df.6 = df.4 |> 
  filter(name == 'Carbon decomposition rate SD (%)') |> 
  left_join(
    data.frame(
      sp = c(2,4,9),
      hetero.block = c(-1056.3684, -1020.3684, -948.3684211),
      hetero.random = c(-496.3684, -180.3684, -0.3684211), 
      hetero.mean = c(-746.3684,-512.3684,-476.3684211)
    )
  ) |>
  filter(round(heterogeneity) == round(hetero.block) |
           round(heterogeneity) == round(hetero.mean) |
           round(heterogeneity) == round(hetero.random))

df.6$hetero = NA
df.6$hetero[round(df.6$heterogeneity) == round(df.6$hetero.block)] = 'Block'
df.6$hetero[round(df.6$heterogeneity) == round(df.6$hetero.mean)] = 'Average'
df.6$hetero[round(df.6$heterogeneity) == round(df.6$hetero.random)] = 'Random'

df.6$hetero = factor(df.6$hetero, levels = c('Block', 'Average', 'Random'))

dd = 
  df.6 |> 
  group_by(hetero, sp) |> 
  summarise(m = mean(value),
            se = sd(value))
p.div.2 = 
  ggplot(data = df.6) + 
  geom_point(data = dd, aes(x = sp, y = m, color = hetero)) +
  geom_line(data = dd, aes(x = sp, y = m, color = hetero)) +
  geom_jitter(data = df.6,
              aes(x = sp, y = value, group = hetero, color = hetero),
              alpha = .1, size = .2, width = .1) +
  annotate(geom = 'segment', x = 9.9, xend = 9.9, y = 10, yend = 4, 
           arrow = arrow(length = unit(.2, 'cm'), type = 'closed'), 
           size = 1) +
  annotate(geom = 'text', x = 10.5, y = 7, angle = -90, label = 'Heterogeneity') +
  labs(x = "Tree species richness", 
       y = "Decomposition rate SD(%)",
       color = "Heterogeneity",
       title = "") +
  scale_x_continuous(trans = 'log2', breaks = c(2,4,9)) +
  scale_color_manual(values = c('gray20', 'darkblue', 'darkred'))+
  theme_bw()

p.div.2

Fig.3 = 
  ggpubr::ggarrange(
    ggpubr::ggarrange(
      p.litter,
      p.decomp,
      nrow = 2,
      align = 'hv'
      ),
    ggpubr::ggarrange(
      p.div,
      p.div.2,
      ncol = 2,
      align = 'hv', 
      common.legend = T
    ),
    nrow = 2,
    heights = c(.6,.3)
    )

Fig.3

ggsave(plot = Fig.3, 
       filename = '3-2_Fig3.png',
       height = 25, width = 20, units = 'cm')

ggsave(plot = Fig.3, 
       filename = '3-2_Fig3.svg',
       height = 25, width = 20, units = 'cm')