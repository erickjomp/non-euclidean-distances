library(tidyverse)
library(lubridate)
library(latex2exp)
library(cowplot)
devtools::load_all("../my_packages/my_themes/")

### inputs 

file_listbeta1 <- 
  "03_raster_interpolation//list_beta1.rds"
base_size_figs <- 7



### PROCESS ####
seasons_labs = c("DJF","MAM", "JJA","SON")
list_beta1 = readRDS(file_listbeta1)


df_beta1 = data.frame(list_beta1) 
names(df_beta1) = 1:12

df_beta1 = df_beta1 %>% 
  mutate(datetime = seq.POSIXt(as.POSIXct("1995-01-01"), as.POSIXct("2019-12-01"),"months"),.before = 1)

df_beta1 = gather(df_beta1,key = region, value = beta1,-1)
df_beta1 = df_beta1 %>% mutate(region = factor(region, levels = 1:12))

# seasons and geog domains
df_beta1 = df_beta1 %>% 
  mutate(domain =ifelse(as.numeric(region) <= 6, "Western\nSlope", 
                        ifelse(as.numeric(region) <= 10, "Eastern\nSlope", "Titicaca\nBasin"))) %>% 
  mutate(domain = factor(domain, levels = unique(domain))) %>% 
  mutate(imonth = lubridate::month(datetime)) %>% 
  mutate(season = case_match(imonth,  # case_when %
                             c(12,1,2)~ seasons_labs[1],
                             c(3,4,5)~ seasons_labs[2],
                             c(6,7,8)~ seasons_labs[3],
                             c(9,10,11)~ seasons_labs[4]))  %>% 
  mutate(season = factor(season, levels = seasons_labs))

(
  plot1 = 
    ggplot(df_beta1, aes(x = region, y = beta1, fill = domain)) + 
    geom_boxplot( linewidth = 0.4, outlier.size = 0.5) +
    # geom_violin(linewidth = 0.4) + 
    # ggh4x::facet_grid2(id_reg~indice,  scales="free_y",
    #                    independent = "y") + 
    my_ggtheme(base_size = base_size_figs)  +
    scale_y_continuous(expand = expansion(mult = 0.05)) +
    labs(x = "Region", y = TeX("$beta_1$ (mm/m)"), fill = "Geographic domain") +
    theme(
      panel.grid.major.y = element_line(color = "grey80", linewidth = 0.5,linetype = "dashed"),
      # panel.grid.minor = element_line(color = "grey90", size = 0.25),
      # panel.background = element_rect(fill = "white")
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.box = "vertical",legend.title.position = "top",
      # legend.title.align = 0.5 
      legend.title = element_text(hjust=0.5)
    ) + 
    # coord_fixed(ratio = 30) +
    # theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    theme(axis.text.x = element_text(size = base_size_figs)) +
    
    theme(strip.text.x = element_text(size = base_size_figs*1.2), #*1.2
          strip.text.y = element_text(size = base_size_figs*1.2)) #*1.2
)


(
  plot2 = 
    ggplot(df_beta1, aes(x = region, y = beta1, fill = domain)) + 
    geom_boxplot( linewidth = 0.4, outlier.size = 0.5) +
    facet_wrap(~season) + 
    # geom_violin(linewidth = 0.4) + 
    # ggh4x::facet_grid2(id_reg~indice,  scales="free_y",
    #                    independent = "y") + 
    my_ggtheme(base_size = base_size_figs)  +
    scale_y_continuous(expand = expansion(mult = 0.025)) +
    labs(x = "Region", y = TeX("$beta_1$ (mm/m)"), fill = "Geographic domain") +
    theme(
      panel.grid.major.y = element_line(color = "grey80", linewidth = 0.5,linetype = "dashed"),
      # panel.grid.minor = element_line(color = "grey90", size = 0.25),
      # panel.background = element_rect(fill = "white")
      legend.position = "bottom"
    ) + 
    # theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    theme(axis.text.x = element_text(size = base_size_figs),
          legend.position = "none") +
    theme(strip.text.x = element_text(size = base_size_figs*1.2), #*1.2
          strip.text.y = element_text(size = base_size_figs*1.2)) #*1.2
)


(
  fig_c <- 
    plot_grid(plot1 + coord_fixed(ratio = 30), plot2,
              rel_heights = c(0.5,1),
              rel_widths = c(1.2,2),
              nrow = 1, labels = c("A","B"))
)


ggsave(plot1, filename = "03_raster_interpolation/fig_boxplot_beta1.pdf",
       width = 90,height = 90, units = "mm")
ggsave(plot1,  filename = "article/figs/fig_boxplot_beta1.pdf",
       width = 90,height = 90, units = "mm")


ggsave(fig_c,  filename = "03_raster_interpolation/fig_boxplot_beta1_seasons.pdf",
       width = 190,height = 100, units = "mm")
ggsave(fig_c,  filename = "article/figs/fig_boxplot_beta1_seasons.pdf",
       width = 190,height = 100, units = "mm")

# ggsave("article/document/figs/fig_boxplot_beta1.pdf",
#        width = 140,height = 180, units = "mm")

