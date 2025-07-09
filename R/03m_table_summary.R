library(tidyverse)
library(terra)
library(sf)

file_ras_avg_yearly_gen <- 
  "03_raster_interpolation//raster_avg_yearly/ras_%s_reg%s.tif"

methods <- c("IED","NED1","NED2","NED3")
labs_methods <- 
  c("IED","IED - NED1","IED - NED2","IED - NED3")
# methods <- c("IDW","OK","IDW-V","KED","IED","NED1","NED2","NED3")
# labs_methods <- 
#   c("IDW","OK","IDW-V","KED","IED","IED - NED1","IED - NED2","IED - NED3")
id_regions <- 1:12

id_reg = 12

sf_regions <- 
  st_read("00_data//shp/regions.shp")

df_res <- 
lapply(methods, function(method){
  values_avg_pry <- sapply(id_regions, function(id_reg){
    file_ras <- 
      sprintf(file_ras_avg_yearly_gen, method, id_reg)
    ras_pry <- rast(file_ras)
    extract(ras_pry,sf_regions[id_reg,],fun = "mean")[,2]
  })
  data.frame(method = method,
             lab_method = labs_methods[methods == method],
             region = id_regions, values_avg_pry)
}) %>% do.call(rbind,.)

df_res$param <- ""
df_res$param_name <- NA

df_lambdas_NED1 <- 
  readRDS("03_raster_interpolation/df_lambdas_NED1.rds")
df_res$param[df_res$method == "NED1"] <- 
  df_lambdas_NED1$lambda
  # sprintf("($\\lambda$ = %s)", df_lambdas_NED1$lambda)
df_res$param_name[df_res$method == "NED1"] <- "$\\lambda$"

df_lambdas_NED2 <- 
  readRDS("03_raster_interpolation//df_lambdas_NED2.rds")
df_res$param[df_res$method == "NED2"] <- 
  df_lambdas_NED2$lambda
  # sprintf("($\\lambda$ = %s)", df_lambdas_NED2$lambda)
df_res$param_name[df_res$method == "NED2"] <- "$\\lambda$"

df_lambdas_NED3 <- 
  readRDS("03_raster_interpolation//df_lambdas_NED3.rds")
df_res$param[df_res$method == "NED3"] <- 
  df_lambdas_NED3$lambda
  # sprintf("($\\lambda$ = %s)", df_lambdas_NED3$lambda)
df_res$param_name[df_res$method == "NED3"] <- "$\\lambda$"
#"(\u03bb = %s)"


library(kableExtra)
df_res <- 
df_res %>% mutate(info = ifelse(param != "", 
                            paste0(round(values_avg_pry),
                                "\n", param ) %>% linebreak(align = "c"),
                            round(values_avg_pry)))


df_res_trans <- 
  df_res %>% 
  mutate(values_avg_pry = round(values_avg_pry)) %>% 
  select(region, lab_method, values_avg_pry) %>%
  mutate(lab_method = factor(lab_method, levels = unique(lab_method))) %>% 
  spread("lab_method", "values_avg_pry") %>% 
  rename(Region = region)

library(kableExtra)
df_table_latex <- 
 df_res_trans %>% 
  kbl("latex",booktabs = T,escape = F,
      # linesep = 
      # col.names = names_cols[-1],
      align = "c") %>% 
  pack_rows(index = auto_index(rep(
    c("Western slope","Eastern slope", "Titicaca basin"), c(6,4,2))
    ))

save_kable(df_table_latex,
           "article/tables/table_prys_raster-inter.tex",
           keep_tex = T)


#### table parameters ####
df_res_param <- 
  df_res %>% filter(!is.na(param_name)) %>% 
  select(region,lab_method, param_name,param) %>% 
  mutate(lab_method = factor(lab_method,
                             levels = unique(lab_method))) %>% 
  rename(Region = region)

df_res_param <- 
df_res_param %>% select(-3) %>% 
  spread("lab_method","param")

names_met <- names(df_res_param)[-1]
names(df_res_param)[-1] <- rep(c("$\\lambda$"),c(3))

names_header <- rep(1,4)
names(names_header) <- c(" ", names_met)

df_res_param_latex <- 
  df_res_param  %>% 
  kbl("latex",booktabs = T,escape = F,
      # linesep = 
      # col.names = names_cols[-1],
      align = "c") %>% 
  pack_rows(index = auto_index(rep(
    c("Western slope","Eastern slope", "Titicaca basin"), c(6,4,2))
  )) %>% 
    add_header_above(names_header)

save_kable(df_res_param_latex,
           "article/tables/table_params_raster-inter.tex",
           keep_tex = T)

