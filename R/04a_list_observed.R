library(tidyverse)
library(sf)
library(gstat)
library(hydroGOF)
Sys.setenv(TZ='GMT')



sf_stas_reg_all <- 
  st_read("01_data_byregion//shp_stations//sf_stas_reg.shp")

sf_stas_regbuf_all <- 
  st_read("01_data_byregion/shp_stations//sf_stas_regbuf.shp")

df_m_all <- readRDS("00_data/data_pr_validated&selected//df_m.rds")

list_df_m_inregs <- lapply(1:12, function(id_reg){
  namcods_sel <-
    sf_stas_reg_all %>%
    filter(region == id_reg) %>% .$namcod
  df_m_all %>% select(dates,namcods_sel)
})


## saving
saveRDS(
  list_df_m_inregs, 
  file = "04_cross_validation//observed/list_df_m_inregs_obs.rds")
