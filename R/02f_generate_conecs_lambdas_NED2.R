library(tidyverse)
library(sf)
library(terra)

lambdas <- read.csv("02_regions_dem_conecs//input/lambdas.csv")[[1]]

file_conec_gen <- 
  "02_regions_dem_conecs/df_conecs_NED2/df_conec_reg%s.rds"
file_conec_gen_new <- 
  "02_regions_dem_conecs/df_conecs_lambdas_NED2//df_conec_reg%s.rds"

file_ras_dem_reg_gen <- 
  "02_regions_dem_conecs/raster/DEM/ras_dem_reg%s.tif"
sf_stas_regbuf_all <- 
  st_read("01_data_byregion/shp_stations//sf_stas_regbuf.shp")
sf_stas_reg_all <- 
  st_read("01_data_byregion/shp_stations//sf_stas_reg.shp")


id_regions <- sort(unique(sf_stas_reg_all$region))


lapply(id_regions, function(id_reg){
  message(paste("REGION",id_reg))
  file_conec <- sprintf(file_conec_gen, id_reg)
  df_conec <- readRDS(file_conec)
  
  for (lambda in lambdas){
    df_conec[[as.character(lambda)]] <- 
      sqrt(df_conec$dif_h ^ 2 + (lambda * df_conec$dif_v) ^ 2)
  }
  
  file_conec_new <- sprintf(file_conec_gen_new, id_reg)
  saveRDS(df_conec,file_conec_new)
  
})
