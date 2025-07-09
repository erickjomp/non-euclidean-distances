# library(raster)
library(terra)
library(tidyverse)
library(sf)
library(hydroGOF)
# library(sp)


source("R/functions/funciones_NED_interpolation.R")
source("R/functions/funciones_NED_interpolation_2.R")

# input -------
lambdas <-
  read.csv("02_regions_dem_conecs//input/lambdas.csv")[[1]]

file_list_mat_dists_gen <-
  "03_raster_interpolation//mat_NED3_raster/list_mat_NED3_lambdas_reg%s.rds"

file_list_mat_dists_loc_gen <-
  "04_cross_validation//mat_NED3_stas/list_mat_NED3_lambdas_reg%s.rds"

sf_stasbuf_all <-
  st_read("02_regions_dem_conecs//shp/sf_stas_regbuf.shp")

df_m_all <-
  readRDS("00_data/data_pr_validated&selected//df_m.rds")

file_dem_gen <-
  "02_regions_dem_conecs//raster/DEM/ras_dem_reg%s.tif"


# output ----
file_save_raster_NED_gen <-
  "03_raster_interpolation/raster_inter/ras_NED3_reg%s.nc"
file_save_raster_avg_NED_gen <-
  "03_raster_interpolation/raster_avg_yearly/ras_NED3_reg%s.tif"
file_df_lambdas <- 
  "03_raster_interpolation/df_lambdas_NED3.rds"
file_list_df_error_lambdas <- 
  "03_raster_interpolation/list_df_error_lambdas_NED3.rds"

df_lambdas_NED3 <- data.frame()
list_df_error_lambdas_NED3 <- list()

for (id_reg in 1:12) {
  message(paste0("REGION ", id_reg))
  sf_stasbuf <- sf_stasbuf_all %>% filter(region == id_reg)
  
  df_m <- df_m_all %>%
    dplyr::select(dates, sf_stasbuf$namcod)
  
  file_list_mat_dists <-
    sprintf(file_list_mat_dists_gen, id_reg)
  list_mat_dists <-
    readRDS(file_list_mat_dists)
  
  file_list_mat_dists_loc <-
    sprintf(file_list_mat_dists_loc_gen, id_reg)
  list_mat_dists_loc <-
    readRDS(file_list_mat_dists_loc)
  
  # list_mat_dists_loc <- list_mat_dists_loc %>%
  #   lapply(function(x) x /10000)
  
  ista <- (sf_stasbuf$in_region == 1)
  
  ras_dem <- rast(sprintf(file_dem_gen, id_reg))
  
  # is_fin_cells <- is.finite(values(ras_dem)[,1])
  list_mat_dists <- lapply(list_mat_dists, function(mat_dists){
    mat_dists[is.na(mat_dists)] <- Inf
    mat_dists
  })
  
  # very old function
  list_out <-
    NED_interpolation(
      formula = pr ~ elev,
      locations = sf_stasbuf,
      data = df_m,
      #1:26
      new_data = ras_dem,
      list_mat_dists = list_mat_dists[1:50],
      list_mat_dists_loc = list_mat_dists_loc[1:50],
      idp = 2,
      # ratio_dif_WS_min = NULL,#0.0001,
      min_0 = T,
      optim_lambda = 1,
      lambdas = lambdas[1:50],
      i_target_stas = ista,
      out_process_data = T
    )
  gc()
  ras_NED <- list_out$newdata_estimates
  names(ras_NED) <- df_m$dates
  # time(ras_NED) <- df_m$dates
  
  ras_avg_yearly_NED <-
    ras_NED %>% (function(x)
      mean(x) * 12)
  
  terra::writeRaster(ras_NED,
                     sprintf(file_save_raster_NED_gen, id_reg),
                     overwrite = TRUE)
  #terra::writeCDF
  
  terra::writeRaster(ras_avg_yearly_NED,
                     sprintf(file_save_raster_avg_NED_gen, id_reg),
                     overwrite = TRUE)
  
  df_lambdas_NED3 <- 
    rbind(df_lambdas_NED3, 
          data.frame(id_reg = id_reg, lambda = list_out$lambda))
  list_df_error_lambdas_NED3[[id_reg]] <- 
    list_out$WS
}

saveRDS(df_lambdas_NED3, file_df_lambdas)
saveRDS(list_df_error_lambdas_NED3,
        file_list_df_error_lambdas)
