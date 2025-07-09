library(tidyverse)
library(sf)
library(gstat)
library(hydroGOF)
Sys.setenv(TZ='GMT')

file_dem_gen <- 
  "02_regions_dem_conecs//raster/DEM/ras_dem_reg%s.tif"

file_save_raster_IED_gen <- 
  "03_raster_interpolation//raster_inter/ras_IED_reg%s.nc"
file_save_raster_avg_IED_gen <- 
  "03_raster_interpolation//raster_avg_yearly/ras_IED_reg%s.tif"

file_save_listbeta1 <- 
  "03_raster_interpolation//list_beta1.rds"

#### READING ####
sf_stas_reg_all <- 
  st_read("01_data_byregion//shp_stations//sf_stas_reg.shp")

sf_stas_regbuf_all <- 
  st_read("01_data_byregion//shp_stations//sf_stas_regbuf.shp")

df_m_all <- readRDS("00_data/data_pr_validated&selected/df_m.rds")


#### PROCESSS ####
list_beta1 = list()    # list of vectors of beta1 values for all months, each vector corresponds to a region
 
for (id_reg in 1:12){
  sf_stas_regbuf <- sf_stas_regbuf_all %>% 
    filter(region == id_reg)
  df_m <- df_m_all %>% 
    dplyr::select(dates, sf_stas_regbuf$namcod)
  
  message(paste("REGION ", id_reg))
  
  ras_dem <- raster::raster(sprintf(file_dem_gen,id_reg))
  # stars_dem <- st_as_stars(ras_dem)
  sp_dem <- as(ras_dem,'SpatialGridDataFrame')
  names(sp_dem) <- "elev"
  
  beta1 = rep(NA,nrow(df_m)) ### linear_regression slope
  
  ras_IED <- 
  lapply(1:nrow(df_m), function(i_mes){
    message(paste("\r",i_mes),appendLF = F)
    values_pr <- df_m[i_mes,-1] %>% as.numeric()
    sf_stas_gstat_0 <- sf_stas_regbuf %>% 
      mutate(pr = values_pr)
    sp_stas_gstat_0 <- as_Spatial(sf_stas_gstat_0)
    
    sp_stas_gstat <- 
      sp_stas_gstat_0[!is.na(sp_stas_gstat_0$pr),]
    
    ########
    lm_model <- lm(pr~elev, sp_stas_gstat@data)
    beta1[i_mes] <<- lm_model$coefficients["elev"]
    
    sp_stas_gstat$pr_res <-
      sp_stas_gstat$pr -
      predict(lm_model, sp_stas_gstat@data)

    spgrid_residuals <-
      idw(
        pr_res~1,
        locations = sp_stas_gstat,
        newdata = sp_dem,
        debug.level = 0
      )[1]

    spgrid_lm <- sp_dem
    names(spgrid_lm) <- "pr_lm"
    spgrid_lm$pr_lm[!is.na(sp_dem$elev)] <-
      predict(lm_model,sp_dem)

    ras_sim <- terra::rast(spgrid_lm) +
      terra::rast(spgrid_residuals)

    terra::values(ras_sim)[,1] <-
      ifelse(terra::values(ras_sim)[,1] < 0, 0,
             terra::values(ras_sim)[,1])

    ras_sim
  }) %>% do.call(c,.)
  
  list_beta1 = append(list_beta1,list(beta1))
  
  names(ras_IED) <- df_m$dates
  
  ras_avg_yearly_IED <-
    ras_IED %>% (function(x) mean(x)*12)

  raster::writeRaster(
    ras_IED,
    sprintf(file_save_raster_IED_gen, id_reg),
    overwrite=TRUE)

  raster::writeRaster(
    ras_avg_yearly_IED,
    sprintf(file_save_raster_avg_IED_gen, id_reg),
    overwrite=TRUE)
}


saveRDS(list_beta1, file_save_listbeta1)
