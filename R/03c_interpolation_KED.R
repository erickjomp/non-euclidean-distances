library(tidyverse)
library(sf)
library(gstat)
library(hydroGOF)
library(sp)
Sys.setenv(TZ='GMT')

source("R/functions/funciones_eje_OK.R")

max_ratio_nugget_sill <- 0.25
min_range <- 30
max_range <- 200
model_sel <- "Sph"
type = "KED" #change for "OK
formula_inter <- as.formula(pr~elev)  #change for as.formula(pr~1)


file_dem_gen <- 
  "02_regions_dem_conecs//raster/DEM/ras_dem_reg%s.tif"

file_save_raster_kriging_gen <- 
  "03_raster_interpolation//raster_inter/ras_%s_reg%s.nc"
file_save_raster_avg_kriging_gen <- 
  "03_raster_interpolation//raster_avg_yearly/ras_%s_reg%s.tif"


sf_stas_reg_all <- 
  st_read("01_data_byregion//shp_stations//sf_stas_reg.shp")

sf_stas_regbuf_all <- 
  st_read("01_data_byregion//shp_stations//sf_stas_regbuf.shp")

df_m_all <- readRDS("00_data/data_pr_validated&selected/df_m.rds")


id_regions <- 1:12

for (id_reg in id_regions){
  sf_stas_regbuf <- sf_stas_regbuf_all %>% 
    filter(region == id_reg)
  df_m <- df_m_all %>% 
    dplyr::select(dates, sf_stas_regbuf$namcod)
  
  message(paste("REGION ", id_reg))
  
  ras_dem <- raster::raster(sprintf(file_dem_gen,id_reg))
  # stars_dem <- st_as_stars(ras_dem)
  sp_dem <- as(ras_dem,'SpatialGridDataFrame')
  names(sp_dem) <- "elev"
  
  ras_kriging <- 
    lapply(1:nrow(df_m), function(i_mes){
      message(paste("\r",i_mes),appendLF = F)
      values_pr <- df_m[i_mes,-1] %>% as.numeric()
      sf_stas_gstat_0 <- sf_stas_regbuf %>% 
        mutate(pr = values_pr)
      sp_stas_gstat_0 <- as_Spatial(sf_stas_gstat_0)
      
      sp_stas_gstat <- 
        sp_stas_gstat_0[!is.na(sp_stas_gstat_0$pr),]
      
      fixed_values <- c(NA,NA,NA)
      
      iter = 0
      
      #### loop to fit variogram
      vario_af <- NULL
      while(TRUE){
        
        try(silent =T,
            vario_af <-
              autofitVariogram_v2(
                formula = formula_inter,
                sp_stas_gstat,
                model = model_sel,
                fix.values = fixed_values,
                diagonal = 314.2857, # thats for 110 km
                #diagonal = 306.2099,#342.8571,#285.7143,
                # 110 km : 314.2857
                # start_vals = c(0,50,NA),
                miscFitOptions = 
                  list(merge.small.bins = F)
              )
            # min.np.bin = 3)),
            
        )
        
        if (!is.null(vario_af)){
          nugget <- vario_af$var_model$psill[1]
          sill <- vario_af$var_model$psill[2] + nugget
          range <- vario_af$var_model$range[2]
          
          ratio_nugget_sill <- (nugget/sill)
          ratio_nugget_sill <-  ifelse(sill == 0, 0, 
                                       ratio_nugget_sill)
          
          condition_nugget <- 
            ratio_nugget_sill  <= max_ratio_nugget_sill
          condition_range_1 <- range >= min_range
          condition_range_2 <- range <= max_range
        }
        else{
          condition_nugget <- F
          condition_range_1 <- T
          condition_range_2 <- T
        }
        if(condition_nugget & condition_range_1 & 
           condition_range_2){
          break  #if conditions are satisfied , the loop ends
        }
        iter <- iter +1
        if(iter > 10){
          message(paste0(" ",i_mes,"-10it"))
          break  #if conditions are satisfied , the loop ends
        }
        
        if ((!condition_range_1 || !condition_range_2) & (
          !is.na(fixed_values[1]) || condition_nugget
        )){ 
          # the 2nd condition is to ensure tryinf nugget first 
          if(!condition_range_1){
            fixed_values[2] <- min_range
          }
          
          if(!condition_range_2){
            fixed_values[2] <- max_range
          }
          if (!is.null(vario_af)){
            vario_af$var_model$range[2] <- fixed_values[2]
          }
          
        } 
        
        if (!condition_nugget){
          if (!is.null(vario_af)){
            fixed_values[1] <- 
              min(sd(sp_stas_gstat$pr),0.25*sill)
            vario_af$var_model$psill[1] <- fixed_values[1] 
          } 
          else{
            fixed_values[1] <- sd(sp_stas_gstat$pr)
          } 
        } 
        
      }
      sp_sim <- 
        krige(
          formula_inter, 
          locations = sp_stas_gstat,
          newdata = sp_dem, 
          model = vario_af$var_model,
          debug.level = 0)
      if (all(sp_stas_gstat$pr == 0)){
        sp_sim$var1.pred[is.finite(ras_dem[])] <- 0
      }
      sp_sim[1] %>% terra::rast()
      
      
    }) %>% do.call(c,.)
  
  ras_avg_yearly_kriging <- 
    ras_kriging %>% (function(x) mean(x)*12)
  
  raster::writeRaster(
    ras_kriging,
    sprintf(file_save_raster_kriging_gen, type, id_reg),
    overwrite=TRUE)
  
  raster::writeRaster(
    ras_avg_yearly_kriging,
    sprintf(file_save_raster_avg_kriging_gen, type, id_reg),
    overwrite=TRUE)
}
