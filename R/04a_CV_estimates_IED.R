### IED method replacing the old version ####
# this version is working with 
library(terra)
library(tidyverse)
library(sf)
library(hydroGOF)
Sys.setenv(TZ='GMT')
# _in from input

source("R/functions/funciones_NED_interpolation_2.R")


lambdas <-
  read.csv("02_regions_dem_conecs//input/lambdas.csv")[[1]]

file_list_mat_dists_loc_gen <-
  "04_cross_validation//mat_NED1_stas/list_mat_NED1_lambdas_reg%s.rds"

sf_stasbuf_all <-
  st_read("02_regions_dem_conecs//shp/sf_stas_regbuf.shp")

df_m_all <-
  readRDS("00_data/data_pr_validated&selected//df_m.rds")

file_dem_gen <- 
  "02_regions_dem_conecs//raster/DEM/ras_dem_reg%s.tif"

list_best_lambdas_CV <- list()
list_df_sims_regs <- list()

# list_df_sims_regs <-
# lapply(1:12, function(id_reg) {
for (id_reg in 1:12){
  message(paste0("REGION ", id_reg))
  sf_stasbuf <- sf_stasbuf_all %>% filter(region == id_reg)
  
  df_m <- df_m_all %>%
    dplyr::select(dates, sf_stasbuf$namcod)
  
  file_list_mat_dists_loc <- 
    sprintf(file_list_mat_dists_loc_gen, id_reg)
  list_mat_dists_loc_in <- 
    readRDS(file_list_mat_dists_loc)
  
  # list_mat_dists_loc <- list_mat_dists_loc %>% 
  #   lapply(function(x) x /10000)
  
  ista <- (sf_stasbuf$in_region == 1)
  
  ras_dem <- rast(
    sprintf(file_dem_gen, id_reg)
  )
  
  # very old function
  list_out <- 
    NED_interpolation(formula = pr~elev, 
                      locations = sf_stasbuf,
                      data = df_m, #1:26
                      new_data = ras_dem,
                      optim_lambda = 1,
                      # list_mat_dists = list_mat_dists_in,
                      list_mat_dists_loc = list_mat_dists_loc_in[1],
                      idp = 2,
                      # ratio_dif_WS_min = NULL,#0.0001,
                      min_0 = T,
                      i_target_stas = ista,
                      out_process_data = T,
                      cross_validation = T,
                      lambdas = lambdas[1])
  gc()
  
  df_sim <- data.frame(list_out$mat_sim)
  names(df_sim) <- names(df_m)[-1]
  df_sim <- cbind(dates = df_m$dates, df_sim)
  

  list_df_sims_regs[[id_reg]] <- 
    df_sim[,c(1, which(ista) + 1)]
  list_best_lambdas_CV[[id_reg]] <- list_out$best_lambdas_CV
  # append(list_df_sims_regs,list(list_out$best_lambdas_CV))
}



## saving
saveRDS(
  list_df_sims_regs, 
  file = "04_cross_validation//estimates/list_df_m_inregs_est_IED.rds")
