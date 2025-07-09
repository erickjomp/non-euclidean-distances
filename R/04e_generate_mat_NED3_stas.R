library(tidyverse)
library(sf)
library(terra)
source("R/functions/funciones_perfiles_NED3.R")
source("R/functions/funciones_perfiles.R") 
Sys.setenv(TZ='GMT')

# input ------
file_sf_nodes_gen <-
  "02_regions_dem_conecs//shp/nodes/sf_nodes_reg%s.shp"

sf_stasbuf_all <- st_read("02_regions_dem_conecs//shp/sf_stas_regbuf.shp")

lambdas <-
  read.csv("02_regions_dem_conecs//input/lambdas.csv")[[1]]

file_dem_gen <- 
  "02_regions_dem_conecs//raster/DEM/ras_dem_reg%s.tif"

# files output -----
file_save_list_dists_gen <- 
  "04_cross_validation//mat_NED3_stas/list_mat_NED3_lambdas_reg%s.rds"


for (id_reg in 1:12){
  message(paste0("REGION ", id_reg))

  ras_dem <- rast(sprintf(file_dem_gen,id_reg))
  
  file_sf_nodes <-
    sprintf(file_sf_nodes_gen, id_reg)
  sf_nodes <-
    st_read(file_sf_nodes,quiet = T)
  
  sf_stasbuf <-
    sf_stasbuf_all %>% filter(region == id_reg)
  
  id_nodes_stas <- sf_stasbuf$id_node
  
  n_stas <- nrow(sf_stasbuf)
  
  mat_dists <- matrix(NA,ncol = n_stas, nrow = n_stas)
  diag(mat_dists) <- 0
  
  list_dists_stas_NED3 <- list()
  for (i in 1:length(lambdas)){
    list_dists_stas_NED3[[i]] <- mat_dists
  }
  names(list_dists_stas_NED3) <- lambdas
  
  for (i in 1:(n_stas-1)){
    message("Working in ", i ,"/", n_stas,"\r",appendLF = F)
    i_other_stas <- (i+1):n_stas
    mat_dists_NED3 <- 
      get_NED3(ras_dem, sf_stasbuf[i,], 
               sf_stasbuf[i_other_stas,],
               ver = F, 
               w_z = lambdas)
    for(i_lambda in 1:length(lambdas)){
      list_dists_stas_NED3[[i_lambda]][i_other_stas,i] <- 
        mat_dists_NED3[i_lambda,]
      list_dists_stas_NED3[[i_lambda]][i,i_other_stas] <- 
        mat_dists_NED3[i_lambda,]
    }
  }
  
  file_save_list_dists <- 
    sprintf(file_save_list_dists_gen, id_reg)
  saveRDS(list_dists_stas_NED3,
          file_save_list_dists)
  message()
}






