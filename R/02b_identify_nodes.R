library(tidyverse)
library(sf)
library(terra)
library(foreach)
library(parallel)
library(doParallel)
source("R/functions/funciones_perfiles.R")


file_log <- "02_regions_dem_conecs/conecs_NED1.log"

path_save_conecs <- "02_regions_dem_conecs/df_conecs_NED1/"
filename_conec_gen <- "df_conec_reg%s.rds"

file_sf_nodes_gen <- 
  "02_regions_dem_conecs//shp/nodes/sf_nodes_reg%s.shp"

file_ras_dem_reg_gen <- 
  "02_regions_dem_conecs/raster/DEM/ras_dem_reg%s.tif"
sf_stas_regbuf_all <- 
  st_read("01_data_byregion//shp_stations//sf_stas_regbuf.shp")
sf_stas_reg_all <- 
  st_read("01_data_byregion//shp_stations//sf_stas_reg.shp")

id_regions <- sort(unique(sf_stas_reg_all$region))


list_id_nodes_stas <- list()

## loop ##

for (id_reg in id_regions){
  message(sprintf("REGION %s",id_reg))

  ras_dem_reg <-
    rast(sprintf(file_ras_dem_reg_gen,id_reg))
  
  sf_stas_regbuf <- sf_stas_regbuf_all %>% filter(region == id_reg)
  sf_stas_reg <- sf_stas_reg_all %>% filter(region == id_reg)
  
  mat_coords_ras <- terra::crds(ras_dem_reg,na.rm = F) %>% 
    as.data.frame() %>% 
    mutate(source = "raster", 
           namcod = paste0("raster_",1:n()))
  mat_coords_stas <- st_coordinates(sf_stas_regbuf) %>% 
    as.data.frame() %>% 
    mutate(source = "stas", 
           namcod = sf_stas_regbuf$namcod) %>% 
    rename(x = "X", y = "Y") 
  mat_coords_nodes <- rbind(mat_coords_ras, mat_coords_stas) 
  
  sf_nodes <- st_as_sf(mat_coords_nodes,
                       coords = c("x","y"),crs = 4326)
  sf_nodes$id <- 1:nrow(sf_nodes)
  
  is_finite_nodes <- c(
    is.finite(ras_dem_reg[][,1]),
    rep(TRUE, nrow(sf_stas_regbuf)))
  
  sf_nodes$is_finite <- is_finite_nodes
  
  file_sf_nodes <- sprintf(file_sf_nodes_gen, id_reg)
  st_write(sf_nodes, file_sf_nodes,delete_dsn = T)
  
  list_id_nodes_stas[[id_reg]] <- sf_nodes %>% 
    filter(source == "stas") %>% st_drop_geometry() %>% 
    mutate(region = id_reg) %>% select(region, id_node = id, namcod)
}


#### updating shapes with id_nodes ####

df_id_nodes_stas <- do.call(rbind, list_id_nodes_stas)

sf_stas_regbuf_all <- sf_stas_regbuf_all %>% 
  left_join(df_id_nodes_stas, by = c("region", "namcod"))

sf_stas_reg_all <- sf_stas_reg_all %>% 
  left_join(df_id_nodes_stas, by = c("region", "namcod"))



st_write(sf_stas_reg_all,
         "02_regions_dem_conecs//shp/sf_stas_reg.shp",
         delete_dsn = T)
st_write(sf_stas_regbuf_all,
         "02_regions_dem_conecs/shp/sf_stas_regbuf.shp",
         delete_dsn = T)
