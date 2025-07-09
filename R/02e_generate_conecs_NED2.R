library(tidyverse)
library(sf)
library(terra)
source("R/functions/funciones_perfiles.R")

path_save_conecs <- "02_regions_dem_conecs//df_conecs_NED2/"
filename_conec_gen <- "df_conec_reg%s.rds"

file_ras_dem_reg_gen <- 
  "02_regions_dem_conecs//raster/DEM/ras_dem_reg%s.tif"
file_sf_nodes_gen <- 
  "02_regions_dem_conecs//shp/nodes/sf_nodes_reg%s.shp"
sf_stas_regbuf_all <- 
  st_read("01_data_byregion//shp_stations//sf_stas_regbuf.shp")
sf_stas_reg_all <- 
  st_read("01_data_byregion//shp_stations//sf_stas_reg.shp")

id_regions <- sort(unique(sf_stas_reg_all$region))


#### verifying whic df_conecs....rds exist ####
files_found <- list.files(path_save_conecs)
filenames_conecs <- sprintf(filename_conec_gen, id_regions)
(
  id_regions_todo <- 
    id_regions[
      which(!(filenames_conecs %in% files_found))
    ]
)


# con_file_log <- file(file_log, open = "w")

#### loop by regions ####
for (id_reg in id_regions_todo){
  message(sprintf("#### WORKING IN REGION %s ####",id_reg))
  #TODO create file to write messages from paralleled operations
  # tic()
  ras_dem_reg <-
    rast(sprintf(file_ras_dem_reg_gen,id_reg)) # moved to foreach
  # toc()
  sf_stas_regbuf <- sf_stas_regbuf_all %>% filter(region == id_reg)
  sf_stas_reg <- sf_stas_reg_all %>% filter(region == id_reg)
  # plot(ras_dem_reg)
  # sf_stas_regbuf$geometry %>% plot(add=T)
  
  file_sf_nodes <- 
    sprintf(file_sf_nodes_gen, id_reg)
  sf_nodes <- st_read(file_sf_nodes)
  
  
  is_finite_nodes <- c(
    is.finite(ras_dem_reg[][,1]),
    rep(TRUE, nrow(sf_stas_regbuf)))
  
  n_nodes <- length(sf_nodes$id)
  
  ids_raster_finite <-  sf_nodes %>% 
    filter(is_finite == 1, source == "raster") %>% .$id
  mat_cons_raster <- terra::adjacent(ras_dem_reg,cells = ids_raster_finite,
                  directions = "8", pairs = T)
  df_conecs_raster <- mat_cons_raster %>% as.data.frame() %>% 
    rename_all(~c("id_start_node", "id_end_node")) %>% 
    filter(id_start_node < id_end_node) %>% 
    filter((id_start_node %in% ids_raster_finite) & 
             (id_end_node %in% ids_raster_finite))
    
  df_cells_stas <- 
    terra::cells(ras_dem_reg, vect(sf_stas_regbuf)) %>% 
    as.data.frame()
  
  mat_cons_stas <- 
    terra::adjacent(ras_dem_reg, df_cells_stas$cell,
                    directions = "8", include = T)
  
  df_conecs_stas <- 
  lapply(1:nrow(mat_cons_stas), function(i){
    data.frame(id_start_node = mat_cons_stas[i,], 
               id_end_node = ncell(ras_dem_reg) + i)
  }) %>% do.call(rbind,.)
  
  df_conecs_stas <- 
    df_conecs_stas %>% 
    filter((id_start_node %in% ids_raster_finite))
  
  df_conecs <- rbind(df_conecs_raster, df_conecs_stas)
  
  ### dif_h and dif_v ###
  df_conecs$dif_h <- 
  st_distance(
    sf_nodes[df_conecs$id_start_node,],
    sf_nodes[df_conecs$id_end_node,],by_element = T
  ) %>% as.numeric()
  
  elevs_nodes <- c(ras_dem_reg[], sf_stas_regbuf$elev)
  
  df_conecs$dif_v <- 
    elevs_nodes[df_conecs$id_end_node] -
    elevs_nodes[df_conecs$id_start_node]
  
  ### 
  file_reg <- sprintf(filename_conec_gen, id_reg)
  file_reg <- file.path(path_save_conecs,file_reg)
  saveRDS(df_conecs,file_reg)
  
}


