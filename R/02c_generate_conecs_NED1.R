library(tidyverse)
library(sf)
library(terra)
library(foreach)
library(parallel)
library(doParallel)
source("R/functions/funciones_perfiles.R")


dist_rad <- 120000

file_log <- "02_regions_dem_conecs/conecs_NED1.log"

path_save_conecs <- "02_regions_dem_conecs/df_conecs_NED1/"
filename_conec_gen <- "df_conec_reg%s.rds"

file_ras_dem_reg_gen <- 
  "02_regions_dem_conecs/raster/DEM/ras_dem_reg%s.tif"
sf_stas_regbuf_all <- 
  st_read("01_data_byregion//shp_stations//sf_stas_regbuf.shp")
sf_stas_reg_all <- 
  st_read("01_data_byregion//shp_stations//sf_stas_reg.shp")

id_regions <- sort(unique(sf_stas_reg_all$region))

#### parallel configuration ####
parallel::detectCores()
n_cores <- 4
my_cluster <- parallel::makeCluster(
  n_cores, 
  type = "PSOCK"
)
print(my_cluster)
doParallel::registerDoParallel(cl = my_cluster)
foreach::getDoParRegistered()

#### verifying whic df_conecs....rds exist ####
files_found <- list.files(path_save_conecs)
filenames_conecs <- sprintf(filename_conec_gen, id_regions)
id_regions_todo <- 
id_regions[
  which(!(filenames_conecs %in% files_found))
]

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
  
  mat_coords_ras <- terra::crds(ras_dem_reg,na.rm = F) 
  mat_coords_stas <- st_coordinates(sf_stas_regbuf)
  mat_coords_nodes <- rbind(mat_coords_ras, mat_coords_stas) 
  
  sf_nodes <- st_as_sf(mat_coords_nodes %>% as.data.frame(),
                       coords = c("x","y"),crs = 4326)
  sf_nodes$id <- 1:nrow(sf_nodes)
  
  is_finite_nodes <- c(
    is.finite(ras_dem_reg[][,1]),
    rep(TRUE, nrow(sf_stas_regbuf)))
  
  rm(ras_dem_reg) # since it can not be send to parallel nodes
  # function wrap was other option, howver its slower
  n_nodes <- length(sf_nodes$id)
  
  df_conecs <- 
    # lapply(1:52, function(id_center_node){
    foreach(id_center_node = sf_nodes$id, #sf_nodes$id,
            .combine = "rbind",
            .packages = c("dplyr","sf","terra")
            ) %dopar% {
      message(paste0(id_center_node,"\r"), appendLF = F)
      # this message is shown in the console after executing
      # sink( type = "message") in the first iteration
      # it was also ok to exexute this line before this message
      # so it would be printed in console and the log file would not be
      # necessary
      if(is_finite_nodes[id_center_node]){
        # log start #
        con_file_log <- file(file_log, open = "w")
        sink(con_file_log, type = "message")
        message(paste0(c("Working in regions",id_regions_todo),
                       collapse = " "))
        message(paste("REGION",id_reg,": node", id_center_node, 
                      "\t(",round(100 * id_center_node/ n_nodes, 2) # change
                      ,"% )"))
        sink(type = "message")
        close(con_file_log)  # add to code in linux
        # log end #
        
        ras_dem_reg <-
          rast(sprintf(file_ras_dem_reg_gen,id_reg))  # only for parallel
        
        dist_loc <- st_distance(sf_nodes[id_center_node,],sf_nodes,) %>% 
          as.numeric()
        id_rad <- which((dist_loc < dist_rad) & is_finite_nodes)
        
        # the next line is to ensure dont repeat the same pairs 
        # in different loops of lapply
        id_rad <- id_rad[id_rad > id_center_node]
        # id_rad <- id_rad[id_rad != id_center_node]
        
        # # tic("starting")
        # sf_lines <- 
        #   lapply(id_rad, function(id_end_node){
        #     sf_line <- mat_coords_nodes[c(id_center_node,
        #                             id_end_node),] %>% st_linestring()
        #   }) %>%  do.call(st_sfc, .) %>% 
        #   st_as_sf(crs = 4326)
        # vec_lines <- vect(sf_lines)
        # # toc()
        
        df_conecs <- NULL
        if(length(id_rad) > 0){
          vec_lines <- 
            lapply(1:length(id_rad), function(i_id_end_node){
              id_end_node <- id_rad[i_id_end_node]
              mat_ex1 <- cbind(geom = rep(i_id_end_node, 2), part = c(1,1))
              mat_coords_nodes[c(id_center_node,id_end_node),] %>% 
                cbind(mat_ex1,.)
            }) %>% do.call(rbind,.) %>% vect(crs = "epsg:4326",type = "lines")
          
          df_profiles <- 
            get_profiles_lines_same_origin(ras_dem_reg, vec_lines, 
                                           sf_nodes[id_center_node,])
          df_info_lines_sel <- which_direct_profile(df_profiles)
          id_end_nodes_sel <- id_rad[df_info_lines_sel$ID]
          
          df_conecs <- NULL
          if(length(id_end_nodes_sel) > 0){
            df_conecs <- data.frame(id_start_node = id_center_node, 
                                    id_end_node = id_end_nodes_sel,
                                    dif_h = df_info_lines_sel$dif_h,
                                    dif_v = df_info_lines_sel$dif_v)
          }
        }
      } else {
        df_conecs <- NULL
      }
      df_conecs
      
    }
    # }) %>% do.call(rbind,.)
  
  file_reg <- sprintf(filename_conec_gen, id_reg)
  file_reg <- file.path(path_save_conecs,file_reg)
  saveRDS(df_conecs,file_reg)
  
}


