library(tidyverse)
library(sf)
library(terra)
library(foreach)
source("R/functions/funciones_perfiles_NED3.R")
source("R/functions/funciones_perfiles.R") 
Sys.setenv(TZ='GMT')

# input ------
file_sf_nodes_gen <-
  "02_regions_dem_conecs//shp/nodes/sf_nodes_reg%s.shp"

sf_stasbuf_all <- st_read("02_regions_dem_conecs//shp/sf_stas_regbuf.shp")

file_ras_dem_gen <- 
  "02_regions_dem_conecs//raster/DEM/ras_dem_reg%s.tif"

lambdas <-
  read.csv("02_regions_dem_conecs//input/lambdas.csv")[[1]]

# files output -----
file_save_list_dists_gen <- 
  "03_raster_interpolation//mat_NED3_raster/list_mat_NED3_lambdas_reg%s.rds"

#### parallel configuration ####
parallel::detectCores()
n_cores <- 3
my_cluster <- parallel::makeCluster(
  n_cores, 
  type = "PSOCK"
)
print(my_cluster)
doParallel::registerDoParallel(cl = my_cluster)
foreach::getDoParRegistered()

#### verifying which rds files exist ####
id_regions <- 1:12

id_regions_todo <- 
  id_regions[
    !file.exists(sprintf(file_save_list_dists_gen, id_regions))
  ]

#log file
file_log <- "03_raster_interpolation//log_mat_NED3_raster.log"

# clearing text file
con_file_log <- file(file_log, open = "w",) #op
close(con_file_log)


# loop for regions ----
# lapply(1:12, function(id_reg) {

foreach(
  id_reg = id_regions_todo, #sf_nodes$id,
  .packages = c("dplyr","sf","terra")) %dopar% {
    
  # con_file_log <- file(file_log, open = "w") #op
  # sink(con_file_log, type = "message") #op
  # message(paste0("REGION ", id_reg))
  # sink(type = "message") #op
  
  file_sf_nodes <-
    sprintf(file_sf_nodes_gen, id_reg)
  sf_nodes <-
    st_read(file_sf_nodes,quiet = T)
  
  ras_dem <- rast(
    sprintf(file_ras_dem_gen,id_reg)
    )
  
  sf_stasbuf <-
    sf_stasbuf_all %>% filter(region == id_reg)
  n_stas <- nrow(sf_stasbuf)
  
  id_nodes_stas <- sf_nodes %>% filter(source == "stas") %>% .$id
  id_nodes_raster <- sf_nodes %>% filter(source == "raster") %>% .$id
  
  sf_nodes_raster <- sf_nodes %>% filter(source == "raster") 
  # added
  mat_dists <- matrix(NA,nrow = n_stas, ncol = nrow(sf_nodes_raster))
  
  list_dists_lambdas <- list()
  for (i in 1:length(lambdas)){
    list_dists_lambdas[[i]] <- mat_dists
  }
  names(list_dists_lambdas) <- lambdas
  
  for (i in 1:(n_stas-1)){
    con_file_log <- file(file_log, open = "a",) #op
    sink(con_file_log, type = "message") #op
    message("REGION ",id_reg,": Working in ", i ,
            "/", n_stas,"\r",appendLF = F)
    sink(type = "message") #op
    close(con_file_log) #op
    
    mat_dists_NED3 <- 
      get_NED3(ras_dem, 
               sf_stasbuf[i,], 
               # sf_nodes_raster,
               ver = F, 
               w_z = lambdas)
    for(i_lambda in 1:length(lambdas)){
      list_dists_lambdas[[i_lambda]][i,] <- 
        mat_dists_NED3[i_lambda,]
    }
  }

  message("\n",appendLF = F)
  
  file_save_list <- 
    sprintf(file_save_list_dists_gen, id_reg)
  saveRDS(list_dists_lambdas,file = file_save_list)
}
