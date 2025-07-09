library(tidyverse)
library(sf)
library(terra)
library(Matrix)
library(igraph)
library(tmap)
source("R/functions/funciones_perfiles_NED3.R")
source("R/functions/funciones_perfiles.R")
id_reg <- 8


df_conec_NED1 <- 
  readRDS("02_regions_dem_conecs//df_conecs_lambdas_NED1/df_conec_reg8.rds")
df_conec_NED2 <- 
  readRDS("02_regions_dem_conecs/df_conecs_lambdas_NED2/df_conec_reg8.rds")

ras_dem_reg <- rast("02_regions_dem_conecs/raster/DEM/ras_dem_reg8.tif")

sf_regions <- st_read("00_data/shp/regions.shp")
sf_nodes <- 
  st_read("02_regions_dem_conecs/shp/nodes/sf_nodes_reg8.shp")

# ext_crop <- ext(c(0,0,0,0))


id_nodes <- c(23327,24424) #23643  # first mountain, second valley
plot(ras_dem_reg)
plot(sf_nodes[id_nodes,]$geometry,add=T)
#23327
#26000

lambdas <- c(20,60,100)
lambdas_fig <- as.character(lambdas)



#### NED 1 #### 
list_ras_dists_NED1 <- list()

for (id_node in id_nodes){
  list_ras_dists_NED1[[as.character(id_node)]] <- 
    lapply(lambdas_fig, function(lambda){
      mat_cons <- 
        sparseMatrix(i = df_conec_NED1$id_start_node,
                     j = df_conec_NED1$id_end_node,
                     x = df_conec_NED1[[lambda]],
                     dims = rep(nrow(sf_nodes), 2),symmetric = T)
      graph_cons <- 
        graph.adjacency(mat_cons, weighted = T, mode = "undirected")
      
      ids_finite_raster <- filter(sf_nodes, is_finite == 1, 
                                  source == "raster")$id
      
      dists_NED1 <- shortest.paths(graph_cons, 
                                   v = id_node,
                                   to = ids_finite_raster) %>% 
        as.numeric()
      
      ras_dists_NED1 <- ras_dem_reg
      values(ras_dists_NED1)[ids_finite_raster] <- dists_NED1/1000
      names(ras_dists_NED1) <- lambda
      ras_dists_NED1
    }) %>% do.call(c,.)
}



# breaks = c(0,10,20,30,40,60,80,100,150,200,300,400,Inf)
# breaks = c(0,10,20,30,40,60,80,100,150,200,Inf)
breaks = c(0,20,40,60,80,100,150,200,Inf)

ext_map <- ext(c(-76,-74.5,-13,-11.5))
# plot(ras_dists_NED1)
tm_shape(list_ras_dists_NED1[[2]],ylim=c(-13.5, -10.5)) + 
  tm_raster(n=20,title = 'Distance (km)',breaks = breaks) +
  # tm_grid(lines = F) +
  # tm_grid() +
  tm_layout(aes.palette = list(seq = "viridis")) + # "-Blues"
  tm_shape(sf_regions[id_reg,]) + tm_borders(col = "black")

file_fig1 <- 
  "02_regions_dem_conecs/raster/fig_dists/ras_dists_NED1_mountain.tif"
writeRaster(list_ras_dists_NED1[[1]],
            filename = file_fig1,
            overwrite = T)
file_fig2 <- 
  "02_regions_dem_conecs/raster/fig_dists/ras_dists_NED1_valley.tif"
writeRaster(list_ras_dists_NED1[[2]],
            filename = file_fig2,
            overwrite = T)

#### NED 2 ####

list_ras_dists_NED2 <- list()
  
for (id_node in id_nodes){
  list_ras_dists_NED2[[as.character(id_node)]] <-
    lapply(lambdas_fig, function(lambda){
      mat_cons <- 
        sparseMatrix(i = df_conec_NED2$id_start_node,
                     j = df_conec_NED2$id_end_node,
                     x = df_conec_NED2[[lambda]],
                     dims = rep(nrow(sf_nodes), 2),symmetric = T)
      graph_cons <- 
        graph.adjacency(mat_cons, weighted = T, mode = "undirected")
      
      ids_finite_raster <- filter(sf_nodes, is_finite == 1, 
                                  source == "raster")$id
      
      dists_NED1 <- shortest.paths(graph_cons, 
                                   v = id_node,
                                   to = ids_finite_raster) %>% 
        as.numeric()
      
      ras_dists_NED1 <- ras_dem_reg
      values(ras_dists_NED1)[ids_finite_raster] <- dists_NED1/1000
      names(ras_dists_NED1) <- lambda
      ras_dists_NED1
    }) %>% do.call(c,.)
}


tm_shape(list_ras_dists_NED2[[2]],ylim=c(-13.5, -10.5)) + 
  tm_raster(n=20,title = 'Distance (km)',breaks = breaks) +
  # tm_grid(lines = F) +
  # tm_grid() +
  tm_layout(aes.palette = list(seq = "viridis")) + # "-Blues"
  tm_shape(sf_regions[id_reg,]) + tm_borders(col = "black")

file_fig1 <- 
  "02_regions_dem_conecs/raster/fig_dists/ras_dists_NED2_mountain.tif"
writeRaster(list_ras_dists_NED2[[1]],
            filename = file_fig1,
            overwrite = T)
file_fig2 <- 
  "02_regions_dem_conecs/raster/fig_dists/ras_dists_NED2_valley.tif"
writeRaster(list_ras_dists_NED2[[2]],
            filename = file_fig2,
            overwrite = T)



#### NED 3 ####

list_mat_ras_dists_NED3 <-
  lapply(id_nodes, function(id_node){
    mat_ras_dists_NED3 <- 
      get_NED3(ras_dem_reg, sf_nodes[id_node,], 
               filter(sf_nodes, source == "raster", is_finite == 1),
               ver = T, 
               w_z = lambdas)
  })
names(list_mat_ras_dists_NED3) <- id_nodes

list_ras_dists_NED3 <- list()
for (id_node in id_nodes){
  list_ras_dists_NED3[[as.character(id_node)]] <-
    lapply(1:length(lambdas_fig), function(i_lambda){
      
      lambda <- lambdas_fig[i_lambda]  
      ras_dists_NED3 <- ras_dem_reg
      mat_ras_dists_NED3 <- 
        list_mat_ras_dists_NED3[[as.character(id_node)]]
      
      ids_finite_raster <- filter(sf_nodes, is_finite == 1, 
                                  source == "raster")$id
      values(ras_dists_NED3)[ids_finite_raster] <- 
        mat_ras_dists_NED3[i_lambda, ] %>% as.numeric()
      
      names(ras_dists_NED3) <- lambda
      ras_dists_NED3
    }) %>% do.call(c,.)
}

tm_shape(list_ras_dists_NED3[[1]],ylim=c(-13.5, -10.5)) + 
  tm_raster(n=20,title = 'Distance (km)',breaks = breaks) +
  # tm_grid(lines = F) +
  # tm_grid() +
  tm_layout(aes.palette = list(seq = "viridis")) + # "-Blues"
  tm_shape(sf_regions[id_reg,]) + tm_borders(col = "black")

file_fig1 <- 
  "02_regions_dem_conecs/raster/fig_dists/ras_dists_NED3_mountain.tif"
writeRaster(list_ras_dists_NED3[[1]],
            filename = file_fig1,
            overwrite = T)
file_fig2 <- 
  "02_regions_dem_conecs/raster/fig_dists/ras_dists_NED3_valley.tif"
writeRaster(list_ras_dists_NED3[[2]],
            filename = file_fig2,
            overwrite = T)









#### EXTRA FOR ZOOM ####
# br_dis <- raster::brick(ras_dists_NED1)
# raster::filledContour(br_dis)
ext_crop <- ext(c(-76,-74.5,-13,-11.5))

ras_dists_NED1_cr <- 
  ras_dists_NED1 %>% crop(ext_crop,snap = "out")

tm_shape(ras_dists_NED1_cr) + 
  tm_raster(n=20,title = 'Distance (km)',breaks = breaks) +
  # tm_facets(by = names) + 
  # tm_grid(lines = F) +
  tm_grid() +
  tm_layout(aes.palette = list(seq = "viridis"))  +  # "-Blues"
  tm_shape(sf_regions[id_reg,]) + tm_borders()
  
