library(tidyverse)
library(igraph)
library(sf)
library(Matrix)

# input ------
file_conecs_lambdas_gen <-
  "02_regions_dem_conecs//df_conecs_lambdas_NED2/df_conec_reg%s.rds"
file_sf_nodes_gen <-
  "02_regions_dem_conecs/shp/nodes/sf_nodes_reg%s.shp"

sf_stasbuf_all <- st_read("02_regions_dem_conecs//shp/sf_stas_regbuf.shp")

lambdas <-
  read.csv("02_regions_dem_conecs//input/lambdas.csv")[[1]]

# files output -----
file_save_list_dists_gen <- 
  "03_raster_interpolation//mat_NED2_raster/list_mat_NED2_lambdas_reg%s.rds"


# loop for regions ----
lapply(1:12, function(id_reg) {
  message(paste0("REGION ", id_reg))
  file_conecs_lambdas <-
    sprintf(file_conecs_lambdas_gen, id_reg)
  df_conecs_lambdas <-
    readRDS(file_conecs_lambdas)
  file_sf_nodes <-
    sprintf(file_sf_nodes_gen, id_reg)
  sf_nodes <-
    st_read(file_sf_nodes,quiet = T)
  
  sf_stasbuf <-
    sf_stasbuf_all %>% filter(region == id_reg)
  
  id_nodes_stas <- sf_nodes %>% filter(source == "stas") %>% .$id
  id_nodes_raster <- sf_nodes %>% filter(source == "raster") %>% .$id
  
  list_dists_lambdas <- 
    lapply(as.character(lambdas), function(lambda){
      message(paste0("lambda = ", lambda, "\r"),appendLF = F)
      mat_cons <-
        sparseMatrix(
          i = df_conecs_lambdas$id_start_node,
          j = df_conecs_lambdas$id_end_node,
          x = df_conecs_lambdas[[lambda]],
          dims = rep(nrow(sf_nodes), 2),
          symmetric = T
        )
      graph_cons <-
        graph.adjacency(mat_cons, weighted = T, mode = "undirected")
      
      
      dists_NED <- shortest.paths(graph_cons,
                                  v = id_nodes_stas,
                                  to = id_nodes_raster)
      return(dists_NED)
    })
  message("\n",appendLF = F)
  names(list_dists_lambdas) <- lambdas
  
  file_save_list <- 
    sprintf(file_save_list_dists_gen, id_reg)
  saveRDS(list_dists_lambdas, file_save_list)
})
