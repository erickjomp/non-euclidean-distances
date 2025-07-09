library(tidyverse)
library(terra)
library(sf)
library(metR)
devtools::load_all("~/projects/my_packages/my_themes/")
list_files_dists <- 
  list(
    c("02_regions_dem_conecs//raster/fig_dists/ras_dists_NED1_mountain.tif",
      "02_regions_dem_conecs/raster/fig_dists/ras_dists_NED2_mountain.tif",
      "02_regions_dem_conecs/raster/fig_dists/ras_dists_NED3_mountain.tif"),
    c("02_regions_dem_conecs/raster/fig_dists/ras_dists_NED1_valley.tif",
      "02_regions_dem_conecs/raster/fig_dists/ras_dists_NED2_valley.tif",
      "02_regions_dem_conecs/raster/fig_dists/ras_dists_NED3_valley.tif")
    )
path_ras_dem <- 
  "02_regions_dem_conecs/raster/DEM/ras_dem_reg8.tif"

files_figs <- 
  c("02_regions_dem_conecs//fig_distances_mountain.pdf",
    "02_regions_dem_conecs//fig_distances_valley.pdf")
files_figs_2 <- 
  c("article/figs/fig_distances_mountain.pdf",
    "article/figs/fig_distances_valley.pdf")
id_nodes <- c(23327,24424)

sf_regions <- st_read("00_data/shp/regions.shp")
sf_nodes <- st_read("02_regions_dem_conecs/shp/nodes/sf_nodes_reg8.shp")
id_reg <- 8


#### for rescaling ####
sf_nodes_raster <- sf_nodes %>% filter(source == "raster",is_finite==1) 
dists_raster_id_1 <- 
  st_distance(sf_nodes[id_nodes[1],], sf_nodes_raster) %>% as.numeric()
dists_raster_id_2 <- 
  st_distance(sf_nodes[id_nodes[2],], sf_nodes_raster) %>% as.numeric()

# just to verify
ids_finite_raster <- filter(sf_nodes, is_finite == 1, 
                            source == "raster")$id
ras_dists_euc_id_1 <- rast(path_ras_dem)    # just to verify
terra::values(ras_dists_euc_id_1)[ids_finite_raster] <- dists_raster_id_1
plot(ras_dists_euc_id_1) # just to verify

# calculations avgs
avg_id_1 <- mean(dists_raster_id_1)
avg_id_2 <- mean(dists_raster_id_2)
avgs <- c(avg_id_1, avg_id_2)/1000  # 1000m = 1km

#### ploting ####
for (i in 1:length(files_figs)){
  files_dists <- list_files_dists[[i]]
  id_node <- id_nodes[i]
  file_fig <- files_figs[i]
  file_fig2 <- files_figs_2[i]
  
  ras_dists_NED1 <- 
    terra::rast(files_dists[1])
  ras_dists_NED2 <- 
    terra::rast(files_dists[2])
  ras_dists_NED3 <- 
    terra::rast(files_dists[3])
  
  # rescaling
  fun_scale <- function(ras0, avg){
    for (j in 1:nlyr(ras0)){
      values <- ras0[[j]][]   # it could be i
      values <- values* avg/mean(values,na.rm = T)
      ras0[[j]][] <- values
    }
    return(ras0)
  }
  
  ras_dists_NED1 <- fun_scale(ras_dists_NED1,avgs[i])
  ras_dists_NED2 <- fun_scale(ras_dists_NED2,avgs[i])
  ras_dists_NED3 <- fun_scale(ras_dists_NED3,avgs[i])
  # edn rescaling
  
  
  gdat_NED1 <- 
    tibble(x = crds(ras_dists_NED1, na.rm = F)[,1],
           y = crds(ras_dists_NED1, na.rm = F)[,2],
           type = "NED1",
           as.data.frame(ras_dists_NED1[]))
  
  gdat_NED2 <-  
    tibble(x = crds(ras_dists_NED2, na.rm = F)[,1],
           y = crds(ras_dists_NED2, na.rm = F)[,2],
           type = "NED2",
           as.data.frame(ras_dists_NED2[]))
  
  gdat_NED3 <-  
    tibble(x = crds(ras_dists_NED3, na.rm = F)[,1],
           y = crds(ras_dists_NED3, na.rm = F)[,2],
           type = "NED3",
           as.data.frame(ras_dists_NED3[]))
  
  gdat <- rbind(gdat_NED1,gdat_NED2,gdat_NED3) %>% 
    gather(lambda, distance, -(1:3))
  
  # breaks = c(0,20,40,60,80,100,150,200,Inf)
  breaks = c(0,20,40,60,80,100,150,Inf)
  
  sf_gdat <- st_as_sf(gdat, coords = c("x","y"), crs = 4326)
  
  
  lambdas <- unique(sf_gdat$lambda)
  lambda_labs <- paste0("\u03bb = ", lambdas)
  
  sf_gdat <- sf_gdat %>%
    mutate(lambda_labs = factor(paste0("\u03bb = ", lambda),
                                levels = lambda_labs))
  
  base_size_figs <- 7
  
  lambdas_labs <- paste0("\u03bb = ", lambdas)
  
  
  # https://stackoverflow.com/questions/61737062/the-command-geom-contour-fill-fails-to-render-in-some-cases-in-ggplot2
  
  (
    plot_distances_mountain <- 
      ggplot(data = sf_gdat) + 
      # geom_tile(aes(fill = z)) +
      # geom_contour_filled(aes(x = gdat$x, y = gdat$y,z = distance)
      #                     ,breaks = breaks,color = NA) + # bins = 18,
      geom_contour_fill(aes(x = gdat$x, y = gdat$y, z = distance,
                            fill = stat(level)),
                        breaks = breaks,color = NA) + # bins = 18,
      # guides(fill = guide_colorsteps(#show.limits = T,
      #   barheight = 15, draw.llim = T)) + # works with default legend.position
      guides(fill = guide_colorsteps(#show.limits = T,
        barheight = 1, draw.llim = T, barwidth = 20
        # title.position = "top"
        )) +
      coord_equal() +
      scale_y_continuous(limits = c(-13.6, -10.6),
                         expand  = c(0,0)) +  
      scale_x_continuous(breaks = c(-77,-76, -75, -74)) +
      # ylim(c(-13.5, -10.5)) + 
      # facet_wrap(~lambda) + 
      facet_grid(lambda_labs ~ type) +
      scale_fill_viridis_d(option = "inferno") +
      my_ggtheme(base_size = base_size_figs) +
      labs(x = NULL, y = NULL, fill = "Distance (km)")  +
      geom_sf(data = sf_regions[id_reg,],fill = NA, col = "black")  +
      geom_sf(data = sf_nodes[id_node,], col = "white", 
              shape = 3, size = 1) + 
      theme(legend.position = "bottom",
            strip.text.x = element_text(size = base_size_figs*1.2),
            strip.text.y = element_text(size = base_size_figs*1.2)) 
  )
  
  # ggsave("02_regions_dem_conecs/fig_distances_mountain.png",
  #        width = 190,units = "mm")
  
  # cairo_pdf(file_fig,
  #           width = 190 * 0.0393701, 
  #           height = 145 * 0.0393101 ) # works for right legend
  cairo_pdf(file_fig,
            # width = 190 * 0.0393701, 
            # height = 180 * 0.0393101 
            width = 140 * 0.0393701, 
            height = 132 * 0.0393101 
            )
  print(plot_distances_mountain)
  dev.off()
  
  cairo_pdf(file_fig2,
            width = 140 * 0.0393701, height = 132 * 0.0393101 )
  print(plot_distances_mountain)
  dev.off()
  
}

