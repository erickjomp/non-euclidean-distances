library(tidyverse)
library(sf)
library(terra)
library(Matrix)
library(igraph)
library(tmap)
source("R/functions/funciones_perfiles_NED3.R")
source("R/functions/funciones_perfiles.R") 
devtools::load_all("~/projects/my_packages/my_themes/")
id_reg <- 8



ras_dem_reg <- rast("02_regions_dem_conecs//raster/DEM/ras_dem_reg8.tif")
sf_nodes <- 
  st_read("02_regions_dem_conecs/shp/nodes/sf_nodes_reg8.shp")

df_conec_NED1 <- 
  readRDS("02_regions_dem_conecs/df_conecs_lambdas_NED1/df_conec_reg8.rds")
df_conec_NED2 <- 
  readRDS("02_regions_dem_conecs/df_conecs_lambdas_NED2/df_conec_reg8.rds")

file_fig <- "02_regions_dem_conecs//fig_profiles_NEDs.pdf"
file_fig_png <- "02_regions_dem_conecs//fig_profiles_NEDs.png"
file_fig_article <- "article/figs/fig_profiles_NEDs.pdf"


#### exploration ####
plot(ras_dem_reg)
ex_box <- ext(-75.8,-74.8, -12.6, -11.6)
plot(ex_box,add=T)
plot(sf_nodes[25373,]$geometry,add=T)
plot(sf_nodes[19721,]$geometry,add=T)
plot(crop(ras_dem_reg, ex_box))
ras_dem_ex <- crop(ras_dem_reg, ex_box)

# 23500
#18308 node esquina superior izquierda
#  + 156
extract(ras_dem_reg, ex_box, cells = T)

##### shortest paths ####
# id_start_node <- 19721 + 156
id_start_node <- 19721
# id_end_node <- 25367
id_end_node <- 25373

#### NED 1 ####
lambda <- "60"
mat_cons <- 
  sparseMatrix(i = df_conec_NED1$id_start_node,
               j = df_conec_NED1$id_end_node,
               x = df_conec_NED1[[lambda]],
               dims = rep(nrow(sf_nodes), 2),symmetric = T)
graph_cons <- 
  graph.adjacency(mat_cons, weighted = T, mode = "undirected")

path_NED1 <- shortest_paths(graph_cons, 
                             from = id_start_node,
                             to = id_end_node)[[1]][[1]]
NED1 <- shortest.paths(graph_cons, 
               v = id_start_node,
               to = id_end_node) %>% as.numeric

sf_path_NED1 <- 
  sf_nodes[path_NED1,] %>% 
  st_coordinates() %>% 
  st_linestring() %>% 
  st_sfc() %>% 
  st_as_sf(crs = 4326) %>% 
  mutate(group = "NED1")

plot(ras_dem_ex)
plot(sf_path_NED1$x, add=T)

df_profile_NED1 <- 
  get_profiles_lines_wsegments_same_origin(
    ras_dem_reg,
    transect_line = sf_path_NED1,
    origin = sf_nodes[id_start_node,],
    cells = T) %>% ungroup() %>% 
  mutate(dist_from_origin = dist_from_origin / 1000) 

df_v_path_NED1 <- 
  df_profile_NED1 %>% filter(cell  %in% path_NED1) %>% 
  select(dist = dist_from_origin, elev = value)

df_profile_NED1 <- df_profile_NED1  %>%
  select(dist = dist_from_origin, elev = value)

plot(df_profile_NED1, type = "b")
lines(df_v_path_NED1, col = "red")


#### NED 2 ####
lambda <- "60"
mat_cons <- 
  sparseMatrix(i = df_conec_NED2$id_start_node,
               j = df_conec_NED2$id_end_node,
               x = df_conec_NED2[[lambda]],
               dims = rep(nrow(sf_nodes), 2),symmetric = T)
graph_cons <- 
  graph.adjacency(mat_cons, weighted = T, mode = "undirected")

path_NED2 <- shortest_paths(graph_cons, 
                            from = id_start_node,
                            to = id_end_node)[[1]][[1]]
NED2 <- shortest.paths(graph_cons, 
                       v = id_start_node,
                       to = id_end_node) %>% as.numeric

sf_path_NED2 <- 
  sf_nodes[path_NED2,] %>% 
  st_coordinates() %>% 
  st_linestring() %>% 
  st_sfc() %>% 
  st_as_sf(crs = 4326) %>% 
  mutate(group = "NED2")

plot(ras_dem_ex)
plot(sf_path_NED1$x, add=T)
plot(sf_path_NED2$x, add=T, col = "blue")

df_profile_NED2 <- 
get_profiles_lines_wsegments_same_origin(
  ras_dem_reg,
  transect_line = sf_path_NED2,
  origin = sf_nodes[id_start_node,],
  cells = T) %>% ungroup() %>% 
  mutate(dist_from_origin = dist_from_origin / 1000) 

df_profile_NED2 <- 
  df_profile_NED2 %>% filter(cell  %in% path_NED2)
df_profile_NED2 <- df_profile_NED2 %>% 
  select(dist = dist_from_origin, elev = value)

df_v_path_NED2 <- df_profile_NED2

plot(df_profile_NED2, type = "b")
lines(df_v_path_NED2, col = "red")

#### NED 3 ####
sf_path_NED3 <-
  sf_nodes[c(id_start_node, id_end_node),] %>% 
  st_coordinates() %>% 
  st_linestring() %>% 
  st_sfc() %>% 
  st_as_sf(crs = 4326) %>% 
  mutate(group = "NED3")

plot(ras_dem_ex)
plot(sf_path_NED1$x, add=T, col = "black",size= 2)
plot(sf_path_NED2$x, add=T, col = "blue", size = 2)
plot(sf_path_NED3, add = T, col = "red", size = 2)
df_profile_NED3 <- 
get_profiles_lines_same_origin(ras_dem_reg,
                               transects_lines = sf_path_NED3,
                               origin = sf_nodes[id_start_node,],
                               cells = T) %>% ungroup() %>% 
  mutate(dist_from_origin = dist_from_origin / 1000) %>% 
  select(dist = dist_from_origin, elev = value)

df_v_path_NED3 <- trazo(df_profile_NED3) 

plot(df_profile_NED3 , type = "b")
lines(df_v_path_NED3, col = "red")


base_size_figs <- 7

sf_paths_NEDs <- rbind(sf_path_NED1, 
                       sf_path_NED2,
                       sf_path_NED3) 

df_ras_dem_ex <- 
crds(ras_dem_ex) %>% as.data.frame() %>% 
  data.frame(elev = values(ras_dem_ex,na.rm = T,mat = F))



colors_paths <- c("blue","red","green4")



#YlOrBr
(
  ggmap_1 <-
    ggplot() +
    geom_raster(data = df_ras_dem_ex, aes(x, y, fill = elev)) +
    geom_sf(data = sf_paths_NEDs, aes(col = group, linetype = group),
            show.legend = F,lwd = base_size_figs/10) +
    scale_linetype_manual(values = c("solid", "longdash", "twodash")) +
    scale_color_manual(values = colors_paths, ) +
    scale_fill_gradientn(
      colours =
        c("#FBCEB1", "#F4BB44", "brown"),
      guide = "colourbar"
    ) +
    #FFAA33
    # low = "yellow", high = "brown",) +
    # scale_fill_binned(
    #   breaks = seq(1000, 6000, 250),
    #   guide = guide_coloursteps(even.steps = FALSE),
    # ) +
    coord_sf(expand = F, ) +
    labs(x = NULL, y = NULL, fill = "Elevation (masl)", color = "Paths") +
    # theme(legend.position="bottom") +
    my_ggtheme(base_size = base_size_figs) +
    theme(
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank(),
      legend.position = "bottom",legend.box = "vertical"
    ) +
    guides(
      color = guide_legend(nrow = 2),
      fill = guide_colourbar(
        ticks.colour = "black",
        reverse = F,
        frame.colour = "black",
        barwidth = 6,barheight = 0.8,
        title.position = "top",
        title.hjust = 0.5,
        # fill = guide_legend(order = 1),
        # colour = guide_legend(order = 2)
      )
    )
)
  




library(tmap)
map_1 <- 
  tm_shape(ras_dem_ex) + 
  tm_raster(title = "Elevation\n(masl)",
            legend.format = list(#text.or.more ='or more',
              fun=function(x) formatC(x, digits=0, format="d"),
              text.to.columns = TRUE,
              text.separator = ' -')) +
tm_shape(sf_paths_NEDs,is.master = T) + 
  tm_lines(col = "group", lty = c("solid", "dashed", "dotdash"),
           palette = colors_paths,
           lwd = 2.5,title.col = "Paths",
           legend.col.show = F) +
  tm_layout(legend.bg.color = "white",
            legend.text.size = base_size_figs/18,
            legend.title.size = 1.2*base_size_figs/18)
  # tm_layout(legend.outside.position = "right",legend.outside = T,
  #           legend.position = c("left","center"),
  #           # legend.just = c("left","center")
  #           )
map_1


df_profile_NED1$group <- "NED1"
df_profile_NED2$group <- "NED2"
df_profile_NED3$group <- "NED3"

df_profiles <- rbind(df_profile_NED1, df_profile_NED2, df_profile_NED3)

df_v_path_NED1$group <- "NED1"
df_v_path_NED2$group <- "NED2"
df_v_path_NED3$group <- "NED3"

df_v_paths <- rbind(df_v_path_NED1, df_v_path_NED2, df_v_path_NED3)

# base_size_figs <- 7
(
fig_profs <- 
df_profiles %>% ggplot(aes(dist, elev)) + geom_line() +
  geom_area(fill = "grey") +
  facet_wrap(~group,scales = "free_x", ncol = 1) + 
  my_ggtheme(base_size = base_size_figs) + 
  geom_line(data = df_v_paths, aes(x = dist, y = elev, 
                                   col = group, linetype = group ), show.legend = T,
            linewidth = base_size_figs/10) +
  scale_linetype_manual(values = c("solid", "longdash", "twodash")) +
  # scale_color_manual(values =  colors_paths, name = "Paths") +
  scale_color_manual(values =  colors_paths) +
  labs(x = NULL, y = "Elevation (masl)", color = "Paths",linetype = "Paths") +
  theme(axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), legend.position = "bottom") + 
  scale_x_continuous(expand = c(0,0)) +
  # scale_y_continuous(limits = range(df_profiles$elev))
  coord_cartesian(y = range(df_profiles$elev)) +
    scale_y_continuous(n.breaks = 4) +
    # guides(color = guide_legend(nrow = 2),)
    theme(legend.key.width = unit(2, "line"))
)

library(cowplot)
(
  fig_c <- 
    plot_grid(ggmap_1, fig_profs,rel_widths = c(1,1.4))
)

# fig_c <- 
#   plot_grid(ggmap_1, fig_profs,rel_widths = c(1,1.5))

# png("article/document/figs/fig5.pdf",
#     width = 190,height = 100,
#     units = "mm",res = 500)
# fig_c
# dev.off()
ggsave(fig_c,filename = file_fig_article,
       width = 140,height = 75, units = "mm")
ggsave(fig_c,
       filename = file_fig,
       width = 140,height = 75, units = "mm")
ggsave(fig_c,
       filename = file_fig_png,
       width = 140,height = 75, units = "mm",dpi = 300)
# ggsave(fig_c,
#        filename = file_fig,
#        width = 90,height = 55, units = "mm")

