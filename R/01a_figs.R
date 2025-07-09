library(tidyverse)
library(sf)
library(tmap)
library(spData)
library(terra)
library(rnaturalearth)
library(rnaturalearthdata)
library(grid)

devtools::load_all("~/projects/my_packages/my_themes/")
Sys.setenv(TZ = 'GMT')

path_figs_article <- "article/figs/"


sf_stas <- 
  st_read("01_data_byregion/shp_stations//sf_stas.shp")

sf_stas_reg <- 
  st_read("01_data_byregion/shp_stations//sf_stas_reg.shp")

sf_stas_regbuf <- 
  st_read("01_data_byregion/shp_stations//sf_stas_regbuf.shp")

sf_regions <- st_read("00_data/shp/regions.shp")
sf_area <- sf_regions %>% st_union() %>% st_as_sf()

ras_peru <- 
  rast("~/projects/global_data/MERIT_DEM/merit_dem_merged.tif")
ras_peru <- ras_peru %>% terra::aggregate(fact = 9, na.rm = T)
values(ras_peru)[values(ras_peru) < 0] <- 0


values_dem <- terra::extract(ras_peru,sf_area)[[2]]
  
df_elev <- rbind(data.frame(elev = sf_stas_reg$elev, 
                            type = "stations"),
                 data.frame(elev = values_dem,
                            type = "srtm"))

#### FIGURE 2 ####
base_size <-  8
bin_w <- 200
sf_stas_reg %>% 
  ggplot(aes(elev)) +
  geom_histogram(mapping = aes(y = 100*bin_w*stat(..density..)),
                 binwidth  = bin_w,position ="identity",
                 alpha = 1,size = base_size/12, col = "white", 
                 fill = "skyblue")  + 
  geom_density(data = df_elev,
               bw = 100,alpha = 0.2,size=base_size/8,
               aes(y = 100*bin_w*stat(density), col = type)) +
  scale_color_manual(values = c("red", "blue"),
                     labels = c("Topography", "Stations")) +
  # labs(x = "Elevation (masl)", y = "Frequency in percent",
  #      col = NULL) + 
  labs(x = "Elevation (masl)", y = "Frequency density (%)",
       col = NULL) +
  # scale_y_continuous() + 
  scale_y_continuous(expand = expansion(mult = c(0.0,0.05)),
                     # labels = as_mapper(~ paste0(.x, " %")),
                     ) + 
  scale_x_continuous(expand = expansion(mult = c(0,0.0))) +
  my_ggtheme(base_size = base_size) +
  # ggpubr::theme_classic2(base_size = 8) + 
  theme(legend.position = c(.85, .8),
        legend.background = element_blank())

ggsave(filename = file.path(path_figs_article, "fig_elevations-hist.pdf"), 
       width = 90,height = 45,units = "mm")
ggsave(filename = 
         file.path("01_data_byregion/fig_elevations-hist.png"), 
       width = 90,height = 45,units = "mm")

#### FIGURE 1 ####
sf_peru <-  GADMTools::gadm_sf_loadCountries("PER")$sf
sf_countries <- 
GADMTools::gadm_sf_loadCountries(c("BOL","PER","CHL",
                                   "BRA","ECU","COL"),)$sf
sf_4countries <- sf_countries %>% filter(ISO %in% 
                                           c("PER","CHL",
                                             "BOL","BRA"))
rm(sf_countries)
# map_nz <- tm_shape(sf_regions) + tm_polygons()

breaks_elev <- c(0:6,Inf)*1000

sf_lakes <- 
  rnaturalearth::ne_download(scale = 10, 
                             type = 'lakes', 
                             category = 'physical') %>% 
  st_as_sf()
sf_lakes <- sf_lakes %>% filter(name_es == "Titicaca")

map_regions1 <- 
  tm_shape(sf_4countries) + tm_fill(col = "white") +
  tm_shape(ras_peru) + 
  tm_raster(alpha = 0.55,breaks = breaks_elev, 
            palette = 
              # terrain.colors(7),
              colorRampPalette(marmap::etopo$colours[7:1])(7),
            legend.format = list(#text.or.more ='or more',
              fun=function(x) formatC(x, digits=0, format="d"),
              text.to.columns = TRUE,
              text.separator = '-'),
            title = 'Elevation (masl)') +
  tm_shape(sf_lakes) + tm_fill(col = "skyblue") + 
  tm_shape(sf_4countries) + tm_borders() +
  tm_shape(sf_regions,is.master = T) + tm_borders(col = 'red') + 
  tm_shape(sf_stas) + tm_dots(col = "black", alpha = 0.5) + 
  tm_shape(sf_regions) +
  tm_text('region',size= 1.4,col= 'blue',fontface="bold") + 
  tm_graticules(lines = F,) +   
  tm_scale_bar(position = c('center','bottom')) + 
  tm_layout(bg.color = "skyblue") +
  tm_legend(legend.position = c("left", "bottom"),
            bg.color = 'white') + 
  tm_compass(position = c("center","top"))

# map_regions2 
sf_world <- ne_countries(scale = "small", returnclass = "sf") %>% 
  filter(region_un == "Americas")
sf_SA <- 
  world %>% 
  filter(continent %in% c("South America")) %>% 
  st_union() %>% st_as_sf()
sf_SA$name = "South\nAmerica"
sf_rect <- st_bbox(sf_regions) %>% st_as_sfc(crs = 4326)
  
map_regions2 <-
tm_shape(world) + tm_borders(col = "grey") + tm_fill(col = "white") +
  tm_shape(sf_SA, is.master=T) + #tm_borders( col = "black") + 
  tm_text(text = "name", size = 0.7,ymod = -1) +
  tm_shape(sf_rect) + tm_borders() + 
  tm_shape(sf_regions) + tm_fill(col = "red") + 
  tm_shape(sf_rect) + tm_borders(col = "blue") + 
  tm_shape(sf_lakes) + tm_fill(col = "skyblue") + 
    tm_layout(bg.color = "skyblue")

map_regions1
# print(map_regions2, 
#       vp = viewport(0.8, 0.8, width = 0.3, height = 0.35))
print(map_regions2, 
      vp = viewport(0.8, 0.78, width = 0.3, height = 0.35))

tmap_save(map_regions1,
          filename = file.path(path_figs_article,"fig_map.jpg"),
          width = 140, units = "mm",dpi = 500,
          insets_tm = map_regions2,
          insets_vp =
            viewport(0.85, 0.79, width = 0.3, height = 0.34))
#viewport(0.85, 0.8, width = 0.3, height = 0.3)
tmap_save(map_regions1,
          filename = file.path(path_figs_article,"fig_map.png"),
          width = 140, units = "mm",dpi = 500,
          insets_tm = map_regions2,
          insets_vp =
            viewport(0.85, 0.79, width = 0.3, height = 0.34))

tmap_save(map_regions1,
          filename = 
            file.path("01_data_byregion//fig_map.png"),
          width = 140, units = "mm",dpi = 500,
          insets_tm = map_regions2,
          # insets_vp =
          #   viewport(0.85, 0.8, width = 0.3, height = 0.3),
          insets_vp =
            viewport(0.85, 0.79, width = 0.3, height = 0.34)
          )
