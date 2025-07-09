library(tidyverse)
library(terra)
# library(tmap)
library(sf)
library(metR)
library(geodata)
devtools::load_all("../my_packages/my_themes/")

#### inputs ####
file_yearly_avg_gen <- 
  "03_raster_interpolation//raster_avg_yearly/ras_%s_reg%s.tif"

method_sel <- "NED1"
id_regs <- 1:12

sf_regions <- st_read("00_data//shp/regions.shp")
sf_stasbuf_all <- 
  st_read("02_regions_dem_conecs//shp/sf_stas_regbuf.shp")


#### outputs ####
file_save_fig <- "03_raster_interpolation//field_all_regs.pdf"
file_save_fig_article <- "article/figs/fig_field_all_regs.pdf"


#### process ####
list_rast <- list()

for (id_reg in id_regs){
  file_to_read <- sprintf(file_yearly_avg_gen, method_sel, id_reg)
  sf_reg_sel <- sf_regions %>% filter(id == id_reg)
  ras_reg <- rast(file_to_read)
  # plot(ras_reg)
  list_rast[[id_reg]] <- terra::mask(ras_reg, sf_reg_sel, touches = TRUE)
}
ras_all <- do.call(mosaic, list_rast) 


#### extra sf objects ####
sf_use_s2(FALSE)
sf_countries <- geodata::gadm(country = c("BRA","BOL", "CHL"), 
                              level = 0, path = tempdir(),
                              resolution = 2) %>% st_as_sf() #%>% 
  # st_crop(expand_bbox(st_bbox(sf_regions)))
sf_peru <- geodata::gadm(country = "PER",path = tempdir(), level = 0,
                         resolution = 2) %>% 
  st_as_sf()


#### plot ####
df_ras_methods <- as.data.frame(ras_all,xy = TRUE)
df_ras_methods <- df_ras_methods %>% 
  gather("method","pry",-c(1:2)) %>% 
  mutate(method = factor(method,levels = names(ras_all)))


base_size_figs <- 7

ggplot(data = df_ras_methods) + 
  geom_sf(data = sf_peru, fill = "white", col = "black") + 
  geom_sf(data = sf_countries, fill = "white", col = "black") +
  geom_contour_fill(aes(x = x, y = y, z = pry,
                        fill = stat(level)),
                    # breaks = breaks,
                    color = NA) +
  guides(fill = guide_colorsteps(#show.limits = T,
    barheight = 1, draw.llim = T, barwidth = 12
    # title.position = "top"
  )) +
  scale_fill_viridis_d(option = "viridis", direction = -1) + 
  # scale_fill_manual(values=wes_palette(n=6, name="Zissou1", type="continuous")) +
  # scale_fill_gradient2(high = "blue",mid = "yellow", low = "orange") +
  # scale_fill_brewer(palette="Dark2") + 
  xlim(st_bbox(sf_regions)[c(1,3)]) +
  ylim(st_bbox(sf_regions)[c(2,4)]) +
  # facet_wrap(~method,ncol = 3) + 
  # geom_sf(data = sf_rec,fill = "white", colour = NA) + #fill = NA
  geom_sf(data = sf_regions, fill = NA, colour = "black", size = 1.1) +
  # geom_sf(data = sf_stasbuf, color = "black", size = 0.8,
  #         shape = 21, fill = "white") +
  my_ggtheme(base_size = base_size_figs) +
  labs(x = NULL, y = NULL, fill = "Mean Annual\nPrecipitation (mm)") + 
  # theme_minimal(base_size = base_size_figs) +
  theme(legend.position = "bottom",
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major = element_line(colour = "transparent"),
        strip.text.x = element_text(size = base_size_figs*1.2),
        strip.text.y = element_text(size = base_size_figs*1.2),
        panel.background = element_rect(fill = 'lightblue'))

ggsave(filename = file_save_fig_article,width = 90, units = "mm")
ggsave(filename = file_save_fig,width = 90, units = "mm")
