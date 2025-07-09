library(tidyverse)
library(sf)
library(terra)


fac_agreg <- 24 # 24 to get 0.02 degrees grid
extra_dist_buf <- 5000
sf_regions <- st_read("00_data/shp/regions.shp")
sf_regions_buf <- 
  st_read("01_data_byregion/shp_stations//sf_regions_buf.shp")
# ras_peru <- rast("~/projects/INTER/raw_data/ras_srtm_peru.tif")
ras_peru <- rast("~/projects/global_data/MERIT_DEM/merit_dem_merged.tif")
path_new_dems <- "02_regions_dem_conecs/raster/DEM/"
sf_stas_regbuf <- 
  st_read("01_data_byregion/shp_stations/sf_stas_regbuf.shp")

ras_peru <- aggregate(ras_peru, 24,na.rm = T)

sf_regions_buf2 <- 
  sf_regions_buf %>% st_buffer(extra_dist_buf,nQuadSegs = 100,
                               endCapStyle = "ROUND")

lapply(1:nrow(sf_regions_buf2), function(i){
  sf_reg_buf <- sf_regions_buf2 %>% filter(region == i) 
  sf_pol <- sf_stas_regbuf %>% filter(region == i) %>% 
    st_union() %>% st_convex_hull()
  sf_mask <- st_union(sf_reg_buf,sf_pol)
  
  ras_dem_reg <- crop(ras_peru, sf_reg_buf, snap = "out")
  ras_dem_reg <- mask(ras_dem_reg,sf_mask, touches = T)

  name_ras_dem_reg <- paste0("ras_dem_reg",i,".tif") 
  
  writeRaster(ras_dem_reg,
              file.path(path_new_dems,name_ras_dem_reg),
              overwrite = T)
})
