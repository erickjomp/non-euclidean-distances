library(tidyverse)
library(sf)
library(terra)
library(kableExtra)



#### IMPUTS ####
sf_stas <- 
  st_read("01_data_byregion/shp_stations//sf_stas.shp")
sf_stas_regbuf <- 
  st_read("01_data_byregion/shp_stations//sf_stas_regbuf.shp")

sf_regions <- st_read("00_data/shp/regions.shp")
sf_area <- sf_regions %>% st_union() %>% st_as_sf()

fac_agg_slope = 24 # 24 to be consistent with grid size used for calculating NED

#### processing ####

ras_peru <- 
  rast("~/projects/global_data/MERIT_DEM/merit_dem_merged.tif")
# rast("~/projects/INTER/raw_data/ras_srtm_peru.tif")
# ras_peru <- ras_peru %>% terra::aggregate(fact = 5, na.rm = T)

df_regs_elev <- 
  lapply(1:nrow(sf_regions), function(i){
    print(i)
    ras_reg <- crop(ras_peru, sf_regions[i,],snap = "out")
    ras_slope <- tan(
      terrain(aggregate(ras_reg, fac_agg_slope), unit = "radians")
    ) * 100
    # values_elev <- terra::extract(ras_peru, sf_regions[i,])[[2]]
    values_elev <- terra::extract(ras_reg, sf_regions[i,])[[2]]
    print("cropped")
    slope <- terra::extract(ras_slope,sf_regions[i,],"mean",
                            na.rm = T)[[2]]
    
    df_reg <- data.frame(region = i)
    df_reg$min_elev <- min(values_elev)
    df_reg$mean_elev <- mean(values_elev)
    df_reg$max_elev<- max(values_elev)
    df_reg$mean_slope <- slope
    df_reg
  }) %>% do.call(rbind, .)

df_regs_elev <- df_regs_elev %>% 
  mutate(min_elev = round(min_elev,),
         mean_elev = round(mean_elev),
         max_elev = round(max_elev))

df_regs_pry_stas <- 
  sf_stas %>% filter(!is.na(region)) %>% 
  st_drop_geometry() %>% 
  group_by(region) %>% 
  summarise(min_pry = min(pry),
            mean_pry = mean(pry),
            max_pry = max(pry),
            n_stas = n())

df_n_stas_waux <- sf_stas_regbuf %>% 
  st_drop_geometry() %>% 
  group_by(region) %>% 
  summarise(n_stas_waux = n())

df_regs_pry_stas <- df_regs_pry_stas %>% 
  left_join(df_n_stas_waux, by = "region") %>% 
  mutate(n_stas_aux = n_stas_waux - n_stas)


df_subareas <- 
  read.csv("01_data_byregion//input/df_subareas.csv")


df_table0 <- sf_regions %>% st_drop_geometry() %>% 
  select(region, area_km2) %>% 
  left_join(df_subareas, by = "region") %>% 
  left_join(df_regs_elev) %>% 
  left_join(df_regs_pry_stas) %>% 
  mutate(density_stas = n_stas/area_km2*1000)

df_table <- df_table0 %>% 
  mutate(range = max_elev - min_elev) %>% 
  select(region, area_km2, min_elev, max_elev, range, mean_slope,
         min_pry, mean_pry, max_pry, n_stas, n_stas_aux,density_stas) %>% 
  mutate_at(vars(c("area_km2", "min_pry", "mean_pry", "max_pry")),
            round, digits = 0) %>% 
  mutate_at(vars("mean_slope"),round, digits = 1) %>% 
  mutate_at(vars("density_stas"),round, digits = 2)



groups_rows <- rle(df_table0$subarea)
groups_rows <- groups_rows$lengths %>% 
  setNames(groups_rows$values)

names_cols <- c("Region", 
                # "Subarea", 
                "Area\n(km\\textsuperscript{2})",
                "Min",
                "Max",
                "Altitudinal\nrange (m)",
                "Mean\nslope\n(\\%)",
                "Min",
                "Mean",
                "Max",
                "Number\nof\nstations",
                "Number\nof\ncomplementary\nstations",
                "Stations\n / 1000 km\\textsuperscript{2}") 
# n_names_cols <- setNames(rep(1,ncol(df_table0)), names_cols)

df_table_latex <- 
  df_table %>% 
  kbl("latex",booktabs = T,escape = F,
      col.names = names_cols %>% 
        linebreak(align = "c"),align = "c") %>%
  add_header_above(c(" " = 2, 
                     "Elevation\n(masl)" = 2,
                     " " = 1,
                     " " = 1,
                     "Mean annual\nprecipitation (mm)" = 3,
                     " " = 2),
                   align = "c") %>%
  # collapse_rows(col_names = T) %>% 
  pack_rows(index = groups_rows)

save_kable(df_table_latex,"article/tables/table_regions-summary.tex",
           keep_tex = T,)

