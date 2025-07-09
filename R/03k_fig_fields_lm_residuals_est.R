library(tidyverse)
library(terra)
# library(tmap)
library(sf)
library(metR)
Sys.setenv(TZ='GMT')

devtools::load_all("../my_packages/my_themes/")
source("R/functions/funciones_NED_interpolation.R")
source("R/functions/funciones_NED_interpolation_2.R")

#### inputs ####
id_reg <-  2
date_sel <- as.POSIXct("2012-02-01")

# df_lambdas <- 
#   readRDS("data/11_raster_interpolation/df_lambdas_NED1.rds") 

# best_lambda <- df_lambdas %>% filter(id_reg == id_reg) %>% .$lambda


#### inouts ####
lambdas <-
  read.csv("02_regions_dem_conecs//input/lambdas.csv")[[1]]

file_list_mat_dists_gen <-
  "03_raster_interpolation//mat_NED1_raster/list_mat_NED1_lambdas_reg%s.rds"

file_list_mat_dists_loc_gen <-
  "04_cross_validation//mat_NED1_stas/list_mat_NED1_lambdas_reg%s.rds"

sf_stasbuf_all <-
  st_read("02_regions_dem_conecs//shp/sf_stas_regbuf.shp")

df_m_all <-
  readRDS("00_data/data_pr_validated&selected//df_m.rds")

file_dem_gen <-
  "02_regions_dem_conecs//raster/DEM/ras_dem_reg%s.tif"

sf_regions <- st_read("00_data//shp/regions.shp")

#### output ####
file_fig <- "03_raster_interpolation//fields_lm_residuals_est.pdf"
file_fig_article <- "article/figs/fig_fields_lm_residuals_est.pdf"


#### process ####
message(paste0("REGION ", id_reg))
sf_stasbuf <- sf_stasbuf_all %>% filter(region == id_reg)

df_m <- df_m_all %>%
  dplyr::select(dates, sf_stasbuf$namcod)

file_list_mat_dists <-
  sprintf(file_list_mat_dists_gen, id_reg)
list_mat_dists <-
  readRDS(file_list_mat_dists)

file_list_mat_dists_loc <-
  sprintf(file_list_mat_dists_loc_gen, id_reg)
list_mat_dists_loc <-
  readRDS(file_list_mat_dists_loc)

# list_mat_dists_loc <- list_mat_dists_loc %>%
#   lapply(function(x) x /10000)

ista <- (sf_stasbuf$in_region == 1)

ras_dem <- rast(sprintf(file_dem_gen, id_reg))

# very old function
list_out <-
  NED_interpolation(
    formula = pr ~ elev,
    locations = sf_stasbuf,
    data = df_m,
    #1:26
    new_data = ras_dem,
    list_mat_dists = list_mat_dists[1:50],
    list_mat_dists_loc = list_mat_dists_loc[1:50],
    idp = 2,
    # ratio_dif_WS_min = NULL,#0.0001,
    min_0 = T,
    optim_lambda = 1,
    lambdas = lambdas[1:50],
    i_target_stas = ista,
    out_process_data = T
  )
gc()

i_date <- which(df_m$dates == date_sel)

#### plot ####
plot(list_out$newdata_estimates_lm[[i_date]])
plot(list_out$newdata_residuals[[i_date]])

# df_est_lm <- as.data.frame(list_out$newdata_estimates_lm[[i_date]] %>% 
#                              mask(sf_regions[id_reg,], touches = TRUE),
#                            xy = TRUE)
ras_est_lm0 <- 
  list_out$newdata_estimates_lm[[i_date]] %>% 
  crop(sf_regions[id_reg,], snap = "out")
ras_est_lm <- ras_est_lm0 %>% mask(sf_regions[id_reg,])
df_est_lm <- as.data.frame(ras_est_lm0,
                           xy = TRUE)
names(df_est_lm)[3] <- "value"

ras_residuals0 <- 
  list_out$newdata_residuals[[i_date]] %>% 
  crop(sf_regions[id_reg,], snap = "out")
ras_residuals <- ras_residuals0 %>% mask(sf_regions[id_reg,])
df_residuals <- as.data.frame(ras_residuals0,
                           xy = TRUE)
names(df_residuals)[3] <- "value"

ras_est0 <- 
  list_out$newdata_estimates[[i_date]] %>% 
  crop(sf_regions[id_reg,], snap = "out")
ras_est <- ras_est0 %>% mask(sf_regions[id_reg,])
df_est <- as.data.frame(ras_est0,
                        xy = TRUE)
names(df_est)[3] <- "value"



sf_reg <- sf_regions %>% filter(id == id_reg)
sf_rec <- 
  # st_bbox(ras_methods[[1]]) %>% 
  st_bbox(ras_est) %>% 
  st_as_sfc(crs = 4326)  %>% st_difference(sf_reg$geometry)

  
base_size_figs <- 7

plot_fun <- function(df_grid, breaks_ras , 
                     palette_viridis = "viridis", 
                     barwidth = 5.5){
  ggplot(data = df_grid) + 
    # geom_tile(aes(fill = z)) +
    # geom_contour_filled(aes(x = gdat$x, y = gdat$y,z = distance)
    #                     ,breaks = breaks,color = NA) + # bins = 18,
    geom_contour_fill(aes(x = df_grid$x, y = df_grid$y, z = value,
                          fill = after_stat(level)),
                      breaks = breaks_ras ,color = NA) +
    guides(fill = guide_colorsteps(#show.limits = T,
      barheight = 0.75, draw.llim = T, barwidth = barwidth
      # title.position = "top"
    )) +
    scale_fill_viridis_d(option = palette_viridis, direction = -1) + 
    # scale_fill_manual(values=wes_palette(n=6, name="Zissou1", type="continuous")) +
    # scale_fill_gradient2(high = "blue",mid = "yellow", low = "orange") +
    # scale_fill_brewer(palette="Dark2") + 
    coord_equal(expand = FALSE) +
    # facet_wrap(~method,ncol = 4) + 
    geom_sf(data = sf_rec,fill = "white", colour = NA) + #fill = NA
    geom_sf(data = sf_reg, fill = NA, colour = "black", size = 0.65) +
    geom_sf(data = sf_stasbuf, color = "black", size = 0.9,
            shape = 21, fill = "white") +
    my_ggtheme(base_size = base_size_figs) +
    labs(x = NULL, y = NULL, fill = "Precipitation\n(mm)") + 
    theme_minimal(base_size = base_size_figs) +
    theme(legend.position = "bottom",
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.grid.major = element_line(colour = "transparent"),
          strip.text.x = element_text(size = base_size_figs*1.2),
          strip.text.y = element_text(size = base_size_figs*1.2),
          plot.title = element_text(hjust = 0.5))
}

max(c(values(ras_est_lm)[,1],values(ras_est)[,1]),na.rm = T)
min(c(values(ras_est_lm)[,1],values(ras_est)[,1]),na.rm = T)

# breaks <- seq(0,200,30)
breaks <- seq(90,290,30)
# breaks[c(1,length(breaks))] <- c(-Inf, Inf)
breaks[c(1,length(breaks))] <- c(-Inf, Inf)

(
  plot_est_lm <- 
    plot_fun(df_grid = df_est_lm,breaks_ras = breaks) +
    # ggtitle("Background field (linear model)")+
    ggtitle("Background field")
  )

(plot_est <- 
    plot_fun(df_grid = df_est,breaks_ras = breaks) + 
    ggtitle("Estimated precipitation field"))

range(values(ras_residuals),na.rm = T)
breaks_residuals <- seq(-80,80,20)

(plot_residuals <- 
    plot_fun(df_grid = df_residuals,breaks_ras = breaks_residuals,
             palette_viridis = "mako", barwidth = 7.5) +
    ggtitle("Residuals field") +
    labs(x = NULL, y = NULL, fill = NULL) + 
    scale_fill_brewer(type = "div", palette = "Spectral"))

library(cowplot)
(
  fig_fields <- 
    plot_grid(plot_est_lm, plot_residuals,plot_est,
              rel_widths = c(1,1,1),nrow = 1)
)


#### figure B - lambda selection ##### 
lambda_sel <- list_out$WS$lambdas[which.min(list_out$WS$error)]


(
  plot_lambda_scheme <- 
    ggplot(list_out$WS, aes(lambdas, error)) + 
    geom_point(data = data.frame(x = lambda_sel, 
                                y =  min(list_out$WS$error)),
               aes(x = x, y = y), col = "red", size = 2.5) +
    geom_point(size= 0.6) +
    labs(y = "Mean RMSE", x = expression(lambda)) +
    # coord_fixed(ratio = 45) + 
    labs(title = expression(lambda~selection~procedure)) +
    my_ggtheme(base_size = base_size_figs) +
    theme(plot.title = element_text(hjust = 0.5),
          strip.text = element_text(size = base_size_figs*1.2))
)

(
  plot_lm <- 
    ggplot(list_out$lm[[i_date]]$model, aes(elev, pr)) + 
    labs(y = "Precipitation (mm)", x = "Elevation (masl)") +
    # coord_fixed(ratio = 45) + 
    labs(title = "Fitted linear model") +
    geom_abline(slope = list_out$lm[[i_date]]$coefficients[2],
                intercept =  list_out$lm[[i_date]]$coefficients[1], col = "blue", lwd = 1) +
    # geom_smooth(method='lm',se = FALSE)
    geom_point(size = 0.8) +
    my_ggtheme(base_size = base_size_figs) +
    theme(plot.title = element_text(hjust = 0.5),
          strip.text = element_text(size = base_size_figs*1.2))
)

(
  fig_lm_lambda <- 
    plot_grid(NULL,plot_lm, NULL,plot_lambda_scheme,NULL,nrow = 1,
              rel_widths = c(0.05,0.45,0.05,0.45,0.05),
              labels = c("A",NA,"B",NA,NA))
)

(
  fig_c <- 
    plot_grid(fig_lm_lambda,fig_fields,rel_heights = c(0.5,1), #rel_widths = c(0.5,1), 
              nrow = 2, labels = c(NA, "C") #,labels = c("A","B")
    )
)



ggsave(fig_c, filename =  file_fig,
       width = 140,height = 128, units = "mm")
ggsave(fig_c, filename =  file_fig_article,
       width = 140,height = 128, units = "mm")


# ggsave(fig_fields, filename =  file_fig,
#        width = 140,height = 85, units = "mm")
# ggsave(fig_fields, filename =  file_fig_article,
#        width = 140,height = 85, units = "mm")

# df_m[df_m$dates == "2012-02-01",] %>% gather(key = "namcod", value = "pr", -1) %>%
#   left_join(sf_stasbuf, .) %>%
#   ggplot(aes(elev, pr)) + geom_point()
# plot(list_out$lm[[i_date]]$model[,2:1])
