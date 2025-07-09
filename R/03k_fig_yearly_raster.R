library(tidyverse)
library(terra)
# library(tmap)
library(sf)
library(metR)
library(wesanderson)
devtools::load_all("../my_packages/my_themes/")
library(ggnewscale)
library(cowplot)
# devtools::load_all("../my_packages/basic_tools/")

id_regs <- c(2,3,4,7,10)

file_yearly_avg_gen <- 
  "03_raster_interpolation//raster_avg_yearly/ras_%s_reg%s.tif"
# list.files("data/11_raster_interpolation/raster_avg_yearly/",
#            pattern = "*_reg3.tif")
file_dem_gen <- 
  "02_regions_dem_conecs//raster/DEM/ras_dem_reg%s.tif"


methods <- c("IED","NED1","NED2","NED3")
labs_methods <-
  c("IED","IED - NED1","IED - NED2","IED - NED3")
# methods <- c("IDW","OK","IDW-V","KED","IED","NED1","NED2","NED3")
# labs_methods <- 
#   c("IDW","OK","IDW-V","KED","IED","IED - NED1","IED - NED2","IED - NED3")

sf_regions <- st_read("00_data//shp//regions.shp")
sf_stasbuf_all <- 
  st_read("02_regions_dem_conecs//shp/sf_stas_regbuf.shp")


list_breaks <- list(
  # c(0,seq(500,1250,by = 250),Inf),  # id_reg = 2
  c(0, seq(600,1400,200), Inf),  # reg 2
  c(0, seq(600,1100,100), Inf),  # reg 3
  c(0, seq(600,1100,100), Inf),  # reg 4
  c(0, seq(600,1400,200), Inf),   # reg 7
  c(0, seq(600,1100,100), Inf)  # reg 10
)

breaks_elev_list = list()
breaks_elev_list[["2"]] = c(seq(0,5000,1000), Inf)
breaks_elev_list[["3"]] = c(-Inf,c(1500,2000,2500,3000,3500,4000,4500,5000), Inf)
breaks_elev_list[["4"]] = c(-Inf,c(1500,2000,2500,3000,3500,4000,4500,5000), Inf)
breaks_elev_list[["7"]] = c(-Inf,seq(2000,4000,500), Inf)
breaks_elev_list[["10"]] = c(-Inf,c(1500,2000,2500,3000,3500,4000,4500,5000), Inf)

colors_elev_list = list()
colors_elev_list[["2"]] = c("#78c679BF", "#c3ea95BF", "#ffff8cBF", "#fdae61BF", "#a6611aBF", "#543005BF")
colors_elev_list[["3"]] = c("#78c679BF", "#9fd36dBF", "#c3ea95BF", "#ffff8cBF", "#fee08bBF", "#fdae61BF", "#f46d43BF", "#a6611aBF", "#543005BF")
# colors_elev_list[["3"]] = c("#78c679BF", "#c3ea95BF", "#ffff8cBF", "#fee08bBF", "#fdae61BF", "#f46d43BF", "#a6611aBF", "#543005BF")
colors_elev_list[["4"]] = c("#78c679BF", "#9fd36dBF", "#c3ea95BF", "#ffff8cBF", "#fee08bBF", "#fdae61BF", "#f46d43BF", "#a6611aBF", "#543005BF")
colors_elev_list[["7"]] = c("#78c679BF", "#c3ea95BF", "#ffff8cBF", "#fdae61BF", "#a6611aBF", "#543005BF")
colors_elev_list[["10"]] = c("#78c679BF", "#9fd36dBF", "#c3ea95BF", "#ffff8cBF", "#fee08bBF", "#fdae61BF", "#f46d43BF", "#a6611aBF", "#543005BF")


#### output
file_save_fig_gen <- 
  "03_raster_interpolation//ras_yearly_avg_methods_reg%s.pdf"
file_save_fig_article_gen <- "article//figs/yearly_avg_methods_reg%s.pdf"
# file_save_fig_gen <- 
#   "plots/11_raster_interpolation/ras_yearly_avg_methods_reg%s.pdf"
# files_save_fig_article <- 
#   c("article/document/figs/fig9.pdf",
#     "article/document/figs/fig10.pdf")



list_plots <- 
  lapply(id_regs, function(id_reg){
    
    sf_reg <- sf_regions %>% filter(region == id_reg)
    sf_stasbuf <- sf_stasbuf_all %>% filter(region == id_reg)
    
    files_ras <- sprintf(file_yearly_avg_gen, methods,id_reg)
    
    ras_methods <- lapply(files_ras, rast) %>% 
      do.call(c, .)
    # ras_methods %>% plot
    names(ras_methods) <- labs_methods
    
    ras_methods0 <- ras_methods %>% crop(sf_reg,snap = "out")
    ras_methods <- ras_methods %>%
      mask(sf_reg,touches = T)
    
    sf_rec <- 
      # st_bbox(ras_methods[[1]]) %>% 
      st_bbox(ras_methods0) %>% 
      st_as_sfc(crs = 4326)  %>% st_difference(sf_reg$geometry)
    
    df_ras_methods <- as.data.frame(ras_methods0,xy = TRUE)
    df_ras_methods <- df_ras_methods %>% 
      gather("method","pry",-c(1:2)) %>% 
      mutate(method = factor(method,levels = names(ras_methods)))
    
    # breaks <-  c(0,seq(500,1250,by = 250),Inf)
    breaks <- list_breaks[[which(id_reg == id_regs)]]
    base_size_figs <- 7
    
    # sf_ras_methods <- st_as_sf(df_ras_methods, coords = c("x","y"),crs = 4326)
    
    
    df_ras_differences = df_ras_methods %>% spread("method","pry") %>% 
      mutate_at(4:6,~ .x - `IED`) %>% select(-IED) %>% gather(method, pry,-c(1,2)) %>% 
      mutate(method = factor(method,levels = unique(method)))
    # mutate(across(`IED - NED1`, `IED - NED2`, `IED - NED3`), ~ .x  -`IED`)
    
    # reading topo
    file_dem = sprintf(file_dem_gen, id_reg)
    ras_dem = terra::rast(file_dem)
    ras_dem0 <- ras_dem #%>% crop(sf_stasbuf,snap = "out")
    
    
    df_ras_dem <- as.data.frame(ras_dem0,xy = TRUE)
    df_ras_dem <- df_ras_dem %>% 
      gather("method","elev",-c(1:2)) #%>% 
    # mutate(method = factor(method,levels = names(ras_methods)))
    
    # names(df_ras_dem)
    
    # plot methods pry
    ( plot_ras_methods <-
        ggplot(data = df_ras_methods) + 
        # geom_tile(aes(fill = z)) +
        # geom_contour_filled(aes(x = gdat$x, y = gdat$y,z = distance)
        #                     ,breaks = breaks,color = NA) + # bins = 18,
        geom_contour_fill(aes(x = df_ras_methods$x, y = df_ras_methods$y, z = pry,
                              fill = stat(level)),
                          breaks = breaks,color = NA) +
        guides(fill = guide_colorsteps(#show.limits = T,
          barheight = 10, draw.lim = T, barwidth = 0.8,
          # barheight = 1, draw.llim = T, barwidth = 12
          # title.position = "top"
        )) +
        scale_fill_viridis_d(option = "viridis", direction = -1) + 
        # scale_fill_manual(values=wes_palette(n=6, name="Zissou1", type="continuous")) +
        # scale_fill_gradient2(high = "blue",mid = "yellow", low = "orange") +
        # scale_fill_brewer(palette="Dark2") + 
        coord_equal() +
        facet_wrap(~method,ncol = 4) + 
        geom_sf(data = sf_rec,fill = "white", colour = NA) + #fill = NA
        geom_sf(data = sf_reg, fill = NA, colour = "black", size = 1.1) +
        geom_sf(data = sf_stasbuf, color = "black", size = 0.8,
                shape = 21, fill = "white") +
        scale_x_continuous(expand = c(0, 0.018)) +
        # scale_y_continuous(expand = c(0, 0)) +
        # my_ggtheme(base_size = base_size_figs) +
        labs(x = NULL, y = NULL,# fill = "Mean Annual\nPrecipitation \n(mm)", 
             title = "Mean annual precipitation (mm)", fill = NULL) + 
        theme_minimal(base_size = base_size_figs) +
        theme(legend.position = "right",
              axis.text = element_blank(),
              axis.ticks = element_blank(),
              plot.title = element_text(size = base_size_figs*1.3,hjust = 0.5),
              panel.grid.major = element_line(colour = "transparent"),
              strip.text.x = element_text(size = base_size_figs*1.2),
              strip.text.y = element_text(size = base_size_figs*1.2))
    )
    
    breaks_diff = c(-Inf, seq(-200,200,50), Inf)
    ( 
      plot_ras_diff <-
        ggplot(data = df_ras_differences) + 
        # geom_tile(aes(fill = z)) +
        # geom_contour_filled(aes(x = gdat$x, y = gdat$y,z = distance)
        #                     ,breaks = breaks,color = NA) + # bins = 18,
        geom_contour_fill(aes(x = df_ras_differences$x, y = df_ras_differences$y, z = pry,
                              fill = stat(level)),
                          breaks = breaks_diff,color = NA) +
        guides(fill = guide_colorsteps(#show.limits = T,
          barheight = 10, draw.llim = T, barwidth = 0.8,
          # barheight = 1, draw.llim = T, barwidth = 12
          # title.position = "top"
        )) +
        scale_fill_brewer(type = "div", palette = "Spectral") +
        # scale_fill_viridis_d(option = "viridis", direction = -1,) + 
        # scale_fill_manual(values=wes_palette(n=6, name="Zissou1", type="continuous")) +
        # scale_fill_gradient2(high = "blue",mid = "yellow", low = "orange") +
        # scale_fill_brewer(palette="Dark2") + 
        coord_equal() +
        facet_wrap(~method,ncol = 4) + 
        geom_sf(data = sf_rec,fill = "white", colour = NA) + #fill = NA
        geom_sf(data = sf_reg, fill = NA, colour = "black", size = 1.1) +
        geom_sf(data = sf_stasbuf, color = "black", size = 0.8,
                shape = 21, fill = "white") +
        # my_ggtheme(base_size = base_size_figs) +
        labs(x = NULL, y = NULL, #, fill = "Precipitation\ndifference\n(mm)") + 
             title = "Precipitation difference:  IED-NED  −  IED (mm)",
             fill = NULL) +  
        scale_x_continuous(expand = c(0, 0.018)) +
        # scale_y_continuous(expand = c(0, 0)) +
        theme_minimal(base_size = base_size_figs) +
        theme(legend.position = "right",
              axis.text = element_blank(),
              axis.ticks = element_blank(),
              plot.title = element_text(size = base_size_figs*1.3,hjust = 0.5),
              panel.grid.major = element_line(colour = "transparent"),
              strip.text.x = element_text(size = base_size_figs*1.2),
              strip.text.y = element_text(size = base_size_figs*1.2))
    )
    
    
    
    # breaks_elev = c(seq(0,5000,1000), Inf)
    breaks_elev = breaks_elev_list[[as.character(id_reg)]]
    df_ras_dem$method = "Elevation"
    
    
    
    
    # topo_colors <- c("#78c679","#c3ea95", "#ffff8c", "#fdae61", "#a6611a", "#543005")
    # topo_colors_alpha <- c("#78c679BF", "#c3ea95BF", "#ffff8cBF", "#fdae61BF", "#a6611aBF", "#543005BF")
    topo_colors_alpha = colors_elev_list[[as.character(id_reg)]]
    
    (
      plot_topo = 
        ggplot(data = df_ras_dem) + 
        # geom_tile(aes(fill = z)) +
        # geom_contour_filled(aes(x = gdat$x, y = gdat$y,z = distance)
        #                     ,breaks = breaks,color = NA) + # bins = 18,
        geom_contour_fill(aes(x = df_ras_dem$x, y = df_ras_dem$y, z = elev,
                              fill = stat(level)),
                          breaks = breaks_elev,color = NA) +
        guides(fill = guide_colorsteps(#show.limits = T,
          barheight = 8, draw.llim = T, barwidth = 0.8,
          # barheight = 0.6, draw.llim = T, barwidth = 7
          # title.position = "top"
        ))  +
        coord_equal() +
        # scale_fill_brewer(type = "div", palette = "Spectral")) +
        scale_fill_manual(values = topo_colors_alpha) + 
        # facet_wrap(~method,ncol = 1)  +
        labs(x = NULL, y = NULL, title = "Elevation (masl)", fill = NULL) +
        # new_scale_fill() +  # <- Reset fill scale here
        geom_sf(data = sf_reg, fill = NA, colour = "black", size = 1.1) +
        geom_sf(data = sf_stasbuf, color = "black", size = 0.8,
                shape = 21, fill = "white") +
        # geom_sf(data = sf_stasbuf, color = "black", size = 0.8,
        #         shape = 21, mapping = aes(fill = pry)) +
        
        theme_minimal(base_size = base_size_figs) +
        scale_x_continuous(expand = c(0, 0)) +
        scale_y_continuous(expand = c(0, 0)) +
        theme(legend.position = "inside",
              legend.position.inside = c(0.0, 0.02),legend.justification = c(0,0),
              axis.text = element_blank(),
              axis.ticks = element_blank(),
              # plot.margin = margin(0, 0, 0, 0),
              panel.grid.major = element_line(colour = "transparent"),
              plot.title = element_text(size = base_size_figs*1.3,hjust = 0.5),
              strip.text.x = element_text(size = base_size_figs*1.2),
              strip.text.y = element_text(size = base_size_figs*1.2),
              legend.margin=margin(0,0,0,0),
              # legend.box.margin=margin(+10,+10,+10,+10)
              )
    )
    
    
    ave2 = cowplot::plot_grid(plot_topo, plot_ras_diff,ncol = 2,rel_widths = c(1,3),
                              labels = c("B","C"))
    ave1 = cowplot::plot_grid(NULL,plot_ras_methods, rel_widths = c(0.03,1))
    plot_final = cowplot::plot_grid(ave1, ave2, nrow = 2,labels = c("A",NA), rel_heights = c(1,1))
    
    file_save_fig = sprintf(file_save_fig_gen,id_reg)
    file_save_fig_article = sprintf(file_save_fig_article_gen,id_reg)
    
    ggsave(plot_final,
           filename = file_save_fig,
           width = 190,height = 175,units = "mm")
    
    ggsave(plot_final,
           filename = file_save_fig_article,
           width = 190,height = 175,units = "mm")
    # ggsave(plot_ras_methods,
    #        filename = sprintf(file_save_fig_gen, id_reg),
    #        width = 90,height = 170,units = "mm")
    # 
    # file_fig_article <- files_save_fig_article[id_reg == id_regs]
    # ggsave(plot_ras_methods,
    #        filename = file_fig_article,
    #        width = 90,height = 170,units = "mm")
  })



