library(tidyverse)
devtools::load_all("~/projects/my_packages/my_themes/")

file_fig <- "03_raster_interpolation//plot_errors_lambdas.pdf"
file_fig_article <- 
  "article/figs/fig_errors_lambdas.pdf"


list_df_errors_lambdas_NED1 <- 
  readRDS("03_raster_interpolation//list_df_error_lambdas_NED1.rds")

df_errors_NED1 <- 
lapply(1:12, function(i){
  list_df_errors_lambdas_NED1[[i]] %>% 
    mutate(method = "IM - NED1",region = i, .before = 1) %>% 
    mutate(min_error = error == min(error))
}) %>% do.call(rbind,.)


list_df_errors_lambdas_NED2 <- 
  readRDS("03_raster_interpolation/list_df_error_lambdas_NED2.rds")

df_errors_NED2 <- 
lapply(1:12, function(i){
  list_df_errors_lambdas_NED2[[i]] %>% 
    mutate(method = "IM - NED2",region = i, .before = 1) %>% 
    mutate(min_error = error == min(error))
}) %>% do.call(rbind,.)


list_df_errors_lambdas_NED3 <- 
  readRDS("03_raster_interpolation/list_df_error_lambdas_NED3.rds")

df_errors_NED3 <- 
  lapply(1:12, function(i){
    list_df_errors_lambdas_NED3[[i]] %>% 
      mutate(method = "IM - NED3",region =  i, .before = 1) %>% 
      mutate(min_error = error == min(error))
  }) %>% do.call(rbind,.)


df_errors <- 
  rbind(df_errors_NED1,
        df_errors_NED2,
        df_errors_NED3)

df_errors <- df_errors %>% 
  mutate(region = paste("Reg.", region),
         region = factor(region,levels = unique(region)))

base_size_fig <- 8

(
  plot_errors <- 
    df_errors %>%
    # filter(region %in% paste("Region", 1:6)) %>% 
    ggplot(aes(lambdas, error, col = min_error)) +
    geom_point(size = 0.5) +
    scale_color_manual(values = c("black","red"),) +
    # facet_grid(region~method)
    geom_point(filter(df_errors,min_error), 
               mapping = aes(lambdas, error),size = 0.9, col = "red") +
    ggh4x::facet_grid2(region ~ method,scales = "free_y") +
    my_ggtheme(base_size = base_size_fig) + 
    labs(x = "\u03bb Coefficient", y = "RMSE") +
    theme(strip.text.x = element_text(size = base_size_fig*1.1),
          strip.text.y = element_text(size = base_size_fig*1.1),
          legend.position = "none" )
)


cairo_pdf(file_fig,
          width = 140 * 0.0393701, 
          height = 170 * 0.0393101 )
print(plot_errors)
dev.off()

cairo_pdf(file_fig_article,
          width = 140 * 0.0393701, 
          height = 170 * 0.0393101 )
print(plot_errors)
dev.off()

