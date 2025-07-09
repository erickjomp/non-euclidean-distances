library(tidyverse)
library(sf)
library(gstat)
# library(geoR)
library(automap)
library(hydroGOF)

Sys.setenv(TZ='GMT')

source("R/functions/funciones_eje_OK.R")

max_ratio_nugget_sill <- 0.25
min_range <- 30
max_range <- 200
model_sel <- "Sph"
type = "KED" #change for "OK
formula_inter <- as.formula(pr~elev)  #change for as.formula(pr~1)

file_save_estimates_gen <- 
  "04_cross_validation//estimates/list_df_m_inregs_est_%s-%s.rds"


sf_stas_reg_all <- 
  st_read("01_data_byregion//shp_stations//sf_stas_reg.shp")

sf_stas_regbuf_all <- 
  st_read("01_data_byregion//shp_stations//sf_stas_regbuf.shp")

df_m_all <- readRDS("00_data/data_pr_validated&selected//df_m.rds")


#### creting input variograms for gstat ####
list_vars_regs_OK_0 <- 
  lapply(1:12, function(id_reg){
    sf_stas_regbuf <- sf_stas_regbuf_all %>% 
      filter(region == id_reg)
    df_m <- df_m_all %>% select(dates, sf_stas_regbuf$namcod)
    
    message(paste("REGION ", id_reg))
    list_vars <- 
      lapply(1:nrow(df_m), function(i_mes){
        message(paste("\r",i_mes),appendLF = F)
        values_pr <- df_m[i_mes,-1] %>% as.numeric()
        sf_stas_gstat_0 <- sf_stas_regbuf %>% 
          mutate(pr = values_pr)
        sp_stas_gstat_0 <- as_Spatial(sf_stas_gstat_0)
        
        lapply(1:nrow(sf_stas_gstat_0), function(i){
          # print(i)
          vario_af <- NULL
          if(is.finite(values_pr[i]) &
             sf_stas_regbuf$in_region[i] == 1){
            sp_stas_gstat <- sp_stas_gstat_0[-i,]
            sp_stas_gstat <- 
              sp_stas_gstat[!is.na(sp_stas_gstat$pr),]
            
            fixed_values <- c(NA,NA,NA)
            
            iter = 0
            
            while(TRUE){
              
              try(silent =T,
                  vario_af <-
                    autofitVariogram_v2(
                      formula = formula_inter,
                      sp_stas_gstat,
                      model = model_sel,
                      fix.values = fixed_values,
                      diagonal = 314.2857, # thats for 110 km
                      #diagonal = 306.2099,#342.8571,#285.7143,
                      # 110 km : 314.2857
                      # start_vals = c(0,50,NA),
                      miscFitOptions = 
                        list(merge.small.bins = F)
                    )
                  # min.np.bin = 3)),
                  
              )
              
              if (!is.null(vario_af)){
                nugget <- vario_af$var_model$psill[1]
                sill <- vario_af$var_model$psill[2] + nugget
                range <- vario_af$var_model$range[2]
                
                ratio_nugget_sill <- (nugget/sill)
                ratio_nugget_sill <-  ifelse(sill == 0, 0, 
                                             ratio_nugget_sill)
                
                condition_nugget <- 
                  ratio_nugget_sill  <= max_ratio_nugget_sill
                condition_range_1 <- range >= min_range
                condition_range_2 <- range <= max_range
              }
              else{
                condition_nugget <- F
                condition_range_1 <- T
                condition_range_2 <- T
              }
              if(condition_nugget & condition_range_1 & 
                 condition_range_2){
                break  #if conditions are satisfied , the loop ends
              }
              iter <- iter +1
              if(iter > 10){
                message(paste0(" ",i_mes,"-10it"))
                break  #if conditions are satisfied , the loop ends
              }
              
              if ((!condition_range_1 || !condition_range_2) & (
                !is.na(fixed_values[1]) || condition_nugget
              )){ 
                # the 2nd condition is to ensure tryinf nugget first 
                if(!condition_range_1){
                  fixed_values[2] <- min_range
                }
                
                if(!condition_range_2){
                  fixed_values[2] <- max_range
                }
                if (!is.null(vario_af)){
                  vario_af$var_model$range[2] <- fixed_values[2]
                }
                
              } 
              
              if (!condition_nugget){
                if (!is.null(vario_af)){
                  fixed_values[1] <- 
                    min(sd(sp_stas_gstat$pr),0.25*sill)
                  vario_af$var_model$psill[1] <- fixed_values[1] 
                } 
                else{
                  fixed_values[1] <- sd(sp_stas_gstat$pr)
                } 
              } 
              
            }
          }
          return(vario_af)
        })
      })
  })


#### observed series by region ####
list_df_m_inregs <- lapply(1:12, function(id_reg){
  namcods_sel <- 
    sf_stas_reg_all %>% 
    filter(region == id_reg) %>% .$namcod
  df_m_all %>% select(dates,namcods_sel)
})

#### estimated series by region ####
list_df_m_inregs_est_OK_0 <- 
  lapply(list_df_m_inregs, function(df_m_inreg){
    df_m_inreg %>% mutate_at(-1, function(x) NA)
  })

for (id_reg in 1:12){
  sf_stas_regbuf <- sf_stas_regbuf_all %>% 
    filter(region == id_reg)
  df_m <- df_m_all %>% select(dates, sf_stas_regbuf$namcod)
  
  message(paste("REGION ", id_reg))
  
  for(i_mes in 1:nrow(df_m)){
    message(paste("\r",i_mes),appendLF = F)
    values_pr <- df_m[i_mes,-1] %>% as.numeric()
    sf_stas_gstat_0 <- sf_stas_regbuf %>% 
      mutate(pr = values_pr)
    sp_stas_gstat_0 <- as_Spatial(sf_stas_gstat_0)
    # sp_stas_gstat_0$pr_trans <- 
    #   fun_trans(sp_stas_gstat_0$pr)
    
    for (i in 1:nrow(sp_stas_gstat_0)){
      if(is.finite(values_pr[i]) &
         sp_stas_gstat_0$in_region[i] == 1){
        sp_stas_gstat <- sp_stas_gstat_0[-i,]
        sp_stas_gstat <- 
          sp_stas_gstat[!is.na(sp_stas_gstat$pr),]
        
        vario <-
          list_vars_regs_OK_0[[id_reg]][[i_mes]][[i]]
        sp_sim <- 
          krige(
            formula_inter, 
            locations = sp_stas_gstat,
            newdata = sp_stas_gstat_0[i,], 
            model = vario$var_model,
            debug.level = 0)
        
        namcod_sel <- sp_stas_gstat_0[i,]$namcod
        
        # sim_val <- sp_sim$var1.pred %>%  fun_trans_rev()
        sim_val <- sp_sim$var1.pred 
        sim_val <- ifelse(sim_val<0,0,sim_val)
        
        list_df_m_inregs_est_OK_0[[id_reg]][[namcod_sel]][i_mes] <-
          sim_val
        
      }
    }
  }
}

### saving ####
saveRDS(list_df_m_inregs_est_OK_0,
        file = sprintf(file_save_estimates_gen, type, model_sel))