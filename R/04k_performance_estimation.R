library(tidyverse)
library(hydroGOF)
library(kableExtra)

methods <- c("KED-Sph","IED","NED1","NED2","NED3")
method_base = "IED" 
methods_NE = c("NED1","NED2","NED3")

labs_methods <- 
  c("KED","IED","IED-NED1","IED-NED2","IED-NED3")

list_obs <- 
  readRDS("04_cross_validation//observed/list_df_m_inregs_obs.rds")
file_list_est_method_gen <- 
  "04_cross_validation//estimates/list_df_m_inregs_est_%s.rds"



#### function to select "best values" ####
# f_isbest = function(values, index){
#   if (index == "bias"){
#     values =  abs(values)
#   }
#   
#   if (index %in% c("RMSE", "MAE")){
#     # is_best = values < (min(values) + frac_sd * sd(values))
#     is_best = values < (min(values) + 0.01 * min(values))
#     
#   } else if (index == "bias"){
#     is_best = values < (min(values) + 0.5)
#   } else if (index %in% c( "cor")){
#     # is_best = values > (max(values) - frac_sd * sd(values))
#     is_best = values > (max(values)  - 0.005)
#   }
#   return (is_best)
# }


#### function to identify if NE distances improves vs base case IED ####
# returns 1 if it is better or -1 if it is worse, 0 otherwise
f_isbetter_worse = function(values_all,value_base,values_asses, index){
  if (index == "bias"){
    values_all =  abs(values_all)
    value_base =  abs(value_base)
    values_asses =  abs(values_asses)
  }
  
  if (index %in% c("RMSE","MAE", "bias")){
    # first criteria
    # is_better = values_asses < (value_base - frac_sd * sd(values_all))
    # is_worse = values_asses > (value_base + frac_sd * sd(values_all))
    # second criteria
    if (index %in% c("RMSE","MAE")){
      is_better = (value_base - values_asses) / value_base >0.02
      is_worse = (value_base - values_asses) / value_base < -0.02
    } else if (index == "bias"){
      is_better =  (value_base - values_asses) > 1
      is_worse =  (value_base - values_asses) < -1
    }
    # is_better = is_better & is_better2
    # is_worse = is_worse & is_worse2
  } else if (index %in% c( "cor")){
    # first criteria
    # is_better = values_asses > (value_base + frac_sd * sd(values_all))
    # is_worse = values_asses < (value_base - frac_sd * sd(values_all))
    # second criteria 
    is_better = (values_asses - value_base) > 0.01
    is_worse = (values_asses - value_base) < -0.01
    
    # is_better = is_better & is_better2
    # is_worse = is_worse & is_worse2
  }
  vec_out = rep(0,length(values_asses))
  vec_out[is_better] = 1L
  vec_out[is_worse] = -1L
  return (vec_out)
}




## loop

df_perf_stas <- 
  lapply(1:length(methods), function(i_met){
    lab_method <- labs_methods[i_met]
    list_est_method <- sprintf(file_list_est_method_gen, methods[i_met]) %>% 
      readRDS()
    
    message(lab_method,"\r", appendLF = F)
    lapply(1:12, function(id_reg){
      df_est <- list_est_method[[id_reg]]
      df_obs <- list_obs[[id_reg]]
      lapply(2:ncol(df_obs), function(i){
        sim <- df_est[[i]]
        obs <- df_obs[[i]]
        data.frame(RMSE = hydroGOF::rmse(sim, obs), 
                   # nse = hydroGOF::NSE(sim, obs), 
                   MAE = hydroGOF::mae(sim, obs),
                   bias = mean(sim - obs,na.rm =T),
                   cor = cor(sim, obs, use = "na.or.complete")
                   # pbias = pbias(sim,obs),
        )
      }) %>% do.call(rbind, .) %>% # summarise_all(mean) %>% 
        mutate(id_reg = id_reg,
               namcod = names(df_est)[-1],
               method = methods[i_met],.before = 1)
    }) %>% do.call(rbind,.)
  }) %>% do.call(rbind,.)


df_perf_stas %>% #setNames(names_w_labs) %>%
  mutate(method = case_match(method, 
                             "KED-Sph" ~ "KED",
                             "IED" ~ "IED", 
                             "NED1" ~ "IED-NED1",
                             "NED2" ~ "IED-NED2",
                             "NED3" ~ "IED-NED3")) %>% 
  write.csv("04_cross_validation//df_perf_stas.csv")


df_perf <- df_perf_stas %>% 
  select(-namcod) %>% 
  group_by(id_reg,method) %>% summarise_all(mean)


df_perf <- 
  df_perf %>% gather(index, value, -c(1,2)) %>% 
  mutate(index = factor(index, levels = unique(index))) %>% 
  spread(method, value) %>% 
  select(id_reg, index, all_of(methods)) %>% 
  arrange(id_reg)

### writing to csv
names(labs_methods) <- methods

names_w_labs <- names(df_perf)
names_w_labs[match(methods,names_w_labs)] <- labs_methods[methods]

df_perf %>% setNames(names_w_labs) %>%
  write.csv("04_cross_validation//df_perf.csv",)


####3 processing for tables ####
# names_cols <- names(df_perf)
# names_cols[1:2] <- c("Region", "Index")
names_w_labs[1:2] <- c("Region", "Index")

#### styling final table ### 
df_perf__table <-  df_perf 

# df_perf__table <- 
#   df_perf %>% 
#   mutate_at(-c(1,2), function(x) ifelse(df_perf$index == "cor", round(x, 3), round(x, 2)))

# groups_rows <- rle(df_perf__table$id_reg)
# groups_rows <- groups_rows$lengths %>% 
#   setNames(groups_rows$values)

# to bold r highlight best values  (to bold,  change bold = TRUE)
df_perf__table <- df_perf__table %>% 
  mutate_at(-c(1,2), 
            function(x) ifelse(df_perf$index == "cor", sprintf("%.3f", x), sprintf("%.2f",x)))

for (i in 1:nrow(df_perf__table)){
  values <- df_perf[i,methods] %>% unlist() # unlist better than as.vector()
  index = df_perf$index[i]
  # is_best = f_isbest(values, index)
  
  #### is better worse for IED
  is_betterworse = f_isbetter_worse(values, values[method_base], values[methods_NE], index)
  is_betterworse_ext = rep(0,length(values)) %>% setNames(names(values))
  is_betterworse_ext[methods_NE] <- is_betterworse
  
  #### is better worse for KED
  is_betterworse_KED = f_isbetter_worse(values, values["KED-Sph"], values[c("IED",methods_NE)], index)
  is_betterworse_ext_KED = rep(0,length(values)) %>% setNames(names(values))
  is_betterworse_ext_KED[c("IED",methods_NE)] <- is_betterworse_KED
  
  
  
  values__table <- df_perf__table[i,-c(1,2)] %>% unlist()
  
  
  # df_perf__table[i, -c(1,2)] <- ifelse(is_best, sprintf("**%s**",values), values)
  values__table <-
    cell_spec(values__table, "latex",
              # bold = is_betterworse_ext_KED != 0,
              background = case_when(is_betterworse_ext == 1 ~ "DodgerBlue!30",
                                     is_betterworse_ext == (-1) ~ "Red!20", .default = "Red!0"),  # 	#FFC1CC
              color = case_when(is_betterworse_ext_KED == 1 ~ "blue",
                                is_betterworse_ext_KED == (-1) ~ "red",
                                .default = "black"),
              # background = ifelse(is_best,"DodgerBlue!30","Blue!0"),   #using in latex:  \usepackage[table, dvipsnames, svgnames]{xcolor}
              latex_background_in_cell = TRUE
    )
  
  df_perf__table[i, methods] <-  values__table
  # ifele(is_better, cell_spec(., ))
}


# it is necessary sepecifically for latex
df_perf__table$index <- sprintf("$%s$", df_perf__table$index)



#### TABLE reg 1 to 6 ####
df_perf__table_1to6 <- df_perf__table %>% filter(id_reg <= 6)

df_table_latex1 <- 
  df_perf__table_1to6 %>% select(-1) %>% 
  kbl("latex",booktabs = T,escape = F,
      col.names = names_w_labs[-1],
      align = "c") %>% 
  pack_rows(index = auto_index(paste("Region", df_perf__table_1to6$id_reg)))
# auto_index()
# pack_rows(index = groups_rows[1:6])

save_kable(df_table_latex1,
           "article/tables/table_indices_reg1-6.tex",
           keep_tex = T)

#### TABLE reg 7 to 12 ####
df_perf__table_7to12 <- df_perf__table %>% filter(id_reg > 6, id_reg <= 12 )

df_table_latex2 <- 
  df_perf__table_7to12 %>% select(-1) %>% 
  kbl("latex",booktabs = T,escape = F,
      col.names = names_w_labs[-1],
      align = "c") %>% 
  pack_rows(index = auto_index(paste("Region", df_perf__table_7to12$id_reg)))
# auto_index()
# pack_rows(index = groups_rows[1:6])

save_kable(df_table_latex2,
           "article/tables/table_indices_reg7-12.tex",
           keep_tex = T)



