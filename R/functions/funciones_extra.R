
rem_stations <- function(data, df_rem_stas){
  missing_stas <- F
  lapply(df_rem_stas$namcod, function(sta){
    if (!(sta %in% names(data)[-1])){
      message(paste(sta, " not found"))
      missing_stas <- T
    }
  })
  
  if(missing_stas){
    stop("There are missing stations")
  }
  
  name_dates <- names(data)[1]
  names_stas <- names(data)[-1]
  stas_chosen <- names_stas[!(names_stas %in% df_rem_stas$namcod)]
  data[c(name_dates, stas_chosen)]
}

join_stations <- function(data, df_pairs , rem_sta2 = T){
  require(tibble)
  missing_stas <- F
  lapply(c(df_pairs$BASE, df_pairs$TO_JOIN), function(sta){
    if (!(sta %in% names(data)[-1])){
      message(paste(sta, " not found"))
      missing_stas <- T
    }
  })
  
  if(missing_stas){
    stop("There are missing stations")
  }
  
  if (!(is.data.frame(data) | tibble::is_tibble(data))){
    stop(paste("data should be a data frame or a tibble"))  
  }
  if (!(is.data.frame(df_pairs) | tibble::is_tibble(df_pairs))){
    stop(paste("df_pairs should be a data frame or a tibble"))  
  }
  if (!inherits(data[[1]],"POSIXct")){
    stop("the first variable of data should be POSIXct")
  }
  if (!("STA2_FROM" %in% names(df_pairs))){
    df_pairs$STA2_FROM = NA
    warning("STA2_FROM not found")
  }
  if (!("STA2_TO" %in% names(df_pairs))){
    df_pairs$STA2_TO = NA
    warning("STA2_TO not found")
  }
  if(!all(is.na(df_pairs$STA2_FROM))){
    if (inherits(df_pairs$STA2_FROM, "Date")){
      df_pairs$STA2_FROM <- as.POSIXct(df_pairs$STA2_FROM)
    }
    else if (!inherits(df_pairs$STA2_FROM, "POSIXct")){
      stop("the variable 'STA2_FROM' of df_pairs should be POSIXct")
    }
  }
  if(!all(is.na(df_pairs$STA2_TO))){
    if (inherits(df_pairs$STA2_TO, "Date")){
      df_pairs$STA2_TO <- as.POSIXct(df_pairs$STA2_TO)
    }
    else if (!inherits(df_pairs$STA2_TO, "POSIXct")){
      stop("the variable 'STA2_TO' of df_pairs should be POSIXct")
    }
  }
  
  
  dates <- data[[1]]
  for (i in 1:nrow(df_pairs)){
    sta_base <-  df_pairs$BASE[i]
    sta_2 <- df_pairs$TO_JOIN[i]
    vec_1 <- data[[sta_base]]
    vec_2 <- data[[sta_2]]

    start_sta2 <- df_pairs$STA2_FROM[i]
    if (!is.na(start_sta2)){
      vec_2 <- ifelse(dates >= start_sta2, vec_2, NA)
    }
    end_sta2 <- df_pairs$STA2_TO[i]
    if (!is.na(end_sta2)){
      vec_2 <- ifelse(dates < end_sta2, vec_2, NA)
    }
    vec <- ifelse(!is.na(vec_1), vec_1, vec_2)
    data[[df_pairs$BASE[i]]] <- vec
    if (rem_sta2){
      data[,sta_2] <- NULL
    }
  }
  return(data)
}



rem_chg_values <- function(data, df_rules, check_value = T){
  namcods_rules <- df_rules$namcod
  dates <- data[[1]]
  
  # check if ther is duplicated dates
  if (any(duplicated(dates))){
    stop("There are duplicated dates")
  }
  
  # check if all namcods of df_rules exist in data
  vec_is_sta <- !(namcods_rules %in% (names(data)[-1]))
  if (any(vec_is_sta)){
    namcods_notfound <- namcods_rules[vec_is_sta]
    text <- as.list(namcods_notfound) %>% c(sep = "\t") %>% 
      paste0("\t not found in data")
    stop(text)
  }
  
  # removing values
  if ("remove" %in% names(df_rules)){
    df_rem <- df_rules[df_rules$remove,]
    for(i in 1:nrow(df_rem)){
      date_i <- df_rem$dates[i]
      namcod_i <- df_rem$namcod[i]
      value_i <- df_rem$value[i]
      i_row <- which(dates == date_i) 
      # check if value match
      if (check_value){
        if (value_i != data[[namcod_i]][i_row]){
          stop(paste(namcod_i, ": Value for",date_i,
                     "does not match"))
        }
      }
      # remove
      data[[namcod_i]][i_row] <- NA
    }
    # updaqting df_rules
    df_rules <- df_rules[!df_rules$remove,]
  }
  if ("change" %in% names(df_rules)){
    df_chg <- df_rules[!is.na(df_rules$change),]
    for(i in 1:nrow(df_chg)){
      date_i <- df_chg$dates[i]
      namcod_i <- df_chg$namcod[i]
      value_i <- df_chg$value[i]
      i_row <- which(dates == date_i) 
      # check if value match
      if (check_value){
        if (value_i != data[[namcod_i]][i_row]){
          stop(paste(namcod_i, ": Value for",date_i,
                     "does not match"))
        }
      }
      # change
      data[[namcod_i]][i_row] <- df_chg$change[i]
    }
  }
  data
}






# df periods should contain some columns: namcod, INICIO, FIN
# namcod is character, INICIO and FIN should be POSIXct
# data should contain as first column called dates
rem_periods <- function(data, df_periods, 
                           ignore_missing_namcod = F){
  # if(!inherits(df_periods$INICIO, "POSIXct")){
  #   
  # }
  if (!ignore_missing_namcod){
    not_in_data <- !(df_periods$namcod %in% names(data)[-1])
    if (any(not_in_data)){
      for (i in which(not_in_data)){
        message(paste(df_periods$namcod[i], "not found!"))
      }
      stop("check namcods")
    }
  }
  
  
  for (fila in 1:nrow(df_periods)){
    if(df_periods$namcod[fila] != '') sta <- df_periods$namcod[fila]
    if(!is.na(df_periods$INICIO[fila]) |  
       !is.na(df_periods$FIN[fila])){
      if(!is.na(df_periods$INICIO[fila])){
        con_ini <- df_periods$INICIO[fila] < data$dates
      } else con_ini <- T
      if(!is.na(df_periods$FIN[fila])){
        con_fin <- df_periods$FIN[fila] > data$dates
      } else con_fin <- T
      data[[sta]][con_ini & con_fin] <- NA
    }
  }
  data
}


### data should have a first column called dataes
sel_stas_by_length <- function(data, n_min,
                               start_date = NULL, 
                               end_date = NULL){
  namcods_0 <- names(data)[-1]
  data_eval <-  data
  
  if (!is.null(start_date)){
    if (!inherits(start_date, "POSIXct")) {
      start_date <- as.POSIXct(start_date)
    }
    data_eval <- data_eval %>% filter(dates >= start_date)
  }
  
  if (!is.null(end_date)){
    if (!inherits(end_date, "POSIXct")) {
      end_date <- as.POSIXct(end_date)
    }
    data_eval <- data_eval %>% filter(dates < end_date)
  }

  
  namcods_sel <- 
    data_eval %>%  
    summarise_at(-1,function(x) sum(!is.na(x))) %>% 
    as.numeric() %>% (function(x) which(x >= n_min)) %>% 
    namcods_0[.]
  
  data[c("dates", namcods_sel)]
}


# this function change coordinates of a sf point considering
# new coordiantes fromn a data frame (containing a namcod variable) 
correct_coords <- function(sf_loc, new_df, names_coords = c("x","y")){
  if (!("namcod"  %in% names(sf_loc))){
    stop("sf_loc should contain namcod")
  }
  if (!("namcod"  %in% names(new_df))){
    stop("new_df should contain namcod")
  }
  try(
    new_df <- new_df[c("namcod", names_coords)]
  ) 
  if (ncol(new_df) == 3){
    
    for (i in 1:nrow(new_df)){
      new_coords <- as.numeric(new_df[i,-1])
      namcod_i <- new_df$namcod[i]
      if (!(namcod_i  %in% sf_loc$namcod)){
        stop(paste("row",i," : ",namcod, " not found in sf_loc"))
      } 
      st_geometry(sf_loc[sf_loc$namcod == namcod_i, ]) <-  
        st_sfc(st_point(new_coords))
    }
    
  } else stop("the format is not correct")
  
  sf_loc
}
