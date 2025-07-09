
# locations should be sf

NED_interpolation_0 <- function(formula, locations, 
                              data = NULL, 
                              new_data = NULL,
                              list_mat_dists = NULL, 
                              list_mat_dists_loc, 
                              only_finite = FALSE, 
                              idp = 2,
                              out_process_data = FALSE, 
                              ratio_dif_WS_min = 0.0001,
                              list_mat_weights = NULL,
                              list_mat_weights_loc = NULL,
                              min_0 = FALSE,
                              i_target_stas = NULL,
                              cross_validation = FALSE){
  gc_after_a_number <- 10
  fac_inf <- 1000 # when a weight is Inf (bc of some dist = 0), ->
  # the Inf is replaced by   (sum of the other weights) x fac_inf 
  # it should be large number
  
  # NA or Nan values were found in som cells of the raster interpolation
  # tthat happened in grid cells where a station was present
  # the estimation of distances (see arguments) should be improved 
  # this estimation should be more accurate for these grid cells
  #         
  # initially it was solved by replacing weights with Inf with 0 and 1
  # but it didnt worked
  #  finally the solutions was just replace Inf with large numbers
  if(cross_validation){
    new_data <- locations
    list_mat_dists <- list_mat_dists_loc
    list_mat_weights <- list_mat_weights_loc
  }
  
  
  if (inherits(i_target_stas, "integer")){
    i_target_stas0 <- i_target_stas
    i_target_stas <- rep(FALSE, nrow(locations))
    i_target_stas[i_target_stas0] <- TRUE
    
  }
  
  ## checking data 
  if(!inherits(new_data, c("SpatRaster", "sf", "numeric"))){
    stop ("new_data must be SpatRaster, sf or numeric")
  }
  if(!inherits(locations, c("sf"))){
    stop ("locations must be sf")
  }
  

  
  ### formula handling 
  c_for <- as.character(formula)
  if ((length(c_for) != 3) || (c_for[1] != "~")){
    stop("formula must be of the the form:  predicted_var ~ explanatory_var")
  } else {
    pred_var <- c_for[2]
    exp_var <- c_for[3]
  }
  
  # lambdas
  lambdas = names(list_mat_dists)
  if(is.null(lambdas)){
    lambdas = paste0("X",1:length(list_mat_dists))
  }
  
  ### data handling 
  if (!is.null(data)){
    is_dates <- FALSE
    is_data_df <- FALSE
    if(is.data.frame(data)){
      is_data_df <- TRUE
      name_data <- names(data)[1]
      if(inherits(data[[1]], c("Date", "POSIXct"))){
        is_dates <- TRUE
        dates <- data[[1]]
        data <- as.matrix(data[,-1])
      }
      else {
        data <- as.matrix(data)
      }
    } # data is alwauys converted to matrix
    n_dates <- nrow(data)
  }

  if(!(exp_var %in% names(locations)) && is.null(data)){
    stop("if data is not provided, locations must contain the explanatory variable")
  } else if (is.null(data)){
    data = matrix(locations[[exp_var]],nrow = 1)
    n_dates = 1
  }
  
  
  if (is.null(list_mat_weights)){
    
    if (!only_finite && inherits(new_data, "SpatRaster")){
      is_finite_cells <- is.finite(terra::values(new_data)[,1])
      list_mat_dists <-
        lapply(list_mat_dists , function(x) x[,is_finite_cells]) 
    }
    
    list_mat_weights <- transform_dists_to_weights(list_mat_dists)
    # rm(list_mat_dists)
  }
  
  if(is.null(list_mat_weights_loc)){
    list_mat_weights_loc <- 
      lapply(list_mat_dists_loc, function(mat_dists){
        diag(mat_dists) <- NA # important for WS (weighting scheme)
        1/(mat_dists ^ idp)
      })
    # rm(list_mat_dists_loc)
  }

  # list_col_has_Inf <- lapply(list_mat_weights, function(mat_weights){
  #   which(apply(mat_weights, 2, function(x) any(is.infinite(x))))
  # })
  
  list_mat_weights <- lapply(list_mat_weights, function(mat_weights){
    apply(mat_weights, 2, function(x){
      is_inf <- is.infinite(x)
      if(any(is_inf)){
        x[is_inf] <- fac_inf * sum(x[!is_inf])
      }
      return(x)
    })
  })
  
  # handli8ng newdata
  if(inherits(new_data, "SpatRaster")){
    names(new_data) <- exp_var
    new_data_exp_values = terra::values(new_data)[,1]
  } else  if (inherits(new_data, "sf")){
    new_data_exp_values <- new_data[[exp_var]]
  } else { # in that case it is numeric
    nes_data_exp_values <- new_data
  } 
    

  # loop
  
  vals_expvar_new_data <- 
    if (inherits(new_data, "SpatRaster")){
      values(new_data)
    } else if (inherits(new_data, "sf")){
      new_data[[exp_var]]
    } else { # numeric
      new_data
    }
  
  
  mat_estimates <- matrix(NA, ncol = n_dates, 
                          nrow = length(vals_expvar_new_data))
  
  if(out_process_data){
    vec_lambdas_sel <- rep(NA,n_dates)
    list_df_WS <- list()
    list_lm <- list()
    mat_estimates_lm <- mat_estimates
    mat_residuals <- mat_estimates
    delta_residuals <- rep(NA,n_dates)
    rel_delta_residuals <- rep(NA,n_dates)
  }
  
  for(i in 1:n_dates){
    if( i %% gc_after_a_number == 0){
      gc()
    }
    
    message(format(round(i/n_dates,3)*100, nsmall = 1)," %\r",appendLF = F)
    values_obs <- data[i,, drop = T]
    is_valid_obs <- !is.na(values_obs)
    sf_sel_locs <- locations[is_valid_obs,]
    sf_sel_locs[[pred_var]] <- values_obs[is_valid_obs]
    
    # list_mat_weights_sel <- 
    #   mapply(
    #     list_mat_weights, list_col_has_Inf,
    #     FUN = function(mat_weights, col_has_Inf){
    #       new_mat <- mat_weights[is_valid_obs,] 
    #       for(i_col in col_has_Inf){
    #         new_mat[,i_col] <- replace_Inf_weights(new_mat[,i_col])
    #       }
    #       new_mat
    #   },SIMPLIFY = F)
    # for some reason this solutions caused outliers in grids of stas
    
    list_mat_weights_sel <- 
      lapply(list_mat_weights, function(mat_weights){
        mat_weights[is_valid_obs,]
      })
    
    list_mat_weights_loc_sel <- 
      lapply(list_mat_weights_loc, function(mat_weights_loc){
        mat_weights_loc[is_valid_obs,is_valid_obs]
      })
    
    list_results <- 
    NED_interpolation__core_0(formula = formula,
                            locations = sf_sel_locs, 
                            new_data = new_data_exp_values,
                            # new_data = new_data,
                            list_mat_weights = list_mat_weights_sel,
                            list_mat_weights_loc = list_mat_weights_loc_sel,
                            fun_WS = hydroGOF::rmse,
                            ratio_dif_WS_min = ratio_dif_WS_min,#0.2/100,
                            # 0.0001 ~ 1% x (dif lambda de 100)
                            return_as_new_data = F,
                            min_0 = min_0,
                            i_target_stas = i_target_stas[is_valid_obs], 
                            cross_validation = cross_validation,
                            out_process_data = out_process_data)
    
    
    if (cross_validation || !out_process_data){
      mat_estimates[,i][is_valid_obs] <- list_results
      out_process_data <- FALSE
    }
    
    
    if(out_process_data){
      mat_estimates[,i] <- list_results[[1]]
      # lambdas
      vec_lambdas_sel[i] <- list_results[[2]]
      
      # WSs
      list_df_WS[[i]] <- list_results[[3]]
      delta_residuals[i] <- 
        abs(list_results[[3]][1,2] - min(list_results[[3]][,2])) 
      rel_delta_residuals[i] <- 
        abs(delta_residuals[i]/ list_results[[3]][1,2] ) 
      
      # lm models
      list_lm[[i]] <- list_results[[4]]
      
      # estimates_lm
      mat_estimates_lm[,i] <- list_results[[5]]
      
      # residuals
      mat_residuals[,i] <- list_results[[6]]
    }
  }
  message("\n",appendLF = F)
  
  #building SpatRaster
  if (inherits(new_data, "SpatRaster")){
    field_estimates <- terra::rast(new_data, nly = n_dates)
    values(field_estimates) <- mat_estimates
    rm(mat_estimates)
    if (out_process_data){
      field_estimates_lm <- terra::rast(new_data, nly = n_dates)
      values(field_estimates_lm) <- mat_estimates_lm
      field_residuals <- terra::rast(new_data, nly = n_dates)
      values(field_residuals) <- mat_residuals
    }
  } else { # sf or numeric
    field_estimates <- t(mat_estimates)
    rm(mat_estimates)
    if (out_process_data){
      field_estimates_lm <- mat_estimates_lm
      field_residuals <- mat_residuals
    }
  }
  
  if(out_process_data){
    return (list(field_estimates = field_estimates, 
                 lambdas_sel = vec_lambdas_sel, 
                 list_df_ws = list_df_WS,
                 list_lm = list_lm,
                 field_estimates_lm = field_estimates_lm,
                 field_residuals = field_residuals,
                 delta_residuals = delta_residuals,
                 rel_delta_residuals = delta_residuals
                 ))
  } else {
    return(field_estimates)
  }
  
}

# locations is sf
# matrices of list_mat_weights_loc must of have diagonal values as NA

# complete_LOOV_weigth_scheme: not implemented yet
# cross_validation: not implemented yet
# list_mat_weights must contain only comlums corresponding to finite values of new_data
NED_interpolation__core_0 <- function(formula,
                                    locations, new_data = NULL, 
                                    list_mat_weights = NULL, 
                                    cross_validation = FALSE,
                                    list_mat_weights_loc,
                                    #complete_LOOV_weigth_scheme= FALSE,
                                    fun_WS = hydroGOF::rmse,
                                    ratio_perc_WS_min = NULL,
                                    ratio_dif_WS_min = NULL,
                                    lambdas = NULL,
                                    i_target_stas = NULL,
                                    return_as_new_data = TRUE,
                                    out_process_data = FALSE,
                                    min_0 = FALSE){
  
  pred_var <- as.character(formula)[2]
  exp_var <- as.character(formula)[3]
  
  #### cross validation -----
  if (cross_validation){
    new_data <-  locations
    list_mat_weights <- list_mat_weights_loc
    ## call recursevely in a loop over all observations 
    #TODO#
    int_estimates <- rep(NA, nrow(locations))
    
    for (i in which(!is.na(locations[[exp_var]]))){
      int_estimates[i] <- 
        NED_interpolation__core_0(
          formula = formula,
          locations = locations[-i,],
          new_data = locations[[exp_var]][i],
          list_mat_weights = 
            lapply(list_mat_weights_loc, function(x) x[-i, i,drop = F]),
          list_mat_weights_loc = 
            lapply(list_mat_weights_loc, function(x) x[-i, -i, drop = F]),
          cross_validation = F, # bc it is crossvalidating  
          i_target_stas = i_target_stas[-i],
          min_0 = min_0,
          ratio_perc_WS_min = ratio_perc_WS_min,
          ratio_dif_WS_min = ratio_dif_WS_min
        )
    }
    return(int_estimates)
  }
  
  #### FORM: lambdas reading ----
  if (is.null(lambdas)){
    is_null_perc_min <- is.null(ratio_perc_WS_min)
    is_null_dif_min <- is.null(ratio_dif_WS_min)
    if ((is_null_perc_min || is_null_dif_min) && 
      (!is_null_perc_min || !is_null_dif_min)){ # strong or
      lambdas <- as.numeric(names(list_mat_weights_loc))
      if(any(is.na(lambdas))){
        stop(paste0(
          "lambdas must be provided or the names of list_mat_loc",
          " shoul be coercible to numeric"))
      }
    } else if (!is_null_perc_min && !is_null_dif_min){
      stop("only ratio_perc_WS_min or ratio_dif_WS_min must be provided")
    }
  }
  
  #### FORM: transforming new_data to numeric----
  # all transformed to numeric since it is faster than spatRaster
  if(inherits(new_data,"numeric")){
    values_exp_new_data <- new_data
  } else if (inherits(new_data, "SpatRaster")){
    values_exp_new_data <- terra::values(new_data)[,1]
  } else if (inherits(new_data,"sf")){ # case sf
    values_exp_new_data <- new_data[[exp_var]]
  } else {
    stop("new_data must be numeric, sf or SpatRaster")
  }
  is_fin_new <- is.finite(values_exp_new_data)
  
  #### linear model, estimates_lm and residuals in locations ----
  model_lm <- lm(formula,locations)
  loc_estimates_lm <- predict(model_lm, locations)
  loc_residuals_lm <- residuals(model_lm)
  
  int_estimates_lm <- 
    predict(model_lm, 
            setNames(data.frame(values_exp_new_data[is_fin_new]), 
                     exp_var))
  
  n_estimates <- length(int_estimates_lm)
  
  mat_loc_residuals_for_int <- # this is for interpolation
    matrix(
      rep(loc_residuals_lm, n_estimates), ncol = n_estimates
    )
  
  mat_loc_residuals <- # this is for WS
    matrix(
      rep(loc_residuals_lm, nrow(locations)), ncol = nrow(locations)
    )
  
  
  # if(complete_LOOV_weigth_scheme){
  #   loc_estimates_lm
  #   lapply(1:nrow(locations), function(i_loc){
  #     partial_lm_model <- lm(formula, locations[-i_loc,])
  #     
  #   })
  # }
  
  #### WEIGHTING SCHEME  (WS)----
  mat_estimates_loc_WS <- 
  lapply(list_mat_weights_loc, function(mat_weights_loc){
    loc_estimates_lm + 
      colSums(mat_weights_loc * mat_loc_residuals, na.rm = T) / 
      colSums(mat_weights_loc, na.rm = T)
  }) %>% do.call(rbind,.) 
  
  # calculatin errors in WS
  if (is.null(i_target_stas)){
    errors_WS <- 
      apply(mat_estimates_loc_WS, 1,
            function(x) fun_WS(x, locations[[pred_var]]))
  } else {
    errors_WS <- 
      apply(mat_estimates_loc_WS[,i_target_stas, drop = F], 1,
            function(x) fun_WS(x, locations[[pred_var]][i_target_stas]))
  }
  
  # choowsing min lambda in WS
  if (is.null(ratio_dif_WS_min) && is.null(ratio_perc_WS_min)){
    i_lambda_sel <- which.min(errors_WS)
  } else {
    i_lambda_sel <- 1
    min_lambda <- lambdas[i_lambda_sel]
    min_error <- errors_WS[i_lambda_sel]
    error0 <- min_error
    
    if (!is.null(ratio_dif_WS_min)) {
      for (i in 2:length(errors_WS)) {
        if (errors_WS[i] < min_error) {
          ratio_change <- (min_error - errors_WS[i]) / error0 /
            (lambdas[i] - min_lambda)
          if (ratio_change > ratio_dif_WS_min) {
            i_lambda_sel <- i
            min_error <- errors_WS[i]
            min_lambda <- lambdas[i]
          }
        }
      }
    } else {
      # case when ratio_perc_WS_min has value
      # test this case
      # TODO (test)
      max_range <- error0 - min(errors_WS)
      for (i in 2:length(errors_WS)) {
        if (errors_WS[i] < min_error) {
          ratio_change <- (min_error - errors_WS[i]) /  max_range /
            (lambdas[i] - min_lambda)
          if (ratio_change > ratio_perc_WS_min) {
            i_lambda_sel <- i
            min_error <- errors_WS[i]
            min_lambda <- lambdas[i]
          }
        }
      }
    }
    
  }
  
  #### CACLCULATING VALUES FIELDS ----
  # adding estimates of lm to residuals interpolations
  int_residuals <- 
    colSums(
      list_mat_weights[[i_lambda_sel]] * mat_loc_residuals_for_int,
      na.rm = T) /
    colSums(list_mat_weights[[i_lambda_sel]], na.rm = T)
  int_estimates <- 
    int_estimates_lm + int_residuals
    
  int_estimates <- 
    replace(values_exp_new_data, is_fin_new, int_estimates)
  int_estimates_lm <- 
    replace(values_exp_new_data, is_fin_new, int_estimates_lm)
  int_residuals <- 
    replace(values_exp_new_data, is_fin_new, int_residuals)
  
  
  #### FORM: formatting for output ---
  if(min_0) int_estimates[int_estimates < 0] <- 0
  
  if (return_as_new_data){
    
    if (inherits(new_data, "SpatRaster")){
      int_estimates <- setValues(new_data, int_estimates)
      if (out_process_data){
        int_estimates_lm <- setValues(new_data, int_estimates_lm)
        int_residuals<- setValues(new_data, int_residuals)
      }
      
    } else if(inherits(new_data, "sf")){
      new_data$var.pred <- int_estimates
      int_estimates <- new_data[,"var.pred"]
      if (out_process_data){
        new_data$var.lm <- int_estimates_lm
        new_data$var.residuals <- int_residuals
        int_estimates_lm <- new_data[var.lm]
        int_residuals <- new_data[var.residuals]
      }
    }
    
  }
  
  
  if(out_process_data){
    if (is.null(lambdas)){
      lambdas <- names(list_mat_weights)
      lambdas_num <- as.numeric(lambdas)
      if(all(!is.na(lambdas_num))){
        lambdas <- lambdas_num
      } 
    }

    
    df_WS <- data.frame(lambdas, errors_WS )
    
    return(list(field_estimates = int_estimates, 
                lambda_selected = lambdas[i_lambda_sel],
                df_WS = df_WS,
                model_lm = model_lm,
                field_lm = int_estimates_lm,
                field_residuals = int_residuals
                ))
  } else {
    return(int_estimates)
  }
  
}
  
# mat_dists:  a matrix or a list of matrices
# idp: power of idw
transform_dists_to_weights <- function(mat_dists, idp = 2){
  if(inherits(list_mat_dists,"list")){
    return(
      lapply(mat_dists, function(mat){
        1/(mat^idp)
      })
    )

  } else if(inherits(list_mat_dists, "matrix")){
    return(
      1/(mat^idp)
    )
  }
}

replace_Inf_weights <- function(x){
  is_inf <- is.infinite(x)
  if(any(is_inf)){
    x[is_inf] <- 1
    x[!is_inf] <- 0
    return(x)
    # return(ifelse(is_inf,1,0))
  } else{
    return(x)
  }
}
  
# ave <- new_data
# 
# microbenchmark::microbenchmark(
#   setValues(ave, int_estimates) 
# )
# 
# microbenchmark::microbenchmark(
#   values(ave)[is.finite(values(ave))] <-  int_estimates
# )
# 
# microbenchmark::microbenchmark(
#   ave[][is.finite(ave[])] <-  int_estimates
# )

# vals_dem <- values(new_data)[,1]
# microbenchmark::microbenchmark(
#   NED_interpolation__core_0(formula = formula,
#                           locations = sf_sel_locs,
#                           new_data = new_data ,
#                           list_mat_weights = list_mat_weights_sel,
#                           list_mat_weights_loc = list_mat_weights_loc_sel,
#                           fun_WS = hydroGOF::rmse,
#                           return_as_new_data = F,
#                           out_WS = T)
# )


