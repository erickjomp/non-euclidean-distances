



# list_mat_weights must contain only comlums corresponding to finite values of new_data

# pred_var and exp_var are character
# locations is sf
# new_data is numeric. this contain the values of exp_var for data to be predicted
# formula must be of the form  name_var_pred ~ name_var_aux
# mat_weights and mat_weights_loc are matrices or a list of matrices
# i_stas_lm is logical or numeric
# param_min_WS must be a list with arguments for choosing best lambda for minimization
#           the arguments of param_min_WS should be fun, mat_weights_loc, ratio_dif_WS_min, ratio_perc_WS_min i_target_stas
#           (see funciones_NED_interpolation.R)


# if mat_weights and mat_weights_loc  are matrices, a vector of pred values is returned
# if they are lists, a matrix is returned

# param_min_WS is provided, the returned object is the vector of pred values corresponding to the minimizing lambda

NED_interpolation__core <- function(formula,
                                    locations,
                                    new_data = NULL,
                                    mat_weights = NULL,
                                    i_stas_lm = NULL,
                                    out_process_data = FALSE,
                                    min_0 = FALSE,
                                    param_min_WS = NULL) {
  
  c_for <- as.character(formula)
  pred_var <- c_for[2]
  exp_var <- c_for[3]
  
  if (!is.null(i_stas_lm)) {
    locations_0 <- locations[i_stas_lm, ]
  } else {
    locations_0 <- locations
  }
  
  if (!is.list(mat_weights)) {
    mat_weights <- list(mat_weights)
    # mat_weights_loc <- list(mat_weights_loc)
  }
  
  
  #### linear model, estimates_lm and residuals in locations ----
  model_lm <- lm(formula, locations_0)
  loc_estimates_lm <- predict(model_lm, locations)
  loc_residuals_lm <- locations[[pred_var]] - loc_estimates_lm
  
  int_estimates_lm <-
    predict(model_lm,
            setNames(data.frame(new_data),
                     exp_var))
  
  n_estimates <- length(int_estimates_lm)
  
  mat_loc_residuals_for_int <- # this is for interpolation
    matrix(rep(loc_residuals_lm, n_estimates), ncol = n_estimates)
  
  
  if (!is.null(param_min_WS)) {
    mat_loc_residuals <- # this is for WS
      matrix(rep(loc_residuals_lm, nrow(locations)), ncol = nrow(locations))
    #TODO  (copy form funciones_NED_interpolation.R)
    # it seems it could be better to create other function
    # mat_weights <- #choose mat_weights
    
  }
  
  #### CACLCULATING VALUES FIELDS ----
  # adding estimates of lm to residuals interpolations
  int_estimates <- list()
    # matrix(NA,
    #        nrow = length(mat_weights),
    #        ncol = length(int_estimates_lm))
  
  for (i in 1:length(mat_weights)) {
    mat_w <- mat_weights[[i]]
    int_residuals <-
      colSums(mat_w * mat_loc_residuals_for_int) /
      colSums(mat_w)
    int_estimates[[i]] <-
      int_estimates_lm + int_residuals
    
    if (min_0) int_estimates[[i]][int_estimates[[i]] < 0] <- 0
  }
  
  
  #### OUTPUT ####
  if (length(int_estimates) == 1)
    int_estimates <- int_estimates[[1]]
  
  if (out_process_data && (length(mat_weights) == 1)) {
    return(
      list(
        estimates = int_estimates,
        estimates_lm = int_estimates_lm,
        estimates_residuals = int_residuals,
        lm= model_lm,
        loc_residuals = loc_residuals_lm
      )
    )
  } else {
    return(int_estimates)
  }
}



NED_interpolation <- function(formula,
                              locations,
                              data = NULL,
                              new_data = NULL,
                              optim_lambda = 0,
                              list_mat_dists = NULL,
                              list_mat_dists_loc = NULL,
                              list_mat_weights = NULL,
                              list_mat_weights_loc = NULL,
                              idp = 2,
                              i_target_stas= NULL,
                              i_stas_lm = NULL,
                              out_process_data = FALSE,
                              min_0 = FALSE,
                              param_min_WS = NULL,
                              lambdas = NULL,
                              cross_validation = FALSE) {
  gc_after_a_number <- 10  # to applt gc()
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
  
  # if (cross_validation) {
  #   new_data <- locations
  #   # list_mat_dists <- list_mat_dists_loc
  #   # list_mat_weights <- list_mat_weights_loc
  # }
  
  if (inherits(i_target_stas, "integer")) {
    i_target_stas0 <- i_target_stas
    i_target_stas <- rep(FALSE, nrow(locations))
    i_target_stas[i_target_stas0] <- TRUE
  }
  
  ## checking data
  if (!inherits(new_data, c("SpatRaster", "sf", "numeric"))) {
    stop ("new_data must be SpatRaster, sf or numeric")
  }
  if (!inherits(locations, c("sf"))) {
    stop ("locations must be sf")
  }
  
  ### formula handling ####
  c_for <- as.character(formula)
  if ((length(c_for) != 3) || (c_for[1] != "~")) {
    stop("formula must be of the the form:  predicted_var ~ explanatory_var")
  } else {
    pred_var <- c_for[2]
    exp_var <- c_for[3]
  }
  
  #### lambdas ####
  # lambdas = names(list_mat_dists)
  # if (is.null(lambdas)) {
  #   lambdas = paste0("X", 1:length(list_mat_dists))
  # }
  
  ### data handling ####
  if (!is.null(data)) {
    is_dates <- FALSE
    is_data_df <- FALSE
    if (is.data.frame(data)) {
      is_data_df <- TRUE
      name_data <- names(data)[1]
      if (inherits(data[[1]], c("Date", "POSIXct"))) {
        is_dates <- TRUE
        dates <- data[[1]]
        data <- as.matrix(data[, -1])
      }
      else {
        data <- as.matrix(data)
      }
    } # data is alwauys converted to matrix
    n_dates <- nrow(data)
  }
  
  #### checking exp_var in locations ###
  if (!(exp_var %in% names(locations)) && is.null(data)) {
    stop("if data is not provided, locations must contain the explanatory variable")
  } else if (is.null(data)) {
    data = matrix(locations[[exp_var]], nrow = 1)
    n_dates = 1
  }
  
  #### calculating list_mat_weights ####
  if (is.null(list_mat_weights) && !is.null(list_mat_dists)) {
    if (inherits(new_data, "SpatRaster")) {
      is_finite_cells <- is.finite(terra::values(new_data)[, 1])
      list_mat_dists <-
        lapply(list_mat_dists , function(x)
          x[, is_finite_cells])
    }
    list_mat_weights <- lapply(list_mat_dists, function(mat_dists) {
      1 / (mat_dists ^ idp)
    })
    # rm(list_mat_dists)
  }
  
  if (is.null(list_mat_weights_loc) && !is.null(list_mat_dists_loc)) {
    list_mat_weights_loc <-
      lapply(list_mat_dists_loc, function(mat_dists) {
        diag(mat_dists) <- NA # important for WS (weighting scheme)
        1 / (mat_dists ^ idp)
      })
    # rm(list_mat_dists_loc)
  }
  
  if (!is.null(list_mat_weights)){
    list_mat_weights <-
      lapply(list_mat_weights, function(mat_weights) {
        apply(mat_weights, 2, function(x) {
          is_inf <- is.infinite(x)
          if (any(is_inf)) {
            x[is_inf] <- fac_inf * sum(x[!is_inf])
          }
          return(x)
        })
      })
  }
  
  if (!is.null(list_mat_weights_loc)){
    list_mat_weights_loc <-
      lapply(list_mat_weights_loc, function(mat_weights_loc) {
        apply(mat_weights_loc, 2, function(x) {
          is_inf <- is.infinite(x)
          if (any(is_inf)) {
            x[is_inf] <- fac_inf * sum(x[!is_inf],na.rm = T)
          }
          return(x)
        })
      })
  }

  
  # handli8ng newdata
    if (inherits(new_data, "SpatRaster")) {
      newdata_exp_values <- terra::values(new_data)[,1]
      is_finite_newdata <- is.finite(newdata_exp_values)
      newdata_exp_values <- newdata_exp_values[is_finite_newdata]
    } else if (inherits(new_data, "sf")) {
      newdata_exp_values <- new_data[[exp_var]]
    } else {
      # numeric
      newdata_exp_values <- new_data
    }
  
  #### preparing strucutes to save data ####
  mat_estimates <- matrix(NA, ncol = n_dates,
                          nrow = length(newdata_exp_values))
  
  #### LOOP - CASE 0: no optimization lambda ####
  if (optim_lambda == 0) {
    if(!cross_validation){
      #TODO
    } else {
      #TODO
    }

    # no optimization, just use the first mat of the list
  }
  
  
  #### LOOP - CASE 1: global optimization lambda ####
  if (optim_lambda == 1) {
    if (!cross_validation){
      output_inter <- 
        calculate_inter_case1(formula = formula,
                              gc_after_number = gc_after_a_number,
                              locations = locations,
                              data = data,
                              new_data = newdata_exp_values,
                              i_target_stas = i_target_stas,
                              i_stas_lm = i_stas_lm,
                              fun_min = hydroGOF::rmse,
                              list_mat_weights_loc = list_mat_weights_loc,
                              list_mat_weights = list_mat_weights,
                              lambdas = lambdas,
                              min_0= min_0,
                              out_process_data = out_process_data)
      
      if(inherits(new_data, "SpatRaster")){
        newdata_estimates <- terra::rast(new_data,nlyrs = nrow(data))
        newdata_estimates[][is_finite_newdata,] <- 
          t(output_inter[["newdata_estimates"]])
        output_inter[["newdata_estimates"]] <- newdata_estimates
        rm(newdata_estimates)
        
        if(out_process_data){
          newdata_estimates_lm <- terra::rast(new_data,nlyrs = nrow(data))
          newdata_residuals <- newdata_estimates_lm
          newdata_estimates_lm[][is_finite_newdata,]  <- 
            t(output_inter[["newdata_estimates_lm"]])
          output_inter[["newdata_estimates_lm"]] <- newdata_estimates_lm
          rm(newdata_estimates_lm)
          newdata_residuals[][is_finite_newdata,]  <- 
            t(output_inter[["newdata_residuals"]])
          output_inter[["newdata_residuals"]] <- newdata_residuals
          rm(newdata_residuals)
        }
      }
      
      return(output_inter)
      
      # if(is.null(list_mat_weights) || is.null(new_data)){
      #   return(mat_estimates)
      # }
    } else{
      mat_estimates <- 
        calculate_inter_case1_CV(
          formula = formula,
          gc_after_number = gc_after_a_number,
          locations = locations,
          data = data,
          # new_data = newdata_exp_values,
          i_target_stas = i_target_stas,
          i_stas_lm = i_stas_lm,
          fun_min = hydroGOF::rmse,
          min_0 = min_0,
          list_mat_weights_loc = list_mat_weights_loc,
          list_mat_weights = list_mat_weights,
          lambdas = lambdas)
      return(mat_estimates)
    }
    
    
    
  }
  
  
  
  #### LOOP - CASE 2: optimization for every time step ####
  #TODO
  
  message()
  
  #building SpatRaster
  if (inherits(new_data, "SpatRaster")) {
    field_estimates <- terra::rast(new_data, nly = n_dates)
    values(field_estimates) <- mat_estimates
    rm(mat_estimates)
    if (out_process_data) {
      field_estimates_lm <- terra::rast(new_data, nly = n_dates)
      values(field_estimates_lm) <- mat_estimates_lm
      field_residuals <- terra::rast(new_data, nly = n_dates)
      values(field_residuals) <- mat_residuals
    }
  } else {
    # sf or numeric
    field_estimates <- t(mat_estimates)
    rm(mat_estimates)
    if (out_process_data) {
      field_estimates_lm <- mat_estimates_lm
      field_residuals <- mat_residuals
    }
  }
  
  if (out_process_data) {
    return (
      list(
        field_estimates = field_estimates,
        lambdas_sel = vec_lambdas_sel,
        list_df_ws = list_df_WS,
        list_lm = list_lm,
        field_estimates_lm = field_estimates_lm,
        field_residuals = field_residuals,
        delta_residuals = delta_residuals,
        rel_delta_residuals = delta_residuals
      )
    )
  } else {
    return(field_estimates)
  }
  
}


calculate_inter_case0 <- function(formula,
                                  gc_after_number,
                                  locations,
                                  data,
                                  new_data,
                                  i_stas_lm = NULL,
                                  mat_weights,
                                  lambdas,
                                  min_0 = FALSE,
                                  pre_message = NULL,
                                  out_process_data = FALSE) {
  n_dates <- nrow(data)
  
  c_for <- as.character(formula)
  pred_var <- c_for[2]
  exp_var <- c_for[3]
  
  if (is.null(i_stas_lm)) {
    locations$is_sta_lm <- TRUE
  }
  else if (is.logical(i_stas_lm)) {
    locations$is_sta_lm <- i_stas_lm
  }
  else {
    #numeric
    locations$is_sta_lm <- rep(NA, nrow(locations))
    locations$is_sta_lm[i_stas_lm] <- TRUE
  }
  
  mat_sim_0 <- matrix(NA, nrow = n_dates,
                      ncol = length(new_data))
  if(out_process_data){
    mat_est_lm_0 <- mat_sim_0
    mat_est_residuals_0 <- mat_sim_0
    list_lm = list()
    mat_loc_residuals <- matrix(NA,nrow = n_dates, 
                                ncol = nrow(locations))
  }
  
  for (i in 1:n_dates) {
    if (i %% gc_after_number == 0) {
      gc()
    }
    message(pre_message,
            format(round(i / n_dates, 3) * 100, nsmall = 1),
            " %\r", appendLF = F)
    values_obs <- data[i, , drop = T]
    is_valid_obs <- !is.na(values_obs)
    sf_sel_locs <- locations[is_valid_obs, ]
    sf_sel_locs[[pred_var]] <- values_obs[is_valid_obs]
    
    output_core <-
      NED_interpolation__core(
        formula,
        locations = sf_sel_locs,
        new_data = new_data,
        mat_weights = mat_weights[is_valid_obs,,drop = FALSE],
        i_stas_lm = sf_sel_locs[["is_sta_lm"]],
        min_0 = min_0,
        out_process_data = out_process_data
      )
    
    if(!out_process_data){
      mat_sim_0[i,] <- output_core
    } else {
      mat_sim_0[i,] <- output_core[[1]]
      mat_est_lm_0[i,] <- output_core[[2]]
      mat_est_residuals_0[i,] <- output_core[[3]]
      list_lm[[i]] <- output_core[["lm"]]
      mat_loc_residuals[i,is_valid_obs] <- 
        output_core[["loc_residuals"]]
    }
    
  }
  
  if(!out_process_data){
    return(mat_sim_0)
  } else {
    return(
      list(newdata_estimates = mat_sim_0,
           newdata_estimates_lm = mat_est_lm_0,
           newdata_residuals = mat_est_residuals_0,
           lm = list_lm,
           loc_residuals = mat_loc_residuals)  
      )
  }
  
  
}


calculate_inter_case1 <- function(formula,
                                  gc_after_number,
                                  locations,
                                  data,
                                  new_data,
                                  i_target_stas = NULL,
                                  i_stas_lm = NULL,
                                  fun_min = hydroGOF::rmse,
                                  list_mat_weights_loc,
                                  list_mat_weights,
                                  lambdas, 
                                  min_0 = FALSE,
                                  pre_message = NULL,
                                  out_process_data = FALSE) {
  n_dates = nrow(data)
  
  c_for <- as.character(formula)
  pred_var <- c_for[2]
  exp_var <- c_for[3]
  
  if (is.null(i_target_stas)) {
    locations$is_target <- TRUE
  }
  else if (is.logical(i_target_stas)) {
    locations$is_target <- i_target_stas
  }
  else {
    #numeric
    locations$is_target <- rep(FALSE, nrow(locations))
    locations$is_target[i_target_stas] <- TRUE
  }
  
  
  if (is.null(i_stas_lm)) {
    locations$is_sta_lm <- TRUE
  }
  else if (is.logical(i_stas_lm)) {
    locations$is_sta_lm <- i_stas_lm
  }
  else {
    #numeric
    locations$is_sta_lm <- rep(NA, nrow(locations))
    locations$is_sta_lm[i_stas_lm] <- TRUE
  }
  
  # creating list to save data
  mat_sim_0 <- matrix(NA, nrow = n_dates,
                      ncol = ncol(data))
  list_sim_loc <- list()
  for (i in 1:length(list_mat_weights_loc)) {
    list_sim_loc[[i]] <- mat_sim_0
  }
  
  #### INTERPOLATION OVER LOCATIONS ####
  for (i in 1:n_dates) {
    # if( i %% gc_after_a_number == 0){
    #   gc()
    # }
    message("\r",pre_message,
            format(round(i / n_dates, 3) * 100, nsmall = 1),
            " %", appendLF = F)
    values_obs <- data[i, , drop = T]
    is_valid_obs <- !is.na(values_obs)
    sf_sel_locs <- locations[is_valid_obs, ]
    sf_sel_locs[[pred_var]] <- values_obs[is_valid_obs]
    
    
    
    for (j in which(sf_sel_locs$is_target)) {
      # is_valid_obs_j <- is_valid_obs
      # is_valid_obs_j[j] <- FALSE
      
      list_mat_weights_loc_sel <-
        lapply(list_mat_weights_loc, function(mat_weights_loc) {
          # mat_weights_loc[is_valid_obs,is_valid_obs]
          mat_weights_loc[is_valid_obs, is_valid_obs][-j,j,drop =FALSE]
        })
      
      # cat("j = ",j)
      list_mat_sim <-
        NED_interpolation__core(
          formula,
          locations = sf_sel_locs[-j, ],
          new_data = sf_sel_locs[[exp_var]][j],
          mat_weights = list_mat_weights_loc_sel,
          i_stas_lm = sf_sel_locs[["is_sta_lm"]],
          min_0 = min_0
        )
      
      for (i_lambda in 1:length(list_mat_sim)) {
        list_sim_loc[[i_lambda]][i, which(is_valid_obs)[j]] <-
          list_mat_sim[[i_lambda]]
      }
      
    }
  }
  
  
  #### LAMBDA SELECTION ####
  
  errors_lambdas <-
    sapply(1:length(list_sim_loc), function(i_lambda) {
      sapply(which(locations$is_target), function(i_target) {
        obs_serie <- data[, i_target]
        sim_serie <- list_sim_loc[[i_lambda]][, i_target]
        fun_min(sim_serie, obs_serie)
      }) |> mean()
    })
  
  df_errors <- data.frame(lambdas = lambdas, error = errors_lambdas)
  i_best_lambda <- which.min(errors_lambdas)
  best_lambda = lambdas[i_best_lambda]
  message("Best_lambda is: ", best_lambda)
  # possible output
  if (is.null(list_mat_weights)) {
    warning(
      "Partial output: There is  nor a list_mat_weights or list_mat_dists")
    return(list(
      lambda = best_lambda,
      df_errors = df_errors,
      mat_sim = list_sim_loc[[i_best_lambda]]
    ))
  }
  
  #### INTERPOLATION ####
  mat_weights <- list_mat_weights[[i_best_lambda]]
  output_inter <- 
  calculate_inter_case0(formula = formula,
                        gc_after_number = gc_after_number,
                        locations = locations,
                        data = data,
                        new_data = new_data,
                        i_stas_lm = i_stas_lm,
                        mat_weights = mat_weights,
                        lambdas = lambdas,
                        pre_message = pre_message,
                        min_0 = min_0,
                        out_process_data = out_process_data) 
  if(!out_process_data){
    list(
      lambda = best_lambda,
      newdata_estimates = output_inter
    )
  } else {
    return(
      append(output_inter, 
             list(lambda = best_lambda,
                  WS = df_errors)
             )
    )
  }

}



calculate_inter_case1_CV <- function(formula,
                                     gc_after_number,
                                     locations,
                                     data,
                                     # new_data,
                                     i_target_stas = NULL,
                                     i_stas_lm = NULL,
                                     min_0 = FALSE,                                     fun_min = hydroGOF::rmse,
                                     list_mat_weights_loc,
                                     list_mat_weights,
                                     lambdas){
  c_for <- as.character(formula)
  pred_var <- c_for[2]
  exp_var <- c_for[3]
  
  if (is.null(i_target_stas)){
    locations$is_target <- TRUE
  } else if (is.numeric(i_target_stas)){
    locations$is_target <- FALSE
    locations$is_target[i_target_stas] <- TRUE
  } else if (is.logical(i_target_stas)){
    locations$is_target <- i_target_stas
  }
  
  if (is.null(i_stas_lm)){
    locations$is_sta_lm <- TRUE
  } else if (is.numeric(i_stas_lm)){
    locations$is_sta_lm <- FALSE
    locations$is_sta_lm[i_stas_lm] <- TRUE
  } else if (is.logical(i_stas_lm)){
    locations$is_sta_lm <- i_stas_lm
  }
  
  mat_sim <- matrix(NA, nrow(data),ncol(data))
  best_lambdas_CV <- rep(NA,ncol(data))
    
  for(i_sta in which(i_target_stas)){
    # message(paste0(i_sta,"/",sum(i_target_stas),"\r"),appendLF = F)
    pre_mes <- paste0(i_sta,"/",sum(i_target_stas)," ")
    # print(pre_mes)
    
    list_mat_weights_loc_sel <- 
      lapply(list_mat_weights_loc, function(mat_weights_loc){
        mat_weights_loc[-i_sta,-i_sta,drop =FALSE]
      })
    list_mat_weights_sel <- 
      lapply(list_mat_weights_loc, function(mat_weights){
        mat_weights[-i_sta,i_sta, drop = FALSE]
      })
    
    list_out <- 
    calculate_inter_case1(formula = formula,
                          gc_after_number = gc_after_number,
                          locations = locations[-i_sta,],
                          data = data[,-i_sta],
                          new_data = locations[i_sta,][[exp_var]],
                          i_target_stas = locations$is_target[-i_sta], ## TRY TO CHANGE
                          i_stas_lm = locations$is_sta_lm[-i_sta],
                          fun_min = fun_min,
                          min_0 = min_0,
                          list_mat_weights_loc = list_mat_weights_loc_sel,
                          list_mat_weights = list_mat_weights_sel,
                          lambdas = lambdas,
                          pre_message = pre_mes)
    
    mat_sim[,i_sta] <- list_out$newdata_estimates[,1]
    best_lambdas_CV[i_sta] <- list_out$lambda
  }
  return(list(mat_sim = mat_sim, best_lambdas_CV = best_lambdas_CV))
}
