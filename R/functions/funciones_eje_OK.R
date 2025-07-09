autofitVariogram_v2 <- function (formula, input_data, model = c("Sph", "Exp", "Gau", 
  "Ste"), kappa = c(0.05, seq(0.2, 2, 0.1), 5, 10), fix.values = c(NA, 
  NA, NA), verbose = FALSE, GLS.model = NA, start_vals = c(NA, 
  NA, NA), miscFitOptions = list(),diagonal = NULL, fit.method = 7, ...) 
{
  if ("alpha" %in% names(list(...))) 
    warning("Anisotropic variogram model fitting not supported, see the documentation of autofitVariogram for more details.")
  miscFitOptionsDefaults = list(merge.small.bins = TRUE, min.np.bin = 5)
  miscFitOptions = modifyList(miscFitOptionsDefaults, miscFitOptions)
  longlat = !is.projected(input_data)
  if (is.na(longlat)) 
    longlat = FALSE
  if(is.null(diagonal)) diagonal = spDists(t(bbox(input_data)), longlat = longlat)[1, 
                                                             2]
  boundaries = c(2, 4, 6, 9, 12, 15, 25, 35, 50, 65, 80, 100) * 
    diagonal * 0.35/100
  if (!is(GLS.model, "variogramModel")) {
    experimental_variogram = variogram(formula, input_data, 
                                       boundaries = boundaries, ...)
  }
  else {
    if (verbose) 
      cat("Calculating GLS sample variogram\n")
    g = gstat(NULL, "bla", formula, input_data, model = GLS.model, 
              set = list(gls = 1))
    experimental_variogram = variogram(g, boundaries = boundaries, 
                                       ...)
  }
  if (miscFitOptions[["merge.small.bins"]]) {
    if (verbose) 
      cat("Checking if any bins have less than 5 points, merging bins when necessary...\n\n")
    while (TRUE) {
      if (length(experimental_variogram$np[experimental_variogram$np < 
                                           miscFitOptions[["min.np.bin"]]]) == 0 | length(boundaries) == 
          1) 
        break
      boundaries = boundaries[2:length(boundaries)]
      if (!is(GLS.model, "variogramModel")) {
        experimental_variogram = variogram(formula, 
                                           input_data, boundaries = boundaries, ...)
      }
      else {
        experimental_variogram = variogram(g, boundaries = boundaries, 
                                           ...)
      }
    }
  }
  if (is.na(start_vals[1])) {
    initial_nugget = min(experimental_variogram$gamma)
  }
  else {
    initial_nugget = start_vals[1]
  }
  if (is.na(start_vals[2])) {
    initial_range = 0.1 * diagonal
  }
  else {
    initial_range = start_vals[2]
  }
  if (is.na(start_vals[3])) {
    initial_sill = mean(c(max(experimental_variogram$gamma), 
                          median(experimental_variogram$gamma)))
  }
  else {
    initial_sill = start_vals[3]
  }
  if (!is.na(fix.values[1])) {
    fit_nugget = FALSE
    initial_nugget = fix.values[1]
  }
  else fit_nugget = TRUE
  if (!is.na(fix.values[2])) {
    fit_range = FALSE
    initial_range = fix.values[2]
  }
  else fit_range = TRUE
  if (!is.na(fix.values[3])) {
    fit_sill = FALSE
    initial_sill = fix.values[3]
  }
  else fit_sill = TRUE
  getModel = function(psill, model, range, kappa, nugget, 
                      fit_range, fit_sill, fit_nugget, verbose) {
    if (verbose) 
      debug.level = 1
    else debug.level = 0
    if (model == "Pow") {
      warning("Using the power model is at your own risk, read the docs of autofitVariogram for more details.")
      if (is.na(start_vals[1])) 
        nugget = 0
      if (is.na(start_vals[2])) 
        range = 1
      if (is.na(start_vals[3])) 
        sill = 1
    }
    obj = try(fit.variogram(experimental_variogram, model = vgm(psill = psill, 
                                                                model = model, range = range, nugget = nugget, kappa = kappa), 
                            fit.ranges = c(fit_range), fit.sills = c(fit_nugget, 
                                                                     fit_sill), debug.level = 0, fit.method = fit.method), TRUE)
    if ("try-error" %in% class(obj)) {
      warning("An error has occured during variogram fitting. Used:\n", 
              "\tnugget:\t", nugget, "\n\tmodel:\t", model, 
              "\n\tpsill:\t", psill, "\n\trange:\t", range, 
              "\n\tkappa:\t", ifelse(kappa == 0, NA, kappa), 
              "\n  as initial guess. This particular variogram fit is not taken into account. \nGstat error:\n", 
              obj)
      return(NULL)
    }
    else return(obj)
  }
  test_models = model
  SSerr_list = c()
  vgm_list = list()
  counter = 1
  for (m in test_models) {
    if (m != "Mat" && m != "Ste") {
      model_fit = getModel(initial_sill - initial_nugget, 
                           m, initial_range, kappa = 0, initial_nugget, 
                           fit_range, fit_sill, fit_nugget, verbose = verbose)
      if (!is.null(model_fit)) {
        vgm_list[[counter]] = model_fit
        SSerr_list = c(SSerr_list, attr(model_fit, "SSErr"))
      }
      counter = counter + 1
    }
    else {
      for (k in kappa) {
        model_fit = getModel(initial_sill - initial_nugget, 
                             m, initial_range, k, initial_nugget, fit_range, 
                             fit_sill, fit_nugget, verbose = verbose)
        if (!is.null(model_fit)) {
          vgm_list[[counter]] = model_fit
          SSerr_list = c(SSerr_list, attr(model_fit, 
                                          "SSErr"))
        }
        counter = counter + 1
      }
    }
  }
  strange_entries = sapply(vgm_list, function(v) any(c(v$psill, 
                                                       v$range) < 0) | is.null(v))
  if (any(strange_entries)) {
    if (verbose) {
      print(vgm_list[strange_entries])
      cat("^^^ ABOVE MODELS WERE REMOVED ^^^\n\n")
    }
    warning("Some models where removed for being either NULL or having a negative sill/range/nugget, \n\tset verbose == TRUE for more information")
    SSerr_list = SSerr_list[!strange_entries]
    vgm_list = vgm_list[!strange_entries]
  }
  if (verbose) {
    cat("Selected:\n")
    print(vgm_list[[which.min(SSerr_list)]])
    cat("\nTested models, best first:\n")
    tested = data.frame(`Tested models` = sapply(vgm_list, 
                                                 function(x) as.character(x[2, 1])), kappa = sapply(vgm_list, 
                                                                                                    function(x) as.character(x[2, 4])), SSerror = SSerr_list)
    tested = tested[order(tested$SSerror), ]
    print(tested)
  }
  result = list(exp_var = experimental_variogram, var_model = vgm_list[[which.min(SSerr_list)]], 
                sserr = min(SSerr_list))
  class(result) = c("autofitVariogram", "list")
  return(result)
}




plotea_variograma <- function(obj_automap,main = NULL){
  afv <- obj_automap
  
  dgt <- function(x) if (x >= 10) 0 else if (x >= 1) 1 else 2
  
  mdl <- afv$var_model
  cls <- as.character(mdl[2, "model"])
  ngt <- sum(mdl[1, "psill"])
  sll <- sum(mdl[, "psill"])
  rng <- sum(mdl[, "range"])
  ## INGLES
  # lbl <- paste("Model:", cls,
  #              "\nNugget:", round(ngt, dgt(ngt)),
  #              "\nSill:", round(sll, dgt(sll)),
  #              "\nRange:", round(rng, dgt(rng)))
  ## ESPANOL
  lbl <- paste("Modelo:", cls,
               "\nPepita:", round(ngt, dgt(ngt)),
               "\nMeseta:", round(sll, dgt(sll)),
               "\nRango:", round(rng, dgt(rng)))
  
  if (cls %in% c("Mat", "Ste")) {
    kpp <- mdl[2, "kappa"]
    lbl <- paste(lbl, "\nKappa:", round(kpp, dgt(kpp)), "")
  }
  
  if (is.null(main))  main = "Experimental variogram and fitted variogram model"
  ## create plot
  xyplot(gamma ~ dist, data = afv$exp_var,
         # main = main, 
         xlab = "Distancia [km]", ylab = "Semivarianza",   # change english
         # xlab = NULL,ylab = NULL,
         panel = function(x, y, ...) {
           gstat::vgm.panel.xyplot(x, y, cex = 1.0, ...)   #cex = 1.2
           ltext(max(x), 0.2 * max(y), lbl, font = 2, cex = .6, adj = c(1, 0.4),   #cex = .9
                 col = "grey30")
         }, 
         # arguments required by gstat::vgm.panel.xyplot()
         labels = NULL, mode = "direct", model = mdl,
         par.settings=list(par.main.text=list(cex=0.8)), # 1.0
         # axes = F,layout = c(2, 4),
         abline=list(h=0, lty = 5),
         xlim = c(0,NA),
         ylim = c(0,NA),
         # scales=list(alternating=3),
         direction = c(afv$exp_var$dir.hor[1], afv$exp_var$dir.ver[1]))
}

plotea_var_panel <- function(list_obj_automap,main = NULL){
  vg_0 <- list_obj_automap
  
  varios <- lapply(1:length(vg_0),function(k){
    x <- vg_0[[k]]
    x$exp_var$id <- as.character(df_sta[names(vg_0),]$ESTACION[k] )#as.factor(k)    #paste0('var',k)
    x$exp_var
  })
  
  varios <- varios[1:8] %>% Reduce(function(x,y) rbind(x,y),.)
  
  
  xyplot(gamma ~ dist | id, data = varios,
         # main = main, 
         xlab = "Distancia [km]", ylab = "Semivarianza",   # change english
         #xlab = NULL,ylab = NULL,
         panel = function(x, y, ...,subscripts) {   # 
           nume <- table(varios$id)[1]
           indi <- (as.numeric(subscripts[1])-1)%/%nume + 1
           mdl <- vg_0[[indi]]$var_model
           
           ## create custom text annotation
           dgt <- function(x) if (x >= 10) 0 else if (x >= 1) 1 else 2
           
           # mdl <- afv$var_model
           cls <- as.character(mdl[2, "model"])
           ngt <- sum(mdl[1, "psill"])
           sll <- sum(mdl[, "psill"])
           rng <- sum(mdl[, "range"])
           
           ## INGLES
           # lbl <- paste("Model:", cls,
           #              "\nNugget:", round(ngt, dgt(ngt)),
           #              "\nSill:", round(sll, dgt(sll)),
           #              "\nRange:", round(rng, dgt(rng)))
           ## ESPANOL
           lbl <- paste("Modelo:", cls,
                        "\nPepita:", round(ngt, dgt(ngt)),
                        "\nMeseta:", round(sll, dgt(sll)),
                        "\nRango:", round(rng, dgt(rng)))
           
           if (cls %in% c("Mat", "Ste")) {
             kpp <- mdl[2, "kappa"]
             lbl <- paste(lbl, "\nKappa:", round(kpp, dgt(kpp)), "")
           }
           
           gstat::vgm.panel.xyplot(x, y, cex = 0.85,model = mdl,...)
           ltext(max(x), 0.2 * max(y), lbl, font = 2, cex = .6, adj = c(1, 0.2),#pos = c(0.8,0.8),#offset=0, #paste(lbl,subscripts[1],indi)
                 col = "grey30")
         }, 
         # arguments required by gstat::vgm.panel.xyplot()
         labels = NULL, mode = "direct", #model = mdl,
         par.settings=list(par.main.text=list(cex=1.0)),
         axes = F,layout = c(2, 4),
         abline=list(h=0, lty = 5),
         scales=list(alternating=1),
         par.strip.text=list(cex=0.8),
         xlim = c(0,NA),
         ylim = c(0,NA),
         direction = c(vg_0[[1]]$exp_var$dir.hor[1], vg_0[[1]]$exp_var$dir.ver[1])
         )
  
}
