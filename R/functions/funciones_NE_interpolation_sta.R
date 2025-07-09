
NE_interpolation_sta_NA_cal_par3 <- function(sta_coor, sta_z, data , var = 'pr' ,ista=NULL,
                                             new=T,niter=3,rasss = NULL,Rmax = 50, par0 = NULL,
                                             cal = c(T,T,T),p_s = c(1,1.5,2,2.5),func = NSE,ret_pot = F){
  spP <- SpatialPoints(sta_coor)
  if (is.na(crs(spP))) crs(spP) <- '+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0' 
  dists_0 <- spDists(spP)
  
  if (is.null(rasss)){
    rasss <- list()
    rasss[[1]] <- dists_0
  }
  
  a <- matrix(rep(sta_z,length(sta_z)),length(sta_z))
  
  if (is.null(ista))   ista <- 1:length(sta_z)
  #sim <- vector(mode='numeric',length = nrow(data))
  # if (var == 'pr'  & new)   del_z <- (t(a)/a)    # recien borrado
  # else  del_z <- (t(a)-a)                        # recien borrado
  
  #del_z_t <- del_z_t[,ista]
  # del_z <- del_z[,ista,drop=F]
  rasss <- lapply(rasss, function(dists) dists[,ista,drop=F])
  dists_0 <- dists_0[,ista,drop=F]
  
  data0 <- as.matrix(data)
  ista_d <- ista
  ista_d <- apply(dists_0 < 40,1,any)
  
  calcula <- function(w_s){   #,funcis = list(RMSE)
    sim <- sapply(1:nrow(data0),function(i){
      # modelo <- lm(p~z, data.frame(p = data0[i,],z = sta_z),subset = which(ista_d))
      modelo <- lm(p~z, data.frame(p = data0[i,ista_d],z = sta_z[ista_d]))
      residuales  <- rep(NA,length(sta_z))
      # residuales[!is.na(data0[i,]) & ista_d] <- modelo$residuals
      
      residuales <- data0[i,] - predict(modelo,newdata = data.frame(z = sta_z)) 
      #   colSums(w*residuales,na.rm=T)/colSums(w)
      # predis
      rmses <- sapply(w_s, function(w){
        predis <- predict(modelo,newdata = data.frame(z = sta_z[ista])) +
          colSums(w*residuales,na.rm=T)/colSums(w[!is.na(residuales),,drop=F])     # [!is.na(residuales),]   ,added 20/12/20
        rmse(predis,data0[i,ista])
      })
      w <- w_s[[which.min(rmses)]]
      predis <- predict(modelo,newdata = data.frame(z = sta_z[ista])) +
        colSums(w*residuales,na.rm=T)/colSums(w[!is.na(residuales),,drop=F])     # [!is.na(residuales),]   ,added 20/12/20
      # rmse(predis,data0[i,ista])
    }) %>% t
    # fun <- rmse #(function (x,y) cor(x,y,use = 'na.or.complete'))
    # sapply(funcis,function(funci) sapply(1:length(ista),function(i) funci(sim[,i],data0[,ista][,i])) %>% mean)
    # sapply(1:length(ista),function(i) func(sim[,i],data0[,ista][,i])) %>% mean
  }

  opti <- sapply(p_s,function(pot){
    w_s <- lapply(1:length(rasss),function(j){
      dists <- rasss[[j]]
      # dists <- dists[,ista,drop=F]
      w <- 1/dists^pot
      w[!is.finite(w)] <- 0
      w[dists_0 > Rmax] <- 0 
      w
    })
    sim <- calcula(w_s)
    sim[sim<0] <- 0
    sapply(1:length(ista),function(i) rmse(sim[,i],data0[,ista][,i])) %>% mean
  })
  
  print(p_s)
  print(opti)
  print(paste0('p = ',p_s[which.min(opti)]))
  print(paste0('opti = ',min(opti)))
  
  pot <- p_s[which.min(opti)]
  w_s <- lapply(1:length(rasss),function(j){
    dists <- rasss[[j]]
    # dists <- dists[,ista,drop=F]
    w <- 1/dists^pot
    w[!is.finite(w)] <- 0
    w[dists_0 > Rmax] <- 0 
    w
  })
  sim <- calcula(w_s)
  sim[sim<0] <- 0
  cat('NSE = ',
      sapply(1:length(ista),function(i) NSE(sim[,i],data0[,ista][,i])) %>% mean ,'\n')
  cat('RMSE = ',
      sapply(1:length(ista),function(i) rmse(sim[,i],data0[,ista][,i])) %>% mean ,'\n') 
  cat('r = ',
      sapply(1:length(ista),function(i) cor(sim[,i],data0[,ista][,i],use = 'na.or.complete')) %>% mean ,'\n')
  
    #cbind(predis,lin = predict(modelo,newdata = data.frame(z = sta_z[ista])),ori = data0[i,ista])
  # dists <- rasss[[80]]
  # dists <- dists[,ista,drop=F]
  #   w <- 1/dists^2
  #   w[!is.finite(w)] <- 0
  #   w[dists_0 > Rmax] <- 0 
  # 
  # ista_dis <- apply(dists_0 < 30,1,any)
  if (ret_pot) p_s[which.min(opti)] else NULL
  # par[c(2)]
}





NE_interpolation_sta_NA_cal_par3_V2 <- function(sta_coor, sta_z, data , var = 'pr' ,ista=NULL,
                                             new=T,niter=3,rasss = NULL,Rmax = 50, par0 = NULL,
                                             cal = c(T,T,T),p_s = c(1,1.5,2,2.5),func = NSE,return_sim = F){
  spP <- SpatialPoints(sta_coor)
  if (is.na(crs(spP))) crs(spP) <- '+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0' 
  dists_0 <- spDists(spP)
  
  if (is.null(rasss)){
    rasss <- list()
    rasss[[1]] <- dists_0
  }
  
  a <- matrix(rep(sta_z,length(sta_z)),length(sta_z))
  
  if (is.null(ista))   ista <- 1:length(sta_z)
  #sim <- vector(mode='numeric',length = nrow(data))
  # if (var == 'pr'  & new)   del_z <- (t(a)/a)    # recien borrado
  # else  del_z <- (t(a)-a)                        # recien borrado
  
  #del_z_t <- del_z_t[,ista]
  # del_z <- del_z[,ista,drop=F]
  rasss <- lapply(rasss, function(dists) dists[,ista,drop=F])
  dists_0 <- dists_0[,ista,drop=F]
  
  data0 <- as.matrix(data)
  if(is.logical(ista)) ista <- which(ista)
  ista_d_0 <- ista
  ista_d_0 <- apply(dists_0 < 40,1,any)               #### CAMBIO
  
  calcula <- function(w_s){   #,funcis = list(RMSE)
    sim <- sapply(1:nrow(data0),function(i){
      # el if es para verificar que al menos haya una estacion con datos
      if(data0[i,ista]  %>%( Negate(is.na)) %>% any()){
        
      predis <- sapply(ista,function(i_ista){
        
        ## NORMAL ##
        ista_d <- ista_d_0
        ista_d[i_ista] <- F
        ####
        distas <- dists_0[,i_ista==ista]   # las distancias de la estacion de analisis
        distas[is.na(data0[i,])] <- 100000   # los q son NA, dist grande
        distas[i_ista] <- 100000
        # ista_d <- distas %>% rank()  %>% (function(x) which(x <= 22))   # posiblemente argymento
        
        
        # modelo <- lm(p~z, data.frame(p = data0[i,],z = sta_z))
        modelo <- lm(p~z, data.frame(p = data0[i,ista_d],z = sta_z[ista_d]))
        residuales  <- rep(NA,length(sta_z))
        # residuales[!is.na(data0[i,]) & ista_d] <- modelo$residuals
        residuales <- data0[i,] - predict(modelo,newdata = data.frame(z = sta_z)) 
        
        ## DIBUJO
        # plot(data.frame(p = data0[i,ista_d],z = sta_z[ista_d]))
        # lines(predict(modelo,newdata = data.frame(z = sta_z)),sta_z)
        rmses <- sapply(w_s, function(w){
          predis <- predict(modelo,newdata = data.frame(z = sta_z[ista])) +
            colSums(w*residuales,na.rm=T)/colSums(w[!is.na(residuales),,drop=F])     # [!is.na(residuales),]   ,added 20/12/20
          rmse(predis,data0[i,ista])
        })
        # print(paste(' i_ista = ',i_ista))    # NECESARIO PUEDE SER
        # print(which.min(rmses))              # NECESARIO PUEDE SER
        # print(dim(w_s[[which.min(rmses)]]))
        w <- w_s[[which.min(rmses)]][,i_ista == ista,drop=F]
        
        predi <- predict(modelo,newdata = data.frame(z = sta_z[i_ista])) + 
          colSums(w*residuales,na.rm=T)/colSums(w[!is.na(residuales),,drop=F])      # [!is.na(residuales),]   ,added 20/12/20
      })
      # predis <- predict(modelo,newdata = data.frame(z = sta_z[ista])) +
      #   colSums(w*residuales,na.rm=T)/colSums(w)
      } else {
        predis <- rep(NA,length(ista))
      }
      # rmse(predis,data0[i,ista])
    }) %>% t
    # fun <- rmse #(function (x,y) cor(x,y,use = 'na.or.complete'))
    # sapply(funcis,function(funci) sapply(1:length(ista),function(i) funci(sim[,i],data0[,ista][,i])) %>% mean)
    # sapply(1:length(ista),function(i) func(sim[,i],data0[,ista][,i])) %>% mean
  }
  
  # pb <- txtProgressBar(min = 1, max = length(p_s), style = 3)
  opti <- sapply(p_s,function(pot){
    # setTxtProgressBar(pb, which(pot == p_s))
    
    w_s <- lapply(1:length(rasss),function(j){
      dists <- rasss[[j]]
      # dists <- dists[,ista,drop=F]
      w <- 1/dists^pot
      w[!is.finite(w)] <- 0
      w[dists_0 > Rmax] <- 0 
      w
    })
    sim <- calcula(w_s)
    sim[sim<0] <- 0
    sapply(1:length(ista),function(i) rmse(sim[,i],data0[,ista][,i])) %>% mean
  })
  
  # print(p_s)
  # print(opti)
  print(paste0('p = ',p_s[which.min(opti)]))
  # print(paste0('opti = ',min(opti)))
  
  pot <- p_s[which.min(opti)]
  w_s <- lapply(1:length(rasss),function(j){
    dists <- rasss[[j]]
    # dists <- dists[,ista,drop=F]
    w <- 1/dists^pot
    w[!is.finite(w)] <- 0
    w[dists_0 > Rmax] <- 0 
    w
    
    # exp #
    # i_border <- apply(dists_0 < 40,1,any) %>% which
    # w[i_border,] <- 0
    # w
    # exp #
  })
  
  ### ---- POSIBLE CAMBIO DE PESOS ----  ###
  # plot(1:131,scalar1(w_s[[1]][,1] ))
  # points(1:131,scalar1(w_s[[79]][,1] ),col='red')
  # 
  # plot(sta_z,scalar1(w_s[[1]][,1] ))
  # points(sta_z,scalar1(w_s[[79]][,1] ),col='red')
  # 
  # plot(dists_0[,1],scalar1(w_s[[1]][,1] ))
  # points(dists_0[,1],scalar1(w_s[[79]][,1] ),col='red')
  
  ### ---------------------- ###
  
  
  sim <- calcula(w_s)
  sim[sim<0] <- 0
  cat('NSE = ',
      sapply(1:length(ista),function(i) NSE(sim[,i],data0[,ista][,i])) %>% mean ,'\n')
  cat('RMSE = ',
      sapply(1:length(ista),function(i) rmse(sim[,i],data0[,ista][,i])) %>% mean ,'\n') 
  cat('r = ',
      sapply(1:length(ista),function(i) cor(sim[,i],data0[,ista][,i],use = 'na.or.complete')) %>% mean ,'\n')
  

  #cbind(predis,lin = predict(modelo,newdata = data.frame(z = sta_z[ista])),ori = data0[i,ista])
  # dists <- rasss[[80]]
  # dists <- dists[,ista,drop=F]
  #   w <- 1/dists^2
  #   w[!is.finite(w)] <- 0
  #   w[dists_0 > Rmax] <- 0 
  # 
  # ista_dis <- apply(dists_0 < 30,1,any)
  
  # par[c(2)]  
  
  # ADDED 08/01/21
  if(return_sim) sim else NULL
}


# scalar1 <- function(x) {x / sqrt(sum(x^2))}
scalar1 <- function(x) {x /sum(x) }






NE_interpolation_sta_NA_cal_par3_valery <- function(sta_coor, sta_z, data , var = 'pr' ,ista=NULL,
                                                      new=T,niter=3,dists = NULL,Rmax = 50, par0 = NULL,
                                                      cal = c(T,T,T),func = NSE,maxi = T,rasss ,p_s){
  spP <- SpatialPoints(sta_coor)
  if (is.na(crs(spP))) crs(spP) <- '+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0' 
  dists_0 <- spDists(spP)
  
  if (is.null(dists)){
    dists <- dists_0
  }
  
  a <- matrix(rep(sta_z,length(sta_z)),length(sta_z))
  
  if (is.null(ista))   ista <- 1:length(sta_z)
  #sim <- vector(mode='numeric',length = nrow(data))
  if (var == 'pr'  & new)   del_z <- (t(a)/a)  
  else  del_z <- (t(a)-a)
  
  #del_z_t <- del_z_t[,ista]
  del_z <- del_z[,ista,drop=F]
  dists <- dists[,ista,drop=F]
  dists_0 <- dists_0[,ista,drop=F]
  
  data0 <- as.matrix(data)
  if (var=='pr'){
    calcula <- function(p,Θ,lambda,fun = func,i_eval = NULL,eval_out=T){    #CAMBIO
      if (!is.null(i_eval)){
        vals_eval <- data0[,i_eval]
        data0[,i_eval] <- NA        
      }
      
      dists_l <- sqrt(dists^2 + (lambda * del_z/1000)^2)     #NUEVO !!!!!!!!!!!
      w <- 1/dists_l^p
      w[!is.finite(w)] <- 0
      w[dists_0 > Rmax] <- 0                            #NUEVO !!!!!!
      #w <- apply(w,2,function(x) x/sum(x))
      
      if (new)   del_z_t <- ((t(a) - Θ )/(a - Θ ))[,ista]  #del_z^Θ     ###AQUI CAMBIO##
      else       del_z_t <- exp(del_z*Θ)
      if (new & Θ == 0)  del_z_t <- del_z^0
      valery_fila <- function(x)   colSums(x*del_z_t*w,na.rm = T)/colSums(w[!is.na(x),,drop=F])
      sim <- t(apply(data0,1,valery_fila))       
      if(length(ista)==1 | sum(ista)==1)  sim <- t(sim)   # esto es solo para el caso de 1sta a analizar
      # fun(sim,data[,ista,drop=F])    #NSEm(sim,data[ista])
      if(eval_out) sapply(1:ncol(dists_0),function(i) fun(sim[,i],data0[,ista][,i])) %>% mean(na.rm=T) %>% return
      else fun(sim[,ista == i_eval],vals_eval) %>% return
    }
    
    ### ------- CALCULA con NEs ------- ###
    calcula_ne <- function(p,Θ,lambda,fun = func,i_eval = NULL,eval_out=T,w_s){    #CAMBIO
      if (!is.null(i_eval)){
        vals_eval <- data0[,i_eval]
        data0[,i_eval] <- NA        
      }
      
      # dists_l <- sqrt(dists^2 + (lambda * del_z/1000)^2)     #NUEVO !!!!!!!!!!!
      # w <- 1/dists_l^p
      # w[!is.finite(w)] <- 0
      # w[dists_0 > Rmax] <- 0                            #NUEVO !!!!!!
      #w <- apply(w,2,function(x) x/sum(x))
      
      #if (new)   del_z_t <- ((t(a) - Θ )/(a - Θ ))[,ista]  #del_z^Θ     ###AQUI CAMBIO##
      # else       del_z_t <- exp(del_z*Θ)
      del_z_t <- exp(del_z*Θ)
      
      # if (new & Θ == 0)  del_z_t <- del_z^0
      if (Θ == 0)  del_z_t <- del_z^0
      resul<- lapply(w_s,function(w){
        valery_fila <- function(x)   colSums(x*del_z_t*w,na.rm = T)/colSums(w[!is.na(x),,drop=F])
        sim <- t(apply(data0,1,valery_fila))       
        if(length(ista)==1 | sum(ista)==1)  sim <- t(sim)
        # sim
        perf <- sapply(1:nrow(df_M),function(i) rmse(sim[i,],data0[,ista][i,]))    # tavlez poner fun
        list(sim,perf)
      })
      rmses <- resul %>% sapply(function(x) x[[2]])
      sims <- resul %>% lapply(function(x) x[[1]])
      rm(resul)
      if (maxi)  iw_sel <- rmses %>% apply(1,which.max) 
      else  iw_sel <- rmses %>% apply(1,which.min)
      
      sim <- lapply(1:length(iw_sel),function(j){
        iw <- iw_sel[j]
        sims[[iw]][j,]
      }) %>% Reduce(function(x,y) rbind(x,y),.)
      names(sim) <- NULL
      # valery_fila <- function(x)   colSums(x*del_z_t*w,na.rm = T)/colSums(w[!is.na(x),,drop=F])
      # sim <- t(apply(data0,1,valery_fila))       
      # if(length(ista)==1 | sum(ista)==1)  sim <- t(sim)   # esto es solo para el caso de 1sta a analizar
      # # fun(sim,data[,ista,drop=F])    #NSEm(sim,data[ista])
      if(eval_out) sapply(1:ncol(dists_0),function(i) fun(sim[,i],data0[,ista][,i])) %>% mean(na.rm=T) %>% return
      else fun(sim[,ista == i_eval],vals_eval) %>% return
    }
  } else if (var=='temp'){
    calcula <- function(p,Θ,lambda,fun = func){   
      dists_l <- sqrt(dists^2 + (lambda * del_z/1000)^2)   #NUEVO !!!!!!!!!
      w <- 1/dists_l^p
      w[!is.finite(w)] <- 0
      w[dists > Rmax] <- 0                            #NUEVO !!!!!!
      
      #w <- apply(w,2,function(x) x/sum(x))
      del_z_t <- del_z*Θ 
      valery_fila <- function(x)   colSums((x+del_z_t)*w,na.rm = T)/colSums(w[!is.na(x),])
      sim <- t(apply(data0,1,valery_fila))       
      # fun(sim,data[ista])  #NSEm(sim,data[ista])
      sapply(1:ncol(dists_0),function(i) fun(sim[,i],data0[,ista][,i])) %>% mean
    }
  }
  if (is.null(par0)){
    if (var=='pr'){
      if (new) par0 <- c(2,1.01,0)    else  par0 <- c(2,4e-04,0)  # par3: 100
    } else if(var=='temp')  par0 <- c(2,-0.0065,10)  
  }
  
  
  cat('idw         ',calcula(par0[1],0,0),'\n')   # solo con new (con exp)  y viejo
  cat('valery def  ',calcula(par0[1],par0[2],0),'\n')
  
  if (niter > 0) {
    calcula2 <- function(Θ,p,lambda,i_eval = NULL,eval_out = T) calcula(p,Θ,lambda,i_eval = i_eval,eval_out = eval_out)
    calcula3 <- function(lambda,Θ,p) calcula(p,Θ,lambda)
    par <-c(p =par0[1],Θ=par0[2],lambda=par0[3],NSE = calcula(par0[1],par0[2],par0[3]))
    if (var=='pr'){
      max_Θ <-ifelse(new,1500,1e-03)  
      min_Θ <-ifelse(new,0,-0.5e-03)  
    }  
    else if (var =='temp')   max_Θ <- -0.012
    ##AQUI CAMBIO##   ifelse(new,3,1e-03) 
    
    # for (i in 1:niter){
    #   if (cal[1]==T) par[c(1,4)] <- as.numeric(optimize(calcula,c(0,4),Θ = par[2],lambda = par[3],
    #                                                     maximum = maxi,tol=0.0001))
    #   if (cal[2]==T) par[c(2,4)] <- as.numeric(optimize(calcula2,c(min_Θ,max_Θ),p = par[1],lambda = par[3],
    #                                                     maximum = maxi,tol=0.0001))
    #   if (cal[3]==T) par[c(3,4)] <- as.numeric(optimize(calcula3,c(0,1500),Θ = par[2],p = par[1],
    #                                                     maximum = maxi,tol=0.0001))
    # }
    
    
    
    potes <- c(1,1.5,2,2.5,3)
    
    
    #### IDW ####    
    
    vales <- sapply(potes,function(pot){
      calcula(pot,0,lambda = par[3],fun = func)    
    })
    if(maxi) i_pot <- which.max(vales_r)
    else i_pot <- which.min(vales)
    pot_idw <- potes[i_pot]
    cat('####   IDW   ####','\n')
    cat(paste0('p = ',pot_idw),'\n')
    ## indexes:
    nse <- calcula(pot_idw,0,lambda = par[3],fun = NSE)
    rmse <- min(vales)
    cor <- calcula(pot_idw,0,lambda = par[3],fun =  (function(x,y) cor(x,y,use = 'na.or.complete')))
    indexes <- c(NSE = nse,RMSE = rmse,cor = cor)
    print(indexes)
    
    
    
    
    #### VALERY ####    
    vales <- lapply(potes,function(pot){
      thetas_2 <- sapply(ista,function(i_val){
        as.numeric(optimize(calcula2,c(min_Θ,max_Θ),p = pot,lambda = par[3],i_eval = i_val,
                            maximum = maxi,tol=0.0001))[1]
      })
      vales_2 <- sapply(1:length(thetas_2),function(k){
        Θ <- thetas_2[k]
        i_eval <- ista[k]
        ## ---- PESOS ---- ##
        w_s <- lapply(1:length(rasss),function(j){
          dists <- rasss[[j]]
          dists <- dists[,ista,drop=F]
          
          w <- 1/dists^pot
          
          w[!is.finite(w)] <- 0
          w[dists_0 > Rmax] <- 0 
          # w <- w[,ista]
          w
          # w <- apply(w,2,function(x) x/sum(x))
          # w
        })
        ## ----       ---- ##
        calcula_ne(pot,Θ,lambda = par[3],fun = func,i_eval = i_eval, eval_out=F,w_s=w_s)
        # calcula(pot,Θ,lambda = par[3],fun = func,i_eval = i_eval, eval_out=F)
        
      })
      rbind(thetas_2,vales_2)
      # mean(vales_2)
    })
    vales_r <- sapply(vales,function(vec) mean(vec[2,])) 
    
    if(maxi) i_pot <- which.max(vales_r)
    else i_pot <- which.min(vales_r)
    
    pot <- potes[i_pot]
    cat('####  VALERY   ####','\n')
    cat(paste0('p = ',pot),'\n')
    # print(paste0('RMSE = ',vales_r[i_pot]))
    # par[1] <- potes[i_pot]
    # par[2] <- vales[1,i_pot]
    # par[4] <- vales[2,i_pot]
    
    # print(par)
    vales <- vales[[i_pot]]
    indexes <- sapply(1:ncol(vales),function(k){
      Θ <- vales[1,k]
      i_eval <- ista[k] 
      w_s <- lapply(1:length(rasss),function(j){
        dists <- rasss[[j]]
        dists <- dists[,ista,drop=F]
        
        w <- 1/dists^pot
        w[!is.finite(w)] <- 0
        w[dists_0 > Rmax] <- 0 
        
        # w <- w[,ista]
        w
        # w <- apply(w,2,function(x) x/sum(x))
        # w
      })
          
      nse = calcula_ne(pot,Θ,lambda = par[3],fun = NSE,i_eval = i_eval, eval_out=F,w_s = w_s)
      rmse = vales[2,k] %>% as.numeric
      cor = calcula_ne(pot,Θ,lambda = par[3],
                    fun = (function(x,y) cor(x,y,use = 'na.or.complete')),
                    i_eval = i_eval, eval_out=F,w_s = w_s)
      
      # nse = calcula(pot,Θ,lambda = par[3],fun = NSE,i_eval = i_eval, eval_out=F)
      # rmse = vales[2,k] %>% as.numeric
      # cor = calcula(pot,Θ,lambda = par[3],
      #                  fun = (function(x,y) cor(x,y,use = 'na.or.complete')),
      #                  i_eval = i_eval, eval_out=F)
      c(nse = nse,rmse = rmse,cor = cor)
    })
    rowMeans(indexes) %>% print()
    # rmses<- sapply(1:ncol(vales),function(k){
    #   Θ <- vales[1,k]
    #   i_eval <- ista[k]
    #   calcula(pot,Θ,lambda = par[3],fun = func,i_eval = i_eval, eval_out=F)
    # })
    # rs <- sapply(1:ncol(vales),function(k){
    #   Θ <- vales[1,k]
    #   i_eval <- ista[k]
    #   calcula(pot,Θ,lambda = par[3],
    #           fun = (function(x,y) cor(x,y,use = 'na.or.complete')),
    #           i_eval = i_eval, eval_out=F)
    # })
    
    #n0 <- c(par[1],par[2],100,calcula(par[1],par[2],10))
    #print(rbind(par,n0))
  }
  # if(niter == 0)  par <- par0
  # if(niter == -1)  par <- c(par0[1],0,0)
  # cat('NSE = ' ,mean(nses),'\n')
  # cat('RMSE = ' ,vales_r[i_pot],'\n')
  # cat('r = ' ,mean(rs),'\n')
  # par[c(2)]
}







