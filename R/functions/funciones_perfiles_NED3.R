

# output: a matrix where rows are lambdas and 
# columns represent different transect lines (from end_coords) #
get_NED3 <- function(srtm,begin_coord,end_coords = NULL,w_z = 1,ver = F){
  

  
  if (inherits(begin_coord,"sf")){
    begin_coord <- st_coordinates(begin_coord)
  }
  
  if (!is.data.frame(begin_coord)){
    df_origin <- data.frame(x = begin_coord[1],y = begin_coord[2])
  } else {
    df_origin <- setNames(begin_coord, c("x","y"))
  }
  sf_origin <- st_as_sf(df_origin, coords = c("x","y"),crs = 4326)
  
  if(!is.finite(terra::extract(srtm,sf_origin)[,2])){
    stop("begin_coord is not in the DEM area with finite values")
  }
  
  if (inherits(end_coords,"sf")){
    end_coords <- st_coordinates(end_coords)
  }
  control_na <- FALSE
  if (is.null(end_coords)){
    control_na <- TRUE
    values_ras_dem <- terra::values(srtm)[,1]
    is_na_cells <- is.na(values_ras_dem)
    is_notna_cells <- !is.na(values_ras_dem)
    end_coords <- terra::crds(srtm) # it only returns crds of not na cells by default
    rm(values_ras_dem)
 }
  
  vec_lines <- 
    lapply(1:nrow(end_coords), function(i){
      mat_ex1 <- cbind(geom = rep(i, 2), part = c(1,1))
      rbind(begin_coord, end_coords[i,]) %>% 
        cbind(mat_ex1,.)
    }) %>% do.call(rbind,.) %>% 
      terra::vect(crs = "epsg:4326",type = "lines")
  
  
  perfiles <- 
    get_profiles_lines_same_origin(srtm, 
                                   transects_lines = vec_lines, 
                                   sf_origin,
                                   all_finite = F,rem_non_finite = T)
  perfiles <- setNames(perfiles, c("ID","dist","elev"))
  
  perfiles <- perfiles %>%
    group_by(ID) %>% group_split()
  perfiles <- lapply(perfiles,function(x) x %>% dplyr::select(-'ID') %>% arrange(dist))
   
  if(!ver){
    lineas <- lapply(perfiles,function(perfil) trazo(perfil,ver = F))    
  } else{
    lineas <- lapply(1:length(perfiles),function(i){
      cat(i,' ')
      trazo(perfiles[[i]],ver = F)
    } )
  }
  # if(dib){
  #   plot(perfil,type = 'l')
  #   lines(linea,col = 'red')
  # }  
  mat_dists <- sapply(lineas,function(linea) dist_trazo(linea,w_z =w_z))/1000
  
  if(control_na){
    mat_dists_new <- matrix(NA, nrow = length(w_z), 
                            ncol = length(is_notna_cells))
    mat_dists_new[,is_notna_cells] <- mat_dists 
    mat_dists <- mat_dists_new
  }
  
  return(mat_dists)
}


trazo <- function(perfil,ver = T){
  linea <- line1(perfil)
  # points(linea,col = 'red')
  maxi <- max_r_elev(perfil,linea)
  dife_maxi <- (perfil$elev-linea$elev)[maxi]
  breaks <- c()
  
  while (dife_maxi > 0 & !(maxi %in% breaks)) {
    breaks <- c(breaks,maxi)
    linea <- line2(perfil ,breaks = breaks)
    maxi <- max_r_elev(perfil,linea)
    dife_maxi <- (perfil$elev-linea$elev)[maxi] %>% round(4)
    if(ver){
      print(breaks)
      cat('I')
    }
  }
  linea
  linea[c(1,sort(breaks),nrow(perfil)),]
}


line1 <-function(perfil,x){
  m <- ( perfil$elev[nrow(perfil)] -perfil$elev[1])/( perfil$dist[nrow(perfil)])
  if(is.nan(m)) m <- 0
  perfil$elev <- perfil$elev[1] + m*perfil$dist
  # perfil_1 <- data.frame(dist = perfil$dist, elev = perfil$elev[1] + m*perfil$dist)
  perfil
}

max_r_elev <- function(perfil,perfil2){
  which.max(perfil$elev-perfil2$elev)
} 

trazo <- function(perfil,ver = T){
  linea <- line1(perfil)
  # points(linea,col = 'red')
  maxi <- max_r_elev(perfil,linea)
  dife_maxi <- (perfil$elev-linea$elev)[maxi]
  breaks <- c()
  
  while (dife_maxi > 0 & !(maxi %in% breaks)) {
    breaks <- c(breaks,maxi)
    linea <- line2(perfil ,breaks = breaks)
    maxi <- max_r_elev(perfil,linea)
    dife_maxi <- (perfil$elev-linea$elev)[maxi] %>% round(4)
    if(ver){
      print(breaks)
      cat('I')
    }
  }
  linea
  linea[c(1,sort(breaks),nrow(perfil)),]
}


#### gives distances in m 
dist_trazo <- function(perfil,w_z = 1){
  perfil$d_2 <- c(perfil$dist[-1],NA)
  perfil$elev_2 <- c(perfil$elev[-1],NA)
  sapply(w_z, function(x){
    perfil$d_gen <- sqrt((perfil$d_2 - perfil$dist)^2 + (x*(perfil$elev_2 - perfil$elev))^2)
    sum(perfil$d_gen,na.rm =T)
  })
}


line2 <- function(perfil,breaks = NA){
  n_fin <- nrow(perfil) 
  # elev0 <- perfil$elev[1]
  perfil_1 <- perfil
  breaks <- c(1,sort(breaks),n_fin)
  
  new_elev_no_breaks <- 
    approx(perfil$dist[breaks], perfil$elev[breaks], perfil$dist[-breaks])$y
  # m <- sapply(1:(length(breaks)-1),function(j){
  #   br_1 <- breaks[j+1]
  #   br_0 <- breaks[j]
  #   (perfil$elev[br_1] - perfil$elev[br_0])/(perfil$dist[br_1] - perfil$dist[br_0])
  # })
  # gr <- findInterval(1:n_fin,breaks)
  # m <- m[gr]
  
  perfil_1$elev[-breaks] <- new_elev_no_breaks
  
  # for(k in 1:(n_fin-1)){
  #   perfil_1$elev[k+1] <- perfil_1$elev[k] + m[k]*(perfil_1$dist[k+1]-perfil_1$dist[k])
  # }
  perfil_1
}
