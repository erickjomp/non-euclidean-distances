# raster is a terra spatRaster
# transct_lines is sf or vect lines
# origin is sf

get_profiles_lines <- function(raster, transects_lines){
  require(terra)
  require(dplyr)
  require(sf)
  df_lines_values <- 
    terra::extract(ras_dem_reg, transects_lines,touches = T, xy = T)
  names(df_lines_values)[2] <- "value"
  df_lines_values <- 
    df_lines_values %>% group_by(ID) %>% 
    filter(all(is.finite(value)))
  
  geom_lines <- geom(transects_lines)
  begin_coords <- geom_lines[!duplicated(geom_lines[,1]),]
  
  #TODO
}

# rem_non_finite just remove intermediate non finite values in profile
get_profiles_lines_same_origin <- function(raster, transects_lines, 
                                           origin, all_finite = T,
                                           rem_non_finite = F,
                                           cells = F){
  if (inherits(transects_lines,"sf")){
    sf_transects_lines <- transects_lines
  } else {
    sf_transects_lines <- st_as_sf(transects_lines)
  }
  
  df_lines_values <- 
    terra::extract(raster, transects_lines,touches = T, xy = T, 
                   cells = cells)
  names(df_lines_values)[2] <- "value"

  if (all_finite){
    df_lines_values <- 
      df_lines_values %>% group_by(ID) %>% 
      filter(all(is.finite(value)))
  } 

  
  sf_cells_profiles <- 
    st_as_sf(df_lines_values[,c("x","y")],coords = c("x","y"),crs = st_crs(sf_transects_lines))
    
  sf_cells_profiles$dist_to_origin <- 
  as.numeric(
    st_distance(origin, sf_cells_profiles)
    )
  
  sf_cells_profiles$id_line <- df_lines_values$ID
  
  sf_transect_lines_rep <- sf_transects_lines[df_lines_values$ID,]

  sf_cells_profiles$dist_to_transect <- 
    as.numeric(
      st_distance(sf_cells_profiles, sf_transect_lines_rep,by_element = T)
    )
  
  sf_cells_profiles <- sf_cells_profiles %>% 
    mutate(dist_in_profile = sqrt(dist_to_origin ^ 2 - dist_to_transect ^ 2))

  
  df_lines_values$dist_from_origin <- sf_cells_profiles$dist_in_profile
  df_lines_values <- arrange(df_lines_values,
                             ID,
                             dist_from_origin)
  
  if(!all_finite & rem_non_finite){
    df_lines_values <- 
      df_lines_values %>% group_by(ID) %>% 
      filter(is.finite(value) | (row_number() %in% c(1,n()))) %>% 
      ungroup(ID)
  }
  
  df_lines_values <- 
    df_lines_values %>% ungroup(ID)
  
  if (cells){
    return(
      df_lines_values[c("ID", "cell", "dist_from_origin",
                        "value")]  #value is elevation
    )
  } else {
    return(
      df_lines_values[c("ID","dist_from_origin","value")]  #value is elevation
    )
  }
  
}


## it is only for plotting (it only works with one line)
get_profiles_lines_wsegments_same_origin <- 
  function(raster, transect_line, 
           origin = NULL, all_finite = T,
           rem_non_finite = F,
           cells = F){
    if (inherits(transect_line,"sf")){
      coords_line <- st_coordinates(transect_line)[,1:2]
      crs_line <- st_crs(transect_line)
      crs_line_char <- as.character(crs_line)[2]
    } else {
      coords_line <- terrra::crds(transect_line)
      crs_line <- terra::crs(transect_line)
      crs_line_char <- crs_line_char
    }
    colnames(coords_line) <- c("x", "y")
    
    sf_points <- 
      st_as_sf(as.data.frame(coords_line), 
               coords = c("x", "y"), crs = crs_line)
    
    dist_points <- 
      st_distance(sf_points[-nrow(sf_points),],
                  sf_points[-1,],by_element = T) %>% 
      as.numeric() %>% cumsum() %>% c(0,.)
    # dist_points <- dist_points 
    # print(dist_points)
    df_line_values <-
      lapply(1:(nrow(coords_line) - 1), function(i) {
        coords_pair <- coords_line[i:(i + 1), ]
        names(coords_pair) <- c("x", "y")
        mat_line <- cbind(geom = 1, part = 1, coords_pair)
        line <- vect(mat_line, "lines",
                     crs = crs_line_char)
        get_profiles_lines_same_origin(
          raster,
          line,
          origin = sf_points[i, ],
          all_finite = all_finite,
          rem_non_finite = rem_non_finite,
          cells = T
        ) %>%
          mutate(dist_from_origin =
                   dist_from_origin + dist_points[i]) #+ dist_points[i]
      }) %>%
      do.call(rbind, .) %>%
      filter(!duplicated(cell))
    if(!cells){
      df_line_values <- df_line_values[names(df_line_values) != "cell"]
    }
    return(df_line_values)
  }

get_profiles_lines_wsegments <- get_profiles_lines_wsegments_same_origin


#this funtion return the IDs (according to variable  ID of df_profiles) 
# of direct profiles (profiles where the straitght lines between the
# start and end doesnt cross any topographic obstacle)

# df_profiles is data frame
# with variables ID, distance_from_origin, value 
# new_names,  vecor of new names of variables

# argumetn:  dif - controls if horizontal and vertical differences 
#                   will be in output 

which_direct_profile <- function(df_profiles, new_names = NULL, 
                                 difs = TRUE){
  require(dplyr)
  # TODO what if new_names != NULL
  df_profiles %>% group_by(ID) %>% group_split() %>% 
    lapply(function(df_profile){

      n_points <- nrow(df_profile)
      dif_h <- df_profile$dist_from_origin[n_points]
      dif_v <- df_profile$value[n_points] -  df_profile$value[1]
      
      if(n_points > 2){

        slope <- dif_v / dif_h
        
        df_profile$value_straight <-  
          df_profile$value[1] + df_profile$dist_from_origin * slope
        
        result <- NULL
        vec_compare <- (df_profile$value <= df_profile$value_straight)
        if (all(vec_compare[-c(1,n_points)])){
          result <- data.frame(ID = df_profile$ID[1], 
                               dif_h = dif_h, dif_v)
        }
      } else {
        result <- 
          data.frame(ID = df_profile$ID[1], 
                     dif_h = dif_h, dif_v)
      }
      result
    }) %>% do.call(rbind,.)
  
}