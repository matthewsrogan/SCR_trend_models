prep_panthera_data <- function(main_data_lists, site, max_density) {
  # Load required libraries
  require(tidyverse)
  require(secr)
  
  # Filter surveys for the specified site
  survey_names <- str_subset(names(main_data_lists), pattern = site)
  site_surveys <- main_data_lists[survey_names]
  nt <- length(site_surveys)
  
  if (nt == 0) {
    stop("No surveys found for site: ", site)
  }
  
  #### Extract data components ####
  cap_hists <- lapply(site_surveys, 
                      function(obj){
                        return(obj$captrs)
                      })
  
  trap_objects <- lapply(site_surveys, 
                         function(obj){
                           return(obj$trap_obj)
                         })
  
  trap_effort <- lapply(trap_objects,
                        function(trap){
                          return(usage(trap))
                        })
  
  masks <- lapply(site_surveys,
                  function(obj){
                    return(obj$msk_obj)
                  })
  
  #### Format traps ####
  # see how many traps there are in each session
  nj <- sapply(trap_objects, nrow)
  
  # get mean coordinates 
  trap_coords <- do.call(bind_rows, trap_objects) %>%
    as_tibble() %>%
    mutate(x_km = x/1000,
           y_km = y/1000)
  
  mean_x <- mean(trap_coords$x_km)
  mean_y <- mean(trap_coords$y_km)
  
  # create empty array of coordinates
  X <- array(0, dim = c(max(nj),
                        nt,
                        2))
  
  # populate array
  for(t in 1:nt){
    trap <- trap_objects[[t]]
    X[1:nrow(trap), t, 1] <- trap$x/1000 - mean_x
    X[1:nrow(trap), t, 2] <- trap$y/1000 - mean_y
  }
  
  #### Format effort ####
  K <- matrix(0, 
              nrow = max(nj),
              ncol = nt)
  
  for(t in 1:nt){
    effort <- as.integer(rowSums(usage(trap_objects[[t]])))
    K[1:length(effort), t] <- effort
  }
  
  space_bounds <- get_space_boundaries(masks)
  
  cntrd_bounds <- (space_bounds - c(mean_x, mean_x, mean_y, mean_y))
  
  xlim <- cntrd_bounds[1:2]
  ylim <- cntrd_bounds[3:4]
  
  A <- as.numeric((xlim[2] - xlim[1]) * (ylim[2] - ylim[1]))
  
  M <- ceiling(max_density * A/100)
  
  #### Format Obs ####
  Y <- array(0, dim = c(M,
                        max(nj),
                        nt))
  
  Z <- matrix(NA, nrow = M, ncol = nt)
  
  Sex2 <- matrix(NA, nrow = M, ncol = nt)
  
  for(t in 1:nt){
    trap_names <- rownames(trap_objects[[t]])
    
    cap <- cap_hists[[t]] %>%
      mutate(id = as.integer(as.factor(.data$ID)),
             trap_id = as.integer(factor(.data$StationID,
                                         levels = trap_names)),
             sex_id = as.integer(factor(.data$Sex,
                                        levels = c("Female", "Male")))) %>%
      count(id, trap_id, sex_id)
    
    for(r in 1:nrow(cap)){
      Y[cap$id[r], cap$trap_id[r], t] <- cap$n[r]
    }
    
    sex_df <- cap %>%
      distinct(id, sex_id)
    
    for(q in 1:nrow(sex_df)){
      Sex2[sex_df$id[q], t] <- sex_df$sex_id[q]
      Z[sex_df$id[q], t] <- 1
    }
  }
  
  # convert categorical sex to binary sex
  Sex <- Sex2 - 1
  
  # calculate maximum individuals in a session
  maxN <- max(colSums(Z, na.rm = T))
  
  data_list <- list(
    Y = Y,
    Z = Z,
    X = X,
    K = K,
    Sex = Sex
  )
  
  constants <- list(
    M = M,
    nt = nt,
    nj = nj,
    A = round(A, 1),
    xlim = xlim,
    ylim = ylim
  )
  
  # Return both data and constants as a list
  return(list(
    data = data_list,
    constants = constants
  ))
}

#### Get state space parameters helper function ####
get_space_boundaries <- function(masks){
  
  minX <- sapply(masks, function(msk)return(min(msk$x/1000, na.rm = T)))
  maxX <- sapply(masks, function(msk)return(max(msk$x/1000, na.rm = T)))
  minY <- sapply(masks, function(msk)return(min(msk$y/1000, na.rm = T)))
  maxY <- sapply(masks, function(msk)return(max(msk$y/1000, na.rm = T)))
  
  return(c(minX = min(minX),
           maxX = max(maxX),
           minY = min(minY),
           maxY = max(maxY)))
}

# Example usage:
# load("../Data/Rogan_data_2025.Rdata")  # This loads secr_cptr_fls
# inputs <- prep_panthera_data(main_data_lists = secr_cptr_fls, 
#                              site = "Mumbwa", 
#                              max_density = 12)
# data_list <- result$data
# constants <- result$constants