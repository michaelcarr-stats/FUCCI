GenerateSummaryStatistics <- function(data, CellTracking = TRUE, Xmax){
  
  if (CellTracking){
    
    # extract variables from list
    red_distance <- data$red_distance
    yellow_distance <- data$yellow_distance
    green_distance <- data$green_distance
    red_time <- data$red_time
    yellow_time <- data$yellow_time
    green_time <- data$green_time
    
    #remove observations where phases werent reached
    #red_distance[red_time == 0] <- NA - by default red is always recorded
    yellow_distance[yellow_time <= 0] <- NA 
    green_distance[green_time <= 0] <- NA
    
    red_distance_avg <- mean(red_distance, na.rm = T) #compute average distance traveled
    yellow_distance_avg <- mean(yellow_distance, na.rm = T) #compute average distance traveled
    green_distance_avg <- mean(green_distance, na.rm = T) #compute average distance traveled
  
    return(c(mean(data$Nred), mean(data$Nyellow), mean(data$Ngreen), red_distance_avg, yellow_distance_avg, green_distance_avg)) 
    
  } else { #spatio temporal data
    
    #extract variables from list
    red_position <- data$red_position
    yellow_position <- data$yellow_position
    green_position <- data$green_position
    
    #calculate median position on left and right sides
    red_distance_med <- apply(matrix(unlist(lapply(red_position, FUN = function(x){c(median(x[x <= Xmax/2], na.rm = T),median(x[x > Xmax/2], na.rm = T))})), ncol = 2, byrow = T), MARGIN = 2, FUN = function(x){mean(x, na.rm = T)})
    yellow_distance_med <- apply(matrix(unlist(lapply(yellow_position, FUN = function(x){c(median(x[x <= Xmax/2], na.rm = T),median(x[x > Xmax/2], na.rm = T))})), ncol = 2, byrow = T), MARGIN = 2, FUN = function(x){mean(x, na.rm = T)})
    green_distance_med <- apply(matrix(unlist(lapply(green_position, FUN = function(x){c(median(x[x <= Xmax/2], na.rm = T),median(x[x > Xmax/2], na.rm = T))})), ncol = 2, byrow = T), MARGIN = 2, FUN = function(x){mean(x, na.rm = T)})
    
    #calculate interquartile range on left and right sides
    red_distance_iqr <- apply(matrix(unlist(lapply(red_position, FUN = function(x){c(IQR(x[x <= Xmax/2], na.rm = T),IQR(x[x > Xmax/2], na.rm = T))})), ncol = 2, byrow = T), MARGIN = 2, FUN = function(x){mean(x, na.rm = T)})
    yellow_distance_iqr <- apply(matrix(unlist(lapply(yellow_position, FUN = function(x){c(IQR(x[x <= Xmax/2], na.rm = T),IQR(x[x > Xmax/2], na.rm = T))})), ncol = 2, byrow = T), MARGIN = 2, FUN = function(x){mean(x, na.rm = T)})
    green_distance_iqr <- apply(matrix(unlist(lapply(green_position, FUN = function(x){c(IQR(x[x <= Xmax/2], na.rm = T),IQR(x[x > Xmax/2], na.rm = T))})), ncol = 2, byrow = T), MARGIN = 2, FUN = function(x){mean(x, na.rm = T)})
    
    return(c(mean(data$Nred), mean(data$Nyellow), mean(data$Ngreen), red_distance_med, yellow_distance_med, green_distance_med, red_distance_iqr, yellow_distance_iqr, green_distance_iqr))
  }
}
