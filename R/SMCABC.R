#set working directory to the directory of this R file
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#load packages
library(Rcpp)
library(RcppArmadillo)
library(readxl)
library(ggplot2)
library(dplyr)
library(tidyr)
library(foreach)
library(MASS)
library(pracma)
library(doParallel)
library(doRNG)
library(openxlsx)
library(stats)

#load functions
library(FUCCIGillespie)
source("SimulationSetup.R") 
source("TransformationFun.R") 
source("GenerateSummaryStatistics.R")

#setup parrallel computing
numCores <- 16
registerDoParallel(numCores)

set.seed(1234)
#########################################
#----------------- Data ----------------#
#########################################

CellTracking = TRUE #TRUE - if using cell trajectory data; FALSE - if using cell density data

#load data
InitPos <- read_excel("../Data/DataProcessing/FUCCI_processed.xlsx", sheet = "InitPos")
FinalPos <- read_excel("../Data/DataProcessing/FUCCI_processed.xlsx", sheet = "FinalPos")
CellTrackingData <- read_excel("../Data/DataProcessing/FUCCI_processed.xlsx", sheet = "CellTracking")
Ntrack <- length(unique(CellTrackingData$ntrack))

Xmax <- 1309.09 #image width
Ymax <- 1745.35 #image height

#summary stats
Nred <- sum(FinalPos$color == 1)
Nyellow <- sum(FinalPos$color == 2)
Ngreen <- sum(FinalPos$color == 3)

if (CellTracking){ # Cell tracking summary statistics

  red_distance <- 
    CellTrackingData %>% 
    group_by(color, ntrack) %>% 
    transmute(Distance = sqrt((x - lag(x))^2+(y - lag(y))^2)) %>% #Radial Distance
    summarise(Distance = sum(Distance, na.rm = T)) %>% 
    ungroup() %>% 
    filter(color == 1) %>% 
    summarise(Distance = sum(Distance)/Ntrack) %>% #average distance
    as.numeric()
  
  yellow_distance <- 
    CellTrackingData %>% 
    group_by(color, ntrack) %>% 
    transmute(Distance = sqrt((x - lag(x))^2+(y - lag(y))^2)) %>% #Radial Distance
    summarise(Distance = sum(Distance, na.rm = T)) %>% 
    ungroup() %>% 
    filter(color == 2) %>%
    summarise(Distance = sum(Distance)/Ntrack) %>% #average distance
    as.numeric()
  
  green_distance <- 
    CellTrackingData %>% 
    group_by(color, ntrack) %>% 
    transmute(Distance = sqrt((x - lag(x))^2+(y - lag(y))^2)) %>% #Radial Distance
    summarise(Distance = sum(Distance, na.rm = T)) %>% 
    ungroup() %>% 
    filter(color == 3) %>% 
    summarise(Distance = sum(Distance)/Ntrack) %>% #average distance
    as.numeric()
  
  sy <- c(Nred, Nyellow, Ngreen, red_distance, yellow_distance, green_distance)  #observed data                 
  
} else { # Cell density summary statistics
   
  red_distance_med <- 
    FinalPos %>% 
    filter(color == 1) %>% #select red cells
    mutate(LHS_x = ifelse(x <= Xmax/2, x, NA), RHS_x = ifelse(x > Xmax/2, x, NA)) %>% #split cells into RHS or LHS based off image width
    summarise(median(LHS_x, na.rm = T), median(RHS_x, na.rm = T)) %>% #compute median position of left and right populations
    as.numeric()
  
  yellow_distance_med <- 
    FinalPos %>% 
    filter(color == 2) %>% #select yellow cells
    mutate(LHS_x = ifelse(x <= Xmax/2, x, NA), RHS_x = ifelse(x > Xmax/2, x, NA)) %>% #split cells into RHS or LHS based off image width
    summarise(median(LHS_x, na.rm = T), median(RHS_x, na.rm = T)) %>% #compute median position of left and right populations
    as.numeric()
  
  green_distance_med <- 
    FinalPos %>% 
    filter(color == 3) %>% #select green cells
    mutate(LHS_x = ifelse(x <= Xmax/2, x, NA), RHS_x = ifelse(x > Xmax/2, x, NA)) %>% #split cells into RHS or LHS based off image width
    summarise(median(LHS_x, na.rm = T), median(RHS_x, na.rm = T)) %>% #compute median position of left and right populations
    as.numeric()
  
  red_distance_iqr <- 
    FinalPos %>% 
    filter(color == 1) %>% #select red cells
    mutate(LHS_x = ifelse(x <= Xmax/2, x, NA), RHS_x = ifelse(x > Xmax/2, x, NA)) %>% #split cells into RHS or LHS based off image width
    summarise(IQR(LHS_x, na.rm = T), IQR(RHS_x, na.rm = T)) %>% #compute interquartile range of left and right populations
    as.numeric()
  
  yellow_distance_iqr <- 
    FinalPos %>% 
    filter(color == 2) %>% #select yellow cells
    mutate(LHS_x = ifelse(x <= Xmax/2, x, NA), RHS_x = ifelse(x > Xmax/2, x, NA)) %>% #split cells into RHS or LHS based off image width
    summarise(IQR(LHS_x, na.rm = T), IQR(RHS_x, na.rm = T)) %>% #compute interquartile range of left and right populations
    as.numeric()
  
  green_distance_iqr <- 
    FinalPos %>% 
    filter(color == 3) %>% #select green cells
    mutate(LHS_x = ifelse(x <= Xmax/2, x, NA), RHS_x = ifelse(x > Xmax/2, x, NA)) %>% #split cells into RHS or LHS based off image width
    summarise(IQR(LHS_x, na.rm = T), IQR(RHS_x, na.rm = T)) %>% #compute interquartile range of left and right populations
    as.numeric()
  
  sy <- c(Nred, Nyellow, Ngreen, red_distance_med, yellow_distance_med, green_distance_med, red_distance_iqr, yellow_distance_iqr, green_distance_iqr) #observed data
  
}

#########################################
#----------------SMC-ABC----------------#
#########################################

## SMC-ABC Inputs
N = 1000; #number of samples from posterior distribution
mcmc_trials = 30; # initial number of MCMC iterations
num_params = 6; #number of parameters to be estimated

# tuning parameters
alpha = 0.5;
c = 0.01;
tol_target = 0;  # target tolerance
pacc_target = 0.01; # target MCMC acceptance ratio

#prior bounds
t_lb = 0; t_ub = 1; #cell cycle transition rate prior
m_lb = 0; m_ub = 10; #cell migration prior

# summary statistic parameters
T_record <- 48; #time summary statistics are recorded

SetupVars <- SimulationSetup(Ntrack, Xmax, Ymax, InitPos, CellTrackingData); # initialise simulation environment

## SMC MAIN
N_a <- ceiling(alpha*N) #number of particles to drop each iteration
#ABC-rejection sampling step 
Samples <- 
  foreach(i = 1:N, .combine = "rbind", .packages = c("Rcpp", "RcppArmadillo", "pracma", "MASS", "FUCCIGillespie")) %dorng% {
    while (TRUE){
      
      theta <- c(runif(3,t_lb,t_ub), runif(3,m_lb,m_ub)) #propose theta
      
      data <- Main_Simulate(theta,SetupVars,T_record, CellTracking) #simulate data
      
      if (data$ExitSimStatus == 1){ #check if simulation finished correctly
        next #reject proposal
      }
      
      sx <- GenerateSummaryStatistics(data, CellTracking, Xmax) #generate summary statistics
      
      if (any(is.na(sx))) { #check integrity of summary statistics
        next #reject proposal
      }
      
      d <- norm(sx - sy, "2") #compute distance between observed and simulated data
      break #accept proposal
    }
    c(theta,d,sx)
  }
colnames(Samples) <- c("Rr", "Ry", "Rg", "Mr", "My", "Mg","d", paste0("sx",seq_along(1:length(sy))))
rownames(Samples) <- NULL  
  
# SMC ABC section

p_accept <- 1; #initialise acceptance rate
tolmax <- max(Samples[,"d"]); #initialise maximum tolerance

Samples_trans <- cbind(TransformationFun(Samples[,1:3],1,FALSE),TransformationFun(Samples[,4:6],10,FALSE),Samples[,-c(1:num_params)]) #transform samples to be unbounded by prior

while (p_accept > pacc_target && tolmax > tol_target) {
  
  Samples_trans <- Samples_trans[order(Samples_trans[,"d"]),] #order samples by distance between observed and simulated data
  tol <- Samples_trans[N-N_a,"d"] #compute tolerance 
  
  rows <- sample(1:(N-N_a),N_a, replace = T) #resample from best N-N_a particles
  Samples_trans[(N-N_a+1):N,] <- Samples_trans[rows,] #replace worst N_a particles with resampled particles
  
  cov_rw <- cov(Samples_trans[1:(N-N_a),1:num_params]) # computing tunning parameter of proposal function
  n_accept <- rep(0, N_a) #initialise acceptance counter
  
  #trial mcmc step
  tmp <-
    foreach (j = (N - N_a + 1):N, .combine = "rbind", .packages = c("Rcpp","RcppArmadillo","pracma", "MASS", "FUCCIGillespie")) %dorng% {
      theta_trans_j <- Samples_trans[j,1:num_params]
      d <-  Samples_trans[j,num_params+1] 
      sx <- Samples_trans[j,(num_params+2):ncol(Samples_trans)] 
      for (k in 1:mcmc_trials) {
        
        theta_trans_prop <- mvrnorm(1, mu = theta_trans_j, Sigma = cov_rw) #propose new theta
        
        #early rejection
        MHr <- exp(sum(log(exp(-theta_trans_prop))) - sum(log((1+exp(-theta_trans_prop))^2)) - (sum(log(exp(-theta_trans_j))) - sum(log((1+exp(-theta_trans_j))^2))))
        if (runif(1) > MHr) {
          next #reject proposal
        }
        
        theta_prop <- c(TransformationFun(theta_trans_prop[1:3],1,TRUE),TransformationFun(theta_trans_prop[4:6],10,TRUE)) #compute inverse transformation  
        
        data <- Main_Simulate(theta_prop, SetupVars, T_record, CellTracking) #simulate data
        
        if (data$ExitSimStatus == 1) { #check if simulation finished correctly
          next
        }
        
        sx_prop <- GenerateSummaryStatistics(data, CellTracking, Xmax) #generate summary statistics
        
        if (any(is.na(sx_prop))) { #check integrity of summary statistics
          next
        }
        
        d_prop <- norm(sx_prop - sy, "2") #calculate euclidean distance between observed and simulated data
        
        if (d_prop < tol){ #accept proposal
          theta_trans_j <- theta_trans_prop 
          n_accept[j - (N-N_a)] <- n_accept[j - (N-N_a)] + 1
          d <- d_prop
          sx <- sx_prop
        }
      }
      c(n_accept[j - (N-N_a)],theta_trans_j,d,sx)
  }
  
  n_accept <- tmp[,1]
  Samples_trans <- rbind(Samples_trans[1:(N-N_a),], tmp[,2:ncol(tmp)])
  rm(tmp)
                  
  p_accept <- sum(n_accept)/(N_a*mcmc_trials) #calculate mcmc acceptance rate
  mcmc_iters <- ceiling(log(c)/log(1-p_accept)) #calculate additional mcmc moves required
  
  #remaining mcmc moves step
  tmp <-
    foreach (j = (N - N_a + 1):N, .combine = "rbind", .packages = c("Rcpp","RcppArmadillo","pracma", "FUCCIGillespie")) %dorng% {
      theta_trans_j <- Samples_trans[j,1:num_params] 
      d <-  Samples_trans[j,num_params+1] 
      sx <- Samples_trans[j,(num_params+2):ncol(Samples_trans)] 
      for (k in 1:(mcmc_iters - mcmc_trials)){
        
        theta_trans_prop <- mvrnorm(1, mu = theta_trans_j, Sigma = cov_rw) #propose new theta
        
        #early rejection
        MHr <- exp(sum(log(exp(-theta_trans_prop))) - sum(log((1+exp(-theta_trans_prop))^2)) - (sum(log(exp(-theta_trans_j))) - sum(log((1+exp(-theta_trans_j))^2))))
        if (runif(1) > MHr) {
          next #reject proposal
        }
        
        theta_prop <- c(TransformationFun(theta_trans_prop[1:3],1,TRUE),TransformationFun(theta_trans_prop[4:6],10,TRUE)) #compute inverse transformation  
        
        data <- Main_Simulate(theta_prop, SetupVars, T_record, CellTracking) #simulate data
        
        if (data$ExitSimStatus == 1) { #check if simulation finished correctly
          next
        }
        
        sx_prop <- GenerateSummaryStatistics(data, CellTracking, Xmax) #generate summary statistics
        
        if (any(is.na(sx_prop))) { #check integrity of summary statistics
          next
        }
        
        d_prop <- norm(sx_prop - sy, "2") #calculate euclidean distance between observed and simulated data
        
        if (d_prop < tol){ #accept proposal
          theta_trans_j <- theta_trans_prop 
          n_accept[j - (N-N_a)] <- n_accept[j - (N-N_a)] + 1
          d <- d_prop
          sx <- sx_prop
        }
      }
      c(n_accept[j - (N-N_a)],theta_trans_j,d,sx)
    }
  n_accept <- tmp[,1]
  Samples_trans <- rbind(Samples_trans[1:(N-N_a),], tmp[,2:ncol(tmp)])
  rm(tmp)
  
  num_mcmc_iters <- max(0,mcmc_iters-mcmc_trials) + mcmc_trials
  p_accept <- sum(n_accept)/(N_a*num_mcmc_iters) #calculate mcmc acceptance rate
  mcmc_iters <- ceiling(log(c)/log(1-p_accept)) #compute number of mcmc iters
  mcmc_trials <- ceiling(mcmc_iters/2) #compute number of pilot mcmc iters for next sequence
  tolmax <- max(Samples_trans[,"d"]) #compute maximum tolerance
  
  #display performance stats
  cat("\n")
  message(paste0("Acceptance Rate: ",round(p_accept,2)))
  message(paste0("Max tolerance: ",round(tolmax,2)))
  message(paste0("% unique: ",round(nrow(unique(Samples_trans))/N,2)))
  cat("\n")
  
  #write performance stats to excel
  write.xlsx(c(p_accept = p_accept, tolmax = tolmax, unique = nrow(unique(Samples_trans))/N), paste0("SMCABC_progress_", Ntrack, ".xlsx"), col.names = TRUE, row.names = FALSE)
}
Samples <- cbind(TransformationFun(Samples_trans[,1:3],1,TRUE),TransformationFun(Samples_trans[,4:6],10,TRUE),Samples_trans[,-c(1:num_params)]) 

# Regression Adjustment
theta_trans <- as_tibble(Samples_trans[,1:num_params])
d_vec <- Samples_trans[,"d"]
summaries <- as_tibble(Samples_trans[,(num_params + 2):ncol(Samples_trans)])
the_weights = as_tibble(0.75*(1-(d_vec/max(d_vec))^2)) #Epanechnikov weighting kernel

theta_trans_regressed <- matrix(NA, nrow = N, ncol = num_params)
for (i in 1:num_params){ 
  data <- bind_cols(theta_trans[,i],summaries,the_weights) %>% setNames(c("y", colnames(summaries), "the_weights"))
  b <- lm(y ~ . - the_weights, data = data, weights = the_weights) #compute coefficients
  theta_trans_regressed[,i] <- as.matrix(theta_trans[,i] - as.matrix(summaries - repmat(sy, N, 1)) %*%  t(t(b$coefficients[-1]))) #apply adjustment
}
theta_regressed <- cbind(TransformationFun(theta_trans_regressed[,1:3],1,TRUE),TransformationFun(theta_trans_regressed[,4:6],10,TRUE)) 
colnames(theta_regressed) <- colnames(Samples)[1:num_params]

#save results
write.xlsx(list(Samples = Samples, theta_regressed = theta_regressed), paste0("SMCABC_Samples_", Ntrack, ".xlsx"), col.names = TRUE, row.names = FALSE)


