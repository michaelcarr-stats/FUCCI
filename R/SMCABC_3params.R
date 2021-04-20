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
source("GenerateSummaryStatistics.R") 
source("TransformationFun.R") 

#setup parrallel computing
numCores = 16
registerDoParallel(numCores)

set.seed(1234)

#########################################
#----------------- Data ----------------#
#########################################

CellTracking <- TRUE #TRUE - if using cell trajectory data; FALSE - if using cell density data

#true value 
theta <- c(0.04,0.17,0.08,4,4,4)
Ntrack <- 20

Xmax <- 1309.09 #image width
Ymax <- 1745.35 #image height

SetupVars <- SimulationSetup(Ntrack, Xmax, Ymax, InitPos = NULL, CellTrackingData = NULL); # initialise simulation environment
data <- Main_Simulate(theta,SetupVars,T_record = 48, CellTracking) #simulate observed data
sy <- GenerateSummaryStatistics(data, CellTracking, Xmax)  #observed data   

#known parameter values
MotilityConstant <- 0 # 0 if transition parameters are held constant, 1 if it is motility
if (MotilityConstant){
  KnownTheta <- theta[4:6]
  sy <- sy[1:3]      
} else {
  KnownTheta <- theta[1:3]
  sy <- sy[-c(1:3)]      
}

#########################################
#----------------SMC-ABC----------------#
#########################################

## SMC-ABC Inputs

N = 1000; #number of samples from posterior distribution
mcmc_trials = 30; # initial R_t value
num_params = 3; #number of parameters to be estimated

# tuning parameters
alpha = 0.5;
c = 0.01;
tol_target = 0;  # target tolerance
pacc_target = 0.01; # target MCMC acceptance ratio

#prior bounds
lb = 0; ub = ifelse(MotilityConstant, 1, 10); 

# summary statistic parameters
T_record <- 48;


## SMC MAIN

N_a <- ceiling(alpha*N) #number of particles to drop each iteration

#ABC-rejection sampling step 
Samples <- 
  foreach(i = 1:N, .combine = "rbind", .packages = c("Rcpp", "RcppArmadillo", "pracma", "MASS", "FUCCIGillespie")) %dorng% {
    while (TRUE){
      
      theta <- runif(3,lb,ub) #propose theta
      
      if (MotilityConstant){
        data <- Main_Simulate(c(theta, KnownTheta),SetupVars,T_record, CellTracking) #simulate data
        
        if (data$ExitSimStatus == 1){ #check if simulation finished correctly
          next #reject proposal
        }
        
        sx <- GenerateSummaryStatistics(data, CellTracking, Xmax)[1:3] #Generate summary statistics
        
      } else {
        data <- Main_Simulate(c(KnownTheta, theta),SetupVars,T_record, CellTracking) #simulate data
        
        if (data$ExitSimStatus == 1){ #check if simulation finished correctly
          next #reject proposal
        }
        
        sx <- GenerateSummaryStatistics(data, CellTracking, Xmax)[-c(1:3)] #Generate summary statistics
      }
      
      if (any(is.na(sx))) { #check integrity of summary statistics
        next #reject proposal
      }
      
      d <- norm(sx - sy, "2") #compute distance between observed and simulated data
      break #accept proposal
    }
    c(theta,d,sx)
  }

if (MotilityConstant){
  colnames(Samples) <- c("Rr", "Ry", "Rg","d", paste0("sx",seq_along(1:length(sy))))
} else {
  colnames(Samples) <- c("Mr", "My", "Mg","d", paste0("sx",seq_along(1:length(sy))))
}
rownames(Samples) <- NULL  

# SMC ABC section

p_accept <- 1; #initialise acceptance rate
tolmax <- max(Samples[,"d"]); #initialise maximum tolerance

Samples_trans <- cbind(TransformationFun(Samples[,1:num_params],ub,FALSE),Samples[,-c(1:num_params)])

while (p_accept > pacc_target && tolmax > tol_target) {
  
  Samples_trans <- Samples_trans[order(Samples_trans[,"d"]),] #order samples by distance between observed and simulated data
  tol <- Samples_trans[N-N_a,"d"] #compute tolerance 
  
  rows <- sample(1:(N-N_a),N_a, replace = T) #resample from best N-N_a particles
  Samples_trans[(N-N_a+1):N,] <- Samples_trans[rows,] #replace worst N_a particles with resampled particles
  
  cov_rw <- cov(Samples_trans[1:(N-N_a),1:num_params]) # computing tunning parameter of proposal function
  n_accept <- rep(0, N - N_a) #initialise acceptance counter
  
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
        
        theta_prop <- TransformationFun(theta_trans_prop,ub,TRUE) #compute inverse transformation  
        
        if (MotilityConstant){
          data <- Main_Simulate(c(theta_prop, KnownTheta),SetupVars,T_record,CellTracking) #simulate data 
          
          if (data$ExitSimStatus == 1){ #check if simulation finished correctly
            next #reject proposal
          }
          
          sx_prop <- GenerateSummaryStatistics(data, CellTracking, Xmax)[1:3] #generate summary statistics
          
        } else {
          data <- Main_Simulate(c(KnownTheta, theta_prop),SetupVars,T_record,CellTracking) #simulate data 
          
          if (data$ExitSimStatus == 1){ #check if simulation finished correctly
            next #reject proposal
          }
          
          sx_prop <- GenerateSummaryStatistics(data, CellTracking, Xmax)[-c(1:3)] #generate summary statistics
        }
        
        if (any(is.na(sx_prop))) { check integrity of summary statistics
          next
        }
        
        d_prop <- norm(sx_prop - sy, "2") #calculate euclidean distance between observed and simulated data
        
        if (d_prop < tol){ #accept proposal
          theta_trans_j <- theta_trans_prop 
          n_accept[j - N_a] <- n_accept[j - N_a] + 1
          d <- d_prop
          sx <- sx_prop
        }
      }
      c(n_accept[j - N_a],theta_trans_j,d,sx)
    }
  
  n_accept <- tmp[,1]
  Samples_trans <- rbind(Samples_trans[1:(N-N_a),], tmp[,2:ncol(tmp)])
  rm(tmp)
  
  p_accept <- sum(n_accept)/(N_a*mcmc_trials) #calculate mcmc acceptance rate
  mcmc_iters <- ceiling(log(c)/log(1-p_accept)) #calculate additional mcmc moves required
  
  # remaining mcmc moves step
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
        
        theta_prop <- TransformationFun(theta_trans_prop,ub,TRUE) #compute inverse transformation  
        
        if (MotilityConstant){
          data <- Main_Simulate(c(theta_prop, KnownTheta),SetupVars,T_record,CellTracking) #simulate data 
          
          if (data$ExitSimStatus == 1){ #check if simulation finished correctly
            next #reject proposal
          }
          
          sx_prop <- GenerateSummaryStatistics(data, CellTracking, Xmax)[1:3] #generate summary statistics
          
        } else {
          data <- Main_Simulate(c(KnownTheta, theta_prop),SetupVars,T_record,CellTracking) #simulate data 
          
          if (data$ExitSimStatus == 1){ #check if simulation finished correctly
            next #reject proposal
          }
          
          sx_prop <- GenerateSummaryStatistics(data, CellTracking, Xmax)[-c(1:3)] #generate summary statistics
        }
        
        if (any(is.na(sx_prop))) { check integrity of summary statistics
          next
        }
        
        d_prop <- norm(sx_prop - sy, "2") #calculate euclidean distance between observed and simulated data
        
        if (d_prop < tol){ #accept proposal
          theta_trans_j <- theta_trans_prop 
          n_accept[j - N_a] <- n_accept[j - N_a] + 1
          d <- d_prop
          sx <- sx_prop
        }
      }
      c(n_accept[j - N_a],theta_trans_j,d,sx)
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
  write.xlsx(c(p_accept = p_accept, tolmax = tolmax, unique = nrow(unique(Samples_trans))/N), paste0("SMCABC_3params_progress_", Ntrack, ".xlsx"), col.names = TRUE, row.names = FALSE)
}
Samples <- cbind(TransformationFun(Samples_trans[,1:num_params],ub,TRUE),Samples_trans[,-c(1:num_params)]) 

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
theta_regressed <- TransformationFun(theta_trans_regressed[,1:3],ub,TRUE) 
colnames(theta_regressed) <- colnames(Samples_trans)[1:num_params]

#save results
write.xlsx(list(Samples = Samples, theta_regressed = theta_regressed), paste0("SMCABC_3params_Samples_", Ntrack, ".xlsx"), col.names = TRUE, row.names = FALSE)
