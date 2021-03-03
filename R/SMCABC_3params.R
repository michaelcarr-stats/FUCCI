library(Rcpp)
library(RcppArmadillo)
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(foreach)
library(MASS)
library(pracma)
library(doParallel)
library(doRNG)
library(openxlsx)
library(stats)
library(FUCCIGillespie)

source("SimulationSetup.R") 
source("GenerateSummaryStatistics.R") 
source("TransformationFun.R") 

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
data <- Main_Simulate(theta,SetupVars,T_record = 48, CellTracking)
sy <- GenerateSummaryStatistics(data, CellTracking, Xmax)

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

N_a <- ceiling(alpha*N)

## SMC MAIN

#ABC-rejection sampling step 
Samples <- 
  foreach(i = 1:N, .combine = "rbind", .packages = c("Rcpp", "RcppArmadillo", "pracma", "MASS", "FUCCIGillespie")) %dorng% {
    while (TRUE){
      
      theta <- runif(3,lb,ub)
      
      if (MotilityConstant){
        data <- Main_Simulate(c(theta, KnownTheta),SetupVars,T_record, CellTracking)
        
        if (data$ExitSimStatus == 1){
          next #reject proposal
        }
        
        sx <- GenerateSummaryStatistics(data, CellTracking, Xmax)[1:3]
        
      } else {
        data <- Main_Simulate(c(KnownTheta, theta),SetupVars,T_record, CellTracking)
        
        if (data$ExitSimStatus == 1){
          next #reject proposal
        }
        
        sx <- GenerateSummaryStatistics(data, CellTracking, Xmax)[-c(1:3)]
      }
      
      if (any(is.na(sx))) {
        next #reject proposal
      }
      
      d <- norm(sx - sy, "2")
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

p_accept <- 1;
tolmax <- max(Samples[,"d"]);

Samples_trans <- cbind(TransformationFun(Samples[,1:num_params],ub,FALSE),Samples[,-c(1:num_params)])

while (p_accept > pacc_target && tolmax > tol_target) {
  
  Samples_trans <- Samples_trans[order(Samples_trans[,"d"]),]
  tol <- Samples_trans[N-N_a,"d"]
  
  rows <- sample(1:(N-N_a),N_a, replace = T)
  Samples_trans[(N-N_a+1):N,] <- Samples_trans[rows,]
  
  cov_rw <- cov(Samples_trans[1:(N-N_a),1:num_params]) # computing tunning parameter
  n_accept <- rep(0, N - N_a)
  
  #trial mcmc step
  tmp <-
    foreach (j = (N - N_a + 1):N, .combine = "rbind", .packages = c("Rcpp","RcppArmadillo","pracma", "MASS", "FUCCIGillespie")) %dorng% {
      theta_trans_j <- Samples_trans[j,1:num_params] 
      d <-  Samples_trans[j,num_params+1] 
      sx <- Samples_trans[j,(num_params+2):ncol(Samples_trans)] 
      for (k in 1:mcmc_trials) {
        
        theta_trans_prop <- mvrnorm(1, mu = theta_trans_j, Sigma = cov_rw)
        MHr <- exp(sum(log(exp(-theta_trans_prop))) - sum(log((1+exp(-theta_trans_prop))^2)) - (sum(log(exp(-theta_trans_j))) - sum(log((1+exp(-theta_trans_j))^2))))
        if (runif(1) > MHr) {
          next #reject proposal
        }
        
        theta_prop <- TransformationFun(theta_trans_prop,ub,TRUE) 
        
        if (MotilityConstant){
          data <- Main_Simulate(c(theta_prop, KnownTheta),SetupVars,T_record,CellTracking)
          
          if (data$ExitSimStatus == 1){
            next #reject proposal
          }
          
          sx_prop <- GenerateSummaryStatistics(data, CellTracking, Xmax)[1:3]
          
        } else {
          data <- Main_Simulate(c(KnownTheta, theta_prop),SetupVars,T_record,CellTracking)
          
          if (data$ExitSimStatus == 1){
            next #reject proposal
          }
          
          sx_prop <- GenerateSummaryStatistics(data, CellTracking, Xmax)[-c(1:3)]
        }
        
        if (any(is.na(sx_prop))) {
          next
        }
        
        d_prop <- norm(sx_prop - sy, "2")
        
        if (d_prop < tol){
          theta_trans_j <- theta_trans_prop #accept proposal
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
  
  p_accept <- sum(n_accept)/((N_a-1)*mcmc_trials)
  mcmc_iters <- ceiling(log(c)/log(1-p_accept))
  
  # remaining mcmc moves step
  tmp <-
    foreach (j = (N - N_a + 1):N, .combine = "rbind", .packages = c("Rcpp","RcppArmadillo","pracma", "FUCCIGillespie")) %dorng% {
      theta_trans_j <- Samples_trans[j,1:num_params] 
      d <-  Samples_trans[j,num_params+1] 
      sx <- Samples_trans[j,(num_params+2):ncol(Samples_trans)] 
      for (k in 1:(mcmc_iters - mcmc_trials)){
        
        theta_trans_prop <- mvrnorm(1, mu = theta_trans_j, Sigma = cov_rw)
        MHr <- exp(sum(log(exp(-theta_trans_prop))) - sum(log((1+exp(-theta_trans_prop))^2)) - (sum(log(exp(-theta_trans_j))) - sum(log((1+exp(-theta_trans_j))^2))))
        if (runif(1) > MHr){
          next #reject proposal
        }
        
        theta_prop <- TransformationFun(theta_trans_prop,ub,TRUE) 
        
        if (MotilityConstant){
          data <- Main_Simulate(c(theta_prop, KnownTheta),SetupVars,T_record, CellTracking)
          
          if (data$ExitSimStatus == 1){
            next #reject proposal
          }
          
          sx_prop <- GenerateSummaryStatistics(data, CellTracking, Xmax)[1:3]
          
        } else {
          data <- Main_Simulate(c(KnownTheta, theta_prop),SetupVars,T_record, CellTracking)
          
          if (data$ExitSimStatus == 1){
            next #reject proposal
          }
          
          sx_prop <- GenerateSummaryStatistics(data, CellTracking, Xmax)[-c(1:3)]
        }
        
        if (any(is.na(sx_prop))) {
          next
        }
        
        d_prop <- norm(sx_prop - sy, "2")
        
        if (d_prop < tol){
          theta_trans_j <- theta_trans_prop #accept proposal
          n_accept[j - N_a] <- n_accept[j - N_a] + 1
          d <- d_prop
          sx_prop <- sx 
        }
      }
      c(n_accept[j - N_a],theta_trans_j,d,sx)
    }
  n_accept <- tmp[,1]
  Samples_trans <- rbind(Samples_trans[1:(N-N_a),], tmp[,2:ncol(tmp)])
  rm(tmp)
  
  num_mcmc_iters <- max(0,mcmc_iters-mcmc_trials) + mcmc_trials
  p_accept <- sum(n_accept)/((N_a-1)*num_mcmc_iters)
  mcmc_iters <- ceiling(log(c)/log(1-p_accept))
  mcmc_trials <- ceiling(mcmc_iters/2)
  tolmax <- max(Samples_trans[,"d"])
  
  cat("\n")
  message(paste0("Acceptance Rate: ",round(p_accept,2)))
  message(paste0("Max tolerance: ",round(tolmax,2)))
  message(paste0("% unique: ",round(nrow(unique(Samples_trans))/N,2)))
  cat("\n")
  write.xlsx(c(p_accept = p_accept, tolmax = tolmax, unique = nrow(unique(Samples_trans))/N), paste0("SMCABC_3params_progress_", Ntrack, ".xlsx"), col.names = TRUE, row.names = FALSE)
}

Samples <- cbind(TransformationFun(Samples_trans[,1:num_params],ub,TRUE),Samples_trans[,-c(1:num_params)]) 

theta_trans <- as_tibble(Samples_trans[,1:num_params])
d_vec <- Samples_trans[,"d"]
summaries <- as_tibble(Samples_trans[,(num_params + 2):ncol(Samples_trans)])
the_weights = as_tibble(0.75*(1-(d_vec/max(d_vec))^2))

theta_trans_regressed <- matrix(NA, nrow = N, ncol = num_params)
for (i in 1:num_params){
  data <- bind_cols(theta_trans[,i],summaries,the_weights) %>% setNames(c("y", colnames(summaries), "the_weights"))
  b <- lm(y~. - the_weights, data = data, weights = the_weights)
  theta_trans_regressed[,i] <- as.matrix(theta_trans[,i] - as.matrix(summaries - repmat(sy, N, 1)) %*%  t(t(b$coefficients[-1])))
}
theta_regressed <- TransformationFun(theta_trans_regressed[,1:3],ub,TRUE) 
colnames(theta_regressed) <- colnames(Samples_trans)[1:num_params]

write.xlsx(list(Samples = Samples, theta_regressed = theta_regressed), paste0("SMCABC_3params_Samples_", Ntrack, ".xlsx"), col.names = TRUE, row.names = FALSE)
