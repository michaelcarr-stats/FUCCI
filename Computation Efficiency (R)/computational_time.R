library(FUCCIGillespie)
library(ggplot2)
library(readxl)
library(openxlsx)
library(tidyverse)
library(Rcpp)
library(foreach)
library(doRNG)
library(MASS)
library(pracma)
library(doParallel)
library(doRNG)

#setup parrallel computing
numCores <- 5
registerDoParallel(numCores)

set.seed(1234)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("../R") #set to location of simulation code
source("SimulationSetup.R") 

InitPos <- read_excel("../Data/DataProcessing/FUCCI_processed.xlsx", sheet = "InitPos")
FinalPos <- read_excel("../Data/DataProcessing/FUCCI_processed.xlsx", sheet = "FinalPos")
CellTrackingData <- read_excel("../Data/DataProcessing/FUCCI_processed.xlsx", sheet = "CellTracking")
Ntrack <- length(unique(CellTrackingData$ntrack))

Xmax <- 1309.09 #image width
Ymax <- 1745.35 #image height
SetupVars <- SimulationSetup(Ntrack, Xmax, Ymax, InitPos, CellTrackingData)
T_record <- 48
CellTracking <- !is.null(CellTrackingData)

N <- 1000
cputime <- 
  foreach(i = 1:N, .combine = "rbind", .packages = c("Rcpp", "RcppArmadillo", "pracma", "MASS", "FUCCIGillespie")) %dorng% {
      
    theta <- c(runif(3,0,1), runif(3,0,10))
    
    s <- Sys.time()
    Main_Simulate(theta, SetupVars = SetupVars, CellTracking = CellTracking, T_record = T_record) #simulate data
    (Sys.time() - s)[[1]]
}


#time in seconds
mean(cputime) #average time
quantile(cputime,.025) #5% empirical quantile
quantile(cputime,.975) # 95% empirical quantile


#time evaluated at the posterior modes
posterior <- read_excel("../Data/SMCABC returned data/SMCABC_DATA.xlsx", sheet = "Experiment_CellTracking")
theta <- apply(posterior[,1:6], MARGIN = 2, FUN = mean)
s <- Sys.time()
Main_Simulate(theta, SetupVars = SetupVars, CellTracking = CellTracking, T_record = T_record) #simulate data
(Sys.time() - s)[[1]]
