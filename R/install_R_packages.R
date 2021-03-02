r_packages_list <- c("Rcpp","RcppArmadillo","readxl","dplyr","tidyr","ggplot2","foreach","MASS","pracma","doParallel","doRNG","openxlsx","devtools","stats")

for (package in r_packages_list){
  if (!(package %in% installed.packages()[,"Package"])){
    install.packages(package, dependencies = T, repos = "https://cran.csiro.au/", lib = "~/R/library")
  }
}

#install simulation model from github
if (!("FUCCIGillespie" %in% installed.packages()[,"Package"])){
  devtools::install_github("michaelcarr-stats/FUCCI/FUCCIGillespie", lib = "~/R/library")  
}
