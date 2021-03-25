
SimulationSetup<- function(ntrack, Xmax, Ymax, InitPos = NULL, CellTrackingData = NULL){

  sourceCpp("InitDomain_x.cpp")
  sourceCpp("InitDomain_y.cpp")
  sourceCpp("InitDomain.cpp")
  sourceCpp("InitialCellCount.cpp")
  source("InitialMotilitySelect.R")
    
  ########################################
  #-------------Input-Values-------------#
  ########################################
  simuNum = 3;
  delta = 20; #Lattice spacing[mum]
  t = 0; #start time
  tstop = 48; #end time
  columnNum = ceiling(Xmax / delta) + 1; #Corresponding total column nodes
  rowNum = ceiling(Ymax / (delta * sqrt(3) / 2)) + 1; #Corresponding total row nodes
  InitialL = 20 * delta; #Length of the initial agent pack
  InitialDensity = 1; #Density of the initial pack
  
  percentR = 0.33; # percentage of the initial domain to be red
  percentY = 0.33; # percentage of the initial domain to be yellow
  percentG = 0.33; # percentage of the initial domain to be green
  
  #initial number of cells in domain - for experiments which dont not use location data
  SetCells = TRUE #whether or not to use this method 
  Nred_0 = 119 #initial number of red cells
  Nyellow_0 = 35 #initial number of yellow cells 
  Ngreen_0 = 121 #initial number of green cells
  
  BC = 2; #Boundary condition, 1.Periodic, 2. No flux
  
  ########################################
  #------------Initialisation------------#
  ########################################
  
  domain_x <- InitDomain_x(rowNum, columnNum, delta) 
  domain_y <- InitDomain_y(rowNum, columnNum, delta)
  
  if (is.null(InitPos) & is.null(CellTrackingData)) { #dont use real data
    percentIS <- matrix(c(percentR,percentY,percentG), nrow = 3, ncol = 1)
    countIS <- matrix(c(Nred_0,Nyellow_0,Ngreen_0), nrow = 3, ncol = 1)
    domain <- InitDomain(domain_x, rowNum, columnNum, Xmax, InitialL, InitialDensity, percentIS, countIS, SetCells)
    
    TrackedCells <- InitialMotilitySelect(ntrack, domain, domain_x, Xmax)
    
    RowPosCell = TrackedCells$RowPosCell
    ColPosCell = TrackedCells$ColPosCell
    CellSelectedIndex = TrackedCells$CellSelectedIndex
    
  } else { #use real data
    i <- round(nrow(domain_x)-(InitPos$y*2)/(delta*sqrt(3)))
    j <- round(InitPos$x/delta + 1)
    
    Indexs <- cbind(i,j, color =InitPos$color) %>% 
      as_tibble() %>% 
      distinct(i,j, .keep_all = T) %>% 
      mutate(i = as.numeric(i), j = as.numeric(j))
    
    domain <- matrix(0, nrow = rowNum, ncol = columnNum)
    for (row in 1:nrow(Indexs)) {
      if (Indexs[row,"color"] == 1) {
        
        domain[Indexs[[row,"i"]],Indexs[[row,"j"]]] <- 1
        
      } else if (Indexs[row,"color"] == 2) {
        
        domain[Indexs[[row,"i"]],Indexs[[row,"j"]]] <- 2
        
      } else if (Indexs[row,"color"] == 3) {
        
        domain[Indexs[[row,"i"]],Indexs[[row,"j"]]] <- 3 
        
      }
    }
    
    # Cell Tracking
    if (is.null(CellTrackingData)) {
      TrackedCells <- InitialMotilitySelect(ntrack, domain, domain_x, Xmax)
      
      RowPosCell = TrackedCells$RowPosCell
      ColPosCell = TrackedCells$ColPosCell
      CellSelectedIndex = TrackedCells$CellSelectedIndex
    } else  {
      RowPosCell <- rep(NA,ntrack)
      ColPosCell <- rep(NA,ntrack)
      CellSelectedIndex <- rep(NA,ntrack)
      for (i in 1:ntrack){
        Initx <- unlist(CellTrackingData[CellTrackingData$frame == 1,"x"])
        Inity <- unlist(CellTrackingData[CellTrackingData$frame == 1,"y"])
        
        tmp <- sqrt((InitPos$x - Initx[i])^2 + (InitPos$y - Inity[i])^2)
        
        index <- match(min(tmp), tmp)
        RowPosCell[i] <- as.numeric(round(rowNum - 2*InitPos[index,"y"]/(delta*sqrt(3))))
        ColPosCell[i] <- as.numeric(round(InitPos[index,"x"]/delta + 1))
        CellSelectedIndex[i] <- (ColPosCell[i] - 1)*rowNum + RowPosCell[i]
      }
    }
    
  }
  NStartCount <- InitialCellCount(domain, domain_x, domain_y, rowNum, columnNum, Xmax, Ymax)
  CellSelected = domain[CellSelectedIndex]

  ########################################
  #------------Store in List-------------#
  ########################################
  
  SetupVars = list(
    ntrack = ntrack,
    simuNum = simuNum,
    delta = delta,
    t = t,
    tstop = tstop,
    Xmax = Xmax,
    Ymax = Ymax,
    columnNum = columnNum, 
    rowNum = rowNum, 
    BC = BC,
    domain = domain,
    domain_x = domain_x,
    domain_y = domain_y,
    NStartCount = NStartCount,
    RowPosCell = RowPosCell - 1, #convert to cpp indexing
    ColPosCell = ColPosCell - 1, #convert to cpp indexing
    CellSelected = CellSelected
  )
  
  return(SetupVars)
  
}


