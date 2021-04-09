
SimulationSetup<- function(ntrack, Xmax, Ymax, InitPos = NULL, CellTrackingData = NULL){
  
  #load functions
  sourceCpp("InitDomain_x.cpp")
  sourceCpp("InitDomain_y.cpp")
  sourceCpp("InitDomain.cpp")
  sourceCpp("InitialCellCount.cpp")
  source("InitialMotilitySelect.R")
    
  ########################################
  #-------------Input-Values-------------#
  ########################################
  simuNum = 3; #number of simulation realisations
  delta = 20; #Lattice spacing[mum]
  t = 0; #start time
  tstop = 48; #end time
  columnNum = ceiling(Xmax / delta) + 1; #Corresponding total column nodes
  rowNum = ceiling(Ymax / (delta * sqrt(3) / 2)) + 1; #Corresponding total row nodes
  InitialL = 20 * delta; #Length of the initial agent pack
  
  #initial number of cells in domain - for experiments which dont not use location data
  #Cells placed based on percentage and density
  InitialDensity = 1; #Density of the initial pack
  percentR = 0.33; # percentage of the initial domain to be red
  percentY = 0.33; # percentage of the initial domain to be yellow
  percentG = 0.33; # percentage of the initial domain to be green
  #Cells placed based on counts
  SetCells = TRUE #whether or not to use this method 
  Nred_0 = 119 #initial number of red cells
  Nyellow_0 = 35 #initial number of yellow cells 
  Ngreen_0 = 121 #initial number of green cells
  
  BC = 2; #Boundary condition, 1.Periodic, 2. No flux
  
  ########################################
  #------------Initialisation------------#
  ########################################
  
  #initialise domain coordinates
  domain_x <- InitDomain_x(rowNum, columnNum, delta) 
  domain_y <- InitDomain_y(rowNum, columnNum, delta)
  
  if (is.null(InitPos) & is.null(CellTrackingData)) { #not using real data
    percentIS <- matrix(c(percentR,percentY,percentG), nrow = 3, ncol = 1) 
    countIS <- matrix(c(Nred_0,Nyellow_0,Ngreen_0), nrow = 3, ncol = 1)
    domain <- InitDomain(domain_x, rowNum, columnNum, Xmax, InitialL, InitialDensity, percentIS, countIS, SetCells) #initialise domain
    
    TrackedCells <- InitialMotilitySelect(ntrack, domain, domain_x, Xmax) #select cells to track
    
    RowPosCell = TrackedCells$RowPosCell # row index in domain of tracked cells
    ColPosCell = TrackedCells$ColPosCell # column index in domain of tracked cells
    CellSelectedIndex = TrackedCells$CellSelectedIndex #cell phase of tracked cells (1,2, or 3)
    
  } else { #use real data
    i <- round(nrow(domain_x)-(InitPos$y*2)/(delta*sqrt(3))) #row index in domain of cells
    j <- round(InitPos$x/delta + 1) #column index in domain of cells
    
    Indexs <- cbind(i,j, color = InitPos$color) %>% #remove duplicate row and column indexs
      as_tibble() %>% 
      distinct(i,j, .keep_all = T) %>% 
      mutate(i = as.numeric(i), j = as.numeric(j))
    
    domain <- matrix(0, nrow = rowNum, ncol = columnNum) #initialise domain
    for (row in 1:nrow(Indexs)) {
      if (Indexs[row,"color"] == 1) {
        
        domain[Indexs[[row,"i"]],Indexs[[row,"j"]]] <- 1 #place red cell on lattice
        
      } else if (Indexs[row,"color"] == 2) {
        
        domain[Indexs[[row,"i"]],Indexs[[row,"j"]]] <- 2 #place yellow cell on lattice
        
      } else if (Indexs[row,"color"] == 3) {
        
        domain[Indexs[[row,"i"]],Indexs[[row,"j"]]] <- 3  #place green cell on lattice
        
      }
    }
    
    # Cell Tracking
    if (is.null(CellTrackingData)) { #default choose cell to track
      TrackedCells <- InitialMotilitySelect(ntrack, domain, domain_x, Xmax)  #select cells to track
      
      RowPosCell = TrackedCells$RowPosCell # row index in domain of tracked cells
      ColPosCell = TrackedCells$ColPosCell # column index in domain of tracked cells
      CellSelectedIndex = TrackedCells$CellSelectedIndex #compute index of cell selected
      
    } else  { #manually choose cells to track
      RowPosCell <- rep(NA,ntrack) #initialise row index in domain of tracked cells
      ColPosCell <- rep(NA,ntrack) #initialise column index in domain of tracked cells
      CellSelectedIndex <- rep(NA,ntrack) #initialise index of cell selected
      
      for (i in 1:ntrack){
        Initx <- unlist(CellTrackingData[CellTrackingData$frame == 1,"x"]) #get initial x coordinate of tracked cell
        Inity <- unlist(CellTrackingData[CellTrackingData$frame == 1,"y"]) #get initial y coordinate of tracked cell
        
        tmp <- sqrt((InitPos$x - Initx[i])^2 + (InitPos$y - Inity[i])^2) #compute radial distance of tracked cell away from all cells
        
        index <- match(min(tmp), tmp) #determine cell closest to tracked cell
        RowPosCell[i] <- as.numeric(round(rowNum - 2*InitPos[index,"y"]/(delta*sqrt(3)))) #assign initial row index based on closest cell
        ColPosCell[i] <- as.numeric(round(InitPos[index,"x"]/delta + 1)) #assign initial column index based on closest cell
        CellSelectedIndex[i] <- (ColPosCell[i] - 1)*rowNum + RowPosCell[i] #compute index of cell selected
      }
    }
    
  }
  NStartCount <- InitialCellCount(domain, domain_x, domain_y, rowNum, columnNum, Xmax, Ymax) #compute initial number of cells
  CellSelected = domain[CellSelectedIndex] #cell phase of tracked cells

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


