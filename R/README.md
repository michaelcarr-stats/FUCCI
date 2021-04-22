# FUCCIGillespie

`FUCCIGillespie` is an `R` package that simulates a Markov proccess of a scratch assay experiment with flurocent cell cycle
indicators using the Gillespie Algorithm. See Simpson et al. (2018) for modelling information.

### How to install
1. In `R` make sure `devtools` is installed. Install with `install.packages("devtools")`.
2. Run `devtools::install_github(""michaelcarr-stats/FUCCI/R/FUCCIGillespie"")`

# SMCABC

`SMCABC.R` is the implementation of the Sequential Monte Carlo Aproximate Bayesian Computation (SMC-ABC) replenishment 
algorithm of Drovandi and Pettitt (2011) which Carr et al. (2021) use to estimate the model parameters in the 
`FUCCIGillespie` package. This script is the main script which should be used to call the subsequent functions within 
this repository, including:

  - `Main_simulate.cpp` from the `FUCCIGillespie` package
  - `TransformationFUN.R` 
  - `GenerateSummaryStatistics.R`
  - `SimulationSetup.R`, which calls:
    - `InitialMotilitySelect.R`
    - `InitialCellCount.cpp`
    - `InitDomain.cpp`
    - `InitDomain_x.cpp`
    - `InitDomain_y.cpp`

To use this function simply download the listed functions into your directory, install the `FUCCIGillespie` package, and 
download and run the `install_R_packages.R` function (details below). Alterations to the statistical variables are made
directly in this folder while alterations to the simulation environment are generally made in `SimulationSetup.R`.

A breif overview of the functions listed above will now be presented:

### `Main_simulate.cpp`

This function simulates a 2D scratch assay experiment with flurocent cell cycle indicators. The function is called with
`Main_Simulate(theta, SetupVars, T_record, CellTracking)` with arguments defined by:
  - theta: a vector of the six model parameters. In order of red, yellow and green cell phase, the first three elements are the cell cycle transition rates and last three elements are the motility rates.
  - SetupVarsList: object returned from SimulationSetup.R specifying simulation environment.
  - T_record: the time at which summary statistics are to be recorded.
  - CellTracking: if TRUE summary statistics recorded will include cell trajectory data. If FALSE summary statistics recorded will include cell density data.

The returned object from this function is a list of summary statistic data.

### `TransformationFUN.R`

This function apply a logistic transformation (or the inverse) between 0 and and the upperbound. The function is called with
`TransformationFUN(theta,b,inverse)` with arguments defined by:
  - `theta`: a matrix or vector of model parameters.
  - `b`: the upperbound of the transformation.
  - `inverse`: if TRUE the function applys the inverse transformation else it applies the logistic transformation between 0 and the specified upperbound.

The returned object from this function is a matrix or vector of model parameters.

### `GenerateSummaryStatistics.R`

This function takes returned data from `Main_Simulate.cpp` and generates summary statistics. The function is called with
`GenerateSummaryStatstics(data, CellTracking)` with arguments defined by
  - `data`: list of summary statistic data returned from `Main_Simulate.cpp`.
  - `CellTracking`: if TRUE (default) then assumes cell trajectory data was returned from `Main_Simulate.cpp` and calculates the average distance travelled through each phase of the cell cycle. If FALSE assumes cell density data was returned from `Main_Simulate.cpp` and calculates the median position and interquartile range of cells in each phase of the cell cycle.

The returned object from this function is a vector of summary statistics.

### `SimulationSetup.R`

This function creates a list of variables used to initalise the simulation model `Main_simulate.cpp`, The function is called with 
`SimulationSetup(ntrack, Xmax, Ymax, InitPos, CellTrackingData)` with arguments defined by
  - `ntrack`: number of cells to track.
  - `Xmax`: width of the experiment. If using experimental data ensure dimensions are the same.
  - `Ymax`: height of the experiment. If using experimental data ensure dimensions are the same.
  - `InitPos`: dataframe of the Initial position of cells with x, y and color named variables representing the 2D cartesian coordinates and cell phase identity, respectively.
  - `CellTrackingData`: dataframe of position of cells throughout experiment with x, y, color, ntrack, frame as named variables representing the 2D cartesian coordinates, cell phase identity, cell trajectory identifier, and image frame  identifier in the sequence of images, respectiely.

Additionally, there are arguments within this function which can be used to customize the simulation environment in differnet ways. 
The returned object from this function is a list of variables which is to be passed into `Main_simulate.cpp`.

# SMCABC_3params

This function is a modified version of `SMCABC.R` which calls the same subfunctions but differs such that either the cell cycle
transition or motility parameters are held constant. Additionally, the data section is formated to use synthetic data sets as
opposed to using the experimental data sets as this function was intended to be used as a testing ground for different summary statistics.


# install_R_packages

This script should be called prior to `SMCABC.R` or `SMCABC_3params.R` to install or load all neccessary packages
including FUCCIGillespie with built in functionality to handle already installed packages.


# References

Carr, M. J., Simpson, M. J. & Drovandi, C. (2021). Estimating parameters of a stochastic cell invasion model with
fluorescent cell cycle labelling using Approximate Bayesian Computation. bioarXiv, doi: 
https://doi.org/10.1101/2021.04.20.440712

Drovandi, C. C. & Pettitt, A. N. (2011). Estimation of parameters for macroparasite population evolution usingapproximate
bayesian computation. Biometrics,67(1), 225–233.

Simpson,  M.  J.,  Jin,  W.,  Vittadello,  S.  T.,  Tambyah,  T.  A.,  Ryan,  J.  M.,  Gunasingh,  G.,  Haass,  N.  K.  
& McCue, S. W. (2018). Stochastic models of cell invasion with fluorescent cell cycle indicators. PhysicaA: Statistical 
Mechanics and its Applications, 510, 375–386.
