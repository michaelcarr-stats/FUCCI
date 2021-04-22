# SMCABC

`SMCABC.m` is the implementation of the Sequential Monte Carlo Aproximate Bayesian Computation (SMC-ABC) replenishment 
algorithm of Drovandi and Pettitt (2011) which Carr et al. (2021) use to estimate the model parameters in the 
`Main_Simulate.m' model. This script is the main script which should be used to call the subsequent functions within 
this repository, including:

  - `Main_Simulate.m`, which calls:
    - `Migration1.m`
    - `Migration2.m`
    - `GenerateSummaryStatData.m`
  - `TransformationFUN.m` 
  - `GenerateSummaryStatistics.m`
  - `SetupStruct.m`, which calls:
    - `InitialMotilitySelect.m`
    - `InitialCellCount.m`
    - `initialise_domain.m`
Alterations to the statistical variables are made directly in this folder while alterations to the simulation environment are 
generally made in `SetupStruct.m`.

A breif overview of the functions listed above will now be presented:

### `Main_simulate.m`

This function simulates a 2D scratch assay experiment with flurocent cell cycle indicators. The function is called with
`Main_simulate(theta, SetupVars, T_record, CellTracking)` with arguments defined by:
  - theta: a vector of the six model parameters. In order of red, yellow and green cell phase, the first three elements are the cell cycle transition rates and last three elements are the motility rates.
  - SetupVarsList: object returned from SimulationSetup.R specifying simulation environment.
  - T_record: the time at which summary statistics are to be recorded.
  - CellTracking: if TRUE summary statistics recorded will include cell trajectory data. If FALSE summary statistics recorded will include cell density data.

The returned object from this function is a list of summary statistic data and an boolean variable which indicates whether the function ran through to completion.

### `TransformationFUN.m`

This function apply a logistic transformation (or the inverse) between 0 and and the upperbound. The function is called with
`TransformationFUN(theta,b,inverse)` with arguments defined by:
  - `theta`: a matrix or vector of model parameters.
  - `b`: the upperbound of the transformation.
  - `inverse`: if TRUE the function applys the inverse transformation else it applies the logistic transformation between 0 and the specified upperbound.

The returned object from this function is a matrix or vector of model parameters.

### `GenerateSummaryStatistics.m`

This function takes returned data from `Main_simulate.m` and generates summary statistics. The function is called with
`GenerateSummaryStatstics(data, CellTracking)` with arguments defined by
  - `data`: list of summary statistic data returned from `Main_simulate.cpp`.
  - `CellTracking`: if TRUE (default) then assumes cell trajectory data was returned from `Main_simulate.cpp` and calculates the average distance travelled through each phase of the cell cycle. If FALSE assumes cell density data was returned from `Main_simulate.cpp` and calculates the median position and interquartile range of cells in each phase of the cell cycle.

The returned object from this function is a vector of summary statistics.

### `SetupStruct.m`

This function creates a list of variables used to initalise the simulation model `Main_simulate.cpp`, The function is called with 
`SimulationSetup(ntrack, Xmax, Ymax, InitPos, CellTrackingData)` with arguments defined by
  - `ntrack`: number of cells to track.
  - `Xmax`: width of the experiment. If using experimental data ensure dimensions are the same.
  - `Ymax`: height of the experiment. If using experimental data ensure dimensions are the same.
  - `InitPos`: dataframe of the Initial position of cells with x, y and color named variables representing the 2D cartesian coordinates and cell phase identity, respectively.
  - `CellTrackingData`: dataframe of position of cells throughout experiment with x, y, color, ntrack, frame as named variables representing the 2D cartesian coordinates, cell phase identity, cell trajectory identifier, and image frame  identifier in the sequence of images, respectiely.

Additionally, there are arguments within this function which can be used to customize the simulation environment in differnet ways. 
The returned object from this function is a list of variables which is to be passed into `Main_simulate.m`.

# References

Carr, M. J., Simpson, M. J. & Drovandi, C. (2021). Estimating parameters of a stochastic cell invasion model with
fluorescent cell cycle labelling using Approximate Bayesian Computation. bioarXiv, doi: 
https://doi.org/10.1101/2021.04.20.440712

Drovandi, C. C. & Pettitt, A. N. (2011). Estimation of parameters for macroparasite population evolution usingapproximate
bayesian computation. Biometrics, 67(1), 225–233.

Simpson,  M.  J.,  Jin,  W.,  Vittadello,  S.  T.,  Tambyah,  T.  A.,  Ryan,  J.  M.,  Gunasingh,  G.,  Haass,  N.  K.  
& McCue, S. W. (2018). Stochastic models of cell invasion with fluorescent cell cycle indicators. PhysicaA: Statistical 
Mechanics and its Applications, 510, 375–386.
