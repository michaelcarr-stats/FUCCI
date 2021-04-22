# Raw Data

This folder contains the mp4 file `FUCCI_WM983c Time Lapse.mp4` of a FUCCI scratch assay experiment taken with permission from
Vittadello et al. (2018). This video file has been converted into individual frames for proccessing in R
and ImageJ. The Data proccessed in ImageJ including cell trajectory data, initial position and final position
data has been collected and saved as FUCCI_DATA_ImageJ.xlsx.

# DataProcessing

This folder contains `DataProcessing.R` use to merge the ImageJ data and the experimental images in
order to retrieve the cell phase identity which was lost when using ImageJ. The output of this function
is a xlsx file named `FUCCI_processed.xlsx` with sheets:

- InitPos, with variables: x, y and color; where x and y are the cartesian coordinates and color is the cell cycle indicator
- FinalPos, with variables: x, y and color; where x and y are the cartesian coordinates and color is the cell cycle indicator
- CellTracking: , with variables: x, y, color, ntrack, frame; where x and y are the cartesian coordinates, 
  color is the cell cycle indicator, ntrack is the cell trajectory identifier and frame is the image in the sequence.

# SMCABC returned data

This folder contains a xlsx file named `SMCABC_DATA.xlsx` which contains the regression adjusted posterior samples used in the 
study by Carr et al. (2021) with sheets:

- Experiment_CellTracking, Samples from the posterior using Experimental data and cell tracking summary statistics with 
variables: Rr, Ry, Rg, Mr, My and Mg which are samples from the parameters of the model; and sx1,...,sx6 which are 
the respective summary statistics
- Experiment_CellDensity, Samples from the posterior using Experimental data and cell density summary statistics with 
variables: Rr, Ry, Rg, Mr, My and Mg which are samples from the parameters of the model; and sx1,...,sx15 which are 
the summary statistics
- ProliferationData, Samples from the posterior using Synthetic data and holding motility parameters constant with 
variables: Rr, Ry and Rg.
- MotilityData, Samples from the posterior using Synthetic data and holding transition parameters constant with 
variables: Mr, My and Mg which are the model motility parameters; Simulation which indicates the synthetic data set; 
ntrack indicates the number of cell tracked if cell tracking data was used; and MotilityData idicates if cell density 
or cell tracking data was used.
- SixParameterData, Samples from the posterior using Synthetic data and holding transition parameters constant with 
variables: Rr, Ry, Rg, Mr, My and Mg which are the model parameters; Simulation which indicates the synthetic data set; 
and MotilityData idicates if cell density or cell tracking data was used.


# References

Carr, M. J., Simpson, M. J. & Drovandi, C. (2021). Estimating parameters of a stochastic cell invasion model with fluorescent cell 
cycle labelling using Approximate Bayesian Computation. bioarXiv, doi: 
https://doi.org/10.1101/2021.04.20.440712

Vittadello, S. T., McCue, S. W., Gunasingh, G., Haass, N. K. & Simpson, M. J. (2018). Mathematical modelsfor cell 
migration with real-time cell cycle dynamics. Biophysical Journal, 114(5), 1241–1253.
