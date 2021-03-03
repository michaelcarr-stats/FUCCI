# Data

This folder contains the mp4 file of a FUCCI scratch assay experiment taken with permission from
Vittadello et al. (2018). This video file has been converted into individual frames for proccessing in R
and ImageJ. The Data proccessed in ImageJ including cell trajectory data, initial position and final position
data has been retrieved and placed in FUCCI_DATA_ImageJ.xlsx.

# DataProcessing

This folder contains `DataProcessing.R` use to merge the ImageJ data and the experimental images in
order to retrieve the cell phase identity which was lost when using ImageJ. The output of this function
is a xlsx file named with sheets:

- InitPos, with variables: x, y and color; where x and y are the cartesian coordinates and color is the cell cycle indicator
- FinalPos, with variables: x, y and color; where x and y are the cartesian coordinates and color is the cell cycle indicator
- CellTracking: , with variables: x, y, color, ntrack, frame; where x and y are the cartesian coordinates, 
  color is the cell cycle indicator, ntrack is the cell trajectory identifier and frame is the image in the sequence.


# References

Vittadello, S. T., McCue, S. W., Gunasingh, G., Haass, N. K. & Simpson, M. J. (2018). Mathematical modelsfor cell 
migration with real-time cell cycle dynamics.Biophysical Journal,114(5), 1241–1253.
