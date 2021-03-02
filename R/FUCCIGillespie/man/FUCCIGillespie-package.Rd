\name{FUCCIGillespie-package}
\alias{FUCCIGillespie-package}
\alias{FUCCIGillespie}
\docType{package}
\title{
  Simulate a 2D scratch assay experiment with flurocent cell cycle indicators.
}
\usage{
Main_Simulate(theta, SetupVars, T_record, CellTracking)
}
\arguments{
  \item{theta}{A vector of the six model parameters. In order of red, yellow and green cell phase, the first three elements are the cell cycle transition rates and last three elements are the motility rates.}
  \item{SetupVars}{List object returned from SimulationSetup.R specifying simulation environment}
  \item{T_record}{The time at which summary statistics are to be recorded}
  \item{CellTracking}{If TRUE summary statistics recorded will include cell trajectory data. If FALSE summary statistics recorded will include cell density data}
}
\value{
A list of summary statistics to be proccessed using GenerateSummaryStatistics.R 
}
\description{
  This package simulates a Markov Process of a scratch assay experiment with fluorescent cell cycle labeling using the Gillespie algorithm.
}
\references{}

