# Initialization.R

## This chunk loads some needed R packages, specifies the file used for the illustration
## of mechanization, specifies the variables needed from the archive file, and
## specifies some time shifts. The data file used here is not part of the standard
## archives. It was produced specially to provide high-rate values of the measured
## variables, and it is therefore preserved in .Rdata format so that file can be
## archived and re-used to ensure reproducibility.

## chunk "initialization.R"
library(knitr)
opts_chunk$set(echo=FALSE, include=FALSE, fig.lp="fig:")
# note that fig.pos="center" gave errors, changed to fig.align
opts_chunk$set(fig.width=6, fig.height=5, fig.align="center", digits=4)
thisFileName <- "KalmanFilterTechNote"
require(Ranadu, quietly = TRUE, warn.conflicts=FALSE)
require(numDeriv)    ## needed for the jacobian() function
library(signal)
library(reshape2)
library(grid)
options(stringsAsFactors=FALSE)

## temporary; remove once Ranadu is updated
source ('~/RStudio/Ranadu/R/ShiftInTime.R')
source ('~/RStudio/Ranadu/R/theme_WAC.R')
source ('~/RStudio/Ranadu/R/XformLA.R')

setwd ('~/RStudio/Ranadu/KalmanFilter')
Directory <- DataDirectory ()
Flight <- "rf15HR" 			
Project = "DWIRU"     		
ProjectDir <- "DEEPWAVE"
fname = sprintf("%s%s/%s%s.nc", Directory, ProjectDir, Project, Flight)
VarList <- c('BLATA', 'BLONGA', 'BNORMA', 'BPITCHR', 'BROLLR', 'BYAWR')
VarList <- c(VarList, 'LAT', 'LON', 'ALT', 'GGLAT', 'GGLON', 'GGALT')
VarList <- c(VarList, 'VEW', 'VNS', 'VSPD', 'GGVEW', 'GGVNS', 'GGVSPD')
VarList <- c(VarList, 'PITCH', 'ROLL', 'THDG')
ReloadData <- FALSE
# ReloadData <- TRUE
SaveRData <- sprintf("%s.Rdata", thisFileName)
if (ReloadData) {
  Data <- getNetCDF (fname, VarList, Start=31000, End=40000)
  ## remove heading adjustments if added during initial processing
  Z <- data.frame(getAttributes(Data$THDG, .print=FALSE))
  if ('CalibrationCoefficients' %in% names(Z)) {
    THDGoffset <- Z$CalibrationCoefficients[1]
  } else {
    THDGoffset <- 0
  }
  ## re signs: assumed lag has been corrected, so use - sign to remove
  if ('TimeLag' %in% names(Z)) {
    Data$THDG <- ShiftInTime(Data$THDG, .rate=25, .shift=-Z$TimeLag[1])
  }
  Data$THDG <- (Data$THDG - THDGoffset) %% 360
  save (Data, file=SaveRData)
} else {
  load (file=SaveRData)
}
Cradeg <- pi/180
.shift = 60
Data$PITCH <- ShiftInTime (Data$PITCH, .rate=25, .shift)
Data$ROLL <- ShiftInTime (Data$ROLL, .rate=25, .shift)
# for (V in c('BYAWR')) {
#   Data[, V] <- ShiftInTime (Data[, V], 25, .shift)
# }
## Get some derivatives used later (but better value are obtained later
## via Savitzgy-Golay polynomials.)
Data$pdot <- c(0, diff (Data$PITCH)) * 25
Data$rdot <- c(0, diff (Data$ROLL)) * 25
Data$hdot <- c(0, diff (Data$THDG))
Data$hdot <- (Data$hdot + 540) %% 360 - 180
Data$hdot <- Data$hdot * 25
