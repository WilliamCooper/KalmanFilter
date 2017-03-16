## ----initialization, echo=FALSE,include=FALSE----------------------------

## This chunk loads some needed R packages,
## specifies the variables needed from the archive file, and
## specifies some time shifts. 

print (sprintf ('Kalman Processor is starting -- %s', Sys.time()))
thisFileName <- "KalmanFilter"
require(Ranadu, quietly = TRUE, warn.conflicts=FALSE)
require(numDeriv, quietly = TRUE, warn.conflicts=FALSE)    ## needed for the jacobian() function
# library(signal)      ## used for filtering -- but ::signal avoids error message
library(reshape2)    ## used with ggplot facet plots
library(grid)
options(stringsAsFactors=FALSE)

setwd ('~/RStudio/KalmanFilter')

##----------------------------------------------------------
## These are the run options to set via command-line or UI:
GeneratePlots <- TRUE
ShowPlots <- FALSE  ## leave FALSE for batch runs or shiny-app runs
UpdateAKRD <- TRUE
UpdateSSRD <- TRUE
SimpleOnly <- FALSE   ## set TRUE to use CorrectPitch/CorrectHeading w/o KF
ALL <- FALSE
NEXT <- FALSE
ReloadData <- FALSE
ReloadData <- TRUE
NSTEP <- 30      ## update interval

Directory <- DataDirectory ()
##----------------------------------------------------------

Project <- 'CSET'
Project <- 'DEEPWAVE'
ProjectDir <- Project
Flight <- 1
Flight <- 17
if (!interactive()) {  ## can run interactively or via Rscript
  run.args <- commandArgs (TRUE)
  if (length (run.args) > 0) {
    if (nchar(run.args[1]) > 1) {
      Project <- run.args[1]
      ProjectDir <- Project
    }
  } else {
    print ("Usage: Rscript KalmanFilter.R Project Flight UpdateAKRD[y/N] UpdateSSRD[y,N] Simple[y/N] Interval[10]")
    print ("Example: Rscript KalmanFilter.R CSET 1 y y n 15")
    stop("exiting...")
  }
  ## Flight
  if (length (run.args) > 1) {
    if (run.args[2] != 'NEXT') {
      Flight <- as.numeric (run.args[2])
    } else {
      ## Find max rf in data directory,
      ## Use as default if none supplied via command line:
      getNext <- function(ProjectDir, Project) {
        Fl <- sort (list.files (sprintf ("%s%s/", DataDirectory (), ProjectDir),
                                sprintf ("%srf..KF.nc", Project)), decreasing = TRUE)[1]
        if (is.na (Fl)) {
          Flight <- 1
        } else {
          Flight <- sub (".*rf", '',  sub ("KF.nc", '', Fl))
          Flight <- as.numeric(Flight)+1
        }
        return (Flight)
      }
      Flight <- getNext(ProjectDir, Project)
    }
  }  
  
  if (length (run.args) > 2) {
    UpdateAKRD <- run.args[3]
    if (UpdateAKRD == 'y' || UpdateAKRD == 'Y' || UpdateAKRD == 'TRUE') {
      UpdateAKRD <- TRUE
    } else {
      UpdateAKRD <- FALSE
    }
  }
  if (length (run.args) > 3) {
    UpdateSSRD <- run.args[4]
    if (UpdateSSRD == 'y' || UpdateSSRD == 'Y' || UpdateSSRD == 'TRUE') {
      UpdateSSRD <- TRUE
    } else {
      UpdateSSRD <- FALSE
    }
  }
  if (length (run.args) > 4) {
    SimpleOnly <- run.args[5]
    if (SimpleOnly == 'y' || SimpleOnly == 'Y' || SimpleOnly == 'TRUE') {
      SimpleOnly <- TRUE
    } else {
      SimpleOnly <- FALSE
    }
  }
  if (length (run.args) > 5) {
    NSTEP <- as.numeric (run.args[6])
  }
  if (length (run.args) > 6) {
    GeneratePlots <- run.args[7]
    if (GeneratePlots == 'y' || GeneratePlots == 'Y' || GeneratePlots == 'TRUE') {
      GeneratePlots <- TRUE
    } else {
      GeneratePlots <- FALSE
    }
  }
} else {  ## This is the interactive section
  x <- readline (sprintf ("Project is %s; CR to accept or enter new project name: ", Project))
  if (nchar(x) > 1) {Project <- x}
  x <- readline (sprintf ("Flight is %d; CR to accept, number or 'ALL' or 'NEXT' for new flight name: ", Flight))
  if (x == 'ALL') {
    ALL <- TRUE
  } else if (x == 'NEXT') {
    Flight <- getNext(ProjectDir, Project)
  } else if (nchar(x) > 0 && !is.na(as.numeric(x))) {
    Flight <- as.numeric(x)
  }
  x <- readline (sprintf ("newAK is %s; CR to accept, Y or T to enable, N or F to disable: ", UpdateAKRD))
  if (nchar(x) > 0) 
    UpdateAKRD <- ifelse ((grepl('^T', x) || grepl('^Y', x)), TRUE, FALSE)
  x <- readline (sprintf ("newSS is %s; CR to accept, Y or T to enable, N or F to disable: ", UpdateSSRD))
  if (nchar(x) > 0) 
    UpdateSSRD <- ifelse ((grepl('^T', x) || grepl('^Y', x)), TRUE, FALSE)
  x <- readline (sprintf ("simple is %s; CR to accept, Y or T to enable, N or F to disable: ", SimpleOnly))
  if (nchar(x) > 0) 
    SimpleOnly <- ifelse ((grepl('^T', x) || grepl('^Y', x)), TRUE, FALSE)
  x <- readline (sprintf ("NSTEP is %d, CR to accept or number to change: ", NSTEP))
  if (nchar(x) > 0 && !is.na(as.numeric(x))) {
    NSTEP <- as.numeric(x)
  }
  x <- readline (sprintf ("Generate Plots is %s; CR to accept, Y or T to enable, N or F to disable: ", GeneratePlots))
  if (nchar(x) > 0) 
    GeneratePlots <- ifelse ((grepl('^T', x) || grepl('^Y', x)), TRUE, FALSE)
}

print (sprintf ('run controls:  Project: %s;  Flight: %d;  UpdateAKRD: %s;  UpdateSSRD: %s;  SimpleOnly: %s;  Time increment: %d',
                Project, Flight, UpdateAKRD, UpdateSSRD, SimpleOnly, NSTEP))

fname = sprintf("%s%s/%srf%02d.nc", Directory, ProjectDir, Project, Flight)
Cradeg <- pi/180
OmegaE <- StandardConstant ('Omega')  ## Earth's rotation rate
OmegaE <- 15*Cradeg/3600              ## better match to INS?
Ree <- 6378137                        ## for radii of curvature
Ecc <- 0.08181919

source ('chunks/RotationCorrection.R')
source ('chunks/STMFV.R')

## adjust GPS velocity components for GPS antenna location
LG <- -4.30 
MaxGap <- 1000
.span <- 25    

## use the calibration determined by comparison to GPS:
# cfa1 <- coef(fa1); cfa2 <- coef(fa2); cfa3 <- coef(fa3)
# Data$BLONGA <- cfa1[1] + cfa1[2] * Data$BLONGA
# Data$BLATA  <- cfa2[1] + cfa2[2] * Data$BLATA
# Data$BNORMA <- cfa3[1] + cfa3[2] * Data$BNORMA
## this loads the cal coefficients cfa1, cfa2, cfa3 determined in Tech Note:
load(file='./BodyAccCal.Rdata')

## These very small adjustments prevent gradual ramping during the flight.
# D1$BROLLR <- D1$BROLLR + 0.00026
# D1$BPITCHR <- D1$BPITCHR + 0.00026
## ----new-data, include=TRUE, cache=FALSE---------------------------------

SaveRData2 <- sprintf("%s2.Rdata", thisFileName)
if (ReloadData) {
  print ('reading netCDF file')
  source ('chunks/AcquireData.R')
  ## add ROC variable
  source ('chunks/ROC.R')
  save (list=c('Data', 'DL', 'dt', 'Rate', 'tau', 'VarList', 'VROC'), file=SaveRData2)
} else {
  load(file=SaveRData2)
}
D1 <- Data  ## make adjustments to a copy; avoid changing original
## adjustments:
source ('chunks/AdjustCal.R')

## GPS shifted vs INS in AdjustCal.R


## transform to the a-frame for comparison to the IRU:
VL <- matrix(c(D1$VEW, D1$VNS, D1$VSPD), ncol=3) 
LA <- matrix (c(D1$vedot, D1$vndot, -D1$vudot - D1$Grav), ncol=3) + RotationCorrection (D1, VL)
AA <- XformLA (D1, LA, .inverse=TRUE)
AA[,3] <- AA[,3] - D1$Grav
fa1 <- lm(D1$BLONGA ~ AA[, 1])
fa2 <- lm(D1$BLATA ~ AA[, 2])
fa3 <- lm(D1$BNORMA ~ AA[, 3])  
fm1 <- lm (D1$vedot ~ D1$LACCX)
fm2 <- lm (D1$vndot ~ D1$LACCY)
fm3 <- lm (D1$vudot ~ D1$LACCZ)

## ----Kalman-setup, include=TRUE, cache=CACHE-----------------------------

## the chunk is 'sourced' here so the same code can be used in KalmanFilter.R
# source ('chunks/Kalman-setup.R')
## Kalman-setup
## initialize matrices needed by the Kalman filter and load the starting-point
## for the error-state vector.

if (!SimpleOnly) {
  source ('chunks/Kalman-setup.R')
} 

## ----Kalman-loop, include=TRUE, cache=CACHE------------------------------

if (!SimpleOnly) {
  print (sprintf ('main loop is starting -- %s', Sys.time()))
  source ('chunks/Kalman-loop.R')
  print (sprintf ('main loop is done -- %s', Sys.time()))
} 

## ----plot-Kalman-position, include=TRUE, 

## Avoid some bad measurements at the start and end
## (e.g., from two-pass filtering of LAT/LON)
minT <- D1$Time[1]- as.integer (D1$Time[1]) %% 60 + 120
maxT <- D1$Time[nrow(D1)]- as.integer (D1$Time[nrow(D1)]) %% 60 - 120
r1 <- max(which(D1$Time == minT)[1], which(D1$TASX > 50)[1])
r2 <- which(D1$TASX > 50); r2 <- r2[length(r2)]
r2 <- min(which(D1$Time == maxT)[1], r2)
r <- r1:r2
DP <- D1[r, ]
if (GeneratePlots && !SimpleOnly) {
  print ("generating plots")
  png (file='KFplots/Position.png', width=900, height=600, res=150)
  print (ggplotWAC (with (DP,
                          data.frame(Time, DLONM, CLONM, DLATM, CLATM, DALT, CALT)),
                    ylab=expression (paste ('difference in position [km]')),
                    panels=3, 
                    labelL=c('KF-GPS', 'correction'),
                    labelP=c('east', 'north', 'up'),
                    legend.position=c(0.2,0.68), theme.version=1))
  invisible(dev.off())
  if (ShowPlots) {
    print (ggplotWAC (with (DP,
                            data.frame(Time, DLONM, CLONM, DLATM, CLATM, DALT, CALT)),
                      ylab=expression (paste ('difference in position [km]')),
                      panels=3, 
                      labelL=c('KF-GPS', 'correction'),
                      labelP=c('east', 'north', 'up'),
                      legend.position=c(0.2,0.68), theme.version=1))
  }
  pct <- as.integer(100/8)
  print (sprintf ('figures %d%% done', pct))
  sdZ <- sd(DP$DALT, na.rm=TRUE)
  sdvew <- sd(DP$GGVEW - DP$VEWKF, na.rm=TRUE)
  sdvns <- sd(DP$GGVNS - DP$VNSKF, na.rm=TRUE)
  meanvz <- mean(DP$ROCKF - DP[, VROC], na.rm=TRUE)
  sdvz <- sd(DP[, VROC] - DP$ROCKF, na.rm=TRUE)
  
  
  ## ----plot-Kalman-velocity
  
  png(file='KFplots/Velocity.png', width=900, height=600, res=150)
  print (ggplotWAC (with (DP,
                          data.frame(Time, DVEW, CVEW, DVNS, CVNS, DROC, CROC)),
                    ylab=expression (paste ('velocity component [m ',s^-1,']')),
                    panels=3, 
                    labelL=c('KF-GPS', 'correction'),
                    labelP=c('east', 'north', 'up'),
                    legend.position=c(0.2, 0.97), theme.version=1))
  invisible(dev.off())
  if (ShowPlots) {
    print (ggplotWAC (with (DP,
                            data.frame(Time, DVEW, CVEW, DVNS, CVNS, DROC, CROC)),
                      ylab=expression (paste ('velocity component [m ',s^-1,']')),
                      panels=3, 
                      labelL=c('KF-GPS', 'correction'),
                      labelP=c('east', 'north', 'up'),
                      legend.position=c(0.2, 0.97), theme.version=1))
  }
  pct <- as.integer(100*2/8)
  print (sprintf ('figures %d%% done', pct))
  
  
  ## ----errors-in-pitch,
  ## must construct using many of the elements of ggplotWAC but in this order:
  d1 <- with(DP, data.frame(Time, CPL, CPLF, CRL, CRLF))
  lines_per_panel <- 2; panels <- 2
  labelL=c('KF value', 'smoothed')
  labelP=c('pitch', 'roll')
  DLP <- nrow(d1)
  VarGroup <- rep (gl (lines_per_panel, DLP, labels=labelL), panels)
  PanelGroup <- gl (panels, lines_per_panel*DLP, labels=labelP)
  dd <- data.frame(reshape2::melt(d1, 1), VarGroup, PanelGroup)
  lvl <- levels(dd$VarGroup)
  d2 <- with(DP, data.frame(Time, "plo"=CPLF-4*SDCPL, "rlo"=CRLF-4*SDCRL,
                            "phi"=CPLF+4*SDCPL, "rhi"=CRLF+4*SDCRL))
  ## note the required order below:
  de <- data.frame (reshape2::melt(d2, 1, c(2,4,3,5)), VarGroup, PanelGroup)
  g7 <- ggplot(dd, aes(x=Time)) + facet_grid(PanelGroup ~ .) 
  g7 <- g7 + geom_ribbon(aes(x=Time, ymin=value, ymax=value), data=de, alpha=0.25)  ##alpha=.25 
  g7 <- g7 + geom_path(aes(x=Time, y=value, colour=VarGroup, linetype=VarGroup, size=VarGroup), data=dd)
  g7 <- g7 + ylim(-0.05,0.05)
  g7 <- g7 + scale_linetype_manual ('', labels=lvl, breaks=lvl, values = c(1,1))
  g7 <- g7 + scale_colour_manual('', labels = lvl, breaks=lvl, values = c('blue', 'red'))
  g7 <- g7 + scale_size_manual('', labels=lvl, breaks=lvl, values = c(0.5, 1.5))
  g7 <- g7 + theme_WAC(1) + theme (legend.position=c(0.5,0.95))
  g7 <- g7 + theme(axis.text.x = element_text (size=11.5, margin=margin(15,0,0,0)))
  g7 <- g7 + theme(axis.title.x = element_text (size=12))
  g7 <- g7 + labs (x='Time [UTC]', y=expression (paste ('l-frame error [',degree,']')))
  ## I'm not sure why this is necessary; printing g7 in the usual way did not leave
  ## the ribbon visible, apparently because the resulting pdf file did not include it.
  ## Now the plot is generated here, but placed in the document in the LaTeX code
  ## preceding this chunk. This out-of-order sequence requires run
  # png (filename='Fig7.png', width=600, height=480, res=150)
  png(file='KFplots/AAlframe.png', width=900, height=600, res=150)
  print (g7)
  invisible(dev.off())
  if (ShowPlots) {print (g7)}
  # invisible (dev.off())
  pct <- as.integer(100*3/8)
  print (sprintf ('figures %d%% done', pct))
  
  
  ## ----a-frame-errors
  
  PC <- CorrectPitch(DP, .span=901)
  DP$PC <- PC[, 1]
  DP$RC <- PC[, 2]
  
  
  d1 <- with(DP, data.frame(Time, CPAF, PC, CRAF, RC))
  lines_per_panel <- 2; panels <- 2
  labelL=c('smoothed KF value', 'PC')
  labelP=c('pitch', 'roll')
  DLP <- nrow(d1)
  VarGroup <- rep (gl (lines_per_panel, DLP, labels=labelL), panels)
  PanelGroup <- gl (panels, lines_per_panel*DLP, labels=labelP)
  dd <- data.frame(reshape2::melt(d1, 1), VarGroup, PanelGroup)
  lvl <- levels(dd$VarGroup)
  d2 <- with(DP, data.frame(Time, "plo"=CPAF-4*SDCPA, "rlo"=CRAF-4*SDCRA,
                            "phi"=CPAF+4*SDCPA, "rhi"=CRAF+4*SDCRA))
  ## note the required order below:
  de <- data.frame (reshape2::melt(d2, 1, c(2,4,3,5)), VarGroup, PanelGroup)
  g9 <- ggplot(dd, aes(x=Time)) + facet_grid(PanelGroup ~ .) 
  g9 <- g9 + geom_ribbon(aes(x=Time, ymin=value, ymax=value), data=de, alpha=0.25) 
  g9 <- g9 + geom_line(aes(x=Time, y=value, colour=VarGroup, linetype=VarGroup), data=dd)
  g9 <- g9 + ylim(-0.05,0.05)
  g9 <- g9 + scale_linetype_manual ('', labels=lvl, breaks=lvl, values = c(1,1))
  g9 <- g9 + scale_colour_manual('', labels = lvl, breaks=lvl, values = c('blue', 'red'))
  g9 <- g9 + theme_WAC(1) + theme (legend.position=c(0.5,0.95))
  g9 <- g9 + theme(axis.text.x = element_text (size=11.5, margin=margin(15,0,0,0)))
  g9 <- g9 + theme(axis.title.x = element_text (size=12))
  g9 <- g9 + labs (x='Time [UTC]', y=expression (paste ('a-frame error [',degree,']')))
  ## I'm not sure why this is necessary; printing g9 in the usual way did not leave
  ## the ribbon visible, apparently because the resulting pdf file did not include it.
  ## Now the plot is generated here, but placed in the document in the LaTeX code
  ## preceding this chunk. This out-of-order sequence requires run
  # png (filename='Fig9.png', width=600, height=480, res=150)
  png(file='KFplots/AAaframe.png', width=900, height=600, res=150)
  print (g9)
  invisible(dev.off())
  if (ShowPlots) {print (g9)}
  # invisible (dev.off())
  pct <- as.integer(100*4/8)
  print (sprintf ('figures %d%% done', pct))
  
  DP$HC <- -CorrectHeading (DP, .plotfile='KFplots/HCPlot.pdf')
  pct <- as.integer(100*5/8)
  print (sprintf ('figures %d%% done', pct))
}

if (SimpleOnly) {
  PC <- CorrectPitch(DP, .span=901)
  DP$PC <- PC[, 1]
  DP$RC <- PC[, 2]
  DP$HC <- -CorrectHeading (DP, .plotfile='KFplots/HCPlot.pdf')
} else {
  DP$dP <- DP$deltaPsi/Cradeg
  HA <- with(DP, sqrt(LACCX^2+LACCY^2))
  DP$dP[HA < 1] <- NA
  sddP <- sd(DP$dP, na.rm=TRUE)
}

## ----plot-Kalman-heading

# plotWAC(subset(DP,,c(Time, dP)))
if (GeneratePlots && !SimpleOnly) {
  png(file='KFplots/HDG.png', width=900, height=600, res=150)
  grid.newpage()
  suppressWarnings (
    with(DP,
         ggplotWAC (data.frame(Time, dP),
                    ylim=c(-0.4,0.4),
                    ylab=expression (paste ('error in heading [',degree,']')),
                    legend.position=NA, position=c(2,2), theme.version=1
         )
    )
  )
  DP$dP <- DP$deltaPsi/Cradeg
  HA <- with(DP, sqrt(LACCX^2+LACCY^2))
  DP$dP[HA < 1] <- NA
  DP$SDH <- sqrt(VCor[r,9]) / (Cradeg * 30)
  DP$CCTHDG <- SmoothInterp (DP$CTHDG, .Length=181)
  DP$CCCTHDG <- DP$CCTHDG
  DP$CCCTHDG[DP$SDH > 0.02] <- NA
  suppressWarnings (
    with (DP,
          ggplotWAC(data.frame(Time, 'Kalman'=CCTHDG, 'Ranadu'=HC, 'spline'=HCS, 'KF'=CCCTHDG), 
                    ylim=c(-0.4,0.4), lwd=c(0.7,1,1,1.5), lty=c(3,1,1,1),
                    ylab=expression (paste ('error in heading [',degree,']')),
                    col=c('blue', 'forestgreen', 'red', 'blue'), position=c(1,2),
                    legend.position=c(0.5,0.95), theme.version=1
          )
    )
  )
  invisible(dev.off())
  DP$dP <- DP$deltaPsi/Cradeg
  HA <- with(DP, sqrt(LACCX^2+LACCY^2))
  DP$dP[HA < 1] <- NA
  sddP <- sd(DP$dP, na.rm=TRUE)
  if (ShowPlots) {
    grid.newpage()
    suppressWarnings (
      with(DP,
           ggplotWAC (data.frame(Time, dP),
                      ylim=c(-0.4,0.4),
                      ylab=expression (paste ('error in heading [',degree,']')),
                      legend.position=NA, position=c(2,2), theme.version=1
           )
      )
    )
    DP$dP <- DP$deltaPsi/Cradeg
    HA <- with(DP, sqrt(LACCX^2+LACCY^2))
    DP$dP[HA < 1] <- NA
    DP$SDH <- sqrt(VCor[r,9]) / (Cradeg * 30)
    DP$CCTHDG <- SmoothInterp (DP$CTHDG, .Length=181)
    DP$CCCTHDG <- DP$CCTHDG
    DP$CCCTHDG[DP$SDH > 0.02] <- NA
    suppressWarnings (
      with (DP,
            ggplotWAC(data.frame(Time, 'Kalman'=CCTHDG, 'Ranadu'=HC, 'spline'=HCS, 'KF'=CCCTHDG), 
                      ylim=c(-0.4,0.4), lwd=c(0.7,1,1,1.5), lty=c(3,1,1,1),
                      ylab=expression (paste ('error in heading [',degree,']')),
                      col=c('blue', 'forestgreen', 'red', 'blue'), position=c(1,2),
                      legend.position=c(0.5,0.95), theme.version=1
            )
      )
    )
  }
  pct <- as.integer(100*6/8)
  print (sprintf ('figures %d%% done', pct))
}

## ----new-AKRD
# load (file='~/RStudio/Reprocessing/AKRD-fit-coef.Rdata')
cffn <- 19.70547
cff <- 21.481
cfs <- c(4.525341674, 19.933222011, -0.001960992)
if (SimpleOnly) {
  PC <- CorrectPitch(D1, .span=901)
  D1$PC <- PC[, 1]
  D1$RC <- PC[, 2]
  D1$HC <- CorrectHeading (D1)
}
D1$QR <- D1$ADIFR / D1$QCF
D1$QR[D1$QCF < 20] <- NA
D1$QR[is.infinite(D1$QR)] <- NA
CutoffFreq <- 600 * Rate
D1$QRS <- zoo::na.approx (as.vector(D1$QR), maxgap=1000*Rate, na.rm = FALSE)
D1$QRS[is.na(D1$QRS)] <- 0
D1$QRS <- signal::filtfilt (signal::butter (3, 2/CutoffFreq), D1$QRS)
D1$QRF <-  D1$QR - D1$QRS
D1$QCFS <- zoo::na.approx (as.vector(D1$QCF), maxgap=1000*Rate, na.rm = FALSE)
D1$QCFS[is.na(D1$QCFS)] <- 0
D1$QCFS <- signal::filtfilt (signal::butter (3, 2/CutoffFreq), D1$QCFS)
D1$AKKF <- cff * D1$QRF + cfs[1] + cfs[2] * D1$QRS + cfs[3] * D1$QCFS
D1$AKKF[D1$QCF < 10] <- NA
D1$SSKF <- 0.008 + 22.301 * D1$BDIFR / D1$QCF
D1$SSKF[D1$QCF < 10] <- NA
if (UpdateAKRD) {
  D1$ATTACK <- D1$AKRD <- D1$AKKF
}
if (UpdateSSRD) {
  D1$SSLIP <- D1$SSRD <- D1$SSKF
}

## ----wind-calculation, include=TRUE, cache=CACHE-------------------------

if (SimpleOnly) {
  ## get the wind variables:
  DataW <- D1
  DataN <- WindProcessor (DataW)
  D1$WICC <- DataN$WIN
  D1$WDCC <- DataN$WDN
  D1$WSCC <- DataN$WSN
  DataW$PITCH <- D1$PITCH - PC[,1]
  DataW$ROLL <- D1$ROLL - PC[,2]
  DataW$THDG <- D1$THDG - HC
  DataW$VEW <- D1$VEWC
  DataW$VNS <- D1$VNSC
  DataW$GGVSPD <- D1$ROC
  DataN <- WindProcessor(DataW, LG=0, CompF=TRUE)    ## suppress GPS lever arm)
  D1$WDSC <- DataN$WDN
  D1$WSSC <- DataN$WSN
  D1$WISC <- DataN$WIN
  if (GeneratePlots) {
    png(file='KFplots/AAaframe.png', width=900, height=600, res=150)
    with (D1[r,], plotWAC (data.frame (Time, PC[r,1], PC[r,2])))
    invisible(dev.off())
    if (ShowPlots) {
      with (D1[r,], plotWAC (data.frame (Time, PC[r,1], PC[r,2])))
    }
  }
} else {
  source ('chunks/wind-calculation.R')
}

if (GeneratePlots) {
  if (SimpleOnly) {
    png(file='KFplots/Wind.png', width=900, height=600, res=150)
    with (D1[r, ], plotWAC (data.frame(Time, WICC, WISC, WICC-WISC), ylim=c(-5,5)))
    invisible (dev.off())
    if (ShowPlots) {
      with (D1[r, ], plotWAC (data.frame(Time, WICC, WISC, WICC-WISC), ylim=c(-5,5)))
    }
    png(file='KFplots/Wind2.png', width=900, height=600, res=150)
    with (D1[r, ], plotWAC (data.frame (Time, WSCC-WSSC, WDCC-WDSC), ylim=c(-10,10)))
    invisible (dev.off())
    if (ShowPlots) {
      with (D1[r, ], plotWAC (data.frame (Time, WSCC-WSSC, WDCC-WDSC), ylim=c(-10,10)))
    }
  } else {
    png(file='KFplots/Wind.png', width=900, height=600, res=150)
    print (ggplotWAC(
      with(D1[r, ], 
           data.frame (Time, 'Kalman'=WIKF, 'original'=WICC, "DWI"=WIKF-WICC, 'SKIP'=WIKF*0)),
      col=c('blue', 'red'), lwd=c(1.4,0.8,1,0), lty=c(1,42),
      ylab=expression(paste('vertical wind [m ', s^-1, ']')),
      panels=2,
      labelL=c('Kalman', 'original'),
      labelP=c('variables', 'difference'),
      legend.position=c(0.8,0.94), theme.version=1)
    )
    invisible(dev.off())
    pct <- as.integer(100*7/8)
    print (sprintf ('figures %d%% done', pct))
    ## ----hw-plot, include=TRUE,
    png(file='KFplots/Wind2.png', width=900, height=600, res=150)
    print (ggplotWAC(
      with(D1[r, ], 
           data.frame (Time, 'direction'=WDKF-WDCC, 'speed'=WSKF-WSCC)),
      ylab=expression(paste('KF correction [', degree, ' (top) or m ', s^-1, ' (bottom)]')),
      panels=2,
      labelL=c('Kalman correction'),
      labelP=c('direction', 'speed'),
      legend.position=c(0.8,0.94), theme.version=1)
    )
    invisible(dev.off())
  }
  if (ShowPlots) {
    print (ggplotWAC(
      with(D1[r, ], 
           data.frame (Time, 'Kalman'=WIKF, 'original'=WICC, "DWI"=WIKF-WICC, 'SKIP'=WIKF*0)),
      col=c('blue', 'red'), lwd=c(1.4,0.8,1,0), lty=c(1,42),
      ylab=expression(paste('vertical wind [m ', s^-1, ']')),
      panels=2,
      labelL=c('Kalman', 'original'),
      labelP=c('variables', 'difference'),
      legend.position=c(0.8,0.94), theme.version=1)
    )
    print (ggplotWAC(
      with(D1[r, ], 
           data.frame (Time, 'direction'=WDKF-WDCC, 'speed'=WSKF-WSCC)),
      ylab=expression(paste('KF correction [', degree, ' (top) or m ', s^-1, ' (bottom)]')),
      panels=2,
      labelL=c('Kalman correction'),
      labelP=c('direction', 'speed'),
      legend.position=c(0.8,0.94), theme.version=1)
    )
  }
}
if (file.exists ('KFplots/HCPlot.pdf')) {system ('convert KFplots/HCPlot.pdf KFplots/HCPlot.png')}
print (sprintf ('plots generated -- %s', Sys.time()))

## ----create-new-netcdf, cache=FALSE--------------------------------------

print (sprintf("making new netCDF file -- %s", Sys.time()))
source ('chunks/create-new-netcdf.R')

## ----modify-new-netcdf, include=TRUE-------------------------------------

source ('chunks/modify-new-netcdf.R')
## modify-new-netcdf.R
print (sprintf ("Kalman Processor is finished -- %s", Sys.time()))

