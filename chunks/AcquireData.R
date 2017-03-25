## AcquireData.R -- read data file and do preliminary processing
## Assumes fname is defined; 

## test if body rotations/accelerations are available
FI <- DataFileInfo (fname)  ## this is slow because it finds lat/long range...
Rate <- FI$Rate
dt <- 1/Rate
## don't use the new AKRD algorithm for flights before 2012:
if (FI$Start < as.POSIXlt ("2012-01-01", tz='UTC')) {updateAKRD <- FALSE}
NeedIRU <- TRUE
if (('BLATA' %in% FI$Variables) && ('BYAWR' %in% FI$Variables)) {
  NeedIRU <- FALSE
}
## SPECIAL TEST ######
# NeedIRU <- TRUE
VROC <- 'GGVSPD'
if (VROC %in% FI$Variables) {
} else {
  VROC <- 'VSPD_G'
  if (VROC %in% FI$Variables) {
  } else {
    VROC <- 'VSPD_A'
    if (VROC %in% FI$Variables) {
    } else {
      print ('required variable for GPS ROC not found; exiting')
      exit ()
    }
  }
}

if (NeedIRU) {
  VarList <- c('LAT', 'LON', 'ALT', 'GGLAT', 'GGLON', 'GGALT',
               'VEW', 'VNS', 'VSPD', 'GGVEW', 'GGVNS', VROC,
               'PITCH', 'ROLL', 'THDG', 'TASX', 'ATTACK', 'SSLIP', 
               'PSXC', 'ATX', 'ACINS', 'ADIFR', 'QCF', 'BDIFR')
  Data <- getNetCDF (fname, VarList)
  Data$Grav <- Gravity (Data$LAT, Data$GGALT)
  Data$Rn <- Ree / (1 - (Ecc*sin(Data$GGLAT*Cradeg))^2)^0.5 + Data$GGALT
  Data$Rm <- Data$Rn * (1-Ecc^2) / (1-(Ecc*sin(Data$GGLAT*Cradeg))^2) + Data$GGALT
  source ('chunks/AddIRUVariables.R')
  VarList <- c(VarList, 'BPITCHR', 'BROLLR', 'BYAWR', 'BLATA', 'BLONGA', 'BNORMA')
} else {
  VarList <- c('LAT', 'LON', 'ALT', 'GGLAT', 'GGLON', 'GGALT',
               'VEW', 'VNS', 'VSPD', 'GGVEW', 'GGVNS', VROC,
               'PITCH', 'ROLL', 'THDG', 'TASX', 'ATTACK', 'SSLIP', 
               'PSXC', 'ATX', 'ACINS', 'ADIFR', 'QCF', 'BDIFR',
               'BLATA', 'BLONGA', 'BNORMA', 'BPITCHR', 'BROLLR', 'BYAWR')
  Data <- getNetCDF (fname, VarList)
  Data$Grav <- Gravity (Data$LAT, Data$GGALT)
  Data$Rn <- Ree / (1 - (Ecc*sin(Data$GGLAT*Cradeg))^2)^0.5 + Data$GGALT
  Data$Rm <- Data$Rn * (1-Ecc^2) / (1-(Ecc*sin(Data$GGLAT*Cradeg))^2) + Data$GGALT
}

Data$Grav <- Gravity (Data$LAT, Data$GGALT)
## remove heading adjustments as added during initial processing
Z <- data.frame(getAttributes(Data$THDG, .print=FALSE))
if ('CalibrationCoefficients' %in% names(Z)) {
  THDGoffset <- Z$CalibrationCoefficients[1]
} else {
  THDGoffset <- 0
}
## re signs: lag has been imposed, so use -ve sign to remove
if ('TimeLag' %in% names(Z)) {
  Data$THDG <- ShiftInTime(Data$THDG, -Z$TimeLag[1])
}
Data$THDG <- (Data$THDG - THDGoffset) %% 360
DL <- nrow(Data)

#interpolate if necessary: otherwise later filters fail
MaxGap <- 1000
for (V in VarList) {
  Data[, V] <- zoo::na.approx (as.vector (Data[, V]), maxgap=MaxGap, na.rm=FALSE)
}

## now add the l-frame accelerations to the data.frame
AB <- matrix(c(Data$BLONGA, Data$BLATA, Data$BNORMA+Data$Grav), ncol=3) #a-frame 
AL <- XformLA (Data, AB)                                    #l-frame
## now corrected for angular effects
## See Noureldin et al, 2013, Eq. (5.55)
VL <- matrix(c(Data$VEW, Data$VNS, Data$VSPD), ncol=3)
AL <- AL - RotationCorrection (Data, VL)  ##### check this sign and prev call above

## the resulting l-frame accelerations
Data$LACCX <- AL[, 1]
Data$LACCY <- AL[, 2]
Data$LACCZ <- AL[, 3] + Data$Grav
Data$LACCZ <- -Data$LACCZ

## smooth to match GPS-velocity derivatives

Data$LACCX <- zoo::na.approx (as.vector(Data$LACCX), maxgap=1000, na.rm=FALSE)
Data$LACCY <- zoo::na.approx (as.vector(Data$LACCY), maxgap=1000, na.rm=FALSE)
Data$LACCZ <- zoo::na.approx (as.vector(Data$LACCZ), maxgap=1000, na.rm=FALSE)
Data$LACCX[is.na(Data$LACCX)] <- 0
Data$LACCY[is.na(Data$LACCY)] <- 0
Data$LACCZ[is.na(Data$LACCZ)] <- 0
.span <- 10*Rate+1
Data$LACCX <- signal::sgolayfilt (Data$LACCX, 3, .span, m=0)
Data$LACCY <- signal::sgolayfilt (Data$LACCY, 3, .span, m=0)
Data$LACCZ <- signal::sgolayfilt (Data$LACCZ, 3, .span, m=0)

## get the pitch and roll in the l-frame:
.thdg <- Data$THDG * Cradeg
Data$PITCHL <- Data$PITCH * cos (.thdg) + Data$ROLL * sin (.thdg)
Data$ROLLL <-  -Data$PITCH * sin (.thdg) + Data$ROLL * cos (.thdg)

