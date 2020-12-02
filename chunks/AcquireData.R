## AcquireData.R -- read data file and do preliminary processing
## Assumes fname is defined; 

## test if body rotations/accelerations are available
FI <- DataFileInfo (fname, LLrange=FALSE)  
Rate <- FI$Rate
dt <- 1/Rate
## don't use the new AKRD algorithm for flights before 2012:
if (FI$Start < as.POSIXlt ("2012-01-01", tz='UTC')) {updateAKRD <- FALSE}
NeedIRU <- TRUE
if (('BLATA' %in% FI$Variables) && ('BYAWR' %in% FI$Variables)) {
  NeedIRU <- FALSE
}

## Allow old-style GPS rates-of-climb
VROC <- 'GGVSPD'
if (VROC %in% FI$Variables) {
} else {
  VROC <- 'GGVSPDB'
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
}

removeOffsets <- function(Data, VarList) {
  ## Remove offsets if imposed during processing:
  for (v in VarList) {
    CC <- attr (Data[, v], 'CalibrationCoefficients')
    if (!is.null(CC[1])) {
      Data[, v] <- (Data[, v] - CC[1]) / CC[2]
      print (sprintf ('adjusted %s to remove cal coeff %.2f', v, CC[1]))
      attr(Data[, v], 'CalibrationCoefficients') <- c(0, 1)  ## Avoid later removing this again
    }
    ## Also list but do not remove time shifts:
    TimeLag <- attr(Data[, v], 'TimeLag')
    if (!is.null(TimeLag)) {
      # Data[, v] <- ShiftInTime (Data[, v], Rate, -as.integer(TimeLag))
      print (sprintf ('time lag is %d ms for %s', as.integer(TimeLag), v))
    }
    if (v == 'THDG') {
      Data$THDG <- Data$THDG %% 360
    }
  }
  return(Data)
}

if (NeedIRU) {
  VarList <- c('LAT', 'LON', 'ALT', 'GGLAT', 'GGLON', 'GGALT',
               'VEW', 'VNS', 'VSPD', 'GGVEW', 'GGVNS', VROC,
               'PITCH', 'ROLL', 'THDG', 'TASX', 'ATTACK', 'SSLIP', 
               'PSXC', 'ATX', 'ACINS', 'ADIFR', 'QCF', 'BDIFR')
  if ('BLONGA' %in% FI$Variables) { ## Allow for case where body acc present, rots absent.
    VarList <- c(VarList, 'BLATA', 'BLONGA', 'BNORMA')
  }
  Data <- getNetCDF (fname, VarList)
  ## Remove offsets if imposed during processing:
  Data <- removeOffsets (Data, VarList)
  Data$Grav <- Gravity (Data$LAT, Data$GGALT)
  Data$Rn <- Ree / (1 - (Ecc*sin(Data$GGLAT*Cradeg))^2)^0.5 + Data$GGALT
  Data$Rm <- Data$Rn * (1-Ecc^2) / (1-(Ecc*sin(Data$GGLAT*Cradeg))^2) + Data$GGALT
  source ('chunks/AddIRUVariables.R')
  VarList <- unique (c (VarList, 'BPITCHR', 'BROLLR', 'BYAWR', 'BLATA', 'BLONGA', 'BNORMA'))
} else {
  VarList <- c('LAT', 'LON', 'ALT', 'GGLAT', 'GGLON', 'GGALT',
               'VEW', 'VNS', 'VSPD', 'GGVEW', 'GGVNS', VROC,
               'PITCH', 'ROLL', 'THDG', 'TASX', 'ATTACK', 'SSLIP', 
               'PSXC', 'ATX', 'ACINS', 'ADIFR', 'QCF', 'BDIFR',
               'BLATA', 'BLONGA', 'BNORMA', 'BPITCHR', 'BROLLR', 'BYAWR')
  Data <- getNetCDF (fname, VarList)
  removeOffsets (Data, VarList)
  Data$Grav <- Gravity (Data$LAT, Data$GGALT)
  Data$Rn <- Ree / (1 - (Ecc*sin(Data$GGLAT*Cradeg))^2)^0.5 + Data$GGALT
  Data$Rm <- Data$Rn * (1-Ecc^2) / (1-(Ecc*sin(Data$GGLAT*Cradeg))^2) + Data$GGALT
}

Data$Grav <- Gravity (Data$LAT, Data$GGALT)
## remove heading adjustments as added during initial processing -- commented now, see above
# Z <- data.frame(getAttributes(Data$THDG, .print=FALSE))
# if ('CalibrationCoefficients' %in% names(Z)) {
#   THDGoffset <- Z$CalibrationCoefficients[1]
#   print (sprintf ('removing %.2f offset from heading', THDGoffset))
# } else {
#   THDGoffset <- 0
# }
# ## re signs: lag has been imposed, so use -ve sign to remove
# if ('TimeLag' %in% names(Z)) {
#   Data$THDG <- ShiftInTime(Data$THDG, .rate=Rate, .shift=-Z$TimeLag[1])
# }
# Data$THDG <- (Data$THDG - THDGoffset) %% 360
DL <- nrow(Data)

#interpolate if necessary: otherwise later filters fail
MaxGap <- 1000
DSave <- Data
for (V in VarList) {
  Data[, V] <- zoo::na.approx (as.vector (Data[, V]), maxgap=MaxGap, na.rm=FALSE, rule=2)
  Data[is.na(Data[, V]), V] <- 0
}
Data <- transferAttributes(DSave, Data)

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

Data$LACCX <- zoo::na.approx (as.vector(Data$LACCX), maxgap=MaxGap, na.rm=FALSE, rule=2)
Data$LACCY <- zoo::na.approx (as.vector(Data$LACCY), maxgap=MaxGap, na.rm=FALSE, rule=2)
Data$LACCZ <- zoo::na.approx (as.vector(Data$LACCZ), maxgap=MaxGap, na.rm=FALSE, rule=2)
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

