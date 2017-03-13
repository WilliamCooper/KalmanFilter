# wind-calculation.R

## get the wind variables:
DataW <- D1
DataN <- WindProcessor (DataW)
D1$WICC <- DataN$WIN
D1$WDCC <- DataN$WDN
D1$WSCC <- DataN$WSN
DataW$PITCH <- D1$PITCHKF
DataW$ROLL <- D1$ROLLKF
DataW$THDG <- D1$THDGKF
DataW$VEW <- D1$VEWKF
DataW$VNS <- D1$VNSKF
DataW$GGVSPD <- D1$ROCKF
DataN <- WindProcessor(DataW, LG=0, CompF=FALSE)    ## suppress comp filter and GPS lever arm)
D1$WDKF <- DataN$WDN
D1$WSKF <- DataN$WSN
D1$WIKF <- DataN$WIN
## add longitudinal and lateral components analogous to UXC and VYC:
.hdg <- D1$THDGKF * Cradeg
.wd <- D1$WDKF * Cradeg + pi
D1$UXKF <- D1$WSKF * (sin(.hdg)*sin(.wd) + cos(.hdg)*cos(.wd))
.hdg <- .hdg - pi/2
D1$VYKF <- D1$WSKF * (sin(.hdg)*sin(.wd) + cos(.hdg)*cos(.wd))
DataW$GGVSPD <- D1[, VROC]
DataN <- WindProcessor(DataW, CompF=FALSE)    ## suppress comp filter and GPS lever arm)
D1$WIKFG <- DataN$WIN
sdWIdif <- sd(D1$WIKF[r]-D1$WICC[r], na.rm=TRUE)

