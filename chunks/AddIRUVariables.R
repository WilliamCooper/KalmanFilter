## AddIRUVariables.R
## Assumes that Ree, Ecc, Rate, Cradeg have been set appropriately and
## data.frame containing archive data is named Data. Adds the IRU pseudo-
## acceleration obtained by differentiating l-frame aircraft velocities and
## transforming to the a-frame. Also differentiates attitude angles and
## transforms the result to the a-frame.

## smooth the measurements when determining derivatives
.span <- 10 * Rate + 1    
## The following are accelerations determined from derivatives of the INS andGPS velocities.
## These should match the measured accelerations after transformation to the l-frame
## and application of the rotation correction:
D1 <- Data
D1$VNS <- zoo::na.approx (as.vector(D1$VNS), maxgap=1000, na.rm=FALSE)
D1$VEW <- zoo::na.approx (as.vector(D1$VEW), maxgap=1000, na.rm=FALSE)
D1$VSPD <- zoo::na.approx (as.vector(D1$VSPD), maxgap=1000, na.rm=FALSE)
D1$vndot <- signal::sgolayfilt (D1$VNS, 3, .span, m=1) * Rate  # m=1 for first deriv.
D1$vedot <- signal::sgolayfilt (D1$VEW, 3, .span, m=1) * Rate
D1$vudot <- signal::sgolayfilt (D1$VSPD, 3, .span, m=1) * Rate
## transform to the a-frame for comparison to the IRU:
G <- D1$Grav
VL <- matrix(c(D1$VEW, D1$VNS, D1$VSPD), ncol=3)
LA <- matrix (c(D1$vedot, D1$vndot, -D1$vudot - G), ncol=3) + RotationCorrection (D1, VL) 
AA <- XformLA (D1, LA, .inverse=TRUE)
AA[,3] <- AA[,3] - G
Data$BLONGA <- AA[, 1]
Data$BLATA  <- AA[, 2]
Data$BNORMA <- AA[, 3]

D1$PITCH <- zoo::na.approx (as.vector(D1$PITCH), maxgap=1000, na.rm=FALSE)
D1$ROLL <- zoo::na.approx (as.vector(D1$ROLL), maxgap=1000, na.rm=FALSE)
D1$pdot <- signal::sgolayfilt (D1$PITCH, 3, .span, m=1) *  Rate
D1$rdot <- signal::sgolayfilt (D1$ROLL, 3, .span, m=1) *  Rate
## handle heading differently to avoid wrap-around problems
D1$hdot <- c(0, diff(D1$THDG))
D1$hdot <- (D1$hdot + 540) %% 360 - 180
D1$hdot <- D1$hdot * Rate
D1$hdot <- SmoothInterp (D1$hdot, .Length=10*Rate+1)
shift <- -500/Rate    ## shift by 1/2 the time interval
D1$hdot <- ShiftInTime(D1$hdot, Rate, shift)
## the derivative of the transformation matrix in terms of hdot etc.
with (D1, {
  ch <- cos(THDG*Cradeg)
  sh <- sin(THDG*Cradeg)
  cp <- cos(PITCH*Cradeg)
  sp <- sin(PITCH*Cradeg)
  cr <- cos(ROLL*Cradeg)
  sr <- sin(ROLL*Cradeg)
  Rd <<- c(hdot*ch*cp-pdot*sh*sp, 
           hdot*(ch*sp*sr-sh*cr)  + pdot*sh*cp*sr    + rdot*(sh*sp*cr-ch*sr),
           hdot*(-sh*sr-ch*sp*cr) + pdot*(-sh*cp*cr) + rdot*(ch*cr+sh*sp*sr),
           hdot*(-sh*cp)          + pdot*(-ch*sp),
           hdot*(-sh*sp*sr-ch*cr) + pdot*ch*cp*sr    + rdot*(ch*sp*cr+sh*sr),
           hdot*(sh*sp*cr-ch*sr)  + pdot*(-ch*cp*cr) + rdot*(ch*sp*sr-sh*cr),
           - pdot*cp,
           pdot*(-sp*sr)    + rdot*cp*cr,
           pdot*sp*cr       + rdot*cp*sr)
})
RdM <- aperm (array (Rd, dim=c(nrow(D1),3,3)))
OmegaA <- array (dim=c(nrow(D1),3,3))
NR <- nrow(D1)
IprogLast <- 0
for (i in 1:nrow(D1)) {
  Iprog <- as.integer(i * 100 / NR)
  if ((Iprog %% 5) == 0 && Iprog != IprogLast) {
    print (sprintf ('IRU loop %d%% done', Iprog))
    IprogLast <- Iprog
  }
  sv <- with(D1[i,], c(LAT, LON, ALT, VEW, VNS, VSPD, PITCH, ROLL, THDG)) 
  rlm <- XformLA (data.frame(PITCH=sv[7], ROLL=sv[8], THDG=sv[9]))
  sv[c(1:2,7:12)] <- sv[c(1:2, 7:12)] * Cradeg
  omega <- as.vector (c(-sv[5] / D1$Rm[i], 
                        OmegaE*cos(sv[1])+sv[4]/(D1$Rn[i]),
                        OmegaE*sin(sv[1])+sv[4]*tan(sv[1])/D1$Rn[i]), mode='numeric')
  ## signs account for b-frame to a-frame: reverse sign of 3, swap 1 and 2
  Oill <- aperm(matrix (c(0, omega[3], omega[2], -omega[3], 0, -omega[1], -omega[2], omega[1], 0), ncol=3))
  Oilb <- Oill %*% rlm
  OmegaA[i,,] <- t(rlm) %*% (RdM[,,i] - Oilb/Cradeg)
}
Data$BPITCHR <- -OmegaA[,1,3]
Data$BROLLR <-   OmegaA[,2,3]
Data$BYAWR  <-  -OmegaA[,1,2]
rm (D1)

