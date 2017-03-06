## AdjustCal.R
## apply adjustments and calibrations to the IRU measurements, and time shifts
## to some other variables. Also add the l-frame accelerations determined by
## differentiation of the GPS measurements

## adjustments:
D1$BYAWR <- (D1$BYAWR - 0.0003) * 1.008
D1$BROLLR <- D1$BROLLR + 0.0005
D1$BPITCHR <- D1$BPITCHR + 0.0003
D1$BNORMA <- D1$BNORMA*1.0086 - 0.0083
D1$BLONGA <- D1$BLONGA * 1.0062 - 0.0002
D1$BLATA <- D1$BLATA * 1.0604 + 0.0030
## consider the possible effects of misalignment of INS wrt AC axes
VROT <- array(c(D1$BPITCHR, D1$BROLLR, D1$BYAWR), dim=c(nrow(D1), 3))
SX <- data.frame(ROLL=rep(0.05,nrow(D1)), PITCH=rep(0.15,nrow(D1)), THDG=rep(0,nrow(D1)))
VROTA <- XformLA(SX, VROT)
# D1$BPITCHR <- VROTA[,2]
# D1$BROLLR <- VROTA[,1]
# D1$BYAWR <- -VROTA[,3]
rm (VROT, SX, VROTA)

## add radii of curvature
D1$Rn <- Ree / (1 - (Ecc*sin(D1$GGLAT*Cradeg))^2)^0.5 + D1$GGALT
D1$Rm <- D1$Rn * (1-Ecc^2) / (1-(Ecc*sin(D1$GGLAT*Cradeg))^2) + D1$GGALT

## time shifts?
# s <- -120
s <- -92
SHIFT <- FALSE
SHIFT <- TRUE
if (SHIFT) {
  D1$GGVEW <- ShiftInTime (D1$GGVEW, Rate, s)
  D1$GGVNS <- ShiftInTime (D1$GGVNS, Rate, s)
  D1[, VROC] <- ShiftInTime (D1[, VROC], Rate, s)
  # D1$BLONGA <- ShiftInTime (D1$BLONGA, Rate, s-80)
  # D1$BLATA <- ShiftInTime (D1$BLATA, Rate, s-80)
  # D1$BNORMA <- ShiftInTime (D1$BNORMA, Rate, s-80)
}
## smooth the measurements when determining derivatives
.span <- 10 * Rate + 1    
## The following are accelerations determined from derivatives of the GPS velocities.
## (Previously used INS velocities to reconstruct IRU values.)
## These should match the measured accelerations after transforming measurements to the l-frame
## and application of the rotation correction:
D1$vndot <- signal::sgolayfilt (D1$GGVNS, 3, .span, m=1) * Rate  # m=1 for first deriv.
D1$vedot <- signal::sgolayfilt (D1$GGVEW, 3, .span, m=1) * Rate
D1$vudot <- signal::sgolayfilt (D1[, VROC], 3, .span, m=1) * Rate
