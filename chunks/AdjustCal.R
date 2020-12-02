## AdjustCal.R
## apply adjustments and calibrations to the IRU measurements, and time shifts
## to some other variables. Also add the l-frame accelerations determined by
## differentiation of the GPS measurements

## this loads the cal coefficients cfa1, cfa2, cfa3 determined in Tech Note, for GV
load(file='./BodyAccCal.Rdata')
if (grepl('130', FI_KF$Platform) && (file.exists('./BodyAccCalC130.Rdata'))) {load(file='./BodyAccCalC130.Rdata')}
## adjustments:
  D1$BLONGA <- cfa1[1] + cfa1[2] * Data$BLONGA
  # D1$BLATA  <- cfa2[1] + cfa2[2] * Data$BLATA
  D1$BNORMA <- cfa3[1] + cfa3[2] * Data$BNORMA
# D1$BYAWR <- (D1$BYAWR - 0.0003) * 1.008
# D1$BROLLR <- D1$BROLLR + 0.0005
# D1$BPITCHR <- D1$BPITCHR + 0.0003
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
  ## Save the shifts (used later in "CorrectHeading" to avoid duplication in shifts)
  attr(D1$GGVEW, 'TimeLag') <- -s
  attr(D1$GGVNS, 'TimeLag') <- -s
  attr(D1[, VROC], 'TimeLag') <- -s
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
