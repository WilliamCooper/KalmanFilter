# Kalman-loop.R

SVEF <- array(dim=c(DL/NSTEP, 15))
CVEF <- array(dim=c(DL/NSTEP, 15))
## initialize error state vector
SVE <- rep (0, 15)  ## respectively: (lat,lon,alt) (vew,vns,vspd),
## (pitch,roll,thdg) (rot. rates) (accel components)
SVE[1:6] <- DZ[1, 1:6]
SVE[7:9] <- 0  ## might initialize using pitch/roll/heading-correction functions here
## start with zero for gyro and accelerometer errors
SVEF[1, ] <- SVE
CVEF[1, ] <- diag (CV)
pctL <- 0
MH <- 1    ## 0 suppresses use of deltaPsi as error in heading
for (i in seq(2*NSTEP, DL, by=NSTEP)) {
  pct <- as.integer(100*i/DL)
  if (pct %% 5 == 0 && pct != pctL) {print (sprintf ('Kalman-loop %d%% done', pct));pctL <- pct}
  SV <- with(D1[i, ], data.frame(LAT, LON, ZROC, VEW, VNS, ROC, PITCH, ROLL, THDG,
                                 BPITCHR, BROLLR, BYAWR, BLATA, BLONGA, BNORMA))
  SV[c(1,2,7:12)] <- SV[c(1,2,7:12)] * Cradeg
  Rn <- D1$Rn[i]
  Rm <- D1$Rm[i]
  Grav <- D1$Grav[i]
  sv <- as.vector (SV, mode='numeric')
  # stmf <- STMFV (sv)
  dcm <- jacobian (STMFV, sv) * dt * NSTEP + diag(15)
  ## modify to include this?
  ## modify this to include decaying error terms for the measurements:
  # dcm[10,10] <- dcm[11,11] <- dcm[12,12] <- -1/tau
  # dcm[13,13] <- dcm[14,14] <- dcm[15,15] <- -1/tau
  ## predict the new error-state vector:
  ## the pitch and roll error derivatives are now those in the a-frame, but
  ##   the pitch/roll error state is in the l-frame. Save the l-frame error state:
  #   SVEL <- SVEA <- SVE
  #   ## Transform l-frame pitch/roll error state to a-frame
  #   SVEA[7] <- cos(sv[9]) * SVEL[7] + sin(sv[9]) * SVEL[8]
  #   SVEA[8] <- -sin(sv[9]) * SVEL[7] + cos(sv[9]) * SVEL[8]
  #   ## apply derivatives to get a-frame change
  #   SVEA <- dcm %*% SVEA
  #   ## transform back to l-frame
  #   SVEL[7] <- cos(sv[9]) * SVEA[7] - sin(sv[9]) * SVEA[8]
  #   SVEL[8] <- sin(sv[9]) * SVEA[7] + cos(sv[9]) * SVEA[8]
  
  SVE <- dcm %*% SVE    ## take time step in error-state vector
  
  ## update the covariance matrix:
  CV <- dcm %*% (CV %*% t(dcm)) + Q
  if (is.na(D1$sdPsi[i]) || (sqrt(D1$LACCX[i]^2+D1$LACCY[i]^2) < 1)) {
    H[7,9] <- 0
    # DZ[i, 7] <- NA
  } else {
    H[7,9] <- MH
    RCV[7, 7] <- 10*D1$sdPsi[i]^2
  }
  ## the Kalman gain:
  Kb <- solve (H %*% CV %*% t(H) + RCV)
  K <- CV %*% t(H) %*% Kb
  DZZ <- DZ[i, ] - H %*% SVE
  SVE <- SVE + K %*% DZZ
  CV <- CV - K %*% H %*% CV
  SVEF[i/NSTEP, ] <- SVE
  CVEF[i/NSTEP, ] <- diag(CV)
}

## interpolate the results over the full data.frame:
IntFilter <- function (X, inRate, outRate) {
  if (inRate == outRate) {return (X)}
  ratio <- as.integer(outRate/inRate)    ## expected to be an integer
  x <- 0:(length(X)-1)
  A <- stats::approx (x, X, n=length(X)*ratio-ratio+1)
  T <- A$y
  T <- signal::filter(signal::sgolay(4,75),T)
  ## now shift to match Rate:
  n <- as.integer (ratio / 2)
  NL = length(T)
  T <- c(rep(T[1],n), T, rep(T[NL],ratio-n-1))  ## OK, even or odd ratio
  return (T)
}
## check for size mis-match
DLSV <- dim(SVEF)[1]
Xtra <- DL - NSTEP * DLSV
Cor <- vector('numeric', DL*15)
dim (Cor) <- c(DL, 15)
VCor <- vector ('numeric', DL*15)
dim (VCor) <- c(DL, 15)
X <- SVEF[, 1]
for (j in 1:15) {
  Cor[, j] <- c(IntFilter (SVEF[, j], 1, NSTEP), rep (SVEF[DLSV, j], Xtra))
  VCor[, j] <- c(IntFilter (CVEF[, j], 1, NSTEP), rep (CVEF[DLSV, j], Xtra))
  VCor[VCor[,j] < 0] <- 0 
  if (j > 6) {next}
  Cor[, j] <- zoo::na.approx (as.vector (Cor[, j]), maxgap=1000, na.rm=FALSE)
  Cor[is.na(Cor[, j]), j] <- 0
  Cor[, j] <- signal::filtfilt (signal::butter (3, 2/600), Cor[, j])
}
# Cor7 <- Cor[, 7]
# VC7 <- VCor[, 7]
# .hdg <- D1$THDG * Cradeg
# Cor[, 7] <- cos(.hdg) * Cor7 - sin(.hdg) * Cor[, 8]
# Cor[, 8] <- sin(.hdg) * Cor7 + cos(.hdg) * Cor[, 8]
# VCor[, 7] <- cos(.hdg) * VC7 - sin(.hdg) * VCor[, 8]
# VCor[, 8] <- sin(.hdg) * VC7 + cos(.hdg) * VCor[, 8]
D1$LATKF <- D1$LAT - Cor[, 1]/Cradeg
D1$LONKF <- D1$LON - Cor[, 2]/Cradeg
## filter the result to smooth the jumps arising from limited INS resolution:
D1$LATKF <- zoo::na.approx (as.vector (D1$LATKF), maxgap=1000, na.rm=FALSE)
D1$LONKF <- zoo::na.approx (as.vector (D1$LONKF), maxgap=1000, na.rm=FALSE)
D1$LATKF[is.na(D1$LATKF)] <- 0
D1$LONKF[is.na(D1$LONKF)] <- 0
D1$LATKF <- signal::filtfilt (signal::butter (3, 2/(10*Rate)), D1$LATKF)
D1$LONKF <- signal::filtfilt (signal::butter (3, 2/(10*Rate)), D1$LONKF)
D1$ALTKF <- D1$ZROC - Cor[, 3]
D1$VEWKF <- D1$VEW - Cor[, 4]
D1$VNSKF <- D1$VNS - Cor[, 5]
D1$ROCKF <- D1$ROC - Cor[, 6]
D1$PITCHKF <- D1$PITCH - Cor[, 7]/Cradeg
D1$ROLLKF <- D1$ROLL - Cor[, 8]/Cradeg
D1$THDGKF <- D1$THDG - Cor[, 9]/Cradeg
D1$BPITCHRKF <- D1$BPITCHR - Cor[, 10]/Cradeg
D1$BROLLRKF <- D1$BROLLR - Cor[, 11]/Cradeg
D1$BYAWRKF <- D1$BYAWR - Cor[, 12]/Cradeg
D1$BLATAKF <- D1$BLATA - Cor[, 13]
D1$BLONGAKF <- D1$BLONGA - Cor[, 14]
D1$BNORMAKF <- D1$BNORMA - Cor[, 15]

D1$DLAT <- D1$LATKF-D1$GGLAT
D1$DLATM <- D1$DLAT * .060 * StandardConstant ('Cmfromnmi')
D1$CLAT <- -Cor[, 1] / Cradeg
D1$CLATM <- D1$CLAT * .060 * StandardConstant ('Cmfromnmi')
D1$DLON <- D1$LONKF-D1$GGLON
D1$DLONM <- D1$DLON * .060 * StandardConstant ('Cmfromnmi') * cos (D1$GGLAT * Cradeg)
D1$CLON <- -Cor[, 2] / Cradeg
D1$CLONM <- D1$CLON * .060 * StandardConstant ('Cmfromnmi') * cos (D1$GGLAT * Cradeg)
D1$DALT <- (D1$ALTKF-D1$GGALT)/1000
D1$CALT <- -Cor[, 3]/1000
D1$DVEW <- D1$VEWKF-D1$GGVEW
D1$CVEW <- -Cor[, 4]
D1$DVNS <- D1$VNSKF-D1$GGVNS
D1$CVNS <- -Cor[, 5]
D1$DROC <- D1$ROCKF-D1[, VROC]
D1$CROC <- -Cor[, 6]
D1$CPITCH <- -Cor[, 7] / Cradeg
D1$CROLL <- -Cor[, 8] / Cradeg
D1$CTHDG <- -Cor[, 9] / Cradeg

XPitch <- function (dP, dR, hdg, p, r, .inverse=FALSE) {
  sh <- sin(hdg); ch <- cos(hdg); sp <- sin(p); cp <- cos(p); sr <- sin(r); cr <- cos(r)
  if (.inverse) {
    dPX <- ch*cp * dP - sh*cp * dR
    dRX <- -(ch*sp*sr - sh*cr) * dP + (sh*sp*sr + ch*cr) * dR
  } else {
    dPX <- ((sh*sp*sr + ch*cr) * dP + sh*cp * dR) / (cp*cr)
    dRX <- ((ch*sp*sr - sh*cr) * dP + ch*cp * dR) / (cp*cr)
  }
  return(matrix(c(dPX, dRX), ncol=2))
}

.hdg <- D1$THDG * Cradeg
.p <- D1$PITCH * Cradeg
.r <- D1$ROLL * Cradeg
dPx <- XPitch (Cor[,7], Cor[,8], .hdg, .p, .r)
D1$CPL <- dPx[,1] / Cradeg
D1$CRL <- dPx[,2] / Cradeg
D1$SDCPL <- sqrt(cos(.hdg)^2*VCor[,7]+sin(.hdg)^2*VCor[,8]) / (30*Cradeg)
D1$SDCRL <- sqrt(sin(.hdg)^2*VCor[,7]+cos(.hdg)^2*VCor[,8]) / (30*Cradeg)
D1$SDCPA <- sqrt(cos(.hdg)^2*D1$SDCPL^2 + sin(.hdg)^2*D1$SDCRL^2)
D1$SDCRA <- sqrt(sin(.hdg)^2*D1$SDCPL^2 + cos(.hdg)^2*D1$SDCRL^2)
D1$CPLF <- signal::filtfilt (signal::butter (3, 1/900), D1$CPL)
D1$CRLF <- signal::filtfilt (signal::butter (3, 1/900), D1$CRL)
dPxa <- XPitch (D1$CPLF, D1$CRLF, .hdg, .p, .r, .inverse=TRUE)
D1$CPAF <- dPxa[,1]
D1$CRAF <- dPxa[,2]
## save corrected values, obtained by subtracting the smoothed a-frame corrections:
D1$PITCHKF <- D1$PITCH - D1$CPAF
D1$ROLLKF <- D1$ROLL - D1$CRAF

HE <- VCor[,9]
HE[HE < 0.00001] <- 0.00001
SS <- smooth.spline(D1$Time, D1$CTHDG, w=1/HE, spar=1.1)
D1$HCS <- predict(SS, as.numeric(D1$Time))$y
D1$THDGKF <- D1$THDG - D1$HCS    ## save the corrected heading

