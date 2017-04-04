# Kalman-setup.R

## the chunk is 'sourced' here so the same code can be used in KalmanFilter.R
# source ('chunks/Kalman-setup.R')
## Kalman-setup
## initialize matrices needed by the Kalman filter and load the starting-point
## for the error-state vector.

# DL <- nrow (D1)
# DLL <- DL - DL %% NSTEP
# if (DLL != DL) {
#   D1 <- D1[1:DLL, ]
#   DL <- DLL
# }
## initial values of the state vector and error-state vector:
SV <- with(D1[1, ], data.frame(LAT, LON, ZROC, VEW, VNS, ROC, PITCH, ROLL, THDG,
                               BPITCHR, BROLLR, BYAWR, BLATA, BLONGA, BNORMA))
SV[c(1,2,7:12)] <- SV[c(1,2,7:12)] * Cradeg
## also need the corresponding noise vector:
tau <- 60
GCF <- function (sv, sp) {    ## sv is the state vector; sp is the data record
  gcf <- vector('numeric', length=15)
  ## For lat and lon, noise evidence indicates 1 m, but when that is used
  ## there are numerical problems in the 'solve' call because of the large range
  ## in values in the matrix. Forcing these 10 times larger avoids that problem.
  gcf[1] <- 10 / sp$Rm
  gcf[2] <- 10 / (sp$Rn * cos (sv$LAT))
  gcf[3] <- 10
  gcf[4:5] <- 0.02    ## 0.3
  gcf[6] <- 0.05
  gcf[7:8] <- 0.02*Cradeg  ## 0.005
  gcf[9] <- 0.02*Cradeg    ## 0.015
  gcf[10:11] <- 0.003*Cradeg    ## 0.015
  gcf[12] <- 0.003*Cradeg       ## 0.015
  gcf[13:15] <- 0.00005 # see text re why this is so small
  # gcf <- as.vector(gcf) * sqrt(2/tau)
  return (gcf)
}
gcf <- as.numeric (GCF (SV, D1[1, ]))

## The measurement model: 
##      Calculate the measurements of acceleration from the GPS to add to the measurement vector.
##      This was explored but not used in the version that produced the the Tech Note
GAEL <- c(D1$LACCX - D1$vedot, D1$LACCY - D1$vndot, D1$LACCZ - D1$vudot)
DL <- nrow(D1)
dim(GAEL) <- c(DL, 3)
## transform to the a-frame
GAE <- XformLA (D1, GAEL, .inverse=TRUE)
## get rotation-rate corrections to apply to GPS measurements
LR <- 4.42; LG <- -4.30
Pdot <- c(0, diff (D1$PITCH*Cradeg)) * Rate  # diff does step-wise differentiation
Hdot <- c(0, diff (D1$THDG*Cradeg))          # see Rate multiplication few lines down
Hdot[is.na(Hdot)] <- 0
Hdot[Hdot > pi] <- Hdot[Hdot > pi] - 2*pi
Hdot[Hdot < -pi] <- Hdot[Hdot < -pi] + 2*pi
Hdot <- Hdot * Rate
cospsi <- cos (D1$THDG*Cradeg)
sinpsi <- sin (D1$THDG*Cradeg)
cosphi <- cos (D1$ROLL*Cradeg)

## note that LG is negative:
DZ <- with(D1, c(Cradeg*(LAT-GGLAT)-LG*cospsi/Rm, Cradeg*(LON-GGLON)-LG*sinpsi/Rn, 
                 ZROC-GGALT, VEW-GGVEW + LG*Hdot*cospsi, 
                 VNS-GGVNS - LG*Hdot*sinpsi, ROC-D1[, VROC] + LG*Pdot*cosphi))
# DZ <- c(DZ, GAE[,1], GAE[,2], GAE[,3]) ## add this later?
## The last three components provide direct feedback to measured acceleration
## in the a-frame but also provide feedback to heading, as developed below
# dim(DZ) <- c(DL, 9)
# dim(DZ) <- c(DL, 6)
# DZ[,1:2] <- DZ[, 1:2] * Cradeg
## now add the pseudo-measurement of heading error found from the accelerations:
D1$deltaPsi <- (atan2 (D1$LACCX, D1$LACCY) - atan2 (D1$vedot, D1$vndot))
## alternate based on (12):
D1$deltaPsi <- (D1$LACCX*(D1$vndot-D1$LACCY)-D1$LACCY*(D1$vedot-D1$LACCX)) /
  (D1$LACCX^2+D1$LACCY^2)
D1$deltaPsi[is.na(D1$deltaPsi)] <- 0
D1$deltaPsi[D1$deltaPsi > pi] <- D1$deltaPsi[D1$deltaPsi > pi] - 2*pi
D1$deltaPsi[D1$deltaPsi < -pi] <- D1$deltaPsi[D1$deltaPsi < -pi] + 2*pi
D1$sdPsi <- zoo::rollapply(D1$deltaPsi, 10, sd, fill=NA)  ## calculate the std dev
## add the heading correction to the measurement vector
DZ <- c(as.vector(DZ), D1$deltaPsi)
dim(DZ) <- c(DL, 7)
## if any are NA, substitute zero:
DZ[is.na(DZ)] <- 0

## The observation matrix: (the first six and last three components of the state error 
## vector are observable, the latter requiring transformation from l-frame to a-frame)
## components 7-9 are connected to IRU-measured rotation rates via the transformation matrix l->a,
## so matrix H must vary with aircraft attitude angles. This is addressed by the GAEL->GAE transform.
# H <- diag(1, nrow=10, ncol=15)
# for (k in 7:9) {
#   H[k,k] <- 0
#   H[k,k+6] <- 1
# }
# H[10,9] <- 1
H <- diag(1, nrow=7, ncol=15)
H[7,7] <- 0
H[7,9] <- 1  ## measurement 7 applies to the heading error, SVE component 9


## at any time step, assume the measurements are contaminated by noise:
RCV <- matrix (rep(0,100), ncol=10)
RCV <- matrix (rep(0,49), ncol=7)
RCV[1,1] <- (1/D1$Rm[1])^2  ## latitude
RCV[2,2] <- (1/(D1$Rn[1]*cos(SV$LAT)))^2
RCV[3,3] <- 5^2 #100^2
RCV[4,4] <- 0.1^2            ## ve
RCV[5,5] <- 0.1^2
RCV[6,6] <- 0.1^2
# RCV[7,7] <- 100.0  # A big value here limits updating of acceleration measurement.
# RCV[8,8] <- 100.0  # The assumption is that there errors arise from an error in heading, not
# RCV[9,9] <- 100.0  # measured acceleration, so the GPS-measured acceleration is used for that.
# RCV[10,10] <- D1$sdPsi[1]^2
# if (is.na(D1$sdPsi[1])) {RCV[10,10] <- 225}  # typical sd is 15 deg.
# # RCV[10, 10] <- 1000  ## suppress effect
RCV[7,7] <- 1000  ## but update this each time step

## initialize covariance matrix with generous variances, because Schuler-oscillation
## errors might be large at the start
CV <- matrix (rep(0,225), ncol=15)
CV[1,1] <- 2000^2 / D1$Rm[1]^2
CV[2,2] <- 2000^2 / (D1$Rn[1]*cos(SV[1]))^2
CV[3,3] <- 500^2
CV[4,4] <- 4
CV[5,5] <- 4
CV[6,6] <- 4
CV[7,7] <- (0.3*Cradeg)^2
CV[8,8] <- CV[7,7]
CV[9,9] <- (1*Cradeg)^2
CV[10,10] <- CV[11,11] <- (0.005*Cradeg)^2
CV[12,12] <- (0.01*Cradeg)^2
CV[13,13] <- CV[14,14] <- CV[15,15] <- 0.0005^2

## Q: (initial estimate):
Q <- diag(gcf^2, 15)

