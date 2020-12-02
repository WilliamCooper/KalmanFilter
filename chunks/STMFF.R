
#' @title State Transformation Matrix (Derivatives)
#' @description Calculates the derivative of a state vector for INS 
#' mechanization.
#' @details The state vector has components representing position,
#' velocity, and attitude angles of the aircraft. Using IRU-provided
#' accelerations and rotation rates, this function calculates the derivative
#' of that state vector, for use e.g. in INS mechanization to propagate 
#' the state vector forward in time. If given a data.frame, it will
#' return the derivative array for each row of the data.frame.
#' @aliases STMFV
#' @importFrom Ranadu (XformLA, Gravity)
#' @export STMFV
#' @param D A data.frame or array with components latitude, longitude, 
#' altitude, eastward groundspeed, northward groundspeed, rate-of-climb,
#' pitch, roll, heading. These must have the names LAT, LON, ZROC, VEW, VNS,
#' ROC, PITCH, ROLL, THDG. The data frame must also include the 
#' body rotation rates about the pitch, roll, and yaw axes (BPITCHR,
#' BROLLR, BYAWR) and the body accelerations in the lateral, longitudinal, 
#' and normal directions (BLATA, BLONGA, BNORMA). At input, units are 
#' degrees for all angles and rotation rates, SI for others. At output, 
#' units are radians/s for the derivatives of the attitude angles. The 
#' data.frame should also contain values (Rn, Rm) for the normal and meridional
#' radii of the Earth (cf. Noureldin p. 48), which are functions of 
#' latitude and altitude. 
#' @param .aaframe (default 'a') If 'l', specifies that the attitude
#' angles should be those in the l-frame and appropriate transformations
#' should be included before and after calculating the derivative to
#' provide results applicable in the l-frame. Anything else skips these
#' transformations.
#' @param .components (default 15) The number of components in the derivative
#' vector to return. This is normally 9 for mechanization, but needs to be
#' 15 when used to calculate the Jacobian of the derivative matrix.
#' @return A 9-component or 15-component array, with one row for each row 
#' in D, which contains the derivatives of each component of the state 
#' vector.  The state vector is assumed to be in the order latitude, 
#' longitude, altitude, eastward groundspeed, northward groundspeed, 
#' rate-of-climb, pitch, roll, heading, with the last three in radians. 

STMFF <- function (D, .components=15, .aaframe='a') {
  ## Version called with a 15-component vector, 
  ## with angles provided in radians
  
  stmf <- vector(mode='numeric', length=15)
  ## If D is provided as a matrix, convert to a data.frame
  ## and assume angles are provided as radians

  lat <- D[1]
#  D[c(1:2, 7:12)] <- D[c(1:2, 7:12)] / Cradeg 
  sinLat <- sin(lat); cosLat <- cos(lat); tanLat <- sinLat / cosLat
  Rn <- Ree / (1 - (Ecc*sinLat)^2)^0.5
  # Rm <- Rn * (1-Ecc^2) / (1-(Ecc*sin(lat))^2) + D$ZROC
  Rm <- Rn * (1 - Ecc^2) / (1-(Ecc*sinLat)^2) + D[3]
  # Rn <- Rn + D$ZROC
  Rn <- Rn + D[3]
  Grav <- Gravity (D[1] / Cradeg, D[3])
  
  ## attitude-angle derivatives: see workflow document
  
  ## This is the inertial correction to rotation rates; 
  ## cf. Noureldin p. 179-180
  
  
  omega1 <- c(0, OmegaE * cosLat, OmegaE * sinLat)
  omega2 <- c (-D[5] / Rm,  D[4] / Rn, D[4] * tanLat / Rn)
  omega <- omega1 + omega2
  VRX <- c(D[11], D[10], -D[12])
  # VRX <- VRX - omegaila
  
  zro <- 0  # c(rep(0, DL))
  SRR <- matrix(c(zro, -VRX[3], VRX[2], 
                  VRX[3], zro, -VRX[1], 
                  -VRX[2], VRX[1], zro), ncol=3)
  
  ## inertial effects. Note sign reversal of component 3:
  SRRI <- matrix(c(zro, omega[3], omega[2],
                   -omega[3], zro, -omega[1],
                   -omega[2], omega[1], zro), ncol=3)
  
  D[7:9] <- D[7:9] / Cradeg   ## convert to degrees for XformLA
  rlm <- XformLA (data.frame(ROLL=D[8], PITCH=D[7], THDG=D[9]))
  ## here "rlm %*% t(rlm) = I" has been used before the 2nd term
  dRLA <- rlm %*% SRR - SRRI %*% rlm
  if (abs (rlm[1,1] / rlm[2,1]) > 10) {  ## avoiding rounding errors; see workflow
    RRT <- c (-dRLA[3,1] / sqrt (1 - rlm[3,1]^2), 
              1 / (1 + (rlm[3,2] / rlm[3,3])^2) * (rlm[3,2] / rlm[3,3]) * 
                (-dRLA[3,2] / rlm[3,2] + dRLA[3,3] / rlm[3,3]), 
              (rlm[2,1] * dRLA[2,1] / rlm[1,1]^2 - dRLA[2,1] / rlm[1,1]))
  } else {
    RRT <- c (-dRLA[3,1] / sqrt (1 - rlm[3,1]^2), 
              1 / (1 + (rlm[3,2] / rlm[3,3])^2) * (rlm[3,2] / rlm[3,3]) * 
                (-dRLA[3,2] / rlm[3,2] + dRLA[3,3] / rlm[3,3]), 
              1 / (1 + (rlm[1,1] / rlm[2,1])^2) * 
                (dRLA[1,1] / rlm[2,1] - dRLA[2,1] * rlm[1,1] / rlm[2,1]^2))
  }
  # }
  
  ## position derivatives
  DR <- c (D[5] / Rm, D[4] / (Rn * cosLat), D[6])
  
  ## accelerations = velocity derivatives
  # AA <- c (D[14], D[13], D[15] + Grav) # a-frame
  AA <- matrix (c(D[14], D[13], D[15] + Grav), ncol=3)
  AL <- XformLA (data.frame(ROLL=D[8], PITCH=D[7], THDG=D[9]), AA)   # l-frame
  ## correct for angular effects
  VL <- c (D[4], D[5], D[6])
  omega <- 2 * omega1 + omega2
  M <- matrix (c (zro, -omega[3], omega[2], omega[3], zro, -omega[1],
                  -omega[2], omega[1], zro), ncol=3)    ## skew-symmetric
  ## R uses column-major representation.
  MP <- aperm (M)
  # for (i in 1:DL) {
  RC <- as.vector(MP %*% VL)
  AL <- AL - RC
  AL[3] <- -AL[3] - Grav
  
  stmf <- c (DR, AL, RRT, rep (0,6))
  
  return (stmf)
}

