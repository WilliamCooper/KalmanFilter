
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
#' @param .components (default 9) The number of components in the derivative
#' vector to return. This is normally 9 for mechanization, but needs to be
#' 15 when used to calculate the Jacobian of the derivative matrix.
#' @return A 9-component or 15-component array, with one row for each row 
#' in D, which contains the derivatives of each component of the state 
#' vector.  The state vector is assumed to be in the order latitude, 
#' longitude, altitude, eastward groundspeed, northward groundspeed, 
#' rate-of-climb, pitch, roll, heading, with the last three in radians. 

STMFV <- function (D, .components=15, .aaframe='a') {
  if (is.vector(D)) {D <- matrix(D, ncol=.components)}
  DL <- nrow(D)
  Cradeg <- pi / 180  ## conversion to radians from degrees
  OmegaE <- 2*pi/86400  ## seems to match INS better than exact value
  OmegaE <- StandardConstant ('Omega')
  stmf <- matrix(0, nrow=DL, ncol=.components)
  ## If D is provided as a matrix, convert to a data.frame
  ## and assume angles are provided as radians
  if (is.matrix(D)) {
    D <- as.data.frame (D)
    D[, c(1:2, 7:12)] <- D[, c(1:2, 7:12)] / Cradeg 
    names(D) <- c('LAT', 'LON', 'ZROC', 'VEW', 'VNS', 'ROC', 
                  'PITCH', 'ROLL', 'THDG', 'BPITCHR', 'BROLLR',
                  'BYAWR', 'BLATA', 'BLONGA', 'BNORMA')
  }

  if (.aaframe == 'l') {  ## if angles provided in l-frame:
    svs7 <- D$PITCH
    svs8 <- D$ROLL  
    thdg <- D$THDG * Cradeg
    D$PITCH <- cos (thdg) * svs7 - sin (thdg) * svs8
    D$ROLL  <- sin (thdg) * svs7 + cos (thdg) * svs8
  }
  
  if (!('ROC' %in% names(D))) {
    D$ROC <- D$VSPD
    D$ZROC <- D$ALT
  }
  ## this assumes that Ree and Ecc are in the environment
  lat <- D$LAT*Cradeg
  Rn <- Ree / (1 - (Ecc*sin(lat))^2)^0.5
  Rm <- Rn * (1-Ecc^2) / (1-(Ecc*sin(lat))^2) + D$ZROC
  Rn <- Rn + D$ZROC
  Grav <- Gravity (D$LAT, D$ZROC)

  ## attitude-angle derivatives: see workflow document
  
  ## This is the inertial correction to rotation rates; 
  ## cf. Noureldin p. 179-180
 
  sinLat <- sin(lat); cosLat <- cos(lat); tanLat <- sinLat / cosLat
  omega1 <- matrix (c (rep(0, DL), OmegaE * cosLat, OmegaE * sinLat), ncol=3)
  omega2 <- with(D, matrix (c (-VNS / Rm,  VEW / Rn, VEW * tanLat / Rn), ncol=3))
  omega <- omega1 + omega2
  # omegaila <- XformLA(D, omega, .inverse=TRUE)
  # omegaila[,3] <- -omegaila[,3]

  ## find the derivative of the transformation matrix:
  ## SRM is skew-symmetric representation of measured rotation rates
  ## Explanation of SRR:
  ##    In the a-frame, the {x,y,z} components of the rotation vector
  ##    are {BROLLR, BPITCHR, -BYAWR}, the last because the a-frame z axis
  ##    is downward. The last is called yaw-rotation rate, but its sign
  ##    gives the heading rotation rate, despite the usual convention that
  ##    yaw angle is opposite to heading angle. These signs have been
  ##    confirmed to give the right rotation directions.
  VRX <- with(D, matrix(c(BROLLR, BPITCHR, -BYAWR), ncol=3)) * Cradeg
  # VRX <- VRX - omegaila

  zro <- c(rep(0, DL))
  SRR <- array(c(zro, -VRX[,3], VRX[,2], 
                 VRX[,3], zro, -VRX[,1], 
                 -VRX[,2], VRX[,1], zro), dim=c(DL,3,3))
  
  ## inertial effects. Note sign reversal of component 3:
  SRRI <- array(c(zro, omega[,3], omega[,2],
                  -omega[,3], zro, -omega[,1],
                  -omega[,2], omega[,1], zro), dim=c(DL,3,3))
  
  dRLA <- array(0, dim=c(DL,3,3))
  RRT <- matrix(0, nrow=DL, ncol=3)
  for (j in 1:DL) {
    rlm <- XformLA (D[j, ])
    ## here "rlm %*% t(rlm) = I" has been used before the 2nd term
    dRLA[j,,] <- rlm %*% SRR[j,,] - SRRI[j,,] %*% rlm
    if (abs (rlm[1,1] / rlm[2,1]) > 10) {  ## avoiding rounding errors; see workflow
      RRT[j,] <- c (-dRLA[j,3,1] / sqrt (1 - rlm[3,1]^2), 
               1 / (1 + (rlm[3,2] / rlm[3,3])^2) * (rlm[3,2] / rlm[3,3]) * 
               (-dRLA[j,3,2] / rlm[3,2] + dRLA[j,3,3] / rlm[3,3]), 
               (rlm[2,1] * dRLA[j,2,1] / rlm[1,1]^2 - dRLA[j,2,1] / rlm[1,1]))
    } else {
      RRT[j,] <- c (-dRLA[j,3,1] / sqrt (1 - rlm[3,1]^2), 
               1 / (1 + (rlm[3,2] / rlm[3,3])^2) * (rlm[3,2] / rlm[3,3]) * 
               (-dRLA[j,3,2] / rlm[3,2] + dRLA[j,3,3] / rlm[3,3]), 
               1 / (1 + (rlm[1,1] / rlm[2,1])^2) * 
               (dRLA[j,1,1] / rlm[2,1] - dRLA[j,2,1] * rlm[1,1] / rlm[2,1]^2))
    }
  }

  ## position derivatives
  DR <- with(D, matrix(c (VNS / Rm, VEW / (Rn * cos (lat)), ROC)), ncol=3)
  
  ## accelerations = velocity derivatives
    AA <- with (D, matrix (c (BLONGA, BLATA, BNORMA + Grav), ncol=3)) # a-frame
    AL <- XformLA (D, AA)   # l-frame
    ## correct for angular effects
    VL <- with (D, matrix (c (VEW, VNS, ROC), ncol=3))
    omega <- 2 * omega1 + omega2
    M <- array (c (zro, -omega[ ,3], omega[ ,2], omega[ ,3], zro, -omega[ ,1],
                  -omega[ ,2], omega[ ,1], zro), dim=c(DL, 3, 3))    ## skew-symmetric
    ## R uses column-major representation.
    MP <- aperm (M)
    RC <- matrix(0, nrow=DL, ncol=3)
    for (i in 1:DL) {
      RC[i, ] <- as.vector(MP[,, i] %*% VL[i, ])
    }
    AL <- AL - RC
    AL[ ,3] <- -AL[ ,3] - Grav

  if (.components != 9) {
    stmf <- matrix (c (DR, AL, RRT, rep (0, DL*(.components-9))), 
                    nrow=DL, ncol=.components)
  } else {
    stmf <- matrix (c (DR, AL, RRT), nrow=DL, ncol=9)
  }
  return (stmf)
}

