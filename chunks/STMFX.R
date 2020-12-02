
#' @title State Transformation Matrix (Derivatives)
#' @description Calculates the derivative of a state vector for INS 
#' mechanization.
#' @details The state vector has components representing position,
#' velocity, and attitude angles of the aircraft. Using IRU-provided
#' accelerations and rotation rates, this function calculates the derivative
#' of that state vector, for use e.g. in INS mechanization to propagate 
#' the state vector forward in time. 
#' @aliases STMFV
#' @importFrom Ranadu (XformLA, Gravity)
#' @export STMFV
#' @param sv A state vector with components latitude, longitude, altitude,
#' eastward groundspeed, northward groundspeed, rate-of-climb,
#' pitch, roll, heading. These 9 components must be followed by the
#' body rotation rates about the pitch, roll, and yaw axes and the
#' body accelerations in the lateral, longitudinal, and normal
#' directions. Units are radians for all angles and rotation rates, SI
#' for others.
#' @param Rn The normal radius of the Earth (cf. Noureldin p. 48).
#' This is a function of latitude and altitude. The default is 6378137, the
#' value at the Equator at sea level.
#' @param Rm The meridion radius of the Earth (cf. Noureldin p. 48). The 
#' default is 6335439, the value at the Equator at sea level.
#' @param .aaframe (default 'a') If 'l', specifies that the attitude
#' angles should be those in the l-frame and appropriate transformations
#' should be included before and after calculating the derivative to
#' provide results applicable in the l-frame. Anything else skips these
#' transformations.
#' @return A 15-component vector containing, in the first 9 components, the 
#' derivatives of each component of the state vector. Zero is returned in 
#' the remaining components. Units as in the input state vector.

STMFW <- function (sv, Rn=6378137, Rm=6335439, .aaframe='a') {
  dm <- dim(sv)
  if (!is.null(dm) && length(dm) > 1) {
    VECT <- TRUE
    DL <- dm[1]
  } else {
    VECT <- FALSE
  }
  Cradeg <- pi / 180  ## conversion to radians from degrees
  OmegaE <- 2*pi/86400  ## seems to match INS better than exact value
  OmegaE <- StandardConstant ('Omega')
  if (VECT) {
    stmf <- matrix(nrow=dm[2], ncol=dm[1])
  } else {
    stmf <- vector('numeric', length=15)
  }
  if (.aaframe == 'l') {  ## if angles provided in l-frame:
    svs7 <- sv[7]
    svs8 <- sv[8]  ## this gets saved for restoration at end of function
    sv[7] <- cos (sv[9]) * svs7 - sin (sv[9]) * sv[8]
    sv[8] <- sin (sv[9]) * svs7 + cos (sv[9]) * sv[8]
  } 
  
  ## Construct a single-row data.frame with angles as needed,
  ## then get the transformation matrix a-frame to l-frame
  rlm <- XformLA (data.frame(PITCH=sv[7]/Cradeg, ROLL=sv[8]/Cradeg, 
                             THDG=sv[9]/Cradeg))
  
  ## This is the inertial correction to rotation rates; 
  ## cf. Noureldin p. 179-180
  sinLat <- sin(sv[1]); cosLat <- cos(sv[1]); tanLat <- tan(sv[1])
  omega1 <- as.vector (c(0, OmegaE * cosLat, OmegaE * sinLat), 
                      mode='numeric')
  omega2 <- as.vector (c(-sv[5] / Rm, 
                         sv[4] / Rn,
                         sv[4] * tanLat / Rn), 
                       mode='numeric')
  omegaila <- t(rlm) %*% (omega1 + omega2)

  ## find the derivative of the transformation matrix:
  ## SRM is skew-symmetric representation of measured rotation rates
  ## Explanation of SRR:
  ##    In the a-frame, the {x,y,z} components of the rotation vector
  ##    are {BROLLR, BPITCHR, -BYAWR}, the last because the a-frame z axis
  ##    is downward. The last is called yaw-rotation rate, but its sign
  ##    gives the heading rotation rate, despite the usual convention that
  ##    yaw angle is opposite to heading angle. These signs have been
  ##    confirmed to give the right rotation directions.
  VRX <- c(sv[11] - omegaila[1], sv[10] - omegaila[2], -sv[12] - omegaila[3])
  SRR <- matrix(c(0, -VRX[3], VRX[2],
                  VRX[3], 0, -VRX[1],
                  -VRX[2], VRX[1], 0), ncol=3)
  ## rlm has been transposed before return from XformLA to make this work:
  dRLA <- rlm %*% (SRR)
  
  # A <- 0.11  ## trial attempt to reduce sin^2 dependence
  # dRKA[1,1] <- (1-A*(1-2*sin(sv[9])^2)) * dRLA[1,1]
  # dRKA[1,2] <- (1-A*(1-2*sin(sv[9])^2)) * dRLA[1,2]
  
  ## position derivatives
  DR <- c(sv[5] / Rm, sv[4] / (Rn * cos (sv[1])), sv[6])
  
  ## accelerations = velocity derivatives
  Grav <- as.numeric (Gravity (sv[1] / Cradeg, sv[3]))
  AA <- as.vector (c (sv[14], sv[13], sv[15] + Grav), mode='numeric') # a-frame
  AL <- as.vector (rlm %*% as.matrix (AA), mode='numeric')   # l-frame
  ## correct for angular effects
  VL <- c(sv[4], sv[5], sv[6])
  ### RotationCorrection
  # RotationCorrection <- function (.data, .V) {
  #   C <- vector ('numeric', 3 * (DL <- nrow (.data))); dim(C) <- c(DL,3)
  #   if (DL == 1) {dim (.V) <- c(1,3)}  ## for consistency, vectorized fn
    # omega <- c (-sv[5] / Rm, 2 * OmegaE * cosLat + sv[4] / Rn,
    #             2 * OmegaE * sinLat + sv[4] * tanLat / Rn)
    omega <- 2*omega1 + omega2
    zro <- 0
    M <- array (c(zro, -omega[3], omega[2], omega[3], zro, -omega[1],
                  -omega[2], omega[1], zro), dim=c(3, 3))    ## skew-symmetric
    ## R uses column-major representation.
    MP <- aperm (M)
    RC <- as.vector(MP %*% VL)
  AL <- as.vector (AL - RC)
  AL[3] <- -AL[3] - Grav
  
  ## attitude-angle derivatives: see workflow document
  if (abs(rlm[1,2]) < 0.001) {  ## avoiding rounding errors; see workflow
    RRT <- c(-dRLA[3,1]/sqrt(1-rlm[3,1]^2), 
             1/(1+(rlm[3,2]/rlm[3,3])^2) * (rlm[3,2]/rlm[3,3]) * 
               (-dRLA[3,2]/rlm[3,2] + dRLA[3,3]/rlm[3,3]), 
             (rlm[1,2]*dRLA[1,1]/rlm[1,1]^2 - dRLA[1,2]/rlm[1,1]))
  } else {
    RRT <- c(-dRLA[3,1]/sqrt(1-rlm[3,1]^2), 
             1/(1+(rlm[3,2]/rlm[3,3])^2) * (rlm[3,2]/rlm[3,3]) * 
               (-dRLA[3,2]/rlm[3,2] + dRLA[3,3]/rlm[3,3]), 
             1/(1+(rlm[1,1]/rlm[1,2])^2) * 
               (dRLA[1,1]/rlm[1,2] - dRLA[1,2]*rlm[1,1]/rlm[1,2]^2))
  }
  ## the following suppressed transformation was used to explore the effect
  ## of possible mis-alignment between IRU and aircraft axes:
  # sx <- sv
  # sx[4] <- RRT[1]
  # sx[5] <- RRT[2]
  # sx[6] <- RRT[3]
  # sx[7] <- 0*Cradeg
  # sx[8] <- 0*Cradeg
  # sx[9] <- 0*Cradeg
  # VROTA <- XformLA (data.frame(PITCH=sx[7]/Cradeg, ROLL=sx[8]/Cradeg, THDG=sx[9]/Cradeg)) %*% RRT
  # RRT[1] <- VROTA[2]
  # RRT[2] <- VROTA[1]
  # RRT[3] <- -VROTA[3]
  # A <- 0.112
  # A <- 0
  # YC <- 1-A*(1-2*sin(sv[9]-7*Cradeg)^2)
  stmf <- c(DR, AL, RRT, rep(0, 6))
  
  ## return sv to original if necessary
  if (.aaframe == 'l') {
    sv[7] <- svs7
    sv[8] <- svs8
    stmf7 <- stmf[7]
    stmf[7] <- cos(sv[9])*stmf7 + sin(sv[9])*stmf[8]
    stmf[8] <- -sin(sv[9])*stmf7 + cos(sv[9])*stmf[8]
  }
  
  return (as.vector (stmf, mode='numeric'))
}

