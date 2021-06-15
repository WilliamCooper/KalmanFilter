## utility-functions.R

## check the rate of the file
ATT <- getAttributes(Data, .print=FALSE)
Rate <- 1
if ('sps25' %in% names (ATT$Dimensions)) {Rate <- 25}  ## only rates 1 and 25 supported
dt <- 1/Rate               
DL <- nrow(Data)
OmegaE <- StandardConstant ('Omega')  ## Earth's rotation rate
Ree <- 6378137                        ## for radii of curvature
Ecc <- 0.08181919
Data$Grav <- Gravity (Data$LAT, Data$GGALT)

source('chunks/RotationCorrection.R')
source('chunks/STMFV.R')

#' @title RotationCorrection
#' @description Correct accelerations for inertial effects, a-to-l-frame 
#' @details Calculates the correction needed to account for the rotation
#' of the Earth and of the l-frame (ENU frame). See Noureldin et al., 2013, 
#' Eqs. 5.55--5.57. Subtract this from the accelerations when transforming
#' to the l-frame. This function assumes that Ree, Ecc, and Cradeg=pi/180 
#' have been defined globally and it references those values.
#' @aliases RotationCorrection
#' @param .data A data frame containing the variables LAT and GGALT 
#' representing latitude (in degrees) and geometric altitude (m). There
#' should be one row for each value in the next argument.
#' @param .V A vector of dimension 3 or nx3 representing aircraft velocity.
#' @return An nx3-dimension correction to be applied to accelerations when
#' transforming from the a-frame to the l-frame, to account for inertial
#' effects.

RotationCorrection <- function (.data, .V) {
  C <- vector ('numeric', 3 * (DL <- nrow (.data))); dim(C) <- c(DL,3)
  if (DL == 1) {dim (.V) <- c(1,3)}  ## for consistency, vectorized fn
  lat <- .data$LAT * Cradeg
  sinLat <- sin(lat); cosLat <- cos(lat); tanLat <- tan(lat)
  Rn <- Ree / (1 - (Ecc * sinLat)^2)^0.5 
  Rm <- Rn * (1 - Ecc^2) / (1 - (Ecc * sinLat)^2) + .data$GGALT
  Rn <- Rn + .data$GGALT
  omega <- c (-.V[,2] / Rm, 2 * OmegaE * cosLat + .V[,1] / Rn,
              2 * OmegaE * sinLat + .V[,1] * tanLat / Rn)
  dim(omega) <- c(DL, 3)
  zro <- rep(0, DL)
  M <- array (c(zro, -omega[3], omega[2], omega[3], zro, -omega[1],
                -omega[2], omega[1], zro), dim=c(DL, 3, 3))    ## skew-symmetric
  ## R uses column-major representation. Matrix multiplication, each row:
  MP <- aperm (M)
  for (i in 1:DL) {C[i,] <- MP[,,i] %*% .V[i,]} 
  return (C)
}

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
#' @param .aaframe (default 'a') If 'l', specifies that the attitude
#' angles should be those in the l-frame and appropriate transformations
#' should be included before and after calculating the derivative to
#' provide results applicable in the l-frame.
#' @return A 15-component vector containing the derivatives of each component
#' of the state vector in the first 9 components and zero in the remaining
#' components.

STMFV <- function (sv, .aaframe='a') { 
  stmf <- vector('numeric', length=15)
  ## transform back to a-frame if angles provided in l-frame
  if (.aaframe == 'l') {  
    svs7 <- sv[7]
    svs8 <- sv[8]
    sv[7] <- cos (sv[9]) * svs7 - sin (sv[9]) * svs8
    sv[8] <- sin (sv[9]) * svs7 + cos (sv[9]) * svs8
  } 
  rlm <- XformLA (data.frame(PITCH=sv[7]/Cradeg, ROLL=sv[8]/Cradeg, 
                             THDG=sv[9]/Cradeg))
  omega <- as.vector (c(-sv[5] / Rm, 
                        OmegaE * cos(sv[1]) + sv[4] / (Rn),
                        OmegaE * sin(sv[1]) + sv[4] * tan(sv[1]) / Rn), 
                      mode='numeric')
  ## signs account for b-frame to a-frame: reverse sign of 3, swap 1 and 2
  Oill <- matrix (c(0, omega[3], omega[2], 
                    -omega[3], 0, -omega[1], 
                    -omega[2], omega[1], 0), ncol=3)
  Oilb <- Oill %*% rlm
  ## find the derivative of the transformation matrix:
  ## SRM is skew-symmetric representation of measured rotation rates
  SRR <- c(0, -sv[12], -sv[10],
           sv[12], 0, sv[11],
           sv[10], -sv[11], 0)
  SRM <- aperm( array (SRR, dim=c(3,3)))
  dRLA <- rlm %*% SRM - Oilb
  # A <- 0.11  ## trial attempt to reduce sin^2 dependence
  # dRKA[1,1] <- (1-A*(1-2*sin(sv[9])^2)) * dRLA[1,1]
  # dRKA[1,2] <- (1-A*(1-2*sin(sv[9])^2)) * dRLA[1,2]
  DR <- c(sv[5] / Rm, sv[4] / (Rn * cos (sv[1])), sv[6])
  Grav <- as.numeric (Gravity (sv[1]/Cradeg, sv[3]))
  AA <- as.vector (c(sv[14], sv[13], sv[15]+Grav), mode='numeric') # a-frame
  AL <- as.vector (rlm %*% as.matrix (AA), mode='numeric')   # l-frame
  ## now correct for angular effects
  VL <- c(sv[4], sv[5], sv[6])
  svdf <- data.frame(LAT=sv[1] / Cradeg, GGALT=sv[3])
  AL <- as.vector (AL - RotationCorrection (svdf, VL), mode='numeric')
  AL[3] <- AL[3] + Grav
  AL[3] <- -AL[3]
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
  stmf <- c(DR, AL, RRT, 0, 0, 0, 0, 0, 0)
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

