
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
  Ree <- 6378137 ## for radii of curvature
  Ecc <- 0.08181919
  Rn <- Ree / (1 - (Ecc * sinLat)^2)^0.5 
  Rm <- Rn * (1 - Ecc^2) / (1 - (Ecc * sinLat)^2) + .data$GGALT
  Rn <- Rn + .data$GGALT
  omega <- c (-.V[,2] / Rm, 2 * OmegaE * cosLat + .V[,1] / Rn,
              2 * OmegaE * sinLat + .V[,1] * tanLat / Rn)
  dim(omega) <- c(DL, 3)
  zro <- rep(0, DL)
  M <- array (c(zro, -omega[, 3], omega[, 2], omega[, 3], zro, -omega[, 1],
                -omega[, 2], omega[, 1], zro), dim=c(DL, 3, 3))    ## skew-symmetric
  ## R uses column-major representation. Matrix multiplication, each row:
  MP <- aperm (M)
  for (i in 1:DL) {
    C[i,] <- MP[,,i] %*% .V[i,]
  } 
  return (C)
}
