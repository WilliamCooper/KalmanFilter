## ROC.R -- chunk to add ROC variable
## needs: data.frame Data containing netCDF variables PSXC, GGLAT, GGALT, ACINS, Grav, ATX
## also assumes Rate is set

PSXC <- zoo::na.approx (as.vector(Data$PSXC), maxgap=1000*Rate, 
                        na.rm=FALSE, rule = 2)
ATX <- zoo::na.approx (as.vector(Data$ATX), maxgap=1000*Rate, 
                        na.rm=FALSE, rule = 2)
DPDT <- c(0, diff(PSXC)) * Rate
g <- Data$Grav
g <- zoo::na.approx (as.vector(g), maxgap=1000*Rate, 
                       na.rm=FALSE, rule = 2)
g[is.na(g)] <- 9.80
WPPRIME <- -StandardConstant('Rd') * (273.15 + ATX) /
              (PSXC * g) * DPDT
ACINS <- zoo::na.approx (as.vector(Data$ACINS), maxgap=1000*Rate, 
                         na.rm=FALSE, rule = 2)
ACINS[is.na(ACINS)] <- 0
WPSTAR <- cumsum(ACINS) / Rate
## Correct for possible offset:
indx <- 1:nrow(Data)
# ff <- lm(WPSTAR ~ indx)
# WPSTAR <- WPSTAR - coef(ff)[2] * indx
WPSTAR <- WPSTAR - WPSTAR[length(WPSTAR)] * indx / nrow(Data)
DIF <- WPPRIME - WPSTAR
DIF <- zoo::na.approx (as.vector(DIF), maxgap=1000*Rate, na.rm=FALSE,
                       rule = 2)
DIF[is.na(DIF)] <- 0
tau <- 600 * Rate
DIF <- signal::filtfilt (signal::butter (3, 2/tau), DIF)
Data$ROC <- WPSTAR + DIF
Data$ZROC <- Data$GGALT[1] + cumsum (Data$ROC) / Rate
rm (DPDT, g, WPPRIME, WPSTAR, DIF)
