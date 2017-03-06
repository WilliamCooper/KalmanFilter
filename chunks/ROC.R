## ROC.R -- chunk to add ROC variable
## needs: data.frame Data containing netCDF variables PSXC, GGLAT, GGALT, ACINS, Grav
## also assumes Rate is set


DPDT <- c(0, diff(Data$PSXC)) * Rate
g <- Data$Grav
g[is.na(g)] <- 9.80
WPPRIME <- -StandardConstant('Rd') * (273.15 + Data$ATX) /
              (Data$PSXC * g) * DPDT
WPSTAR <- cumsum(Data$ACINS)
DIF <- WPPRIME - WPSTAR
DIF <- zoo::na.approx (as.vector(DIF), maxgap=1000, na.rm=FALSE)
tau <- 300
DIF <- signal::filtfilt (signal::butter (3, 2/tau), DIF)
Data$ROC <- WPSTAR + DIF
Data$ZROC <- Data$GGALT[1] + cumsum (Data$ROC)
rm (DPDT, g, WPPRIME, WPSTAR, DIF)
