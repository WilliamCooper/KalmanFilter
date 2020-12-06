library(Ranadu)

KFVARS <- c("AKKF", "ALTKF",  "LATKF", "LONKF", "PITCHKF", "ROCKF",
            "ROLLKF",  "SSKF",  "THDGKF",  "UXKF",  "VEWKF",
            "VNSKF", "VYKF",   "WDKF",   "WIKF",
            "WSKF")
STDVARS <- c('AKRD', 'GGALT', 'LATC', 'LONC', 'PITCH', 'GGVSPD', 'ACINS',
             'LAT', 'LON',
             'ROLL', 'SSRD', 'THDG', 'UXC', 'VYC', 'WIC', 'WSC', 'WDC',
             'VEW', 'VNS', 'GGLAT', 'GGLON', 'GGVEW', 'GGVNS')
Data <- getNetCDF(setFileName('WCR-TEST', 'tf01KF'), standardVariables(unique(c(STDVARS, KFVARS))))
if(any(Data$TASX > 80)) {
  itas <- which(Data$TASX > 80)
  Data <- Data[itas[1]:itas[length(itas)], ]
}
minT <- Data$Time[1]
maxT <- Data$Time[nrow(Data)]
formatTime <- function (time) {
  t <- as.POSIXlt (time)
  tt <- sprintf ("%d%02d%02d", t$hour, t$min, t$sec)
  tt <- as.integer(tt)
  return (tt)
}
panelylim <- NA
lastData <- ' '
lastFLT <- ' '
lastType <- ' '
PJ <- list.dirs(path = DataDirectory(), full.names = FALSE, recursive = FALSE)
PJ <- PJ[-which('lost+found' == PJ)]
## Leave in alphabetical order, except for the first which is the latest modified.
## If the Kalman Filter has added to the directory in a previous run, that will
## change the modified time for the directory and place that directory to the top
## for the next run.
FullPJ <- paste0(DataDirectory(), PJ)
iw <- which.max(file.mtime(FullPJ))
PJ <- c(PJ[iw], PJ[-iw])
## Keep only directories with a rf01 or tf01
for (P in PJ) {
  if (grepl('HIPPO', P)) {
    fn <- sprintf ('%sHIPPO/%srf01.nc', DataDirectory (), P)
  } else {
    fn <- sprintf ('%s%s/%srf01.nc', DataDirectory (), P, P)
    if (!file.exists (fn)) {
      fn <- sub ('\\.nc', '.Rdata', fn)
    }
    if (!file.exists (fn)) {
      fn <- sprintf ('%s%s/%stf01.nc', DataDirectory (), P, P)
    }
    if (!file.exists (fn)) {
      fn <- sub ('\\.nc', '.Rdata', fn)
    }
  }
  if (!file.exists (fn)) {PJ[PJ==P] <- NA}
}
PJ <- PJ[!is.na(PJ)]

