## modify-new-netcdf.R

## variables needed for attributes and new wind calculation:
VarList <- c('TASX', 'ATTACK', 'SSLIP', 'GGVEW', 'GGVNS', VROC, 'VEW', 'VNS', 'THDG', 'ROLL', 'PITCH',
             'LAT', 'LON', 'VSPD')
VarListRef <- VarList
FI <- DataFileInfo (fname)
VarList <- VarListRef

if (!('GGVSPD' %in% FI$Variables)) {
  if ('GGVSPDB' %in% FI$Variables) {
    VarList [which (VarList == 'GGVSPD')] <- 'GGVSPDB'
  } else if ('VSPD_A' %in% FI$Variables) {
    VarList [which (VarList == 'GGVSPD')] <- 'VSPD_A'
  } else if ('VSPD_G' %in% FI$Variables) {
    VarList [which (VarList == 'GGVSPD')] <- 'VSPD_G'
  } else {
    print ('ERROR: no VSPD variable found')
    exit()
  }
}
for (Var in VarList) {
  if (!(Var %in% FI$Variables)) {
    print (sprintf (' required variable %s not found in file %s; skipping...', Var, fname))
    exit()
  }
}

netCDFfile <- nc_open (fnew, write=TRUE) 
Rate <- 1
Dimensions <- attr (D1, "Dimensions")
Dim <- Dimensions[["Time"]]
if ("sps25" %in% names (Dimensions)) {
  Rate <- 25
  Dim <- list(Dimensions[["sps25"]], Dimensions[["Time"]])
}
if ("sps50" %in% names (Dimensions)) {
  Rate <- 50
  Dim <- list(Dimensions[["sps50"]], Dimensions[["Time"]])
}
DATT <- D1  ## save to ensure that attributes are preserved

## variables to add to the netCDF file:
VarNew <- c('LATKF', 'LONKF', 'ALTKF', 'VEWKF', 'VNSKF', 'ROCKF', 'PITCHKF', 'ROLLKF', 'THDGKF',
            'WDKF', 'WSKF', 'WIKF', 'WICC', 'WDCC', 'WSCC', 'UXKF', 'VYKF', 
            'AKKF', 'SSKF')
VarOld <- c('LAT', 'LON', 'ALT', 'VEW', 'VNS', 'VSPD', 'PITCH', 'ROLL', 'THDG', 
            'WD', 'WS', 'WI', 'WIC', 'WDC', 'WSC', 'UXC', 'VYC', 'AKRD', 'SSRD')
VarUnits <- c('degrees', 'degrees', 'm', 'm/s', 'm/s', 'm/s', 'degrees', 'degrees', 'degrees',
              'degrees', 'm/s', 'm/s', 'm/s', 'degrees', 'm/s', 'm/s', 'm/s',
              'degrees', 'degrees')
VarStdName <- c('latitude, KF', 'longitude, KF', 'altitude MSL, KF',
                 'eastward groundspeed, KF', 'northward groundspeed, KF', 'rate of climb, KF',
                 'pitch, KF', 'roll, KF', 'heading, KF',
                 'wind direction, KF', 'wind speed, KF', 'vertical wind, KF', 
                 'vertical wind, KF and GGVSPD', 'WIC recalc', 'WDC recalc', 'WSC recalc',
                 'longitudinal wind, KF', 'lateral wind, KF', 
                 'angle of attack, KF', 'sideslip, KF')
VarLongName <- c('INS latitude, Kalman-filter-corrected',
                'INS longitude, Kalman-filter-corrected',
                'INS altitude, Kalman-filter-corrected',
                'INS eastward ground speed, Kalman-filter-corrected',
                'INS northward ground speed, Kalman-filter-corrected',
                'INS rate of climb, Kalman-filter-corrected',
                'INS aircraft pitch angle, Kalman-filter-corrected',
                'INS aircraft roll angle, Kalman-filter-corrected',
                'INS aircraft true heading angle, Kalman-filter-corrected',
                'horizontal wind direction, Kalman-filter-corrected',
                'horizontal wind speed, Kalman-filter-corrected',
                'vertical wind speed, Kalman-filter-corrected',
                'original WIC, recalculated',
                'original WDC, recalculated',
                'original WSC, recalculated',
                'longitudinal component of the horizontal wind, Kalman-filter-corrected',
                'lateral component of the horizontal wind, Kalman-filter-corrected',
                'angle of attack, from Kalman-filter routine', 'sideslip, from Kalman-filter routine')

## create the new variables
varCDF <- list ()
for (i in 1:length(VarNew)) {
  if (VarNew[i] == 'AKKF' && !UpdateAKRD) {next}
  if (VarNew[i] == 'SSKF' && !UpdateSSRD) {next}
  varCDF[[i]] <- ncvar_def (VarNew[i],  
                            units=VarUnits[i], 
                            dim=Dim, 
                            missval=as.single(-32767.), prec='float', 
                            longname=VarLongName[i])
  if (i == 1) {
    newfile <- ncvar_add (netCDFfile, varCDF[[i]])
  } else {
    newfile <- ncvar_add (newfile, varCDF[[i]])
  }
  ATV <- ncatt_get (netCDFfile, VarOld[i])
  copy_attributes (ATV, VarNew[i], newfile)
  ncatt_put (newfile, VarNew[i], attname="standard_name", 
             attval=VarStdName[i])
  if (Rate == 1) {
    ncvar_put (newfile, varCDF[[i]], D1[, VarNew[i]])
  } else if (Rate == 25) {
    ncvar_put (newfile, varCDF[[i]], D1[, VarNew[i]], count=c(25, nrow(D1)/25))
  }
}
nc_close (newfile)

