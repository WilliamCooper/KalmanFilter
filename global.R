#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

suppressMessages (suppressWarnings (
  library(Ranadu, quietly=TRUE, warn.conflicts=FALSE))
)
require(numDeriv)
## global stuff here
Project <- 'WCR-TEST'
Flight <- 1
ALL <- FALSE
ALLHR <- FALSE
NEXT <- FALSE
NEXTHR <- FALSE
newAK <- TRUE
newSS <- TRUE
simple <- FALSE
NSTEP = 10
HighRate <- FALSE
showPlots <- TRUE
viewPlot <- 1
genPlot <- TRUE
firstRun <- TRUE
## See if this is being run as user "shiny":
SHINY <- ifelse(Sys.info()['user'] == 'shiny', TRUE, FALSE)

PJ <- c('WCR-TEST', 'MethaneAIR', 'ACCLIP-TEST', 'OTREC', 'ECLIPSE2019', 
        'OTREC-TEST', 'WECAN', 'SOCRATES', 'WECAN-TEST', 'ECLIPSE', 
        'ARISTO2017', 'ORCAS', 'CSET', 'NOREASTER', 'HCRTEST',
        'DEEPWAVE', 'CONTRAST', 'SPRITE-II', 'MPEX', 'DC3', 'RICO',
        'TORERO', 'HIPPO-5', 'HIPPO-4', 'HIPPO-3', 'HIPPO-2',
        'HIPPO-1','PREDICT', 'START08', 'PACDEX', 'TREX')
## Replace this by constructing a list of available projects
## by searching the DataDirectory(). That way a new project will
## be incorporated automatically.
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

messg <- 'waiting for run'
progressExists <- FALSE

getNext <- function(Project) {
  if (HighRate) {
    Fl <- sort (list.files (sprintf ("%s%s/", DataDirectory (), Project),
                            sprintf ("%srf..hKF.nc", Project)), decreasing = TRUE)[1]
    ## Consider also if a processed version exists in KFoutput:
    Fl <- sort (c(Fl, list.files ("KFoutput", sprintf ("%srf..hKF.nc", Project))),
                decreasing = TRUE)[1]
  } else {
    Fl <- sort (list.files (sprintf ("%s%s/", DataDirectory (), Project),
                            sprintf ("%srf..KF.nc", Project)), decreasing = TRUE)[1]
    ## Consider also if a processed version exists in KFoutput:
    Fl <- sort (c(Fl, list.files ("KFoutput", sprintf ("%srf..KF.nc", Project))),
                decreasing = TRUE)[1]
  }
  if (is.na (Fl)) {
    Flight <- 1
  } else {
    Flight <- sub (".*rf", '',  sub ("KF.nc", '', Fl))
    Flight <- sub('h$', '', Flight)
    Flight <- as.numeric(Flight)+1
    ## If the next flight doesn't exist, return NA
    checkFile <- paste0(DataDirectory(), Project, '/', Project, 'rf', 
                        sprintf ('%02d', Flight), '.nc')
    if (HighRate) {
      checkFile <- sub('.nc$', 'h.nc', checkFile)
    }
    if (!file.exists(checkFile)) {
      return (NA)
    }
  }
  return (Flight)
}

ShowProgress <- function(NSTEP, progress, Flight) {
  flt <- Flight
  if(is.numeric(Flight)) {
    flt <- sprintf('rf%02d', Flight)
  }
  PLOOP <- 1
  TimeEstimate <- 30 * 9 * 10 / NSTEP  ## for 9-h flight
  while (PLOOP) {
    Sys.sleep (1)
    PLOOP <- PLOOP + 1
    if (PLOOP > TimeEstimate) {break}
    M <- system('tail -n 1 KFlog', intern=TRUE)
    if (grepl ('Kalman-loop', M)) {
      PLOOP <- FALSE
      progress$set(message = 'main-loop progress', 
                   detail = sprintf('flight %s', flt), value=1)
      break
    }
    if (grepl('IRU loop', M)) {
      P <- sub ('.*loop ', '', sub ('% do.*', '', M))
      # print (sprintf ('progress is %d', as.integer(P)))
      progress$set(message = 'retrieve IRU msrmts', 
                   detail = sprintf('flight %s', flt), value=as.integer (P))    
    }
  }
  PLOOP <- 1
  while (PLOOP) {
    Sys.sleep (1)
    PLOOP <- PLOOP + 1
    if (PLOOP > TimeEstimate) {break}
    M <- system('tail -n 1 KFlog', intern=TRUE)
    if (grepl ('main loop is done', M) || grepl ('generating plots', M) || 
        grepl ('making new', M) || grepl('figures', M) || grepl('new-netcdf', M)) {
      PLOOP <- FALSE
      break
    }
    if (grepl('Kalman-loop', M)) {
      P <- sub ('.*loop ', '', sub ('% do.*', '', M))
      # print (sprintf ('progress is %d', as.integer(P)))
      progress$set(message = 'main-loop progress',
                   detail = sprintf('flight %s', flt), value=as.integer (P))
    }
  }
  PLOOP <- 1
  while (PLOOP) {
    Sys.sleep (0.5)
    PLOOP <- PLOOP + 1
    if (PLOOP > 120) {break}
    M <- system('tail -n 1 KFlog', intern=TRUE)
    if (grepl ('plots generated', M) || grepl ('making new', M) || grepl ('new-netcdf', M)) {
      PLOOP <- FALSE
      break
    }
    if (grepl('figures', M)) {
      P <- sub ('.*figures ', '', sub ('% do.*', '', M))
      # print (sprintf ('progress is %d', as.integer(P)))
      progress$set(message = 'generating figures', 
                   detail = sprintf('flight %s', flt), value=as.integer (P))
    }
  }
  progress$set(message = 'creating new netCDF file', value=0)
  PLOOP <- 1
  while (PLOOP) {
    Sys.sleep (1)
    PLOOP <- PLOOP + 1
    if (PLOOP > TimeEstimate) {break}
    M <- system('tail -n 1 KFlog', intern=TRUE)
    if (grepl ('finished', M)) {
      PLOOP <- FALSE
      break
    }
    if (grepl ('new-netcdf', M)) {
      P <- sub ('.*new-netcdf ', '', sub ('% do.*', '', M))
      # print (sprintf ('progress is %d', as.integer(P)))
      progress$set(message = 'creating new netCDF file', 
                   detail = sprintf('flight %s', flt), value=as.integer (P))
    }
  }
  PLOOP <- 1
  while (PLOOP) {
    Sys.sleep (1)
    PLOOP <- PLOOP + 1
    if (PLOOP > TimeEstimate) {break}
    M <- system('tail -n 1 KFlog', intern=TRUE)
    if (grepl ('finished', M)) {
      PLOOP <- FALSE
      break
    }
  }
}

runScript <- function (ssn) {
  if (file.exists ('../KFplots/Position.png')) {
    system('rm ../KFplots/*png ../KFplots/*pdf')
  }
  if (file.exists ('KFlog')) {
    system('rm KFlog')
  }
  
  if (ALL) {
    ## get list of files to process:
    if (HighRate) {
      Fl <- sort (list.files (sprintf ("%s%s/", DataDirectory (), Project),
                              sprintf ("%srf..h.nc", Project)))
    } else {
      Fl <- sort (list.files (sprintf ("%s%s/", DataDirectory (), Project),
                              sprintf ("%srf...nc", Project)))
    }
    if (!is.na (Fl[1])) {
      print('in runScript, with ALL TRUE, files to be processed:')
      print(Fl)
      for (Flt in Fl) {
        FltKF <- sub ('.nc$', 'KF.nc', Flt)
        if (file.exists (sprintf ("%s%s/%s", 
                                  DataDirectory (), Project, FltKF))) {next}
        if (file.exists ('KFplots/Position.png')) {
          system('rm KFplots/*png KFplots/*pdf')
        }
        Flight <- sub('h', '', sub('.*rf', '', sub ('.nc$', '', Flt)))
        Flight <- as.numeric (Flight)
        updateNumericInput (ssn, 'Flight', value=Flight)
        if (HighRate) {
          progress$set(message = 'read data, initialize', 
                       detail = sprintf('flight %dh', Flight),
                       value=0)
          cmd <- sprintf('Rscript KalmanFilter.R %s %s %s %s %s %d %s | tee -a KFlog', 
                         Project, sprintf('rf%02dh', Flight), newAK, newSS, simple, NSTEP, genPlot)
        } else {
          progress$set(message = 'read data, initialize', 
                       detail = sprintf('flight %d', Flight),
                       value=0)
          cmd <- sprintf('Rscript KalmanFilter.R %s %d %s %s %s %d %s | tee -a KFlog', 
                         Project, Flight, newAK, newSS, simple, NSTEP, genPlot)
        }
        print (sprintf ('run command: %s', cmd))
        system (cmd, wait=FALSE)
        ShowProgress (NSTEP, progress, Flight)
      }
    }
  } else if (NEXT) {
    Flight <- getNext(Project)
    updateNumericInput (ssn, 'Flight', value=Flight)
    progress$set(message = 'read data, initialize', 
                 detail = sprintf('flight %d', Flight),
                 value=0)
    if (HighRate) {
      ## LowRate processed file must be present:
      chkFile <- paste0(DataDirectory(), Project, '/', Project, sprintf('rf%02dKF.nc', Flight))
      if (!file.exists(chkFile)) {
        print (sprintf ('unable to process next HR file; need LR file %s', chkfile))
        return()
      } else {
        cmd <- sprintf('Rscript KalmanFilter.R %s %s %s %s %s %d %s | tee -a KFlog', 
                       Project, sprintf('rf%02dh', Flight), newAK, newSS, simple, 
                       NSTEP, genPlot)
      }
    } else {
      cmd <- sprintf('Rscript KalmanFilter.R %s %d %s %s %s %d %s | tee -a KFlog', 
                     Project, Flight, newAK, newSS, simple, NSTEP, genPlot)
    }
    system (cmd, wait=FALSE)
    ShowProgress (NSTEP, progress, Flight)
  } else {
    if (HighRate) {
      progress$set(message = 'read data, initialize', 
                   detail = sprintf('flight rf%02dh', Flight),
                   value=0)
      cmd <- sprintf('nice -10 Rscript KalmanFilter.R %s %s %s %s %s %d %s | tee -a KFlog', 
                     Project, sprintf('rf%02dh', Flight), newAK, newSS, simple, NSTEP, genPlot)
    } else {
      progress$set(message = 'read data, initialize', 
                   detail = sprintf('flight %d', Flight),
                   value=0)
      cmd <- sprintf('nice -10 Rscript KalmanFilter.R %s %d %s %s %s %d %s | tee -a KFlog', 
                     Project, Flight, newAK, newSS, simple, NSTEP, genPlot)
    }
    system (cmd, wait=FALSE)
    ShowProgress (NSTEP, progress, Flight)
  }
  ## When the high-rate processor finishes, unpack the appropriate plots to KFplots:
  plotFile <- paste0(Project, 'Plotsrf', sprintf('%02d', Flight), '.tgz')
  system(sprintf('cd KFplots;tar xfz %s', plotFile))
  return()
}

