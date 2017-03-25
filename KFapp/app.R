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
## global stuff here
Project <- 'DEEPWAVE'
Flight <- 17
ALL <- FALSE
NEXT <- FALSE
newAK <- TRUE
newSS <- TRUE
simple <- FALSE
NSTEP=10
showPlots <- TRUE
viewPlot <- 1
genPlot <- TRUE
firstRun <- TRUE
PJ <- c('ARISTO2017', 'ORCAS', 'CSET', 'NOREASTER', 'HCRTEST',
        'DEEPWAVE', 'CONTRAST', 'SPRITE-II', 'MPEX', 'DC3', 'RICO',
        'TORERO', 'HIPPO-5', 'HIPPO-4', 'HIPPO-3', 'HIPPO-2',
        'HIPPO-1','PREDICT', 'START08', 'PACDEX', 'TREX')
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
  Fl <- sort (list.files (sprintf ("%s%s/", DataDirectory (), Project),
                          sprintf ("%srf..KF.nc", Project)), decreasing = TRUE)[1]
  if (is.na (Fl)) {
    Flight <- 1
  } else {
    Flight <- sub (".*rf", '',  sub ("KF.nc", '', Fl))
    Flight <- as.numeric(Flight)+1
  }
  return (Flight)
}

ShowProgress <- function(NSTEP, progress, Flight) {
  PLOOP <- 1
  TimeEstimate <- 1.5 * 60 * 9 * 10 / NSTEP  ## for 9-h flight
  while (PLOOP) {
    Sys.sleep (1)
    PLOOP <- PLOOP + 1
    if (PLOOP > TimeEstimate) {break}
    M <- system('tail -n 1 KFlog', intern=TRUE)
    if (grepl ('main loop', M)) {
      PLOOP <- FALSE
      progress$set(message = 'main-loop progress', 
                   detail = sprintf('flight %d', Flight), value=1)
      next
    }
    if (grepl('IRU loop', M)) {
      P <- sub ('.*loop ', '', sub ('% do.*', '', M))
      # print (sprintf ('progress is %d', as.integer(P)))
      progress$set(message = 'retrieve IRU msrmts', 
                   detail = sprintf('flight %d', Flight), value=as.integer (P))
    }
  }
  PLOOP <- 1
  while (PLOOP) {
    Sys.sleep (1)
    PLOOP <- PLOOP + 1
    if (PLOOP > TimeEstimate) {break}
    M <- system('tail -n 1 KFlog', intern=TRUE)
    if (grepl ('main loop is done', M) || grepl ('generating plots', M) || grepl ('making new', M)) {
      PLOOP <- FALSE
      next
    }
    if (grepl('Kalman-loop', M)) {
      P <- sub ('.*loop ', '', sub ('% do.*', '', M))
      # print (sprintf ('progress is %d', as.integer(P)))
      progress$set(message = 'main-loop progress', 
                   detail = sprintf('flight %d', Flight), value=as.integer (P))
    }
  }
  PLOOP <- 1
  while (PLOOP) {
    Sys.sleep (0.5)
    PLOOP <- PLOOP + 1
    if (PLOOP > 120) {break}
    M <- system('tail -n 1 KFlog', intern=TRUE)
    if (grepl ('plots generated', M) || grepl ('making new', M)) {
      PLOOP <- FALSE
      next
    }
    if (grepl('figures', M)) {
      P <- sub ('.*figures ', '', sub ('% do.*', '', M))
      # print (sprintf ('progress is %d', as.integer(P)))
      progress$set(message = 'generating figures', 
                   detail = sprintf('flight %d', Flight), value=as.integer (P))
    }
  }
  progress$set(message = 'creating new netCDF file', value=0)
  PLOOP <- 1
  while (PLOOP) {
    Sys.sleep (1)
    PLOOP <- PLOOP + 1
    if (PLOOP > 250) {break}
    M <- system('tail -n 1 KFlog', intern=TRUE)
    if (grepl ('finished', M)) {
      PLOOP <- FALSE
      next
    }
    if (grepl ('new-netcdf', M)) {
      P <- sub ('.*new-netcdf ', '', sub ('% do.*', '', M))
      # print (sprintf ('progress is %d', as.integer(P)))
      progress$set(message = 'creating new netCDF file', 
                 detail = sprintf('flight %d', Flight), value=as.integer (P))
    }
  }
  # PLOOP <- 1
  # while (PLOOP) {
  #   Sys.sleep (1)
  #   PLOOP <- PLOOP + 1
  #   if (PLOOP > 150) {break}
  #   M <- system('tail -n 1 KFlog', intern=TRUE)
  #   if (grepl ('finished', M)) {
  #     PLOOP <- FALSE
  #     next
  #   }
  # }
}

runScript <- function (ssn) {
  if (file.exists ('../KFplots/Position.png')) {
    system('rm ../KFplots/*png ../KFplots/*pdf')
  }
 
  if (ALL) {
    ## get list of files to process:
    Fl <- sort (list.files (sprintf ("%s%s/", DataDirectory (), Project),
                            sprintf ("%srf...nc", Project)))
    if (!is.na (Fl[1])) {
      for (Flt in Fl) {
        FltKF <- sub ('.nc$', 'KF.nc', Flt)
        if (file.exists (sprintf ("%s%s/%s", 
                                  DataDirectory (), Project, FltKF))) {next}
        if (file.exists ('../KFplots/Position.png')) {
          system('rm ../KFplots/*png ../KFplots/*pdf')
        }
        Flight <- sub('.*rf', '', sub ('.nc$', '', Flt))
        Flight <- as.numeric (Flight)
        updateNumericInput (ssn, 'Flight', value=Flight)
        progress$set(message = 'read data, initialize', 
                     detail = sprintf('flight %d', Flight),
                     value=0)
        cmd <- sprintf('Rscript ../KalmanFilter.R %s %d %s %s %s %d %s | tee -a KFlog', 
                       Project, Flight, newAK, newSS, simple, NSTEP, genPlot)
        print (sprintf ('run commnad: %s', cmd))
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
    cmd <- sprintf('Rscript ../KalmanFilter.R %s %d %s %s %s %d %s | tee -a KFlog', 
                   Project, Flight, newAK, newSS, simple, NSTEP, genPlot)
    system (cmd, wait=FALSE)
    ShowProgress (NSTEP, progress, Flight)
  } else {
    progress$set(message = 'read data, initialize', 
                 detail = sprintf('flight %d', Flight),
                 value=0)
    cmd <- sprintf('nice -10 Rscript ../KalmanFilter.R %s %d %s %s %s %d %s | tee -a KFlog', 
                   Project, Flight, newAK, newSS, simple, NSTEP, genPlot)
    system (cmd, wait=FALSE)
    ShowProgress (NSTEP, progress, Flight)
  }
  return()
}

# Define UI for application that controls KalmanFilter
ui <- fluidPage(
  includeCSS("../www/styles.css"),
  tags$head(tags$script(HTML('Shiny.addCustomMessageHandler("jsCode",function(message) {eval(message.value);});'))),
  # Application title
  # titlePanel("Kalman-Filter Processor"),
  fluidRow (
    column (6, actionButton ('Run', h3("Click Here to Run the Kalman Processor"),
                             style="color: #fff; background-color: #337ab7; border-color: #2e6da4"))
  ),
  sidebarLayout(
    sidebarPanel(h4('Run Arguments:'),
                 fluidRow (
                   column (7, selectInput (inputId='Project', label=NULL,
                                           choices=PJ, selected=Project, width='100px')),
                   column (5, checkboxInput ('simple', label='Only Simple?', value=FALSE))),
                 fluidRow (
                   column (5, numericInput (inputId='Flight', label='Flight', value=Flight,
                                            min=1, max=99, step=1, width='80px')),
                   column (3, checkboxInput ('ALL', label='ALL?',
                                             value=FALSE)),
                   column (3, checkboxInput ('NEXT', label='Next',
                                             value=FALSE))
                 ),
                 fluidRow (
                   column (3, checkboxInput ('newAK', label='AK?', value=newAK)),
                   column (3, checkboxInput ('newSS', label='SS?', value=newSS)),
                   column (6, numericInput ('NSTEP', label='step (s)', value=NSTEP, min=5,
                                            max=60, step=1, width='80pc'))
                 ),
                 fluidRow (
                   column (4, checkboxInput ('genPlot', label='plots?', value=genPlot)),
                   column (4, numericInput ('viewPlot', label='view', value=1,
                                            min=1, max=8, step=1, width='80'))
                 )
    ),
    
    mainPanel(
      textOutput('runPar'),
      plotOutput("resultPlot")
    )
  )
)

server <- function(input, output, session) {
  
  observeEvent (input$Run, {
    runScript(session)
    progress$set(message = 'processing is complete', detail=sprintf ('Flight %d', Flight), value=99)
    updateNumericInput(session, 'viewPlot', value=3)
  })
  
  exprProject <- quote ({
    if (input$Project != Project) {
      Project <<- input$Project
    }
  })
  obsProject <- observe (exprProject, quoted=TRUE)
  
  exprFlight <- quote ({
    if (input$Flight != Flight) {
      Flight <<- input$Flight
      progress$set(message = 'ready to run', detail = sprintf('flight %d', Flight))
    }
  })
  obsFlight <- observe (exprFlight, quoted=TRUE) 
  
  exprNSTEP <- quote ({
    if (input$NSTEP != NSTEP) {
      NSTEP <<- input$NSTEP
    }
  })
  obsNSTEP <- observe (exprNSTEP, quoted=TRUE)  
  
  exprNEXT <- quote ({
    if (input$NEXT != NEXT) {
      NEXT <<- input$NEXT
      if (NEXT) {Flight <- 'NEXT'}
    }
  })
  obsNEXT <- observe (exprNEXT, quoted=TRUE)
  
  exprALL <- quote ({
    if (input$ALL != ALL) {
      ALL <<- input$ALL
      if (ALL) {
        js_string <- 'alert("For ALL, see instructions. Previously generated files are skipped; delete to reprocess. This run may take a very long time.")'
        session$sendCustomMessage(type='jsCode', list(value = js_string))
      }
    }
  })
  obsALL <- observe (exprALL, quoted=TRUE)
  
  exprNewAK <- quote ({
    if (input$newAK != newAK) {
      newAK <<- input$newAK
    }
  })
  obsNewAK <- observe (exprNewAK, quoted=TRUE)  
  
  exprNewSS <- quote ({
    if (input$newSS != newSS) {
      newSS <<- input$newSS
    }
  })
  obsNewSS <- observe (exprNewSS, quoted=TRUE)  
  
  exprSimple <- quote ({
    if (input$simple != simple) {
      simple <<- input$simple
    }
  })
  obsSimple <- observe (exprSimple, quoted=TRUE)
  
  exprGenPlot <- quote ({
    if (input$genPlot != genPlot) {
      genPlot <<- input$genPlot
    }
  })
  obsGenPlot <- observe (exprGenPlot, quoted=TRUE)## get the wind variables:
  
  exprViewPlot <- quote ({
    if (input$viewPlot != viewPlot) {
      viewPlot <<- input$viewPlot
    }
  })
  obsViewPlot <- observe (exprViewPlot, quoted=TRUE)
  
  # dat <- data.frame(x = numeric(0), y = numeric(0))
  # 
  # withProgress(message = 'Making plot', value = 0, {
  #   # Number of times we'll go through the loop
  #   n <- 10
  #   
  #   for (i in 1:n) {
  #     # Each time through the loop, add another row of data. This is
  #     # a stand-in for a long-running computation.
  #     dat <- rbind(dat, data.frame(x = rnorm(1), y = rnorm(1)))
  #     
  #     # Increment the progress bar, and update the detail text.
  #     incProgress(1/n, detail = paste("Doing part", i))
  #     
  #     # Pause for 0.1 seconds to simulate a long computation.
  #     Sys.sleep(0.1)
  #   }
  # })
  output$runPar <- renderText({
    if (!progressExists) {
      progress <- Progress$new(session, min=0, max=100)
    } 
    # on.exit(progress$close())
    
    progress$set(message = 'ready to run')
    progress$set (value=1)
    progress <<- progress
    progressExists <<- TRUE
    messg <<- sprintf ("%srf%02d.nc dt=%d AK=%s SS=%s Simple=%s",
                       input$Project, input$Flight, input$NSTEP, input$newAK, input$newSS, input$simple)
    messg <- NULL
  })
  
  output$resultPlot <- renderImage({
    
    plotNo <- input$viewPlot
    pname <- c('../KFplots/Position.png',
               '../KFplots/Velocity.png',
               '../KFplots/AAlframe.png',
               '../KFplots/AAaframe.png',
               '../KFplots/HDG.png',
               '../KFplots/Wind.png',
               '../KFplots/Wind2.png',
               '../KFplots/HCPlot.png')
    # Return a list containing the filename
    list(src = pname[plotNo],
         contentType = 'image/png',
         width = 900,
         height = 600,
         alt = "Waiting for Plots")
  }, deleteFile=FALSE)
}

# Run the application 
shinyApp(ui, server)

