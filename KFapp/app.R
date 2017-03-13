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
Project <- 'CSET'
Flight <- 1
ALL <- FALSE
Next <- FALSE
newAK <- TRUE
newSS <- TRUE
simple <- FALSE
NSTEP=10
showPlots <- TRUE
viewPlot <- 1
genPlot <- TRUE
firstRun <- TRUE
PJ <- c('ORCAS', 'CSET', 'NOREASTER', 'HCRTEST',
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

runScript <- function () {
  system('rm ../KFplots/*png')
  cmd <- sprintf('Rscript ../KalmanFilter.R %s %d %s %s %s %d', Project, Flight, newAK, newSS, simple, NSTEP)
  system(cmd)
  return()
}

# Define UI for application that controls KalmanFilter
ui <- fluidPage(
  
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
	        column (5, numericInput (inputId='Flight', label='Flight', value=1,
	                                        min=1, max=99, step=1, width='80px')),
	        column (3, checkboxInput ('ALL', label='ALL?',
	                                  value=FALSE)),
	        column (3, checkboxInput ('next', label='Next',
	                                  value=FALSE))
	      ),
	      fluidRow (
	        column (3, checkboxInput ('newAK', label='AK?', value=TRUE)),
	        column (3, checkboxInput ('newSS', label='SS?', value=TRUE)),
	        column (6, numericInput ('NSTEP', label='step (s)', value=10, width='80pc'))
	      ),
	      fluidRow (
	        column (4, checkboxInput ('genPlot', label='plots?', value=TRUE)),
	        column (4, numericInput ('viewPlot', label='view', value=1,
	                                 min=1, max=7, step=1, width='80'))
	      )
	),
    
    mainPanel(
      textOutput('runPar'),
      plotOutput("resultPlot")
    )
  )
)

server <- function(input, output) {
  
  observeEvent (input$Run, {
    runScript()
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
    }
  })
  obsFlight <- observe (exprFlight, quoted=TRUE) 
  
  exprNSTEP <- quote ({
    if (input$NSTEP != NSTEP) {
      NSTEP <<- input$NSTEP
    }
  })
  obsNSTEP <- observe (exprNSTEP, quoted=TRUE)  

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
  obsGenPlot <- observe (exprGenPlot, quoted=TRUE)
  
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
    messg <<- sprintf ("%srf%02d.nc dt=%d AK=%s SS=%s Simple=%s", 
                       input$Project, input$Flight, input$NSTEP, input$newAK, input$newSS, input$simple)
    messg
  })
  
  output$resultPlot <- renderImage({
      plotNo <- input$viewPlot
      pname <- c('../KFplots/Position.png',
                   '../KFplots/Velocity.png',
                   '../KFplots/AAlframe.png',
                   '../KFplots/AAaframe.png',
                   '../KFplots/HDG.png',
                   '../KFplots/Wind.png',
                   '../KFplots/Wind2.png')
      # Return a list containing the filename
      list(src = pname[plotNo],
           contentType = 'image/png',
           width = 600,
           height = 450,
           alt = "Waiting for Plots")
    
  }, deleteFile=FALSE)
}

# Run the application 
shinyApp(ui = ui, server = server)

