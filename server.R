
server <- function(input, output, session) {
  
  observeEvent (input$Run, {
    runScript(session)
    progress$set(message = 'processing is complete', detail=sprintf ('Flight %d', Flight), value=99)
    updateNumericInput(session, 'viewPlot', value=3)
  })
  
  observeEvent (input$info, {
    showModal(modalDialog(
      includeHTML('HTML/guide.html'),
      title = "User Info for the Kalman Filter",
      size='m',
      easyClose = TRUE
    ))
  })
  
  exprProject <- quote ({
    if (input$Project != Project) {
      Project <<- input$Project
    }
  })
  obsProject <- observe (exprProject, quoted=TRUE)
  
  exprOutDir <- quote ({
    if (input$dir == 'KF') {
      outputDirectory <<- paste0(DataDirectory(), Project, '/KF/')
    } else if (input$dir == 'standard') {
      outputDirectory <<- paste0(DataDirectory(), Project, '/')
    } else {
      outputDirectory <<- input$dir
    }
  })
  obsOutDir <- observe (exprOutDir, quoted=TRUE)
  
  exprtypeFlight <- quote ({
    if (input$type != typeFlight) {
      typeFlight <<- input$type
      flt <<- sprintf('%s%02d', typeFlight, Flight)
      progress$set(message = 'ready to run', detail = sprintf('flight %s%02d', 
                                                              typeFlight,Flight))
    }
  })
  obstypeFlight <- observe (exprtypeFlight, quoted=TRUE) 
  
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
      if (NEXT) {
        Flight <- 'NEXT'
        ## ensure that ALL is FALSE
        updateCheckboxInput(session, 'ALL', value = FALSE)
        ## Check that the next file is present:
        if (is.na(getNext(input$Project))) {
          js_string <- 'alert("The NEXT file is not present. Resetting.")'
          session$sendCustomMessage(type='jsCode', list(value = js_string))
          NEXT <<- FALSE
          updateCheckboxInput(session, 'NEXT', value = FALSE)
        }
      }
    }
  })
  obsNEXT <- observe (exprNEXT, quoted=TRUE)
  
  exprALL <- quote ({
    if (input$ALL != ALL) {
      ALL <<- input$ALL
      if (ALL) {
        Flight <- 'ALL'
        updateCheckboxInput(session, 'NEXT', value = FALSE)
        js_string <- 'alert("For ALL, see instructions. Previously generated files are skipped; delete to reprocess. This run may take about 5 min per flight.")'
        session$sendCustomMessage(type='jsCode', list(value = js_string))
      }
    }
  })
  obsALL <- observe (exprALL, quoted=TRUE)
  
  exprHR <- quote ({
    if (input$HR != HighRate) {
      HighRate <<- input$HR
      if (HighRate) {
        js_string <- 'alert("For 25-Hz processing, see instructions. A 1-Hz file with the analogous name must be processed first, and that processing should use the AK and SS checkboxes.")'
        session$sendCustomMessage(type='jsCode', list(value = js_string))
      }
    }
  })
  obsHR <- observe (exprHR, quoted=TRUE)
  
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
    messg <<- sprintf ("%s%s%02d.nc dt=%d AK=%s SS=%s Simple=%s",
                       input$Project, input$type, input$Flight, input$NSTEP, 
                       input$newAK, input$newSS, input$simple)
    messg <- NULL
  })
  
  output$resultPlot <- renderImage({
    
    plotNo <- input$viewPlot
    if (plotNo == 8) {
      if (!file.exists('KFplots/HCPlot.png')) {
        png(file = "KFplots/HCPlot.png")
        plot(c(0,10), c(0,10), type = 'n', xlab='', ylab='')
        text(5,5, labels='no CorrectHeading() plot:\nlikely reason is too few qualifying turns')
        dev.off()
      }
    }
    pname <- c('KFplots/Position.png',
               'KFplots/Velocity.png',
               'KFplots/AAlframe.png',
               'KFplots/AAaframe.png',
               'KFplots/HDG.png',
               'KFplots/Wind.png',
               'KFplots/Wind2.png',
               'KFplots/HCPlot.png')
    # Return a list containing the filename
    list(src = pname[plotNo],
         contentType = 'image/png',
         width = 900,
         height = 600,
         alt = "Waiting for Plots")
  }, deleteFile=FALSE)
}

