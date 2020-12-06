#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

shinyServer(function(input, output, session) {
    
    HRD <- reactiveVal(data.frame())
    LRD <- reactiveVal(data.frame())
    observeEvent (input$yrst, {
        panelylim <<- NA
        updateSliderInput (session, 'times', value=c(minT, maxT))
    })
    exprProject <- quote ({
        if (input$Project != Project) {
            Project <<- input$Project
        }
    })
    obsProject <- observe (exprProject, quoted=TRUE)
    observeEvent (input$panel_brush, {
        Brush <<- input$panel_brush
        xmin <- as.integer(input$panel_brush$xmin)
        xmax <- as.integer(input$panel_brush$xmax)
        ymin <- input$panel_brush$ymin
        ymax <- input$panel_brush$ymax
        panelylim <<- c(ymin, ymax)
        # if (input$ybrush) {
        #     if (input$panel_brush$outputId == 'display') {
        #         panel1ylim <<- c(ymin, ymax)
        #     }
        # }
        T1 <- as.POSIXlt(xmin, origin='1970-01-01', tz='UTC')
        T2 <- as.POSIXlt(xmax, origin='1970-01-01', tz='UTC')
        TB1 <- T1$hour*10000 + T1$min*100 + T1$sec
        TB2 <- T2$hour*10000 + T2$min*100 + T2$sec
            print (sprintf ('display brush times are %d %d', TB1, TB2))
            print (sprintf ('display is %s, ordinate limits are %f %f',
                            input$panel_brush$outputId, ymin, ymax))
        updateSliderInput (session, 'times', value=c(T1, T2))
        times <<- c(T1, T2)
    })
    
    output$resultPlot <- renderPlot({
        
        if(input$HR) {
            if((lastData == 'HR') && (input$FLT == lastFLT) && (input$type == lastType)) {
                Data <- HRD()
            } else {
                lastFLT <<- input$FLT
                lastType <<- input$type
                fname <- setFileName('WCR-TEST', sprintf('%s%02dhKF', input$type,  
                                                         as.numeric(input$FLT)))
                Data <- getNetCDF(fname,
                                  standardVariables(unique(c(STDVARS, KFVARS))))
                if(any(Data$TASX > 80)) {
                    itas <- which(Data$TASX > 80)
                    Data <- Data[itas[1]:itas[length(itas)], ]
                }
                lastData <<- 'HR'
                HRD(Data)
            } 
        } else {
            if((lastData == 'LR') && (input$FLT == lastFLT) && (input$type == lastType)) {
                Data <- LRD()
            } else {
                lastFLT <<- input$FLT
                lastType <<- input$type
                fname <- setFileName('WCR-TEST', sprintf('%s%02dKF', input$type,
                                                         as.numeric(input$FLT)))
                ## modify to access the requested directory:
                if (input$dir != 'standard') {
                    if (input$dir == 'KF') {
                        fname <- sub(paste0(Project, '/'), paste0(Project, '/KF/'), fname)
                    } else if (input$dir == 'KFoutput') {
                        fname <- paste0('../KalmanFilter/KFoutput/', sub('.*/', '', fname))
                    }
                }
                print (sprintf ('loading data from %s', fname))
                Data <- getNetCDF(fname,
                                  standardVariables(unique(c(STDVARS, KFVARS))))
                if(any(Data$TASX > 80)) {
                    itas <- which(Data$TASX > 80)
                    Data <- Data[itas[1]:itas[length(itas)], ]
                }
                lastData <<- 'LR'
                LRD(Data)
            }
        }
        
        Vars1 <- c('LAT', 'LATKF', 'GGLAT', 'LATdiff')
        Vars2 <- c('LON', 'LONKF', 'GGLON', 'LONdiff')
        Vars3 <- c('PITCH', 'PITCHKF', 'PITCHdiff')
        Vars4 <- c('ROLL', 'ROLLKF', 'ROLLdiff')
        Vars5 <- c('THDG', 'THDGKF', 'THDGdiff')
        Vars6 <- c('WDC', 'WDKF', 'WDdiff')
        Vars7 <- c('WSC', 'WSKF', 'WSdiff')
        Vars8 <- c('WIC', 'WIKF', 'WIdiff')
        Vars9 <- c('AKRD', 'AKKF', 'AKdiff')
        Vars10 <- c('SSRD', 'SSKF', 'SSdiff')
        Data$PITCHdiff <- with(Data, PITCHKF - PITCH)
        Data$ROLLdiff <- with(Data, ROLLKF - ROLL)
        Data$THDGdiff <- with(Data, THDGKF - THDG)
        Data$AKdiff <- with(Data, AKKF - AKRD)
        Data$SSdiff <- with(Data, SSKF - SSRD)
        Data$LATdiff <- with(Data, LATKF - LAT)
        Data$LONdiff <- with(Data, LONKF - LON)
        Data$WDdiff <- with(Data, WDKF - WDC)
        Data$WSdiff <- with(Data, WSKF - WSC)
        Data$WIdiff <- with(Data, WIKF - WIC)
        
        ylim <- NULL
        Vars <- switch(input$pn, 
                       '1'=Vars1,
                       '2'=Vars2,
                       '3'=Vars3,
                       '4'=Vars4,
                       '5'=Vars5,
                       '6'=Vars6,
                       '7'=Vars7,
                       '8'=Vars8,
                       '9'=Vars9,
                       '10'=Vars10)
        Rline <- mean(Data[, Vars[1]], na.rm=TRUE)
        iw <- which(grepl('diff', Vars))
        Data[, Vars[iw]] <- 10.^input$mag * Data[, Vars[iw]] + Rline
        lwd <- c(2,2,2,2)
        ltyp <- c(1,1,1,1)
        col <- c('blue', 'red', 'forestgreen', 'forestgreen')
        lwd[iw] <- 3
        ltyp[iw] <- 4
        col[iw] <- 'magenta'
        if (any(grepl('WDKF', Vars)) || any(grepl('THDG', Vars))) {
            ylim <- c(0, 360)
        } else {
            ylim <- with(Data, c(min(get(Vars[1]), na.rm=TRUE), max(get(Vars[1]), na.rm=TRUE)))
            yldiff <- with(Data, c(min(get(Vars[iw]), na.rm=TRUE), max(get(Vars[iw]), na.rm=TRUE)))
            if (yldiff[1] < ylim[1]) {ylim[1] <- yldiff[1]}
            if (yldiff[2] > ylim[2]) {ylim[2] <- yldiff[2]}
        }
        if (!is.na(panelylim[1])) {ylim <- panelylim}
        SE <- getStartEnd(Data)
        print(input$times)
        se <- c(formatTime(input$times[1]), formatTime(input$times[2]))
        se[2] <- se[2] - 1
        print(se)
        Data %>% selectTime(se[1], se[2]) %>% 
            select(Time, all_of(Vars)) %>% 
            plotWAC(col=col, ylim = ylim)
                # plotWAC(lwd=lwd, lty=ltyp, col = col)
        title(sprintf('diff variable: multiplied by %.1f, plotted wrt dotted line', 10^input$mag),
              line = 1)
        abline(h = Rline, lty=3)
    })
    
    output$KFplot <- renderImage({
        ## Does the plot file exist?
        pfile <- paste0('../..//KalmanFilter/KFplots/', Project, 'Plots',input$type,
                        sprintf('%02d', as.numeric(input$FLT)), '.tgz')
        print (sprintf ('with KF, pfile is %s', pfile))
        print (sprintf ('wd is %s', getwd()))
        print (sprintf ('exists? %s', file.exists(pfile)))
        if (file.exists(pfile)) {
            if (!dir.exists('KFplts')) {
                dir.create('KFplts')   
            }
            setwd('./KFplts')
            system(sprintf ('tar xvfz ../%s', pfile))
            setwd('..')
        }
        pf <- switch(input$op,
                     '1'='Position.png',
                     '2'='Velocity.png',
                     '3'='AAlframe.png',
                     '4'='AAaframe.png',
                     '5'='HDG.png',
                     '6'='Wind.png',
                     '7'='Wind2.png',
                     '8'='HCPlot.png')
        
        list(src = sprintf('KFplts/%s', pf), contentType = 'image/png', width=500)
    }, deleteFile = FALSE)
    
    
    
    
})

