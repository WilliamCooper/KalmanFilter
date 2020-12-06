#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(Ranadu)

shinyUI(fluidPage(

    # Application title
    titlePanel("Output from the Kalman Filter"),

     sidebarLayout(
        sidebarPanel(
            fluidRow (
                column (7, selectInput (inputId='Project', label=NULL,
                                        choices=PJ, width='100px'))
            ),
            fluidRow(
                column(5,
                       numericInput(inputId='pn', label='plot #', value = 1, min = 1, max = 10, step = 1)),
                column(3,
                       checkboxInput('HR', '25-Hz'))
            ),
            fluidRow(
                column(4, selectInput('type', 'flt type', choices = c('rf', 'tf', 'ff'))),
                column(4, numericInput('FLT', 'flight:', value = 1)),
                column(4, numericInput(inputId='op', label = 'orig plot #', value = 1, min = 1,
                                       max = 8, step = 1))
            ),
            fluidRow(
                sliderInput("mag", "log magnification factor for difference variables:",
                            min=0, max=4, value=2, step=0.1)
            ),
            fluidRow(
                column(6,
                       selectInput('dir', 'data dir:', choices = c('KF', 'standard', 'KFoutput'))),
                column(6,
                       actionButton('yrst', 'reset limits'))
            ),
            fluidRow(
                sliderInput("times", label='time range', min=minT, max=maxT,
                            value=c(minT, maxT),
                            animate=TRUE,
                            step=60,
                            timeFormat='%T', dragRange=TRUE,
                            timezone='+0000')
            )
            
        ),

        mainPanel(
            # On some PCs, must use Fn F4
            tags$script(HTML("$(function() {
                                            $(document).keyup(function(e) {
                                            if (e.which == 115) {
                                            $('#resetT').click()
                                            }
                                            });
                                            })")),
            
            plotOutput("resultPlot", dblclick = "plot_dblclick",
                       brush=brushOpts(id='panel_brush', delay=1000, delayType='debounce',
                                       resetOnNew=TRUE)),
            imageOutput('KFplot', width=400, height=300)
        )
    )
))
