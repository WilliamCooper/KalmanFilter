
# Define UI for the application that controls the KalmanFilter:
ui <- fluidPage(
  includeCSS("./www/styles.css"),
  tags$head(tags$script(HTML('Shiny.addCustomMessageHandler("jsCode",function(message) {eval(message.value);});'))),
  # Application title
  # titlePanel("Kalman-Filter Processor"),
  fluidRow (
    column (1, actionButton ('info', label='help')),
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
                   column (3, selectInput('type', 'type flt', choices=c('rf', 'tf', 'ff'))),
                   column (5, numericInput (inputId='Flight', label='Flight', value=Flight,
                                            min=1, max=99, step=1, width='80px')),
                   column (2, checkboxInput (inputId='HR', label='25Hz', value=FALSE))),
                 fluidRow (
                   column (4, checkboxInput ('ALL', label='ALL?',
                                             value=FALSE)),
                   column (4, checkboxInput ('NEXT', label='Next',
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
                 ), 
                 fluidRow (
                   selectInput('dir', 'output', choices=c('KF', 'standard', 'KFoutput/'))
                 )
    ),
   
    mainPanel(
      textOutput('runPar'),
      plotOutput("resultPlot")
    )
  )
)

