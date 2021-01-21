# Load libraries, assuming these have been installed
# First, gslea

library(ccca) # includes dependency: datatable, but nothing else
library(shinyjs)
library(shiny)

ui <- fluidPage(
  titlePanel('Climate change conditioned advice (CCCA) for fisheries'),
  navlistPanel(
    widths=c(2,10),
    id='tabs',
    tabPanel(title = 'Step 1. Run the P/B~E model',
      sidebarPanel(
        p('First, `ccca` calculates the annual P/B of the population from a survey index time series. This is then modelled as a function of an `E` variable, either ecological or environmental.'),
        # h3('Additional datasets here? Upload via CSV? Or only turbot?'),
        hr(), 
        checkboxInput('turbot', 'Use the default turbot dataset', value=T),
        p('Or upload a new dataset as a CSV file, making sure the file includes columns "Year", "Catch", "Index" and "E" at least.'),
        # p('Any uploaded dataset must have columns, at least,'),
        # tags$b('Year, Catch, Index, E'),
        fileInput('paramfile', label='Upload CSV dataset', multiple=F, 
          accept=c("text/csv",
                   "text/comma-separated-values,text/plain",
                    ".csv"), 
          placeholder = "No file uploaded, using default dataset."),
        hr(),
        uiOutput('input.refyears'),
        p('Define additional parameters'),
        numericInput('q', 'Catchability (q)', value=1, min=0),
        # hr(),
        p('Model the relationship between stock productivity and E'),
        selectInput('modtype', 'Form of the model', 
          choices = c(
            'Polynomial' = 'poly', 
            'GAM' = 'gam',
            'GAM (adaptive)' = 'gam.adaptive',
            'SCAM (MPI)' = 'mpi',
            'SCAM (MPD)' = 'mpd',
            'SCAM (CX)' = 'cx',
            'SCAM (CV)' = 'cv',
            'SCAM (MICX)' = 'micx',
            'SCAM (MDCX)' = 'mdxc',
            'SCAM (MDCV)' = 'mdcv',
            'Linear model with slope=0, resampled from P/B' = 'avg'
        ), selected='poly'),
        numericInput('polydegree', 'Degree of the polynomial', min=1, value=1, step=1),
        numericInput('knots', 'Number of knots for the adaptive GAM', value=1),
        actionButton('runPBf', 'Initialize the P/B data and model as a function of E', width='100%')
      ),
      mainPanel(
        column(5, align = "center",
          plotOutput("timeSeries", width='100%'), 
          plotOutput("PBE_dataplot", width='100%')
        ),
        column(5, align="center",
          plotOutput("kobe", width='100%'),
          plotOutput("normalE", width='100%')
        )
      )
    ),
    tabPanel(title = 'Step 2. Generate more complex models',
      sidebarPanel(),
      mainPanel()
    )
  )
)

server <- function(input, output, session){
  # Pre-run parameters if turbot is selected
  # Default disabled inputs
  shinyjs::disable('uploadedDataset')
  shinyjs::disable('knots')

  dataset <<- reactive({
    if(input$turbot==T || is.null(input$paramfile)){
      data(params)
      mapply(assign, names(params), params, MoreArgs=list(envir=globalenv()))
      out <- turbot
    } else {
      infile <- input$paramfile
      out <- read.csv(infile$datapath)
    }
    out
  })
  # Year span from the dataset
  output$input.refyears <- renderUI({
    sliderInput('refyear.slider', 
      label='Reference years for P/B', 
      min=min(dataset()$Year), 
      max=max(dataset()$Year), 
      step=1, 
      round=T,
      sep="", 
      value=c(min(dataset()$Year), max(dataset()$Year))
    )
  })

  # Reactive model inputs
  observeEvent(input$input.modtype, {
    if(input$modtype!='poly'){
      shinyjs::enable('polydegree')
    }
    if(input$modtype=='gam.adaptive'){
      shinyjs::enable('knots')
    }
  })

  # Run the initial model
  observeEvent(input$runPBf, {
    # paste(head(dataset()))
    PB <<- ccca::PB.f(dataset=dataset(), ref.years=seq(input$refyear.slider[1], input$refyear.slider[2]), q=input$q)
    PvsE <<- PBE.fit.f(PB, model.type=input$modtype, poly.degree=input$polydegree, knots=input$knots)
    kobePlot <- ccca::Kobe.f(PB=PB, E=PB$E)
      # ccca::colramp.legend(col1="red", col2="blue", ncol=length(PB$E), 2.5, 3.5, 2.7, 4.5)
    output$kobe <- renderPlot({
    ### PLOT
    })
  })
  
  # Graphical outputs
  ### TAB 2: 
}

shiny::shinyApp(ui = ui, server = server)