# Load libraries, assuming these have been installed
# First, gslea

library(ccca) # includes dependency: datatable, but nothing else
library(shinyjs)
library(shiny)
library(shinybusy)

ui <- fluidPage(
  shinyjs::useShinyjs(),
  titlePanel('Climate change conditioned advice (CCCA) for fisheries'),
  tabsetPanel(
#    widths=c(2,10),
    id='tabs',
    tabPanel(title = 'Step 1. Run the P/B~E model',
      sidebarPanel(
        p('First, `ccca` calculates the annual P/B of the population from a survey index time series. This is then modelled as a function of an `E` variable, either ecological or environmental.'),
        p('In the following steps, you can specify fishing strategies (using the mean exploitation rate from the dataset) and project climate projections for that fishing scenario.'),
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
        h4('Define additional parameters'),
        p("If the built-in turbot data are selected, the below parameters will be pre-specified and auto-filled. If you have uploaded a .csv file, specify the below parameters."),
        # uiOutput('input.step1params'),
        uiOutput('input.q'),
        uiOutput('input.Emean.shift.warm'),
        uiOutput('input.Emean.shift.cold'),
        uiOutput('input.E.var.inc'),
        numericInput('E.var.inc.extra', 'Additional variance for E (var) exploration', min=0, value=1.5),
        # uiOutput('input.Evar'),
        hr(),
        uiOutput('input.modtype'), 
        uiOutput('input.poly.degree'),
        uiOutput('input.knots'),        
        actionButton('runPBf', 'Initialize the P/B data and model as a function of E', width='100%', class="btn btn-primary")
      ),
      mainPanel(
        fluidRow(
          column(4, align="center", 
            plotOutput("timeSeries", width='100%')
          ), 
          column(4, align="center", 
            plotOutput("PBE_dataplot", width='100%')
          ), 
          column(4, align="center", 
            plotOutput("normalE", width='100%')
          )
        ), 
        fluidRow(
          column(10, align="center",  
            plotOutput("kobe", width='100%')
          ),
          column(2, align="center",
            plotOutput("kobe_legend", width='100%')
          )
        ),
        fluidRow(
          column(10,
            helpText("If these plots look accurate, proceed to the next step")
          ),
          column(2, align="right", 
            actionButton("next_tab1", "Next", width='100%', class="btn btn-secondary")
          )
        ) 
      )
    ),
    tabPanel(title = 'Step 2. Create E and PB projections',
      sidebarPanel(
        p('Next, project E values based on the choice of projection parameters, then use this to develop P/B values (based on the P/B vs E fit from Step 1).'),
        p('(If choosing the turbot example, those values are auto-filled below).'),
        hr(),
        uiOutput('input.proj.years'),
        uiOutput('input.N'),
        uiOutput('input.add.resids'),
        add_busy_spinner(spin = "fading-circle"),
        ### uiOutput('step2params'),
        actionButton('runEPBprojections', 'Calculate projection with these parameters', width='100%', class="btn btn-primary")
      ),
      mainPanel(
        p("The distribution of the P/B values for the future projections. It is the P/B that results by sampling E distribution function to simulate a future climate and while altering the mean and variance if desired. The P/B distribution is then deermined by running sampled E value through the fitted P/B vs E relationship."),
        fluidRow(
          column(1),
          column(10, align="center",  
            plotOutput("pbFuturePlot", width='100%')
          ),
          column(1)
        ),
        fluidRow(
          column(10,
            helpText("If this plot looks accurate, proceed to the next step")
          ),
          column(2, align="right", 
            actionButton("next_tab2", "Next", width='100%')
          )
        )
      )
    ),
    tabPanel(title = 'Step 3. Generate more complex models with exploitation',
      sidebarPanel(
        p('Next, project biomass over time given some fishing exploitation (using both density dependent and independent models).'),
        hr(),
        p('Define fishing exploitation parameters.'),
        p('(If choosing the turbot example, those values are auto-filled below).'),
        uiOutput('input.moratorium'),
        hr(),
        uiOutput('fs.length'), 
        uiOutput('input.fs'),
        hr(),
        uiOutput('input.fyear.slider'),
        uiOutput('input.ref.pt'),
        uiOutput('input.risk'),
        uiOutput('input.time.frame'),
        hr(),
        uiOutput('input.Bstart.mult'),
        uiOutput('input.K'),
        uiOutput('input.theta'),
        add_busy_spinner(spin = "fading-circle"),
        actionButton('runProjectionF', 'Calculate projection with these parameters', width='100%', class="btn btn-primary")
      ),
      mainPanel(
        p("Projections of biomass given the E scenario (for both density independent and dependent models). The biomass reference level is depicted as the horizontal dashed line. Confidence intervals are 90%."),
        plotOutput("fishFuturePlot", width='100%'),
        p("The probability that the biomass is greater than the reference level each year of the projection. The horizontal dashed line represents the risk tolerance of not achieving the objective. i.e. For the objective to be met, the biomass line should be above the risk line at the end of the time.frame period specified in the parameters list."),
        plotOutput("pBiomass10", width='100%'),
        p("The maximum exploitation rate that will achieve the objective in the specified time period given the E scenario projected."),
        plotOutput("maxExpRate", width='100%'),
        fluidRow(
          column(10,
            helpText("If this plot looks accurate, proceed to the next step")
          ),
          column(2, align="right", 
            actionButton("next_tab3", "Next", width='100%')
          )
        )
      )
    ),
    tabPanel(title = 'Step 4. Generate a contour plot of safe operating space',
      sidebarPanel(
        p("Calculate projections for each combination of E x F scenario (F values taken from Step 3, designate a vector of E shifts below). Then, generate a contour plot of safe operating space."),
        uiOutput('input.Ems.length'),
        uiOutput('input.Emean.shifts'),
        hr(),
        uiOutput('input.N.CCF'),
        add_busy_spinner(spin = "fading-circle"),
        actionButton('runCCF', 'Calculate projections', width='100%', class="btn btn-primary")
      ), 
      mainPanel(
        plotOutput("contour")
      )
    )
  )
)

server <- function(input, output, session){
  # Modified Kobe plotting function:
  Kobe.f.noL= function(PB, E, Bref.multiplier=1, col1="blue", col2="red", ...){

  	Year= PB$Year
  	F.rel= PB$F.rel
  	Index.q= PB$Index.q
  	base= PB$refererence.years==1
  	E.kobe= E/mean(E[base])
  	F.kobe= F.rel/mean(F.rel[base])
  	B.kobe= Index.q*Bref.multiplier/mean(Index.q[base])
  
    plot(B.kobe,F.kobe,type="b",pch="    ",xlab=expression("B/B"["base"]),ylab=expression("F/F"["base"]),...)
    E.categ= floor(E*4)/4 #quarter degree C categories
    tempcol=colorRampPalette(c(col1, col2))(length(E.kobe))
    temperaturecolours= tempcol[order(E.categ)]
    last.year= length(F.kobe)
    points(B.kobe,F.kobe,pch=21,bg=temperaturecolours,col=temperaturecolours,cex=1)
    points(B.kobe[1],F.kobe[1],pch=21,bg=temperaturecolours[1],col=temperaturecolours[1],cex=3)
    text(B.kobe[1],F.kobe[1],PB$Year[1],col="white",cex=0.55,font=2)
    points(B.kobe[last.year],F.kobe[last.year],pch=21,bg=temperaturecolours[last.year],col=temperaturecolours[last.year],cex=3)
    
    text(B.kobe[last.year],F.kobe[last.year],PB$Year[last.year],col="white",cex=0.55,font=2)
    abline(h=1,col="grey")
    abline(v=1,col="grey")
#    legend("topright",bty="n",cex=0.7,legend=c(paste0("Bbase=",round(B.kobe),
#                                       paste0("Fbase=",round(F.kobe,3)))))
  }
  colramp.legend.mod= function(col1="red", col2="blue", ncol, xleft, ybottom, xright, ytop, ...){
    tempcol=colorRampPalette(c(col1, col2))(ncol)
    legend_image <- as.raster(matrix(tempcol, ncol=1))
    rasterImage(legend_image, xleft=xleft, ybottom=ybottom, xright=xright, ytop=ytop)
    rasterImage(legend_image, 12,2,13,8)
    text(xright*1.5,ytop,labels=round(max(PB$E,na.rm=T),1),cex=1)
    text(xright*1.5,ybottom,labels=round(min(PB$E,na.rm=T),1),cex=1)
  }
  # Pre-run parameters if turbot is selected
  # Default disabled inputs
  shinyjs::disable('uploadedDataset')
  # shinyjs::disable('input.poly.degree')
  # shinyjs::disable('input.knots')
  # Disable "next" buttons unless that step's action button has been pushed
  shinyjs::disable("next_tab1")
  shinyjs::disable("next_tab2")
  shinyjs::disable("next_tab3")

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

  ### IF TURBOT SELECTED, AUTO-FILL UI ELEMENTS
  observeEvent(input$turbot, {
    # If not using the pre-loaded dataset:
    if(input$turbot == FALSE){

      # STEP 1 UI COMPONENTS
      
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

      output$input.q <- renderUI({
        numericInput('q', 'Catchability (q)', value=1, min=0)
      })
      output$input.Emean.shift.warm <- renderUI({
        numericInput('Emean.shift.warm', 'Mean upper shift of environmental variable (e.g. temperature is E, a value of 0.5 = 0.5 degrees WARMER)', value = 0.5, min=0)
      })
      output$input.Emean.shift.cold <- renderUI({
        numericInput('Emean.shift.cold', 'Mean lower shift environmental variable (e.g. temperature is E, a value of -0.5 = 0.5 degrees COLDER)', value = -0.5, max=0)
      })
      output$input.E.var.inc <- renderUI({
        numericInput('E.var.inc', 'The change in the variance of the E distribution for a future climate scenario (scalar, 1 for as is, >1 decreases variance)', min=0, step=1, value=1)
      })

      output$input.modtype <- renderUI({
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
        ), selected='poly')
      })
      output$input.poly.degree <- renderUI({
        numericInput('polydegree', 'Degree of the polynomial', min=1, value=1, step=1)
      })
      output$input.knots <- renderUI({
        numericInput('knots', 'Number of knots for the adaptive GAM', min=1, value=1, step=1)
      })
      # STEP 2 UI COMPONENTS
      output$input.proj.years <- renderUI({
        numericInput('proj.years', 'Number of years for projection into the future', min=1, value=10) #, post=' years')
      })
      output$input.N <- renderUI({
        numericInput('N', 'Number of different realizations of the future to create', value=2000, min=1)
      })
      output$input.add.resids <- renderUI({
        checkboxInput('add.resids', 'Add residuals to P/B values at each time step of projection', value=1)
      })

      ## STEP 3 UI ELEMENTS
      output$input.moratorium <- renderUI({
        checkboxInput('moratorium', 'Put under moratorium (i.e. exploitation rate = 0 going forward)', value=F)
      })

      output$fs.length <- renderUI({
        numericInput('input.fs.length', 'Number of F sequences to test (input different F values below)', min=1, step=1, value=1)
      })

      output$input.fs <- renderUI({
        fluidRow(
          p("Alternative F values to test"),
            column(12,
            lapply(1:max(1,input$input.fs.length), function(i){
              column(4,
                numericInput(paste0("f",i),paste0("F[",i,"]"),min=0,max=1,value=(i/input$input.fs.length))
              )
            })
          )
        )
      })

      output$input.fyear.slider <- renderUI({
        sliderInput('fyear.slider', 
          label='Reference years to calculate mean exploitation rate, F (the mean exploitation rate of this years will be used in projections)', 
          min=min(dataset()$Year), 
          max=max(dataset()$Year), 
          step=1, 
          round=T,
          sep="", 
          value=c(min(dataset()$Year), max(dataset()$Year))
        )
      })

      output$input.ref.pt <- renderUI({
        numericInput('ref.pt', 'Assuming your reference period represents the population at Bmsy, set a reference point as a proportion of B over these years (e.g. 1 indicates Bmsy, 0.4 indicates 40% of Bmsy; like Blim).', value=0.4, min=0)
      })

      output$input.risk <- renderUI({
        numericInput('risk', 'Acceptable risk of not achieving the reference point', value=0.5, min=0, max=1)
      })

      output$input.time.frame <- renderUI({
        numericInput('time.frame', 'Desired time (years) to reach reference point', value = 5, min=1)
      })      
      output$input.Bstart.mult <- renderUI({
        sliderInput('Bstart.mult', "Proportion of last year's biomass to initialize projection", value=1, min=0, max=1)
      })
      output$input.K <- renderUI({
        numericInput('K', "The multiplier of maximum observed biomass to carrying capacity", value=2)
      })
      output$input.theta <- renderUI({
        numericInput('theta', "Skew of the density dependence factor (1 = Schaeffer)", min=1e-6, value=1)
      })

      # STEP 4
      output$input.Ems.length <- renderUI({
        numericInput('Ems.length', 'Number of mean E shifts to test (input different F values below)', min=1, step=1, value=1)
      })

      output$input.Emean.shifts <- renderUI({
        fluidRow(
          p("Mean shifts in E (indicate cooling scenarios with negative numbers, warming with positive)"),
            column(12,
            lapply(1:max(1,input$Ems.length), function(i){
              column(4,
                numericInput(paste0("e",i),paste0("E[",i,"]"),value=0)
              )
            })
          )
        )
      })

      output$input.N.CCF <- renderUI({
        numericInput('N.CCF', 'Number of MC runs for determining climate change fishery scenarios', min=1, value=1000)
      })
    } else {

      # STEP 1 
      # If using the turbot dataset, auto-load these parameters
            # Year span from the dataset
      output$input.refyears <- renderUI({
        sliderInput('refyear.slider', 
          label='Reference years for P/B', 
          min=min(dataset()$Year), 
          max=max(dataset()$Year), 
          step=1, 
          round=T,
          sep="", 
          value=c(min(params$ref.years), max(params$ref.years))
        )
      })

      output$input.q <- renderUI({
        numericInput('q', 'Catchability (q)', value=as.numeric(params$q), min=0)
      })
      output$input.Emean.shift.warm <- renderUI({
        numericInput('Emean.shift.warm', 'Mean upper shift of environmental variable (e.g. temperature is E, a value of 0.5 = 0.5 degrees WARMER)', value = as.numeric(params$Emean.shift), min=0)
      })
      output$input.Emean.shift.cold <- renderUI({
        numericInput('Emean.shift.cold', 'Mean lower shift environmental variable (e.g. temperature is E, a value of -0.5 = 0.5 degrees COLDER)', value = -as.numeric(params$Emean.shift), max=0)
      })
      output$input.modtype <- renderUI({
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
        ), selected=params$model.type)
      })
      output$input.E.var.inc <- renderUI({
        numericInput('E.var.inc', 'The change in the variance of the E distribution for a future climate scenario (scalar, 1 for as is, <1 decreases variance)', min=0, step=1, value=params$E.var.inc)
      })
      output$input.poly.degree <- renderUI({
        numericInput('polydegree', 'Degree of the polynomial', min=1, value=params$poly.degree, step=1)
      })
      output$input.knots <- renderUI({
        numericInput('knots', 'Number of knots for the adaptive GAM', min=1, value=params$knots, step=1)
      })
      # STEP 2 UI COMPONENTS
      output$input.proj.years <- renderUI({
        numericInput('proj.years', 'Number of years for projection into the future', min=1, value=params$proj.years) #, post=' years')
      })
      output$input.N <- renderUI({
        numericInput('N', 'Number of different realizations of the future to create', value=params$N, min=1)
      })
      output$input.add.resids <- renderUI({
        checkboxInput('add.resids', 'Add residuals to P/B values at each time step of projection', value=params$add.resids)
      })

      ## STEP 3 UI ELEMENTS
      output$input.moratorium <- renderUI({
        checkboxInput('moratorium', 'Put under moratorium (i.e. exploitation rate = 0 going forward)', value=params$moratorium)
      })

      output$fs.length <- renderUI({
        numericInput('input.fs.length', 'Number of F sequences to test for climate-conditioned advice (inputs for these F appear below)', min=1, step=1, value=length(params$fs))
      })

      output$input.fs <- renderUI({
        fluidRow(
          p("Alternative F values to test"),
          column(12,
            lapply(1:input$input.fs.length, function(i){
              column(4,
                numericInput(paste0("f",i),paste0("F[",i,"]"),min=0,max=1,value=params$fs[i])
              )
            })
          )
        )
      })

      output$input.fyear.slider <- renderUI({
        sliderInput('fyear.slider', 
          label='Reference years to calculate mean exploitation rate, F (the mean exploitation rate of this years will be used in projections)', 
          min=min(dataset()$Year), 
          max=max(dataset()$Year), 
          step=1, 
          round=T,
          sep="", 
          value=c(2014,2018)
        )
      })

      output$input.ref.pt <- renderUI({
        numericInput('ref.pt', 'Assuming your reference period represents the population at Bmsy, set a reference point as a proportion of B over these years (e.g. 1 indicates Bmsy, 0.4 indicates 40% of Bmsy; like Blim).', min=0, value=params$ref.pt)
      })

      output$input.risk <- renderUI({
        numericInput('risk', 'Acceptable risk of not achieving the reference point', value=params$risk, min=0, max=1)
      })

      output$input.time.frame <- renderUI({
        numericInput('time.frame', 'Desired time (years) to reach reference point', value = params$time.frame, min=1)
      })
      output$input.Bstart.mult <- renderUI({
        sliderInput('Bstart.mult', "Proportion of last year's biomass to initialize projection", value=as.numeric(params$Bstart.mult), min=0, max=1)
      })
      output$input.K <- renderUI({
        numericInput('K', "The multiplier of maximum observed biomass to carrying capacity", value=as.numeric(params$K))
      })
      output$input.theta <- renderUI({
        numericInput('theta', "Skew of the density dependence factor (1 = Schaeffer)", min=1e-6, value=as.numeric(params$theta))
      })

      # STEP 4
      output$input.Ems.length <- renderUI({
        numericInput('Ems.length', 'Number of mean E shifts to test (input different F values below)', min=1, step=1, value=length(params$Emean.shifts))
      })

      output$input.Emean.shifts <- renderUI({
        fluidRow(
          p("Mean shifts in E (indicate cooling scenarios with negative numbers, warming with positive)"),
            column(12,
            lapply(1:max(1,input$Ems.length), function(i){
              column(4,
                numericInput(paste0("e",i),paste0("E[",i,"]"),value=params$Emean.shifts[i])
              )
            })
          )
        )
      })

      output$input.N.CCF <- renderUI({
        numericInput('N.CCF', 'Number of MC runs for determining climate change fishery scenarios', min=1, value=as.numeric(params$N.CCF))
      })

    }
  })

  # In step 2, print out the type of distribution being used: 
  output$Eproj.modtype <- renderText({
    paste0("FUTURE FEATURE!!! MODIFY THIS \n\n E is distributed according to a ", input$modtype, " dist")
  })
  
  ### Reactive ui - not related to default dataset
  # STEP 1
###  observeEvent(input$modtype, {
###    if(input$modtype=='poly'){
###      shinyjs::enable('input.poly.degree')
###    } else {
###      shinyjs::disable('input.poly.degree')
###    }
###    if(input$modtype=='gam.adaptive'){
###      shinyjs::enable('input.knots')
###    } else {
###      shinyjs::disable('input.knots')
###    }
###  })

  # REACTIONS/ACTIONS
  # STEP 1:
  # Run the initial model
  observeEvent(input$runPBf, {
    shinyjs::enable("next_tab1")
    # paste(head(dataset()))
    PB <<- ccca::PB.f(dataset=dataset(), ref.years=seq(input$refyear.slider[1], input$refyear.slider[2]), q=input$q)
      # ccca::colramp.legend(col1="red", col2="blue", ncol=length(PB$E), 2.5, 3.5, 2.7, 4.5)
    PvsE <<- ccca::PBE.fit.f(PB, model.type=input$modtype, poly.degree=input$polydegree, knots=input$knots)
    PvsE.null <<- ccca::PBE.fit.f(PB, model.type="avg", knots=input$knots, poly.degree=input$polydegree)
    output$timeSeries <- renderPlot({
      matplot(dataset()$Year, cbind(dataset()$Index, dataset()$Catch),type="l",lty=c(1,1),lwd=2,xlab="Year",ylab="Survey biomass and catch",col=c("black","green"))
      yaxis2.f(dataset()$Year, dataset()$E,ylabel=expression('Env. variable [often temperature ('^o*C*')]'),type="l",cex=1.1,lwd=2,lty=1,col="red")
      legend("topleft",bty="n",legend=c("Survey","Catch","Env. variable"),lwd=2,lty=c(1),cex=0.7,col=c("black","green","red"))
      title("Time series of survey, catch, and E")
    })

    output$PBE_dataplot <- renderPlot({
      na.year= nrow(PB)
      plot(na.omit(cbind(PB$E,PB$PB)),pch=20,xlab="E",ylab="P/B",col="darkgrey",type="n")
      text(PB$E[-na.year],PB$PB[-na.year],PB$Year[-na.year],cex=.7)
      pred.x= seq(min(PB$E)*.90,max(PB$E)*1.05,length=1000)
      lines(pred.x,predict(PvsE.null,newdata=data.frame(E=pred.x)),lwd=2,col="grey")
      lines(pred.x,predict(PvsE,newdata=data.frame(E=pred.x)),lwd=2)
      title("P/B vs E")
    })

    output$normalE <- renderPlot({
      Enorm <<- norm.fit.f(E=PB$E)
      Nrand <<- 1000000
      Edist.a <<- Enorm$estimate[1]
      Edist.b <<- Enorm$estimate[2]
      Ebase <<- density(norm.plot.f(Nrand=Nrand, Edist.a=Edist.a, Edist.b=Edist.b,Emean.shift=0,E.var.inc=input$E.var.inc))
      Ewarm <<- density(norm.plot.f(Nrand=Nrand, Edist.a=Edist.a, Edist.b=Edist.b,Emean.shift=input$Emean.shift.warm,E.var.inc=input$E.var.inc))
      Ecold <<- density(norm.plot.f(Nrand=Nrand, Edist.a=Edist.a, Edist.b=Edist.b,Emean.shift=input$Emean.shift.cold,E.var.inc=input$E.var.inc))
      Evar <<- density(norm.plot.f(Nrand=Nrand, Edist.a=Edist.a, Edist.b=Edist.b,Emean.shift=0,E.var.inc=input$E.var.inc.extra))
      plot(Ebase, xlab=expression('Env. variable'), ylab="Density",xlim=c(0,6),lwd=2,main="")
      lines(Ewarm, lwd=2,col="red")
      lines(Ecold, lwd=2,col="blue")
      lines(Evar, lwd=2,col="green")
      title("Normal distributions for E \n(with upper and lower shifts)")
      legend("topright",                    # Add legend to plot
       legend = c("E (base)", "E (lower)", "E (upper)", "E (var)"),
       col = c("black","red", "blue", "green"),
       pch = 16)
    })

    output$kobe <- renderPlot({
      Kobe.f.noL(PB=PB, E=PB$E)
      title("Kobe plot")
    })
    output$kobe_legend <- renderPlot({
      plot.new()
      colramp.legend.mod(PB, ncol=length(PB$E), xleft=0, ybottom=0, xright=0.3, ytop=1, col1="red", col2="blue")
    })
  })
  ### TAB 2: STEP 2 MORE COMPLEX MODELS
#  observeEvent(input$moratorium, {
#    if(input$moratorium){
#      shinyjs::disable('input.fishyears')
#    } else {
#      shinyjs::enable('input.fishyears')
#    }
#  })
  # uiInputs
  observeEvent(input$runEPBprojections, {
    shinyjs::enable("next_tab2")
    ###### Enorm <<- norm.fit.f(E=PB$E)
    ###### Nrand <<- 1000000
    ###### Edist.a <<- Enorm$estimate[1]
    ###### Edist.b <<- Enorm$estimate[2]
    Eproj <<- Eprojnorm.f(Edist.a=Edist.a, Edist.b=Edist.b, Emean.shift=input$Emean.shift.warm, proj.years=input$proj.years, input$N)
    Eproj.warm <<- Eprojnorm.f(Edist.a=Edist.a, Edist.b=Edist.b, Emean.shift=input$Emean.shift.warm, proj.years=input$proj.years, input$N)
    Eproj.cold <<- Eprojnorm.f(Edist.a=Edist.a, Edist.b=Edist.b, Emean.shift=input$Emean.shift.cold, proj.years=input$proj.years, input$N)

    PBproj <<- PB.for.projection.f(PvsE=PvsE, Eproj, add.residuals=input$add.resids)
    PBproj.null <<- PB.for.projection.f(PvsE=PvsE.null,Eproj, add.residuals=input$add.resids)
    PBproj.warm <<- PB.for.projection.f(PvsE=PvsE,Eproj.warm, add.residuals=input$add.resids)
    PBproj.cold <<- PB.for.projection.f(PvsE=PvsE,Eproj.cold, add.residuals=input$add.resids)

    Emean.shift.var <<- 0.0
    Edist.b.var <<- Edist.b*input$E.var.inc.extra
    Eproj.var <<- Eprojnorm.f(Edist.a=Edist.a, Edist.b=Edist.b, Emean.shift=Emean.shift.var, proj.years=input$proj.years, input$N)
    PBproj.var <<- PB.for.projection.f(PvsE=PvsE, Eproj.var, add.residuals=input$add.resids)

    ylim_max <- max(c(
      max(density(PBproj)$y),
      max(density(PBproj.null)$y),
      max(density(PBproj.cold)$y),
      max(density(PBproj.warm)$y),
      max(density(PBproj.var)$y)
    ))
    output$pbFuturePlot <- renderPlot({
      plot(density(PBproj.null),xlab="P/B",ylab="Density",lwd=2,main="",col="grey", ylim=c(0,ylim_max))
      lines(density(PBproj.cold),lwd=2,col="blue")
      lines(density(PBproj.warm),lwd=2,col="red")
      lines(density(PBproj.var),lwd=2,col="green")
      lines(density(PBproj),lwd=2,col="black")
    })
  })

  # STEP 3
  observeEvent(input$runProjectionF, {
    shinyjs::enable('next_tab3')
    # Convert series of F inputs into a vector of Fs to test:
    fs <<- vector()
    for(i in 1:input$input.fs.length){
      fs[i] <<- eval(parse(text=paste0("input$f",i)))
      # print(eval(parse(text=paste0("input$f",i))))
    }

    Fstrat <<- F.strategy(PB, input$fyear.slider[1]:input$fyear.slider[2], moratorium=input$moratorium)

    # Messy, not yet bug fixed
    Bproj <<- projection.f(PB=PB, Bstart.mult=input$Bstart.mult, PBproj=PBproj, Fstrat, K=input$K, theta=input$theta)
    Fout <<- Fseq.f(PB,PBproj=PBproj,Fseq=fs,time.frame=input$time.frame, N=input$N, K=input$K)
    PofF <<- PofF.f(PB,Fout,ref.pt=input$ref.pt)
    Bproj.summary <<- Bproj.summary.f(PB,Bproj,PBproj,Eproj)
    P <<- rankprob.f(Bproj,PB,input$ref.pt)

    Bproj.null <<- projection.f(PB=PB, Bstart.mult=input$Bstart.mult, PBproj=PBproj.null, Fstrat, K=input$K, theta=input$theta)
    Fout.null <<- Fseq.f(PB,PBproj=PBproj.null,Fseq=fs,time.frame=input$time.frame, N=input$N, K=input$K)
    PofF.null <<- PofF.f(PB,Fout.null,ref.pt=input$ref.pt)
    Bproj.summary.null <<- Bproj.summary.f(PB,Bproj.null,PBproj.null,Eproj)
    P.null <<- rankprob.f(Bproj.null,PB,input$ref.pt)

    Bproj.warm <<- projection.f(PB=PB, Bstart.mult=input$Bstart.mult, PBproj=PBproj.warm, Fstrat, K=input$K, theta=input$theta)
    Fout.warm <<- Fseq.f(PB,PBproj=PBproj.warm,Fseq=fs,time.frame=input$time.frame, N=input$N, K=input$K)
    PofF.warm <<- PofF.f(PB,Fout.warm,ref.pt=input$ref.pt)
    Bproj.summary.warm <<- Bproj.summary.f(PB,Bproj.warm,PBproj.warm,Eproj.warm)
    P.warm <<- rankprob.f(Bproj.warm,PB,input$ref.pt)

    Bproj.cold <<- projection.f(PB=PB, Bstart.mult=input$Bstart.mult, PBproj=PBproj.cold, Fstrat, K=input$K, theta=input$theta)
    Fout.cold <<- Fseq.f(PB, PBproj=PBproj.cold, Fseq=fs, time.frame=input$time.frame, N=input$N, K=input$K)
    PofF.cold <<- PofF.f(PB, Fout.cold, ref.pt=input$ref.pt)
    Bproj.summary.cold <<- Bproj.summary.f(PB, Bproj.cold, PBproj.cold, Eproj.cold)
    P.cold <<- rankprob.f(Bproj.cold, PB, input$ref.pt)

    Bproj.var <<- projection.f(PB=PB, Bstart.mult=input$Bstart.mult, PBproj=PBproj.var, Fstrat, K=input$K, theta=input$theta)
    Fout.var <<- Fseq.f(PB,Fseq=fs, PBproj=PBproj.var, time.frame=input$time.frame, N=input$N, K=input$K)
    PofF.var <<- PofF.f(PB, Fout.var, ref.pt=input$ref.pt)
    Bproj.summary.var <<- Bproj.summary.f(PB, Bproj.var, PBproj.var, Eproj.var)
    P.var <<- rankprob.f(Bproj.var, PB, input$ref.pt)

    output$fishFuturePlot <- renderPlot({
      par(mfcol=c(5,1),mar=c(1,4,2,2),omi=c(.6,2.2,.1,2.2))
      Bref= ref.pt*sum(PB$Index.q * PB$refererence.years)/sum(PB$refererence.years)
      plot(Bproj.summary.null$year,Bproj.summary.null$B.di.CI.med,type="n",ylim=c(0,max(Bproj.summary.null$B.di.CI.high)),xlab="",ylab="")
      confint(Bproj.summary.null$year,Bproj.summary.null$B.di.CI.low,Bproj.summary.null$B.di.CI.high,col="grey")
      confint(Bproj.summary.null$year,Bproj.summary.null$B.dd.CI.low,Bproj.summary.null$B.dd.CI.high,col="lightblue")
      lines(Bproj.summary.null$year,Bproj.summary.null$B.dd.CI.med,lwd=2,col="blue")
      lines(Bproj.summary.null$year,Bproj.summary.null$B.dd.CI.low,lwd=1,lty=2,col="blue")
      lines(Bproj.summary.null$year,Bproj.summary.null$B.dd.CI.high,lwd=1,lty=2,col="blue")
      lines(Bproj.summary.null$year,Bproj.summary.null$B.di.CI.med,lwd=2)
      lines(Bproj.summary.null$year,Bproj.summary.null$B.di.CI.low,lwd=1,lty=2,col="black")
      lines(Bproj.summary.null$year,Bproj.summary.null$B.di.CI.high,lwd=1,lty=2,col="black")
      lines(Bproj.summary.null$year,rep(Bref,length(Bproj.summary.null$year)),lty=2)
      legend("topright",legend=c("Density independent","Density dependent"),lwd=2,col=c("black","blue"),bty="n",cex=0.75)
      legend("topleft",legend="Null",bty="n",cex=0.75)

      Bref= ref.pt*sum(PB$Index.q * PB$refererence.years)/sum(PB$refererence.years)
      plot(Bproj.summary$year,Bproj.summary$B.di.CI.med,type="n",ylim=c(0,max(Bproj.summary$B.di.CI.high)),xlab="",ylab="")
      confint(Bproj.summary$year,Bproj.summary$B.di.CI.low,Bproj.summary$B.di.CI.high,col="grey")
      confint(Bproj.summary$year,Bproj.summary$B.dd.CI.low,Bproj.summary$B.dd.CI.high,col="lightblue")
      lines(Bproj.summary$year,Bproj.summary$B.dd.CI.med,lwd=2,col="blue")
      lines(Bproj.summary$year,Bproj.summary$B.dd.CI.low,lwd=1,lty=2,col="blue")
      lines(Bproj.summary$year,Bproj.summary$B.dd.CI.high,lwd=1,lty=2,col="blue")
      lines(Bproj.summary$year,Bproj.summary$B.di.CI.med,lwd=2)
      lines(Bproj.summary$year,Bproj.summary$B.di.CI.low,lwd=1,lty=2,col="black")
      lines(Bproj.summary$year,Bproj.summary$B.di.CI.high,lwd=1,lty=2,col="black")
      lines(Bproj.summary$year,rep(Bref,length(Bproj.summary$year)),lty=2)
      legend("topleft",legend="Mean temperature",bty="n",cex=0.75)
      yaxis2.f(Bproj.summary$year,Bproj.summary$ E.CI.med,ylabel="",type="l",cex=1,,lwd=2,lty=1,col="red")

      Bref= ref.pt*sum(PB$Index.q * PB$refererence.years)/sum(PB$refererence.years)
      plot(Bproj.summary.warm$year,Bproj.summary.warm$B.di.CI.med,type="n",ylim=c(0,max(Bproj.summary.warm$B.di.CI.high)),xlab="",ylab="")
      confint(Bproj.summary.warm$year,Bproj.summary.warm$B.di.CI.low,Bproj.summary.warm$B.di.CI.high,col="grey")
      confint(Bproj.summary.warm$year,Bproj.summary.warm$B.dd.CI.low,Bproj.summary.warm$B.dd.CI.high,col="lightblue")
      lines(Bproj.summary.warm$year,Bproj.summary.warm$B.dd.CI.med,lwd=2,col="blue")
      lines(Bproj.summary.warm$year,Bproj.summary.warm$B.dd.CI.low,lwd=1,lty=2,col="blue")
      lines(Bproj.summary.warm$year,Bproj.summary.warm$B.dd.CI.high,lwd=1,lty=2,col="blue")
      lines(Bproj.summary.warm$year,Bproj.summary.warm$B.di.CI.med,lwd=2)
      lines(Bproj.summary.warm$year,Bproj.summary.warm$B.di.CI.low,lwd=1,lty=2,col="black")
      lines(Bproj.summary.warm$year,Bproj.summary.warm$B.di.CI.high,lwd=1,lty=2,col="black")
      lines(Bproj.summary.warm$year,rep(Bref,length(Bproj.summary.warm$year)),lty=2)
      legend("topleft",legend=paste0(input$Emean.shift.warm, " °C warmer"),bty="n",cex=0.75)
      yaxis2.f(Bproj.summary.warm$year,Bproj.summary.warm$ E.CI.med,ylabel="",type="l",cex=1,,lwd=2,lty=1,col="red")

      Bref= ref.pt*sum(PB$Index.q * PB$refererence.years)/sum(PB$refererence.years)
      plot(Bproj.summary.cold$year,Bproj.summary.cold$B.di.CI.med,type="n",ylim=c(0,max(Bproj.summary.cold$B.di.CI.high)),xlab="",ylab="")
      confint(Bproj.summary.cold$year,Bproj.summary.cold$B.di.CI.low,Bproj.summary.cold$B.di.CI.high,col="grey")
      confint(Bproj.summary.cold$year,Bproj.summary.cold$B.dd.CI.low,Bproj.summary.cold$B.dd.CI.high,col="lightblue")
      lines(Bproj.summary.cold$year,Bproj.summary.cold$B.dd.CI.med,lwd=2,col="blue")
      lines(Bproj.summary.cold$year,Bproj.summary.cold$B.dd.CI.low,lwd=1,lty=2,col="blue")
      lines(Bproj.summary.cold$year,Bproj.summary.cold$B.dd.CI.high,lwd=1,lty=2,col="blue")
      lines(Bproj.summary.cold$year,Bproj.summary.cold$B.di.CI.med,lwd=2)
      lines(Bproj.summary.cold$year,Bproj.summary.cold$B.di.CI.low,lwd=1,lty=2,col="black")
      lines(Bproj.summary.cold$year,Bproj.summary.cold$B.di.CI.high,lwd=1,lty=2,col="black")
      lines(Bproj.summary.cold$year,rep(Bref,length(Bproj.summary.cold$year)),lty=2)
      legend("topleft",legend=paste0(input$Emean.shift.cold, " °C colder"),bty="n",cex=0.75)
      yaxis2.f(Bproj.summary.cold$year,Bproj.summary.cold$ E.CI.med,ylabel="",type="l",cex=1,,lwd=2,lty=1,col="red")

      Bref= ref.pt*sum(PB$Index.q * PB$refererence.years)/sum(PB$refererence.years)
      plot(Bproj.summary.var$year,Bproj.summary.var$B.di.CI.med,type="n",ylim=c(0,max(Bproj.summary.var$B.di.CI.high)),xlab="",ylab="")
      confint(Bproj.summary.var$year,Bproj.summary.var$B.di.CI.low,Bproj.summary.var$B.di.CI.high,col="grey")
      confint(Bproj.summary.var$year,Bproj.summary.var$B.dd.CI.low,Bproj.summary.var$B.dd.CI.high,col="lightblue")
      lines(Bproj.summary.var$year,Bproj.summary.var$B.dd.CI.med,lwd=2,col="blue")
      lines(Bproj.summary.var$year,Bproj.summary.var$B.dd.CI.low,lwd=1,lty=2,col="blue")
      lines(Bproj.summary.var$year,Bproj.summary.var$B.dd.CI.high,lwd=1,lty=2,col="blue")
      lines(Bproj.summary.var$year,Bproj.summary.var$B.di.CI.med,lwd=2)
      lines(Bproj.summary.var$year,Bproj.summary.var$B.di.CI.low,lwd=1,lty=2,col="black")
      lines(Bproj.summary.var$year,Bproj.summary.var$B.di.CI.high,lwd=1,lty=2,col="black")
      lines(Bproj.summary.var$year,rep(Bref,length(Bproj.summary.var$year)),lty=2)
      legend("topleft",legend=paste0("sd x ", input$E.var.inc),bty="n",cex=0.75)
      yaxis2.f(Bproj.summary.var$year,Bproj.summary.var$ E.CI.med,ylabel="",type="l",cex=1,,lwd=2,lty=1,col="red")

      mtext(side=1,outer=F,text="Year",line=4)
      mtext(outer=T,side=2,text="Biomass (t)",line=-1)
      mtext(outer=T,side=4,text="Temperature (°C)",line=1)
    })

    output$pBiomass10 <- renderPlot({
      par(mfcol=c(5,1),mar=c(1,4,2,2),omi=c(.6,2.2,.1,2.2))
      matplot(P.null[,1],P.null[,-1],type='l',xlab="",ylab="",lwd=3,ylim=c(0,1),lty=1,col=c("black","blue"))
      legend("topright",legend=c("Density independent","Density dependent"),lwd=2,col=c("black","blue"),bty="n",cex=0.75)
      legend("topleft",legend="Null",bty="n",cex=0.75)
      abline(h=1-input$risk,lty=2,col="grey")
      box()

      matplot(P[,1],P[,-1],type='l', xlab="", ylab="", lwd=3, ylim=c(0,1), lty=1, col=c("black","blue"))
      legend("topleft",legend="Mean temperature",bty="n",cex=0.75)
      abline(h=1-input$risk,lty=2,col="grey")
      box()

      matplot(P.warm[,1],P.warm[,-1],type='l',xlab="",ylab="",lwd=3,ylim=c(0,1),lty=1,col=c("black","blue"))
      legend("topleft",legend=paste0(input$Emean.shift.warm, " °C warmer"), bty="n",cex=0.75)
      abline(h=1-input$risk,lty=2,col="grey")
        box()

      matplot(P.cold[,1],P.cold[,-1],type='l',xlab="",ylab="",lwd=3,ylim=c(0,1),lty=1,col=c("black","blue"))
      legend("topleft",legend=paste0(input$Emean.shift.cold, " °C colder"), bty="n",cex=0.75)
      abline(h=1-input$risk,lty=2,col="grey")
      box()

      matplot(P.var[,1],P.var[,-1],type='l',xlab="",ylab="",lwd=3,ylim=c(0,1),lty=1,col=c("black","blue"))
      legend("topleft",legend=paste0("sd x ", input$E.var.inc), bty="n", cex=0.75)
      abline(h=1-input$risk,lty=2,col="grey")
      box()

      mtext(side=1,outer=F,text="Year",line=4)
      mtext(outer=T,side=2,text="Probability of being at or above biomass objective",line=-1)
    })

    output$maxExpRate <- renderPlot({

      par(mfcol=c(5,1),mar=c(1,4,2,2),omi=c(.6,2.2,.1,2.2))
      matplot(PofF.null[,1],PofF.null[,-1],xlab="", ylab="" ,ylim=c(0,1),
      type="l",lwd=3,xaxs="i",yaxs="i",lty=1,col=c("black","blue"))
      legend("topright",legend=c("Density independent","Density dependent"),lwd=2,col=c("black","blue"),bty="n",cex=0.75)
      legend("topleft",legend="Null",bty="n",cex=0.75)
      # print(PofF.null$P.di)
      # print(PofF.null$P.dd)
      di.intersection <- switch(input$modtype,
        poly = predict(lm(f~poly(P.di, degree=input$polydegree),data=PofF.null),newdata=data.frame(P.di=1-input$risk)),
        gam = predict(gam(f~s(P.di),data=PofF.null),newdata=data.frame(P.di=1-input$risk)),
        gam.adaptive=predict( gam(f~s(P.di, k=input$knots,bs="ad"), data=PofF.null),newdata=data.frame(P.di=1-input$risk)),
        mpi= predict(scam(f~s(P.di, bs="mpi"),data=PofF.null),newdata=data.frame(P.di=1-input$risk)),
        mpd= predict(scam(f~s(P.di, bs="mpd"),data=PofF.null),newdata=data.frame(P.di=1-input$risk)),
        cx= predict(scam(f~s(P.di, bs="cx"),data=PofF.null),newdata=data.frame(P.di=1-input$risk)),
        cv= predict(scam(f~s(P.di, bs="cv"),data=PofF.null),newdata=data.frame(P.di=1-input$risk)),
        micx= predict(scam(f~s(P.di, bs="micx"),data=PofF.null),newdata=data.frame(P.di=1-input$risk)),
        micv= predict(scam(f~s(P.di, bs="micv"),data=PofF.null),newdata=data.frame(P.di=1-input$risk)),
        mdcx= predict(scam(f~s(P.di, bs="mdcx"),data=PofF.null),newdata=data.frame(P.di=1-input$risk)),
        mdcv= predict(scam(f~s(P.di, bs="mdcv"),data=PofF.null),newdata=data.frame(P.di=1-input$risk)),
        avg= predict(lm(f - 0*P.di ~ 1, data=PofF.null), newdata=data.frame(P.di=1-input$risk))
      )
      # print("dd:")
      # print(PofF.null$P.dd)
      dd.intersection <- switch(input$modtype,
        poly = predict(lm(f~poly(P.dd, degree=input$polydegree),data=PofF.null),newdata=data.frame(P.dd=1-input$risk)),
        gam = predict(gam(f~s(P.dd),data=PofF.null),newdata=data.frame(P.dd=1-input$risk)),
        gam.adaptive=predict( gam(f~s(P.dd, k=input$knots,bs="ad"), data=PofF.null),newdata=data.frame(P.dd=1-input$risk)),
        mpi= predict(scam(f~s(P.dd, bs="mpi"),data=PofF.null),newdata=data.frame(P.dd=1-input$risk)),
        mpd= predict(scam(f~s(P.dd, bs="mpd"),data=PofF.null),newdata=data.frame(P.dd=1-input$risk)),
        cx= predict(scam(f~s(P.dd, bs="cx"),data=PofF.null),newdata=data.frame(P.dd=1-input$risk)),
        cv= predict(scam(f~s(P.dd, bs="cv"),data=PofF.null),newdata=data.frame(P.dd=1-input$risk)),
        micx= predict(scam(f~s(P.dd, bs="micx"),data=PofF.null),newdata=data.frame(P.dd=1-input$risk)),
        micv= predict(scam(f~s(P.dd, bs="micv"),data=PofF.null),newdata=data.frame(P.dd=1-input$risk)),
        mdcx= predict(scam(f~s(P.dd, bs="mdcx"),data=PofF.null),newdata=data.frame(P.dd=1-input$risk)),
        mdcv= predict(scam(f~s(P.dd, bs="mdcv"),data=PofF.null),newdata=data.frame(P.dd=1-input$risk)),
        avg= predict(lm(f - 0*P.dd ~ 1, data=PofF.null), newdata=data.frame(P.dd=1-input$risk))
      )

#      dd.intersection <- predict(gam(f~s(P.dd),data=PofF.null),newdata=data.frame(P.dd=1-input$risk))
      rect(0,0,di.intersection,1-input$risk,lty=2,border="darkgrey")
      rect(0,0,dd.intersection,1-input$risk,lty=2,border="darkgrey")
      box()

      matplot(PofF[,1],PofF[,-1],xlab="", ylab="" ,ylim=c(0,1),
      type="l",lwd=3,xaxs="i",yaxs="i",lty=1,col=c("black","blue"))
      legend("topleft",legend="Mean temperature",bty="n",cex=0.75)

      ### di.intersection= predict(gam(f~s(P.di),data=PofF),newdata=data.frame(P.di=1-input$risk))
      di.intersection <- switch(input$modtype,
        poly = predict(lm(f~poly(P.di, degree=input$polydegree),data=PofF),newdata=data.frame(P.di=1-input$risk)),
        gam = predict(gam(f~s(P.di),data=PofF),newdata=data.frame(P.di=1-input$risk)),
        gam.adaptive=predict( gam(f~s(P.di, k=input$knots,bs="ad"), data=PofF),newdata=data.frame(P.di=1-input$risk)),
        mpi= predict(scam(f~s(P.di, bs="mpi"),data=PofF),newdata=data.frame(P.di=1-input$risk)),
        mpd= predict(scam(f~s(P.di, bs="mpd"),data=PofF),newdata=data.frame(P.di=1-input$risk)),
        cx= predict(scam(f~s(P.di, bs="cx"),data=PofF),newdata=data.frame(P.di=1-input$risk)),
        cv= predict(scam(f~s(P.di, bs="cv"),data=PofF),newdata=data.frame(P.di=1-input$risk)),
        micx= predict(scam(f~s(P.di, bs="micx"),data=PofF),newdata=data.frame(P.di=1-input$risk)),
        micv= predict(scam(f~s(P.di, bs="micv"),data=PofF),newdata=data.frame(P.di=1-input$risk)),
        mdcx= predict(scam(f~s(P.di, bs="mdcx"),data=PofF),newdata=data.frame(P.di=1-input$risk)),
        mdcv= predict(scam(f~s(P.di, bs="mdcv"),data=PofF),newdata=data.frame(P.di=1-input$risk)),
        avg= predict(lm(f - 0*P.di ~ 1, data=PofF), newdata=data.frame(P.di=1-input$risk))
      )

      ### dd.intersection= predict(gam(f~s(P.dd),data=PofF),newdata=data.frame(P.dd=1-input$risk))
      dd.intersection <- switch(input$modtype,
        poly = predict(lm(f~poly(P.dd, degree=input$polydegree),data=PofF),newdata=data.frame(P.dd=1-input$risk)),
        gam = predict(gam(f~s(P.dd),data=PofF),newdata=data.frame(P.dd=1-input$risk)),
        gam.adaptive=predict( gam(f~s(P.dd, k=input$knots,bs="ad"), data=PofF),newdata=data.frame(P.dd=1-input$risk)),
        mpi= predict(scam(f~s(P.dd, bs="mpi"),data=PofF),newdata=data.frame(P.dd=1-input$risk)),
        mpd= predict(scam(f~s(P.dd, bs="mpd"),data=PofF),newdata=data.frame(P.dd=1-input$risk)),
        cx= predict(scam(f~s(P.dd, bs="cx"),data=PofF),newdata=data.frame(P.dd=1-input$risk)),
        cv= predict(scam(f~s(P.dd, bs="cv"),data=PofF),newdata=data.frame(P.dd=1-input$risk)),
        micx= predict(scam(f~s(P.dd, bs="micx"),data=PofF),newdata=data.frame(P.dd=1-input$risk)),
        micv= predict(scam(f~s(P.dd, bs="micv"),data=PofF),newdata=data.frame(P.dd=1-input$risk)),
        mdcx= predict(scam(f~s(P.dd, bs="mdcx"),data=PofF),newdata=data.frame(P.dd=1-input$risk)),
        mdcv= predict(scam(f~s(P.dd, bs="mdcv"),data=PofF),newdata=data.frame(P.dd=1-input$risk)),
        avg= predict(lm(f - 0*P.dd ~ 1, data=PofF), newdata=data.frame(P.dd=1-input$risk))
      )
      rect(0,0,di.intersection,1-input$risk,lty=2,border="darkgrey")
      rect(0,0,dd.intersection,1-input$risk,lty=2,border="darkgrey")
      box()

      # warm
      matplot(PofF.warm[,1],PofF.warm[,-1],xlab="", ylab= "",ylim=c(0,1),
      type="l",lwd=3,xaxs="i",yaxs="i",lty=1,col=c("black","blue"))
      legend("topleft",legend="0.5 °C warmer",bty="n",cex=0.75)
      ### di.intersection= predict(gam(f~s(P.di),data=PofF.warm),newdata=data.frame(P.di=1-input$risk))
      di.intersection <- switch(input$modtype,
        poly = predict(lm(f~poly(P.di, degree=input$polydegree),data=PofF.warm),newdata=data.frame(P.di=1-input$risk)),
        gam = predict(gam(f~s(P.di),data=PofF.warm),newdata=data.frame(P.di=1-input$risk)),
        gam.adaptive=predict( gam(f~s(P.di, k=input$knots,bs="ad"), data=PofF.warm),newdata=data.frame(P.di=1-input$risk)),
        mpi= predict(scam(f~s(P.di, bs="mpi"),data=PofF.warm),newdata=data.frame(P.di=1-input$risk)),
        mpd= predict(scam(f~s(P.di, bs="mpd"),data=PofF.warm),newdata=data.frame(P.di=1-input$risk)),
        cx= predict(scam(f~s(P.di, bs="cx"),data=PofF.warm),newdata=data.frame(P.di=1-input$risk)),
        cv= predict(scam(f~s(P.di, bs="cv"),data=PofF.warm),newdata=data.frame(P.di=1-input$risk)),
        micx= predict(scam(f~s(P.di, bs="micx"),data=PofF.warm),newdata=data.frame(P.di=1-input$risk)),
        micv= predict(scam(f~s(P.di, bs="micv"),data=PofF.warm),newdata=data.frame(P.di=1-input$risk)),
        mdcx= predict(scam(f~s(P.di, bs="mdcx"),data=PofF.warm),newdata=data.frame(P.di=1-input$risk)),
        mdcv= predict(scam(f~s(P.di, bs="mdcv"),data=PofF.warm),newdata=data.frame(P.di=1-input$risk)),
        avg= predict(lm(f - 0*P.di ~ 1, data=PofF.warm), newdata=data.frame(P.di=1-input$risk))
      )

      dd.intersection <- switch(input$modtype,
        poly = predict(lm(f~poly(P.dd, degree=input$polydegree),data=PofF.warm),newdata=data.frame(P.dd=1-input$risk)),
        gam = predict(gam(f~s(P.dd),data=PofF.warm),newdata=data.frame(P.dd=1-input$risk)),
        gam.adaptive=predict( gam(f~s(P.dd, k=input$knots,bs="ad"), data=PofF.warm),newdata=data.frame(P.dd=1-input$risk)),
        mpi= predict(scam(f~s(P.dd, bs="mpi"),data=PofF.warm),newdata=data.frame(P.dd=1-input$risk)),
        mpd= predict(scam(f~s(P.dd, bs="mpd"),data=PofF.warm),newdata=data.frame(P.dd=1-input$risk)),
        cx= predict(scam(f~s(P.dd, bs="cx"),data=PofF.warm),newdata=data.frame(P.dd=1-input$risk)),
        cv= predict(scam(f~s(P.dd, bs="cv"),data=PofF.warm),newdata=data.frame(P.dd=1-input$risk)),
        micx= predict(scam(f~s(P.dd, bs="micx"),data=PofF.warm),newdata=data.frame(P.dd=1-input$risk)),
        micv= predict(scam(f~s(P.dd, bs="micv"),data=PofF.warm),newdata=data.frame(P.dd=1-input$risk)),
        mdcx= predict(scam(f~s(P.dd, bs="mdcx"),data=PofF.warm),newdata=data.frame(P.dd=1-input$risk)),
        mdcv= predict(scam(f~s(P.dd, bs="mdcv"),data=PofF.warm),newdata=data.frame(P.dd=1-input$risk)),
        avg= predict(lm(f - 0*P.dd ~ 1, data=PofF.warm), newdata=data.frame(P.dd=1-input$risk))
      )
      ### dd.intersection= predict(gam(f~s(P.dd),data=PofF.warm),newdata=data.frame(P.dd=1-input$risk))
      rect(0,0,di.intersection,1-input$risk,lty=2,border="darkgrey")
      rect(0,0,dd.intersection,1-input$risk,lty=2,border="darkgrey")
      box()

      # cold
      matplot(PofF.cold[,1],PofF.cold[,-1],xlab="", ylab= "",ylim=c(0,1),
      type="l",lwd=3,xaxs="i",yaxs="i",lty=1,col=c("black","blue"))
      legend("topleft",legend="0.5 °C colder",bty="n",cex=0.75)
      ### di.intersection= predict(gam(f~s(P.di),data=PofF.cold),newdata=data.frame(P.di=1-input$risk))
      ### dd.intersection= predict(gam(f~s(P.dd),data=PofF.cold),newdata=data.frame(P.dd=1-input$risk))
      di.intersection <- switch(input$modtype,
        poly = predict(lm(f~poly(P.di, degree=input$polydegree),data=PofF.cold),newdata=data.frame(P.di=1-input$risk)),
        gam = predict(gam(f~s(P.di),data=PofF.cold),newdata=data.frame(P.di=1-input$risk)),
        gam.adaptive=predict( gam(f~s(P.di, k=input$knots,bs="ad"), data=PofF.cold),newdata=data.frame(P.di=1-input$risk)),
        mpi= predict(scam(f~s(P.di, bs="mpi"),data=PofF.cold),newdata=data.frame(P.di=1-input$risk)),
        mpd= predict(scam(f~s(P.di, bs="mpd"),data=PofF.cold),newdata=data.frame(P.di=1-input$risk)),
        cx= predict(scam(f~s(P.di, bs="cx"),data=PofF.cold),newdata=data.frame(P.di=1-input$risk)),
        cv= predict(scam(f~s(P.di, bs="cv"),data=PofF.cold),newdata=data.frame(P.di=1-input$risk)),
        micx= predict(scam(f~s(P.di, bs="micx"),data=PofF.cold),newdata=data.frame(P.di=1-input$risk)),
        micv= predict(scam(f~s(P.di, bs="micv"),data=PofF.cold),newdata=data.frame(P.di=1-input$risk)),
        mdcx= predict(scam(f~s(P.di, bs="mdcx"),data=PofF.cold),newdata=data.frame(P.di=1-input$risk)),
        mdcv= predict(scam(f~s(P.di, bs="mdcv"),data=PofF.cold),newdata=data.frame(P.di=1-input$risk)),
        avg= predict(lm(f - 0*P.di ~ 1, data=PofF.cold), newdata=data.frame(P.di=1-input$risk))
      )

      dd.intersection <- switch(input$modtype,
        poly = predict(lm(f~poly(P.dd, degree=input$polydegree),data=PofF.cold),newdata=data.frame(P.dd=1-input$risk)),
        gam = predict(gam(f~s(P.dd),data=PofF.cold),newdata=data.frame(P.dd=1-input$risk)),
        gam.adaptive=predict( gam(f~s(P.dd, k=input$knots,bs="ad"), data=PofF.cold),newdata=data.frame(P.dd=1-input$risk)),
        mpi= predict(scam(f~s(P.dd, bs="mpi"),data=PofF.cold),newdata=data.frame(P.dd=1-input$risk)),
        mpd= predict(scam(f~s(P.dd, bs="mpd"),data=PofF.cold),newdata=data.frame(P.dd=1-input$risk)),
        cx= predict(scam(f~s(P.dd, bs="cx"),data=PofF.cold),newdata=data.frame(P.dd=1-input$risk)),
        cv= predict(scam(f~s(P.dd, bs="cv"),data=PofF.cold),newdata=data.frame(P.dd=1-input$risk)),
        micx= predict(scam(f~s(P.dd, bs="micx"),data=PofF.cold),newdata=data.frame(P.dd=1-input$risk)),
        micv= predict(scam(f~s(P.dd, bs="micv"),data=PofF.cold),newdata=data.frame(P.dd=1-input$risk)),
        mdcx= predict(scam(f~s(P.dd, bs="mdcx"),data=PofF.cold),newdata=data.frame(P.dd=1-input$risk)),
        mdcv= predict(scam(f~s(P.dd, bs="mdcv"),data=PofF.cold),newdata=data.frame(P.dd=1-input$risk)),
        avg= predict(lm(f - 0*P.dd ~ 1, data=PofF.cold), newdata=data.frame(P.dd=1-input$risk))
      )

      rect(0,0,di.intersection,1-input$risk,lty=2,border="darkgrey")
      rect(0,0,dd.intersection,1-input$risk,lty=2,border="darkgrey")
      box()

      # increased variance
      matplot(PofF.var[,1],PofF.var[,-1],xlab="", ylab= "",ylim=c(0,1),
      type="l",lwd=3,xaxs="i",yaxs="i",lty=1,col=c("black","blue"))
      legend("topleft",legend="sd x 1.5",bty="n",cex=0.75)
      ### di.intersection= predict(gam(f~s(P.di),data=PofF.var),newdata=data.frame(P.di=1-input$risk))
      ### dd.intersection= predict(gam(f~s(P.dd),data=PofF.var),newdata=data.frame(P.dd=1-input$risk))
      di.intersection <- switch(input$modtype,
        poly = predict(lm(f~poly(P.di, degree=input$polydegree),data=PofF.var),newdata=data.frame(P.di=1-input$risk)),
        gam = predict(gam(f~s(P.di),data=PofF.var),newdata=data.frame(P.di=1-input$risk)),
        gam.adaptive=predict( gam(f~s(P.di, k=input$knots,bs="ad"), data=PofF.var),newdata=data.frame(P.di=1-input$risk)),
        mpi= predict(scam(f~s(P.di, bs="mpi"),data=PofF.var),newdata=data.frame(P.di=1-input$risk)),
        mpd= predict(scam(f~s(P.di, bs="mpd"),data=PofF.var),newdata=data.frame(P.di=1-input$risk)),
        cx= predict(scam(f~s(P.di, bs="cx"),data=PofF.var),newdata=data.frame(P.di=1-input$risk)),
        cv= predict(scam(f~s(P.di, bs="cv"),data=PofF.var),newdata=data.frame(P.di=1-input$risk)),
        micx= predict(scam(f~s(P.di, bs="micx"),data=PofF.var),newdata=data.frame(P.di=1-input$risk)),
        micv= predict(scam(f~s(P.di, bs="micv"),data=PofF.var),newdata=data.frame(P.di=1-input$risk)),
        mdcx= predict(scam(f~s(P.di, bs="mdcx"),data=PofF.var),newdata=data.frame(P.di=1-input$risk)),
        mdcv= predict(scam(f~s(P.di, bs="mdcv"),data=PofF.var),newdata=data.frame(P.di=1-input$risk)),
        avg= predict(lm(f - 0*P.di ~ 1, data=PofF.var), newdata=data.frame(P.di=1-input$risk))
      )

      dd.intersection <- switch(input$modtype,
        poly = predict(lm(f~poly(P.dd, degree=input$polydegree),data=PofF.var),newdata=data.frame(P.dd=1-input$risk)),
        gam = predict(gam(f~s(P.dd),data=PofF.var),newdata=data.frame(P.dd=1-input$risk)),
        gam.adaptive=predict( gam(f~s(P.dd, k=input$knots,bs="ad"), data=PofF.var),newdata=data.frame(P.dd=1-input$risk)),
        mpi= predict(scam(f~s(P.dd, bs="mpi"),data=PofF.var),newdata=data.frame(P.dd=1-input$risk)),
        mpd= predict(scam(f~s(P.dd, bs="mpd"),data=PofF.var),newdata=data.frame(P.dd=1-input$risk)),
        cx= predict(scam(f~s(P.dd, bs="cx"),data=PofF.var),newdata=data.frame(P.dd=1-input$risk)),
        cv= predict(scam(f~s(P.dd, bs="cv"),data=PofF.var),newdata=data.frame(P.dd=1-input$risk)),
        micx= predict(scam(f~s(P.dd, bs="micx"),data=PofF.var),newdata=data.frame(P.dd=1-input$risk)),
        micv= predict(scam(f~s(P.dd, bs="micv"),data=PofF.var),newdata=data.frame(P.dd=1-input$risk)),
        mdcx= predict(scam(f~s(P.dd, bs="mdcx"),data=PofF.var),newdata=data.frame(P.dd=1-input$risk)),
        mdcv= predict(scam(f~s(P.dd, bs="mdcv"),data=PofF.var),newdata=data.frame(P.dd=1-input$risk)),
        avg= predict(lm(f - 0*P.dd ~ 1, data=PofF.var), newdata=data.frame(P.dd=1-input$risk))
      )

      rect(0,0,di.intersection,1-input$risk,lty=2,border="darkgrey")
      rect(0,0,dd.intersection,1-input$risk,lty=2,border="darkgrey")
      box()

      mtext(side=1,outer=F,text="Exploitation rate",line=4)
      mtext(outer=T,side=2,text="Probability of being at or above biomass objective in 10 years",line=-1)
    })
  })

  # STEP 4 - plot contour plot
  observeEvent(input$runCCF, { 
  ems <<- vector()
    for(i in 1:input$Ems.length){
      ems[i] <<- eval(parse(text=paste0("input$e",i)))
      # print(eval(parse(text=paste0("input$f",i))))
    }
  output$contour <- renderPlot({
    ECCF <<- Eproj.list.f(Emean.shifts=ems, N=input$N.CCF, proj.years=input$proj.years, 
      Edist.a=Edist.a, Edist.b=Edist.b)
    PBCCF <<- PBproj.list.f(PvsE=PvsE, Eprojection=ECCF)
    CCF.raw <<- P.R.for.EF.f(E.CCF=ECCF, PB.CCF=PBCCF, 
      Fs=fs, PB=PB, ref.pt=input$ref.pt, Bstart.mult=input$Bstart.mult, K=input$K, theta=input$theta)

    CCF.contour <<- interp(x=CCF.raw$E.med,y=CCF.raw$Fval,z=CCF.raw$P.di)
    contour(CCF.contour,xlab="Median temperature (°C)",ylab="Exploitation rate",xaxs="i",yaxs="i")
    risk.equi.exp.rate.di <- contourLines(CCF.contour$x,CCF.contour$y,CCF.contour$z,nlevels=1,levels=1-input$risk)
    confint(risk.equi.exp.rate.di[[1]]$x,risk.equi.exp.rate.di[[1]]$y*0,risk.equi.exp.rate.di[[1]]$y,col=rgb(0, 1, 0,0.5))
    lines(PB$E,PB$F.rel,col="slateblue",lwd=1)
    points(PB$E,PB$F.rel,col="slateblue",pch=20)
    year.endpoints= match(range(PB$Year),PB$Year)
    points(PB$E[year.endpoints],PB$F.rel[year.endpoints],pch=22,cex=3,bg="white",col="slateblue")
    text(PB$E[year.endpoints],PB$F.rel[year.endpoints],PB$Year[year.endpoints],cex=.5)
  })
})

  ### NAVIGATION BETWEEN TABS
  # Observe hitting "next" button:
  observeEvent(input$next_tab1, {
    updateTabsetPanel(session, "tabs", "Step 2. Create E and PB projections")
  })
  observeEvent(input$next_tab2, {
    updateTabsetPanel(session, "tabs", "Step 3. Generate more complex models with exploitation")
  })
  observeEvent(input$next_tab3, {
    updateTabsetPanel(session, "tabs", "Step 4. Generate a contour plot of safe operating space")
  })
}

shiny::shinyApp(ui = ui, server = server)