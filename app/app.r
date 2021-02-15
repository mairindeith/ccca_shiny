# Load libraries, assuming these have been installed
# First, gslea

library(ccca) # includes dependency: datatable, but nothing else
library(shinyjs)
library(shiny)

ui <- fluidPage(
  titlePanel('Climate change conditioned advice (CCCA) for fisheries'),
  tabsetPanel(
#    widths=c(2,10),
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
        numericInput('Emean.shift.warm', 'Mean upper shift of environmental variable (e.g. temperature is E, a value of 0.5 = 0.5 degrees WARMER)', value = 0.5, min=0),
        numericInput('Emean.shift.cold', 'Mean lower shift environmental variable (e.g. temperature is E, a value of -0.5 = 0.5 degrees COLDER)', value = -0.5, max=0),
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
            actionButton("next_tab1", "Next", width='100%')
          )
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
      # ccca::colramp.legend(col1="red", col2="blue", ncol=length(PB$E), 2.5, 3.5, 2.7, 4.5)
    PvsE <<- PBE.fit.f(PB, model.type=input$modtype, poly.degree=input$polydegree, knots=input$knots)
    PvsE.null <<- PBE.fit.f(PB, model.type="avg", knots=knots, poly.degree=poly.degree)
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
      Enorm= norm.fit.f(E=PB$E)
      Nrand=1000000
      Edist.a=Enorm$estimate[1]
      Edist.b=Enorm$estimate[2]
      Ebase=density(norm.plot.f(Nrand=Nrand, Edist.a=Edist.a, Edist.b=Edist.b,Emean.shift=0,E.var.inc=1))
      Ewarm=density(norm.plot.f(Nrand=Nrand, Edist.a=Edist.a, Edist.b=Edist.b,Emean.shift=input$Emean.shift.warm,E.var.inc=1))
      Ecold=density(norm.plot.f(Nrand=Nrand, Edist.a=Edist.a, Edist.b=Edist.b,Emean.shift=input$Emean.shift.cold,E.var.inc=1))
      Evar=density(norm.plot.f(Nrand=Nrand, Edist.a=Edist.a, Edist.b=Edist.b,Emean.shift=0,E.var.inc=1.5))
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

  # Observe hitting "next" button:
    observeEvent(input$next_tab1, {
      updateTabsetPanel(session, "tabs", "Step 2. Generate more complex models")
    })
  
  ### TAB 2: STEP 2 MORE COMPLEX MODELS

}

shiny::shinyApp(ui = ui, server = server)