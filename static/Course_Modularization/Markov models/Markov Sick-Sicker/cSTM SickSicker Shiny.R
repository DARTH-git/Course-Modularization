library(shiny)

fun_SSS <- function(minage, maxage ,p.HS1 ,p.S1H ,p.S1S2,p.HD  ,rr.S1 ,rr.S2 ,d.r,c.H,c.S1,c.S2,c.Trt,c.D,u.S1,  u.S2 , u.Trt ){
  v.n  <- c("H", "S1", "S2", "D")               # the 4 states of the model: Healthy (H), Sick (S1), Sicker (S2), Dead (D)
  n.s  <- length(v.n)     
  Strategies <- c("No Treatment", "Treatment")  # strategy names 
  v.n  <- c("H", "S1", "S2", "D")               # the 4 states of the model: Healthy (H), Sick (S1), Sicker (S2), Dead (D)
  n.t <- maxage - minage 
  r.HD    <- - log(1 - p.HD)  
  r.S1D   <- rr.S1 * r.HD  	  
  r.S2D   <- rr.S2 * r.HD  	  
  p.S1D   <- 1 - exp(- r.S1D) 
  p.S2D   <- 1 - exp(- r.S2D) 
  u.D <- 0
  u.H <- 1
  v.dwe <- v.dwc <- 1 / ((1 + d.r) ^ (0:n.t))  # discount weight 
  # (equal discounting is assumed for costs and effects)
  
  #### 04 Model initialization ####
  # create transition probability matrices for both treatment and no treatment
  m.P_no_trt <- m.P_trt <- matrix(NA, nrow = n.s, ncol = n.s, 
                                  byrow = TRUE, dimnames = list(v.n, v.n))
  
  # fill in the transition probability matrix
  ### from H
  m.P_no_trt["H", "H"]  <- 1 - (p.HS1 + p.HD)
  m.P_no_trt["H", "S1"] <- p.HS1
  m.P_no_trt["H", "S2"] <- 0
  m.P_no_trt["H", "D"]  <- p.HD
  
  ### from S1
  m.P_no_trt["S1", "H"]  <- p.S1H
  m.P_no_trt["S1", "S1"] <- 1 - (p.S1H + p.S1S2 + p.S1D)
  m.P_no_trt["S1", "S2"] <- p.S1S2
  m.P_no_trt["S1", "D"]  <- p.S1D
  
  ### from S2
  m.P_no_trt["S2", "H"]  <- 0
  m.P_no_trt["S2", "S1"] <- 0
  m.P_no_trt["S2", "S2"] <- 1 - p.S2D
  m.P_no_trt["S2", "D"]  <- p.S2D 
  
  ### from D
  m.P_no_trt["D", "H"]   <- 0 
  m.P_no_trt["D", "S1"]  <- 0
  m.P_no_trt["D", "S2"]  <- 0 
  m.P_no_trt["D", "D"]   <- 1
  
  # the two transition probability matrices are the same
  m.P_trt <- m.P_no_trt
  m.P_no_trt # print the matrix to check the value 
  
  # create the markov trace 
  # matrix M is capturing the proportion of the cohort in each state at each cycle
  m.M_no_trt <- m.M_trt <- matrix(NA, nrow = n.t + 1, ncol = n.s,
                                  dimnames = list(paste("cycle", 0:n.t, sep = " "), v.n))
  
  # The cohort starts as healthy
  m.M_no_trt[1, ] <- m.M_trt[1, ] <- c(1, 0, 0, 0) # initiate the Markov trace 
  
  #### 05 Process of the model ####
  for (t in 1:n.t){
    ######### using transition matrices ###########
    # calculate the proportion of the cohort in each state at time t
    m.M_no_trt[t + 1, ] <- t(m.M_no_trt[t, ]) %*% m.P_no_trt
    m.M_trt[t + 1, ]   <- t(m.M_trt[t, ])    %*% m.P_trt
  } # close the loop
  
  v.u_trt    <- c(u.H, u.Trt, u.S2, u.D)
  v.u_no_trt <- c(u.H,  u.S1, u.S2, u.D)
  
  v.c_trt    <- c(c.H, c.S1 + c.Trt, c.S2 + c.Trt, c.D)
  v.c_no_trt <- c(c.H, c.S1,         c.S2,         c.D)
  
  # estimate mean QALYs and costs
  v.E_no_trt <- m.M_no_trt %*% v.u_no_trt
  v.E_trt    <- m.M_trt    %*% v.u_trt
  
  v.C_no_trt <- m.M_no_trt %*% v.c_no_trt
  v.C_trt    <- m.M_trt    %*% v.c_trt
  
  ### discount costs and QALYs
  te_no_trt <- t(v.E_no_trt) %*% v.dwe  # 1x31 %*% 31x1 -> 1x1
  te_trt    <- t(v.E_trt)    %*% v.dwe
  
  tc_no_trt <- t(v.C_no_trt) %*% v.dwc
  tc_trt    <- t(v.C_trt)    %*% v.dwc
  
  
  
  DC          <- tc_trt - tc_no_trt      # calculate the difference in discounted costs between the two strategies 
  names(DC)   <- "Incremental costs"
  DE          <- te_trt - te_no_trt      # calculate the difference in discounted effects between the two strategies 
  names(DE)   <- "QALYs gained"
  ICER        <- DC / DE                 # calculate the ICER
  names(ICER) <- "ICER" 
  results     <- c(DC, DE, ICER)         # combine the results 
  
  # create full incremental cost-effectiveness analysis table 
  C <- round(c(tc_no_trt, tc_trt), 2)  # bind and round the total costs of the two strategies
  E <- round(c(te_no_trt, te_trt), 2)  # bind and round the total effects of the two strategies
  
  DC   <- c("", as.character(round(DC, 2)))   # round the delta of the costs (No Treatment is reference)
  DE   <- c("", as.character(round(DE, 2)))   # round the delta of the effects (No Treatment is reference)
  ICER <- c("", as.character(round(ICER, 2))) # round the ICER 
  
  table_Markov_SickSicker <- cbind(Strategies, C, E, DC, DE, ICER)     # combine all data in a table
  table_Markov_SickSicker <- as.data.frame(table_Markov_SickSicker)    # create a data frame 
  table_Markov_SickSicker                                              
  
  return(list(trace = m.M_no_trt,
              table = table_Markov_SickSicker,
              n.t = n.t))
  
}

###Function
ui <- fluidPage(
  # App title ----
  titlePanel("Sick-Sicker Markov Model"),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(width=5,
                 tabsetPanel(type="tabs",
                             tabPanel("Model Structure",
                                      # Input: Age ----
                                      sliderInput(inputId = "age",
                                                  label = "Cohort Age",value=c(25,55),
                                                  min = 0,max=100,step=1),
                                      # Input: Discount Rate ----
                                      numericInput(inputId = "Dis",
                                                   label = "Discount Rate",
                                                   value = 0.03,
                                                   min=0,
                                                   step=0.0001)
                                      ),
                             tabPanel("Probabilities", 
                                      # Input: Probabilities----
                                      numericInput(inputId = "HS",
                                                   label = "Probability of Transition Healthy to Sick",
                                                   value = 0.15,
                                                   min=0,max=1),
                                      
                                      numericInput(inputId = "SH",
                                                   label = "Probability of Transition Sick to Healthy",
                                                   value = 0.5,
                                                   min=0,max=1),
                                      
                                      numericInput(inputId = "SSR",
                                                   label = "Probability of Transition Sick to Sicker",
                                                   value = 0.105,
                                                   min=0,max=1),
                                      
                                      numericInput(inputId = "HD",
                                                   label = "Probability of Transition Healthy to Dead",
                                                   value = 0.005,
                                                   min=0,max=1),
                                      # Input: Rate Ratio ----
                                      numericInput(inputId = "SD",
                                                   label = "Rate Ratio of Sick to Dead",
                                                   value = 3,
                                                   min=0),
                                      numericInput(inputId = "SRD",
                                                   label = "Probability of Transition Sicker to Dead",
                                                   value = 10,
                                                   min=0)
                                      ),
                             tabPanel("Utilities",                # Input: Utilities ----
                                      
                                      numericInput(inputId = "uS",
                                                   label = "Utility of Sick",
                                                   value = 0.75,
                                                   max=1,
                                                   min=0),
                                      
                                      numericInput(inputId = "uSR",
                                                   label = "Utility of Sicker",
                                                   value = 0.5,
                                                   max=1,
                                                   min=0),
                                      
                                      numericInput(inputId = "uT",
                                                   label = "Utility of Sick Patients on Treatment",
                                                   value = 0.95,
                                                   max=1,
                                                   min=0)
                                      ),
                             tabPanel("Costs",                # Input: Costs ----
                                      numericInput(inputId = "cH",
                                                   label = "Cost of Healthy",
                                                   value = 2000,
                                                   min=0),
                                      
                                      numericInput(inputId = "cS",
                                                   label = "Cost of Sick",
                                                   value = 4000,
                                                   min=0),
                                      
                                      numericInput(inputId = "cSR",
                                                   label = "Cost of Sicker",
                                                   value = 15000,
                                                   min=0),
                                      
                                      numericInput(inputId = "cT",
                                                   label = "Cost of Treatment",
                                                   value = 12000,
                                                   min=0)
                                      )
                             
                   
                 ),
                 actionButton("button", "Run")
    ),
    # Main panel for displaying outputs ----
    mainPanel(width=7,
              p("Evalulating the Cost-Effectiveness of a Treatment to Improve Quality of Life for Sick Patients using a Sick-Sicker Markov Model"),
              tabsetPanel(type="tabs",
                          tabPanel("Trace Plot",  
                                   #Output: Matplot ----
                                   plotOutput(outputId = "Plot",width="500px",height="500px")),
                          tabPanel("Cost-Effectiveness Analysis",
                                   #Output: Print CE Analysis
                                   tableOutput(outputId="Stopping"))
                                   )
                                   )
                          ),
  #BANTER!!!
  p("Copyright: ",tags$b("David Rios, Anna Heath")," and Petros Pechlivanoglou",Sys.time()))
                          

server <- function(input, output) {
  
  # Histogram of the Old Faithful Geyser Data ----
  # with requested number of bins
  # This expression that generates a histogram is wrapped in a call
  # to renderPlot to indicate that:
  #
  # 1. It is "reactive" and therefore should be automatically
  #    re-executed when inputs (input$bins) change
  # 2. Its output type is a plot
  
  observeEvent(input$button, {
    withProgress(message = 'Performing Health Economic Analysis', value = 0, {
      ages<-as.numeric(input$age)
      CEAnalysis<-fun_SSS(minage = ages[1], maxage = ages[2] ,
                   input$HS ,input$SH , input$SSR,input$HD  , input$SD, input$SRD ,
                   input$Dis,
                   input$cH,input$cS,input$cSR,input$cT,c.D=0,
                   input$uS,input$uSR , input$uT )
    })
    
    output$Stopping<-renderTable({
      CEAnalysis$table
    },digits=6)
    
    output$Plot <- renderPlot({
      matplot(0:CEAnalysis$n.t,CEAnalysis$trace, type = 'l',
              ylab = "Probability of state occupancy",
              xlab = "Cycle",
              main = "Markov Trace")              # create a plot of the Markov trace
      legend("topright", c("H", "S1", "S2", "D")  , col = 1:4, lty = 1:4, bty = "n")
    })
    
  })
  
}

shinyApp(ui = ui, server = server)
debug(fun_SSS)
