#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
rm(list=ls())
library(shiny)
library(dplyr)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # set the working directory 

# load functions that are useful  for
source("Function_MicroSim.R")  # execution of Microsimulation 
source("Functions_Input.R")    # handling input
source("Functions.R")          # general functions 

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Three-state microsimulaiton model"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            sliderInput(inputId = "p_HS",
                        label = "Probability of transitioning from \"Healthy\" to \"Sick\" :",
                        min = 0,
                        max = 1,
                        value = 0.2),
        numericInput(inputId = "n_i",
                     label   =  "Number of Individuals to simulate",
                     min     = 2,
                     max     = 10000,
                     value   = 1000 ),
        numericInput(inputId = "n_t",
                     label   =  "Number of cycles",
                     min     = 2,
                     max     = 100,
                     value   = 30 )),

        # Show a plot of the generated distribution
        mainPanel(
            
        plotOutput( outputId = "Trace"),
        tableOutput(outputId = "results")
        )
    )  
)
#input = list(p_HS = 0.5)
# Define server logic required to draw a histogram
server <- function(input, output) {
 output$Trace <- renderPlot({
    # browser()
     
        params  <- generate_params(input) 
        run_sim <- calculate_ce_out(params)
        # generate bins based on input$bins from ui.R
        plot_m_TR(run_sim$Trace, params)
        
    })
 output$results <- renderTable(expr ={
     params  <- generate_params(input) 
     run_sim <- calculate_ce_out(params)
      run_sim$df_ce
 })

   # output$Res <- renderPlot({
        
    #    df_ce <- data.frame(Cost   = run_sim$tc_hat,
   #                         Effect = run_sim$te_hat)
     #   return(df_ce)
        
    #})
}

# Run the application 
shinyApp(ui = ui, server = server)
