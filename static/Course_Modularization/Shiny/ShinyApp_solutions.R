#---------------------------------------------------------------------------#
#### R functions to develop the shiny app                               ####
#---------------------------------------------------------------------------#

# Developed by the Decision Analysis in R for Technologies in Health (DARTH) group
# Fernando Alarid-Escudero, PhD (1) 
# Eva A. Enns, MS, PhD (2)	
# M.G. Myriam Hunink, MD, PhD (3,4)
# Hawre J. Jalal, MD, PhD (5) 
# Eline M. Krijkamp, MSc (3)
# Petros Pechlivanoglou, PhD (6) 

# In collaboration of: 		
# 1 Center for Research and Teaching in Economics (CIDE), Drug Policy Program, Mexico
# 2 University of Minnesota School of Public Health, Minneapolis, MN, USA
# 3 Erasmus MC, Rotterdam, The Netherlands
# 4 Harvard T.H. Chan School of Public Health, Boston, USA
# 5 University of Pittsburgh Graduate School of Public Health, Pittsburgh, PA, USA
# 6 The Hospital for Sick Children, Toronto and University of Toronto, Toronto ON, Canada

# clear memory of any R objects 
rm(list=ls())

# load libraries

if (!require('pacman')) install.packages('pacman'); library(pacman) # use this package to conveniently install other packages
# load (install if required) packages from CRAN
p_load("shiny", "dplyr", "devtools") 
# install_github("DARTH-git/darthtools", force = TRUE) #Uncomment if there is a newer version
p_load_gh("DARTH-git/darthtools")

# load functions that are useful  for:
source("Function_MicroSim_3state.R")  # execution of Microsimulation 
source("Functions_Input.R")           # handling input

# Define UI for application that presents the results of the 3-state model
# this is the "layout of the 
ui <- fluidPage(

# Application title and subtite - note the use of HTML language ("h3")
titlePanel("Three-state microsimulaiton model"),
h3("A Simple Example of R Shiny capabilities"),

# Sidebar with :
#     a slider input for p_HS and 
#     a numeric input for number of individuals in simulation
#     a numeric input for time horizon
    sidebarLayout(
        sidebarPanel(
          h2("Input parameters", align = "center"),
          br(),
          br(),
          sliderInput(inputId    = "p_HS",
                        label    = "Probability of transitioning from \"Healthy\" to \"Sick\" :",
                        min      = 0,
                        max      = 1,
                        value    = 0.2),
           numericInput( inputId = "n_i",
                         label   =  "Number of Individuals to simulate",
                         min     = 2,
                         max     = 10000,
                         value   = 1000 ),
           numericInput( inputId = "n_t",
                         label   =  "Number of cycles",
                         min     = 2,
                         max     = 100,
                         value   = 30 ),
           sliderInput( inputId  = "c_S",
                        label    =  "Cost in the Sick State",
                        min      = 100,
                        max      = 1000,
                        value    = 500),
           sliderInput( inputId  = "u_S",
                        label    =  "Utility in the Sick State",
                        min      = 0,
                        max      = 1,
                        value    = 0.85)),
    # Presetn a trace plot in the main panel together with the table of results
        mainPanel(
          br(),
        h2("Main Results",align= "center"),

        h3("Trace Plot"),
        plotOutput( outputId = "Trace"),
        
        h3("Table of Results"),
        tableOutput(outputId = "results")
        )
    )  
)

# in the server section is where the model running is happening. 
# the server function loads the user input from the user interface above and creates an output list with objects that have same  names as the outputId values

server <- function(input, output) {
# object named "Trace": the Trace plot stored throuht the renderPlot function 
 output$Trace <- renderPlot({

# we load all input parameters that will be used in the model. For more details on how we do that please see file "Functions_Input.R" 
        params  <- generate_params(input) 

        
# we run the model using the function we have specified during our Microsimulation day      
        run_sim <- calculate_ce_out(params)

# generate Markov Trace plot based on a (modified) plot_m_TR
        plot_trace_microsim_shiny(run_sim$Trace, params)
    })
 
# the second object in the output list is the table of results. We use the function renderTable to grab the results table.
 
 output$results <- renderTable(expr ={
     # we load all input parameters that will be used in the model. For more details on how we do that please see file "Functions_Input.R" 
   
     params  <- generate_params(input) 
     
     # we run the model using the function we have specified during our Microsimulation day      
     run_sim <- calculate_ce_out(params)
     
     # here we print the data frame with the results .
     run_sim$df_ce
 })
}

# Run the application 
shinyApp(ui = ui, server = server)
