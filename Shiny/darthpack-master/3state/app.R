#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

#source("Functions_PSA.R")
# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Three-state microsimulaiton model"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            sliderInput("p_HS",
                        "Probability of transitioning from \"Healthy\" to \"Sick\" :",
                        min = 0,
                        max = 1,
                        value = 0.2)
        ),

        # Show a plot of the generated distribution
        mainPanel(
           plotOutput(outputId = "Trace"),
           tableOutput(outputId = "CEA")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
    
    model_inputs <- 
    model_inputs$p_HS <-input$p_HS
    
    run_sim <- Microsim(model_input,seed = 1)
    
    
    
    output$Trace <- renderPlot({
        # generate bins based on input$bins from ui.R
        run_sim$m.M
    })
    
    output$Res <- renderPlot({
        # generate bins based on input$bins from ui.R
        run_sim$dc
        
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
