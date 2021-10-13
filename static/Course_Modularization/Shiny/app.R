# Load required packages
library(shiny)
library(ggplot2)

# Import and subset data
framingham <- read.csv("framingham.csv")
# select columns with quantitative data of interest
quantvars <- c("AGE", "TOTCHOL", "HDLC", "LDLC", "SYSBP", 
               "DIABP", "BMI", "HEARTRTE", "GLUCOSE")
framingham <- framingham[, quantvars] 

# Sets up appearance of app and fields user input
ui <- fluidPage(
  
  # Specify title of app
  titlePanel("Framingham Summary"),
  
  # Create an input variable with inputID "var"
  selectInput(inputId = "var", label = "Variable", choices = names(framingham)),
  
  # Add output$plot ('plot') to the app
  plotOutput('plot1'),
  
  # Add output$table ('table') to the app
  tableOutput('table1')
)

# The server function takes the inputs and creates the outputs
server <- function(input, output) {
  
  # Create plot1 output object that is reactive to changes in input$var
  output$plot1 <- renderPlot({
    
    # Rows with missing values in the variable of interest are removed from the data.frame
    framingham_plot <- framingham[which(!is.na(framingham[, input$var])),]
    
    # Create and return ggplot2 object
    ggplot(data = framingham_plot, aes(x = .data[[input$var]])) +
      geom_density()
  })
  
  # Create table1 output object that is reactive to changes in input$var
  output$table1 <- renderTable({
    data.frame(Median = median(framingham[, input$var], na.rm = TRUE),
               Mean = mean(framingham[, input$var], na.rm = TRUE),
               StandardDeviation = sd(framingham[, input$var], na.rm = TRUE),
               N = sum(!is.na(framingham[, input$var])),
               Missing = sum(is.na(framingham[, input$var])))
  })
}

# The final function brings together the ui and the server function to create the final Shiny app object
shinyApp(ui = ui, server = server)

# To run the completed app, left click the "Run App" button 
# in the top right of the app.R script window 



