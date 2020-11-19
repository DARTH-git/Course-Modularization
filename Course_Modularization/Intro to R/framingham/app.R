#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Framingham Data - Stratified Boxplots"),
    
    # fill in this part with continuous variables to show on the y-axis
    selectInput("variable_y", "Variable:",),
  
    selectInput("variable_x", "Variable:",
                c("Sex" = "SEX",
                  "Current smoker" = "CURSMOKE",
                  "Diabetic" = "DIABETES", 
                  "Prevalent Myocardial Infarction" = "PREVMI",
                  "Prevalent stroke" = "PREVSTRK")),
    
        # Show a plot of the generated distribution
        mainPanel(
           plotOutput("distPlot")
        )
    )

# Define server logic required to draw a histogram
server <- function(input, output) {

    output$distPlot <- renderPlot({
      
        # read in the framingham dataset
        data <- read.csv('framingham.csv', header = TRUE)
        
        # clean code (same code from before)
        library(dplyr)
        data <- data %>% 
          mutate(SEX = ifelse(!is.na(SEX), ifelse(SEX == 1, 'men', 'women'), NA)) %>%
          mutate(CURSMOKE = ifelse(!is.na(CURSMOKE), ifelse(CURSMOKE== 1, 'current smoker', 'not a current smoker'), NA)) %>%
          mutate(DIABETES = ifelse(!is.na(DIABETES), ifelse(DIABETES== 1, 'diabetic', 'non-diabetic'), NA)) %>%
          mutate(PREVMI= ifelse(!is.na(PREVMI), ifelse(PREVMI == 1, 'prevalent myocardian infarction', 'no prevalent myocardian infarction'), NA)) %>%
          mutate(PREVSTRK = ifelse(!is.na(PREVSTRK), ifelse(PREVSTRK == 1, 'prevalent stroke', 'no prevalent stroke'), NA))
        
        # extract the selected variables from UI
        y <- data[, input$variable_y]
        x <- data[, input$variable_x]

        # draw the boxplot of y stratified by x

    })
}

# Run the application 
shinyApp(ui = ui, server = server)
