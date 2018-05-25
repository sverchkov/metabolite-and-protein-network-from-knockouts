#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)

# Define UI for application that draws a histogram
shinyUI(fluidPage(
  
  # Application title
  titlePanel("H3K Data Explorer"),
  
  # Sidebar with a slider input for number of bins 
  sidebarLayout(
    sidebarPanel(
      selectInput( "knockoutSelect", "Knockout", choices = NULL ),
      checkboxGroupInput( "moleculeInput",
                          "Molecules to show:",
                          c( Proteins="Protein", Metabolites="Metabolite", Lipids="Lipid" ) )
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
       plotOutput("volcanoPlot")
       #textOutput("info")
    )
  )
))
