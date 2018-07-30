#
# This is the server logic of a Shiny web application. You can run the 
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)
library(ggplot2)
library(dplyr)

# Load dataset (hard-coded)
full_omics_table <- readRDS( "data/full-omics-table.rds" ) %>%
  mutate( Significance = ifelse( `p-Value` < 0.05,
                                 ifelse( abs( `Mean log2FC` ) > 1, "p<0.5 and |log_2 FC|>1", "p<0.5" ),
                                 "p>=0.05" ),
          call = Significance == "p<0.5 and |log_2 FC|>1", "p<0.5" )

knockouts <- unique( full_omics_table$Knockout )

# Server!
shinyServer(function(input, output, session) {
  
  updateSelectInput( session, "knockoutSelect", choices = knockouts )
  
  #output$info <- renderText( input$moleculeInput )
  
  filteredData <- reactive({
    filtered_table <- full_omics_table %>%
      filter( Knockout %in% input$knockoutSelect, `Molecule Type` %in% input$moleculeInput )
    
    if( input$filterCommon ){
      select_molecules <- filtered_table %>%
        group_by( `Molecule ID` ) %>%
        summarize( common = all( call ) ) %>%
        ungroup() %>%
        filter( common )
      
      filtered_table <- inner_join( filtered_table, select_molecules, by = "Molecule ID" )
    }
    
    filtered_table
  })
  
  output$volcanoPlot <- renderPlot({
    ggplot( filteredData(), aes( x = `Mean log2FC`, y = -log10(`p-Value`), color = Significance, shape = `Molecule Type` ) ) +
      xlab( expression(mean ~ log[2] ~ fold ~ change) ) +
      ylab( expression(- log[10] ~ p ~ "-value") ) +
      geom_point()
  })

  output$selectedData <- renderDataTable( filteredData() )  
})
