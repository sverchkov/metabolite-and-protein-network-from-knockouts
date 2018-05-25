#' Convert table to long-form
#' 
#' @param wide_table wide-format table where columns are correspond to replicates in a KO-CRISPR#-Bio# format
#' @param knockout_columns the list of knockout columns to gather
#' @param value what to name the value column being gathered
#' @return
#' @import dplyr
#' @import tidyr
#' @import rlang
#' @author Yuriy Sverchkov
makeLongTable <- function( wide_table, knockout_columns, value = sym("log2FC") ){
  
  # Derive knockout breakdown table
  ko_breakdown <- Reduce( rbind, strsplit( knockout_columns, "-", fixed = T ) )
  ko_breakdown <- tibble(
    `Knockout ID` = knockout_columns,
    `Knockout` = ko_breakdown[,1],
    `CRISPR Replicate` = ko_breakdown[,2],
    `Biological Replicate` = ko_breakdown[,3] )

  wide_table %>%
    gather( key = "Knockout ID", value = !! value, knockout_columns ) %>%
    left_join( ko_breakdown, by = "Knockout ID" )
}