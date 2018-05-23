#######################################
## Prepare data for small experiment ##
#######################################

library("dplyr")
library("tidyr")
library("ggplot2")

# Knockouts of interest:
knockouts <- c( "COQ5", "COQ7", "COQ8A", "COQ8B", "COQ9" )

knockouts_KEGG <- bind_rows( Map( function ( gene_name ) {
  result <- keggFind( "hsa", gene_name )
  id <- names( result )
  tibble( `Gene Name` = gene_name, `KEGG ID` = id, Description = result )
}, knockouts ) )

#enzyme_list <- names( genes.by.enzymes )
#enzyme_table <- tibble( Enzyme = enzyme_list, `KEGG gene ID` = genes.by.enzymes )

knockout_reference_table <- left_join( knockouts_KEGG, master_knockout_table, by = c( `Gene Name` = "KO" ) )

# Let's work with the unimputed proteomics table for this.
proteomicsKOs <- knockout_reference_table %>% filter( Proteomics ) %>% select( Knockout )
proteomics_subset <- proteomics_unimputed %>% select( c( identifying_columns, proteomicsKOs$Knockout ) )
long_table <- proteomics_subset %>% gather( key = "Knockout", value = "log2fc", proteomicsKOs$Knockout )

long_table <- left_join( long_table, knockout_reference_table )

# Error-catching for t-test
safe.t.test <- function( x ){
  tryCatch( t.test(x)$p.value, error = function( anything ) NA )
}

# T-test to make calls
call_table <- long_table %>% group_by( `Fasta headers`, `Gene Name` ) %>%
  summarize( p.value = safe.t.test( log2fc ), `Mean log2fc` = mean( log2fc ), `Number of replicates` = n() )

hit_counts <- call_table %>% group_by( `Fasta headers` ) %>% filter( p.value < 0.05 ) %>% summarize( n = n() )

# Pick out a sample of proteins
sample_of_proteins <- Reduce( c, Map(
  function ( num ) sample( ( hit_counts %>% filter( n == num ) )$`Fasta headers`, size = 5 ),
  1:5 ) )

small_table <- call_table %>% filter( `Fasta headers` %in% sample_of_proteins )

ggplot( small_table, aes( x = `Fasta headers`, y = `Gene Name`, color = `Mean log2fc` ) ) + geom_point() + scale_color_distiller( palette="Spectral" )

####################################
### Getting matching metabolomics ##
####################################
metabolomics_subset <- metabolomics_unimputed %>% select( c( identifying_columns_m, proteomicsKOs$Knockout ) )
metabolomics_long_table <- metabolomics_subset %>% gather( key = "Knockout", value = "log2fc", proteomicsKOs$Knockout )
metabolomics_long_table <- left_join( metabolomics_long_table, knockout_reference_table )
metabolomics_call_table <- metabolomics_long_table %>% group_by( `Metabolites`, `Gene Name` ) %>%
  summarize( p.value = safe.t.test( log2fc ), `Mean log2fc` = mean( log2fc ), `Number of replicates` = n() )