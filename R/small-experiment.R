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

ggplot( long_table, aes( `Protein names`, log2fc ) ) + geom_point()
        