library( dplyr )
library( tidyr )

source( "R/load-data.R")
source( "R/makeLongTable.R" )

m_df <- makeLongTable( select( metabolomics, -`KEGG IDs` ),
                       setdiff( colnames( metabolomics ), c("Metabolites", "KEGG IDs") ) )

l_df <- makeLongTable( select( lipidomics, -`KEGG ID` ),
                       setdiff( colnames( lipidomics ), c("Lipid", "KEGG ID" ) ) )

p_df <- bind_rows( lapply(
  list( proteomics1, proteomics2, proteomics3, proteomics4 ),
  function ( df ) {
    df2 <- select( df, -`Protein names`, -`Gene names`, -`Fasta headers` )
    makeLongTable( df2, colnames( df2 )[-1] )
  }
) )

long_table <- bind_rows( lapply(
  list( m_df, l_df, p_df ),
  . %>% rename( `Molecule ID` = 1 )
) )

saveRDS( long_table, file = "processed-data/long_table.rds" )
