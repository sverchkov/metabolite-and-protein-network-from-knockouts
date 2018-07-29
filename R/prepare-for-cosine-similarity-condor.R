# What I tried to do for presentation July 30, 2018
library("dplyr")

# Load table
computed_omics_table <- readRDS( "H3KExplorer/data/full-omics-table.rds" )

# Make table shape
spread_table <- computed_omics_table %>%
  mutate( t = `Mean log2FC`/`SD of log2FC` ) %>%
  tidyr::spread( key = Knockout, value = t, fill = 0 ) %>%
  select( -`Mean log2FC`, -`SD of log2FC`, -`p-Value` ) %>%
  group_by( `Molecule ID`, `Molecule Name`, `Molecule Type` ) %>%
  summarize_all( sum ) %>%
  ungroup() %>%
  arrange( desc(`Molecule Type`), `Molecule ID` ) %>%
  mutate( id = dplyr::row_number() )

saveRDS( spread_table, "processed-data/spread-table-for-similarities.rds" )

# Make matrix-shaped df
matrix_df <- spread_table %>% select( -`Molecule ID`, -`Molecule Name`, -`Molecule Type`, -id )

for ( i in 1:nrow( matrix_df ) )
  saveRDS( as.vector( matrix_df[i,] ), sprintf( "processed-data/for-condor/row-%d.rds", i ) )