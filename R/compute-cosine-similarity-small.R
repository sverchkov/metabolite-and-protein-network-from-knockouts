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
  sample_n( 10 ) %>%
  mutate( id = dplyr::row_number() )

# Make matrix-shaped df
matrix_df <- spread_table %>% select( -`Molecule ID`, -`Molecule Name`, -`Molecule Type`, -id )

# Get dot products
n_molecules <- nrow( matrix_df )
dot_products <-
  bind_rows( Map( function ( a ) {
    v1 <- as.vector( matrix_df[a,] )
    bind_rows( Map( function ( b ) {
      prd <- sum( as.vector( matrix_df[b,] ) * v1 )
      tibble( a = a, b = b, product = prd )
    }, 1:a ) )
  }, 1:n_molecules ) )

# Get norms
norms <- dot_products %>% filter( a == b ) %>%
  transmute( v = a, norm = sqrt( product ) )

# Get cosine simiarity = a.b/|a||b|
similarities <- dot_products %>% filter( a != b ) %>%
  left_join( norms, by = c( a = "v" ) ) %>% rename( a_norm = norm ) %>%
  left_join( norms, by = c( b = "v" ) ) %>% rename( b_norm = norm ) %>%
  mutate( similarity = product/(a_norm*b_norm) ) %>%
  arrange( desc( abs( similarity ) ) )

saveRDS( similarities, "processed-data/cosine-similarities.rds" )
saveRDS( spread_table, "processed-data/spread-table-for-similarities.rds" )

