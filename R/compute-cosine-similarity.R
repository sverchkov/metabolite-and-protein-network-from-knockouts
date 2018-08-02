# What I tried to do for presentation July 30, 2018
library("dplyr")
library("futile.logger")

flog.threshold(TRACE)

# Load table
tested_omics_table <- readRDS( "processed-data/tested-omics-table.rds" )

# Make table shape
spread_table <- tested_omics_table %>%
  mutate( t = `Mean log2FC`/`SD of log2FC` ) %>%
  tidyr::spread( key = Knockout, value = t, fill = 0 ) %>%
  select( -`Mean log2FC`, -`SD of log2FC`, -`p-Value` ) %>%
  group_by( `Molecule ID`, `Molecule Name`, `Molecule Type` ) %>%
  summarize_all( sum ) %>%
  ungroup() %>%
  arrange( desc(`Molecule Type`), `Molecule ID` ) %>%
  mutate( id = dplyr::row_number() )

# Figure out how many proteins there are
n_proteins <- ( spread_table %>% filter( "Protein" == `Molecule Type` ) %>% summarize( count=n() ) )$count

# Make matrix-shaped df
matrix_df <- spread_table %>% select( -`Molecule ID`, -`Molecule Name`, -`Molecule Type`, -id )

# Get number of molecules for convenience
n_molecules <- nrow( matrix_df )

# Estimate the number of dot products we will compute
n_dots <- n_molecules * (n_molecules + 1)/2

# Get dot products, progressive computation
dot_products <- NULL
result_rows <- 0
added_rows <- 0

for ( a in 1:n_molecules ){
  v1 <- as.vector( matrix_df[a,] )
  block <- bind_rows( Map( function ( b ) {
      prd <- sum( as.vector( matrix_df[b,] ) * v1 )
      tibble( a = a, b = b, product = prd )
    }, 1:a ) )
  
  dot_products <- rbind( dot_products, block )
  added_rows <- added_rows + nrow( block )
  if ( added_rows > 10000 ){
    result_rows <- result_rows + added_rows
    added_rows <- 0
    flog.trace("Computed %s of %s dot products (%s%%)", result_rows, n_dots, 100*result_rows/n_dots )
    saveRDS( dot_products, "processed-data/dot-products.rds" )
  }
}

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
