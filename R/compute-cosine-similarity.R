# What I tried to do for presentation July 30, 2018
library("dplyr")
library("futile.logger")
library("parallel")

flog.threshold(TRACE)

# Load table
feature_df <- readRDS( "processed-data/feature_t_df.rds" )

# Get the id list for future reference
id_list <- pull( feature_df, `Molecule ID` )

# Remove id column to get matrix
matrix_df <- feature_df %>% select( -`Molecule ID` )

# Get number of molecules for convenience
n_molecules <- nrow( matrix_df )

# Get a list of row vectors
vectors <- lapply( 1:n_molecules, function( n ) as.vector( matrix_df[n,] ) )

# Estimate the number of dot products we will compute
n_dots <- n_molecules * (n_molecules + 1)/2

# Define pattern for intermediate files
filepattern <- "temporary-data/dot-products/dots-%d.rds"

# Get dot products, parallel computation on the lower-level
cl <- makeCluster( detectCores()-1 )
added_rows <- 0
result_rows <- 0
for ( n in 1:n_molecules ){
  
  outfile <- sprintf( filepattern, n )
  
  if( !file.exists( outfile ) ){
    
    block <- tibble(
      a = id_list[n],
      b = id_list[1:n],
      product = parSapply(
        cl,
        vectors[1:n],
        function ( v1, v2 ) sum( v1 * v2, na.rm = T ),
        vectors[n] )
      )
    
    saveRDS( block, outfile )
  }
  
  added_rows <- added_rows + n
  if ( added_rows > 10000 ){
    result_rows <- result_rows + added_rows
    added_rows <- 0
    flog.trace("Computed %s of %s dot products (%7.3f%%)", result_rows, n_dots, 100*result_rows/n_dots )
  }
}

# Read files into one table
flog.trace( "Reading temp files into one table..." )
dot_products <- bind_rows( lapply( 1:n_molecules, function( a ) readRDS( sprintf( filepattern, a ) ) ) )
flog.trace( "...table ready." )

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
