library("dplyr")

input_matrix <- readRDS( "processed-data/matrix-for-cosim.rds" )

dot_prod_df <- NULL

n <- nrow( input_matrix )

for ( a in n:1 ){
  block <- bind_rows( Map( function( b ){
    tibble( a=a, b=b, product=( sum( input_matrix[a,]*input_matrix[b,] ) ) )
  }, n:a ) )
  
  dot_prod_df <- rbind( dot_prod_df, block )
  saveRDS( "processed-data/dot-products.rds" )
}