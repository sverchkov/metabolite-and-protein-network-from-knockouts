# Command line arguments:
# 1. RDS file containing feature matrix, with IDs in first column
# 2. Number n at which to compute the dot products of the nth row with rows 1:n
# 3. Output file name

cli_args <- commandArgs( trailingOnly = T )

getArg <- function ( i, default ){
  value <- cli_args[i]
  if ( is.na( value ) ) default
  else value
}

in_file <- getArg( 1, "feature_t_df.rds" )
n <- getArg( 2, 1 )
out_file <- getArg( 3, sprintf( "dot-products-%d.rds", n ) )

feature_df <- readRDS( in_file )

v1 <- as.vector( feature_df[n,-1] )

block <- do.call( rbind, lapply( 1:n, function ( b ) {
  prd <- sum( as.vector( feature_df[b,-1] ) * v1, na.rm = T )
  data.frame( a = feature_df[n,1], b = feature_df[b,1], product = prd )
} ) )

saveRDS( block, file = out_file )