## Fgnem stability evaluation

# Load things
library("dplyr")
# TODO: Load things

# Function for converting {-1,0,1} matrix to 1-hot m x n x 3
convertToOneHot <- function ( m ) {
  array(
    data = c( as.vector( m == 1 ), as.vector( m == 0 ), as.vector( m == -1 ) ),
    dim = c( dim( m ), 3 ),
    dimnames = c( dimnames( m ), list( c( "1", "0", "-1" ) ) )
  )
}

# Get mean counts
average_acc <- (1/b) * Reduce( `+`, Map( function ( obj ) convertToOneHot( obj$acc ), nem_bootstrap_results ) )

# Turn edges into DF
row_names <- rownames( average_acc )
col_names <- colnames( average_acc )
n_rows <- length( row_names )
n_cols <- length( col_names )

edge_list <- tibble(
  source = rep( row_names, times = n_cols*2 ),
  target = rep( rep( col_names, each = n_rows ), times = 2 ),
  interaction = rep( c( "positive", "negative" ), times = n_rows*n_cols ),
  proportion = as.vector( average_acc[,,c("1","-1")] ) )

edge_list <- filter( edge_list, proportion != 0 )
