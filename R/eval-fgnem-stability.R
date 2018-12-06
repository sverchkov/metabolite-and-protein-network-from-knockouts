## Fgnem stability evaluation

# Load things
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