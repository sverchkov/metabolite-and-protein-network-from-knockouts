# Profiling matrix construction

if ( !exists( "limited_edges" ) ) stop()

usingLists <- function( edges ){
  n <- length( unique( c( edges$a, edges$b ) ) )
  a_list = rep( list(rep( list(0), n) ), n )
  for ( row in 1:nrow( edges ) ){
    i <- edges$a[row]
    j <- edges$b[row]
    a_list[[i]][[j]] <- edges$weight[row]
    a_list[[j]][[i]] <- edges$weight[row]
  }
  Reduce( rbind, Map( function( l ) Reduce( cbind, l ), a_list ) )
}

buildingMatrix <- function( edges ){
  n <- length( unique( c( edges$a, edges$b ) ) )
  mapply( function(i){
    mask_a <- edges$a == i
    mask_b <- edges$b == i
    jays <- c( edges$b[mask_a], i, edges$a[mask_b] )
    w <- c( edges$weight[mask_a], 0, edges$weight[mask_b] )
    w[order(jays)]
  }, 1:n, SIMPLIFY = TRUE )
}

buildingMatrix2 <- function( edges ){
  n <- length( unique( c( edges$a, edges$b ) ) )
  mapply( function(i){
    mask_a <- edges$a == i
    mask_b <- edges$b == i
    jays <- c( edges$b[mask_a], edges$a[mask_b] )
    w <- numeric(n)
    w[jays] <- c( edges$weight[mask_a], edges$weight[mask_b] )
    w
  }, 1:n, SIMPLIFY = TRUE )
}

usingParallel <- function ( edges ) {
  n <- length( unique( c( edges$a, edges$b ) ) )
  parallel::parSapply( cl, X = 1:n, FUN = function(i){
    mask_a <- edges$a == i
    mask_b <- edges$b == i
    jays <- c( edges$b[mask_a], edges$a[mask_b] )
    w <- numeric(n)
    w[jays] <- c( edges$weight[mask_a], edges$weight[mask_b] )
    w
  }, simplify = TRUE )
}

profvis::profvis( {
  mb_am <- buildingMatrix( limited_edges )
  lists_am <- usingLists( limited_edges )
  parallel_am <- usingParallel( limited_edges )
  mb2_am <- buildingMatrix2( limited_edges )
  pkg_am <- CommunityInference::adjacencyMatFromDF( limited_edges )
} )



# molecule_amat <- usingLists( molecule_edges )
