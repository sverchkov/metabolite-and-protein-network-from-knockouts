transitively_reduce <- function ( adj ) {
  adj_prev <- diag( nrow( adj ) )
  diag( adj ) <- 0
  while( !isTRUE( all.equal( adj, adj_prev ) ) ) {
    adj_prev <- adj
    for ( a in 1:nrow( adj ) )
      for( b in 1:ncol( adj ) ) if ( a != b )
        for( c in 1:ncol( adj ) ) if ( a != c && b != c )
          if ( adj[a,b] == adj[a,c]*adj[c,b] )
            adj[a,b] <- 0
  }
  return ( adj )
}