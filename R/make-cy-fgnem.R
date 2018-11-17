# Draw FGNEM in cytoscape

source( "R/transitively_reduce.R" )

stopifnot( exists("fgnem") )

library( futile.logger )
flog.threshold( TRACE )

adj <- fgnem$acc
the_cols <- colnames( adj )

adj <- transitively_reduce( adj )

edges_df <- do.call( rbind, lapply(
  rownames( adj ),
  function ( the_row ) {
    pos <- the_cols[ adj[ the_row, ] ==  1 ]
    neg <- the_cols[ adj[ the_row, ] == -1 ]
    if( length( pos ) + length( neg ) <= 0 )
      NULL
    else
      data.frame( source = the_row,
                  target = c( pos, neg ),
                  interaction = c( rep( "positive", length(pos) ), rep( "negative", length(neg) ) ),
                  stringsAsFactors = F )
  } ) )

RCy3::createNetworkFromDataFrames(
  nodes = NULL,
  edges = edges_df,
  title = paste("FGNEM",date()),
  collection = "H3K Networks" )

# Let's try to merge some nodes
for ( the_col in colnames( adj ) ){
  for( other_col in colnames( adj ) ) if( the_col != other_col ){
    if( all( adj[,the_col] == adj[,other_col] ) ) flog.info( "%s = %s", the_col, other_col )
  }
} # We found out that no nodes are identical