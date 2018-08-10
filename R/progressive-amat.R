progressiveADJM <- function( edges ){
  n <- length( unique( c( edges$a, edges$b ) ) )
  rows <- list()
  for( i in 1:n ){
    mask_a <- edges$a == i
    mask_b <- edges$b == i
    jays <- c( edges$b[mask_a], edges$a[mask_b] )
    w <- numeric(n)
    w[jays] <- c( edges$weight[mask_a], edges$weight[mask_b] )
    rows[[i]] <- w
    
    if( ( i %% 100 ) == 0 ){
      saveRDS(rows, "processed-data/partial_amat.rds" )
      flog.info("Wrote %s/%s rows", i, n )
    }
  }
  
  do.call( rbind, rows )
}

molecule_adj <- progressiveADJM( molecule_edges )

