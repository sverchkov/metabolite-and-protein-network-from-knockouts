# Write a DAG based on missing output files

to_write <- character( length = 0 )

for ( n in 1:n_molecules ){
  f <- sprintf( filepattern, n )
  if ( !file.exists( f ) ){
    # Write to file
    to_write <- c( to_write,
                   sprintf( "JOB DOT-%d assume-r-dot.condor", n ),
                   sprintf( "VARS DOT-%d n=\"%d\"", n, n ),
                   sprintf( "VARS DOT-%d result=\"dot-products/dots-%d.rds\"", n, n ) )
  }
}

writeLines( text = to_write, con = "condor-bundles/dot-product/assume-r-dots.dag" )