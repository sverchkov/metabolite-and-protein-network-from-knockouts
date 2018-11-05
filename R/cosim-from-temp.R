# Get cosine similarities from (potentially partial) dot product files

flog.trace( "Reading temp files into one table..." )
dot_products <- bind_rows(
  lapply(
    1:n_molecules,
    function( a ) {
      f <- sprintf( filepattern, a )
      if ( file.exists( f ) )
        readRDS( f )
      else
        NULL
    } ) )
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