# What I tried to do for presentation July 30, 2018
library("dplyr")
library("futile.logger")

flog.threshold(TRACE)

# Load table
tested_omics_table <- readRDS( "processed-data/tested-omics-table.rds" )

# Make table shape
spread_table <- tested_omics_table %>%
  mutate( t = `Mean log2FC`/`SD of log2FC` ) %>%
  tidyr::spread( key = Knockout, value = t, fill = 0 ) %>%
  select( -`Mean log2FC`, -`SD of log2FC`, -`p-Value` ) %>%
  group_by( `Molecule ID`, `Molecule Name`, `Molecule Type` ) %>%
  summarize_all( sum ) %>%
  ungroup() %>%
  arrange( desc(`Molecule Type`), `Molecule ID` ) %>%
  mutate( id = dplyr::row_number() )

# Figure out how many proteins there are
n_proteins <- ( spread_table %>% filter( "Protein" == `Molecule Type` ) %>% summarize( count=n() ) )$count

# Make matrix-shaped df
matrix_df <- spread_table %>% select( -`Molecule ID`, -`Molecule Name`, -`Molecule Type`, -id )

# We use the number of dimensions (KOs) later
n_kos <- ncol( matrix_df )

# Get number of molecules for convenience
n_molecules <- nrow( matrix_df )

# Estimate the number of dot products we will compute
n_dots <- n_molecules * (n_molecules + 1)/2

# Get dot products, progressive computation
dot_products <- NULL
result_rows <- 0
added_rows <- 0

for ( a in 1:n_molecules ){
  v1 <- as.vector( matrix_df[a,] )
  block <- bind_rows( Map( function ( b ) {
      prd <- sum( as.vector( matrix_df[b,] ) * v1 )
      tibble( a = a, b = b, product = prd )
    }, 1:a ) )
  
  dot_products <- rbind( dot_products, block )
  added_rows <- added_rows + nrow( block )
  if ( added_rows > 10000 ){
    result_rows <- result_rows + added_rows
    added_rows <- 0
    flog.trace("Computed %s of %s dot products (%s%%)", result_rows, n_dots, 100*result_rows/n_dots )
    saveRDS( dot_products, "processed-data/dot-products.rds" )
  }
}

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
saveRDS( spread_table, "processed-data/spread-table-for-similarities.rds" )

# notes: see https://stats.stackexchange.com/questions/85916/distribution-of-scalar-products-of-two-random-unit-vectors-in-d-dimensions
# for null distribution of cosine similarities.
# Having a null allows us to obtain p-values
# For network we probably need to do multiple testing correction?

# Bonferroni correction:
alpha = 0.05 / nrow( similarities )

# Cosine similarity p-Value:
cutoff <- 1 - qbeta( alpha, (n_kos-1)/2, (n_kos-1)/2 ) * 2

similarity_edges <- similarities %>% filter( abs(similarity ) > cutoff )

similarity_nodes <-
  left_join( union( similarity_edges %>% select( id = a ), similarity_edges %>% select( id = b ) ),
             spread_table %>% select( id, Name = `Molecule Name`, Type = `Molecule Type` ),
             by = "id" ) %>%
             mutate( id = paste0( "M", id ) )

# Cytoscape doesn't like integer IDs for some reason
similarity_edges <- similarity_edges %>%
  transmute( source = paste0( "M", a ), target = paste0( "M", b ), similarity )

# Create cytoscape network with RCy3
net_id <- RCy3::createNetworkFromDataFrames(
  nodes = as.data.frame( similarity_nodes ),
  edges = as.data.frame( similarity_edges ),
  title = paste( "Cosine Similarity Network", date() ),
  collection = "H3K Networks" )

# Visual style
RCy3::setNodeColorMapping( # Color nodes by molecule type
  table.column = "Type",
  table.column.values = c("Protein", "Metabolite", "Lipid"),
  colors = c("#8888AA", "#AA4444", "#AAAA44"),
  mapping.type = "d",
  network = net_id )

RCy3::setEdgeColorMapping( # Map edge color to similarity
  table.column = "similarity",
  table.column.values = c( -1, -cutoff, cutoff, 1 ),
  colors = c( "#FF0000", "#F01010", "#908080", "#808090", "#1010F0", "#0000FF" ),
  mapping.type = "c",
  network = net_id )