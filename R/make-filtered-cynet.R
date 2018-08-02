#' Make p-value-filtered cytoscape graph from similarities
#' 
#' @author Yuriy Sverchkov

# Libraries
library("dplyr")

# Load necessary data
similarities <- readRDS( "processed-data/cosine-similarities.rds" )
spread_table <- readRDS( "processed-data/spread-table-for-similarities.rds" )

# notes: see https://stats.stackexchange.com/questions/85916/distribution-of-scalar-products-of-two-random-unit-vectors-in-d-dimensions
# for null distribution of cosine similarities.
# Having a null allows us to obtain p-values
# For network we probably need to do multiple testing correction?

# We are assuming that the non-KO columns in spread_table are Molecule ID, Molecule Name, Molecule Type, and id
n_kos <- ncol( spread_tabe ) - 4

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

