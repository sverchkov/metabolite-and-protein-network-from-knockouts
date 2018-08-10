#' Make p-value-filtered cytoscape graph from similarities
#' 
#' @author Yuriy Sverchkov

# Libraries
library("dplyr")

# Load filtered similarities
filtered_similarity <- readRDS( "processed-data/filtered-cosine-similarities.rds" )

similarity_nodes <-
  left_join( union( similarity_edges %>% select( id = a ), similarity_edges %>% select( id = b ) ),
             spread_table %>% select( id, Name = `Molecule Name`, Type = `Molecule Type` ),
             by = "id" ) %>%
  mutate( id = paste0( "M", id ) )

# Cytoscape doesn't like integer IDs for some reason
similarity_edges <- filtered_similarity %>%
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

