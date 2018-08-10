#' Make p-value-filtered cytoscape graph from similarities
#' 
#' @author Yuriy Sverchkov

# Libraries
library("dplyr")

# Load necessary data
similarities <- readRDS( "processed-data/cosine-similarities.rds" )
spread_table <- readRDS( "processed-data/spread-table-for-similarities.rds" )

# Node data
cy_nodes <-
  spread_table %>%
  select( id, Name = `Molecule Name`, Type = `Molecule Type` ) %>%
  mutate( id = paste0( "M", id ) )

# Cytoscape doesn't like integer IDs for some reason
similarity_edges <- similarities %>%
  transmute( source = paste0( "M", a ), target = paste0( "M", b ), similarity )

# Create cytoscape network with RCy3
net_id <- RCy3::createNetworkFromDataFrames(
  nodes = as.data.frame( cy_nodes ),
  edges = as.data.frame( similarity_edges ),
  title = paste( "Full Cosine Similarity Network", date() ),
  collection = "H3K Networks" )

# Visual style
RCy3::setNodeColorMapping( # Color nodes by molecule type
  table.column = "Type",
  table.column.values = c("Protein", "Metabolite", "Lipid"),
  colors = c("#8888AA", "#AA4444", "#AAAA44"),
  mapping.type = "d",
  network = net_id )

max_line_width = 3
min_line_width = 0.5

RCy3::setEdgeLineWidthMapping( # Map edge width to similarity
  table.column = "similarity",
  table.column.values = c(-1, 0, 1),
  widths = c( max_line_width, min_line_width, max_line_width ),
  mapping.type = "c",
  network = net_id )

RCy3::setEdgeColorMapping( # Map edge color to similarity
  table.column = "similarity",
  table.column.values = c( -1, 0, 1 ),
  colors = c( "#F01010", "#808080", "#0000FF" ),
  mapping.type = "c",
  network = net_id )

