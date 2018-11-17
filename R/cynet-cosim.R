#' Make p-value-filtered cytoscape graph from similarities
#' 
#' 2018, November 16
#' 
#' @author Yuriy Sverchkov

# Libraries
library("dplyr")

# Load filtered similarities
filtered_similarity <- readRDS( "processed-data/filtered-cosine-similarities.rds" )

# Load molecule ID table
molecule_ids <- readRDS("processed-data/molecule-ids.rds" )

# Get molecule names and types
#indeces <- match( molecule_communities$node, molecule_id_table$id )
molecule_names <- molecule_ids %>% pull( `Molecule Name` )
molecule_types <- molecule_ids %>% pull( `Molecule Type` )
molecule_id_array <- molecule_ids %>% pull( `Molecule ID` )

# Make nodes df
similarity_nodes <- molecule_ids %>%
  select( -id ) %>%
  select(
    id = `Molecule ID`,
    type = `Molecule Type`,
    name = `Molecule Name` ) %>%
  as.data.frame()

# Cytoscape doesn't like integer IDs for some reason
similarity_edges <- filtered_similarity %>%
  transmute( source = a, target = b, similarity ) %>%
  as.data.frame()

# Create cytoscape network with RCy3
net_id <- RCy3::createNetworkFromDataFrames(
  nodes = similarity_nodes,
  edges = similarity_edges,
  title = paste( "Cosine Similarity Network", date() ),
  collection = "H3K Networks" )

# Visual style
style_name <- "H3K Molecule Network"

RCy3::createVisualStyle(
  style.name = style_name,
  defaults = list(NODE_SHAPE="round rectangle"),
  mappings = list(
    RCy3::mapVisualProperty( 'node label', 'name', 'p' ),
    RCy3::mapVisualProperty( 'node shape', 'type', 'd',
                             c("Protein", "Metabolite", "Lipid"),
                             c("round rectangle", "circle", "triangle") ),
    RCy3::mapVisualProperty( 'edge unselected paint', 'similarity', 'c',
                             c( -1, -0.85, 0.85, 1 ),
                             c( "#FF0000", "#F01010", "#908080", "#808090", "#1010F0", "#0000FF" ) )
  ) )

RCy3::setVisualStyle( style_name )
