#' Make p-value-filtered cytoscape graph from similarities
#' 
#' @author Yuriy Sverchkov

# Libraries
library("dplyr")

# Load filtered similarities
filtered_similarity <- readRDS( "processed-data/filtered-cosine-similarities.rds" )

# Load grouped nodes
molecule_communities <- readRDS( "processed-data/molecule-communities-aslp.rds" )

# Load molecule ID table
molecule_id_table <- readRDS("processed-data/molecule-ids.rds" )

# Get molecule names and types
indeces <- match( molecule_communities$node, molecule_id_table$id )
molecule_names <- molecule_id_table$`Molecule Name`[indeces]
molecule_types <- molecule_id_table$`Molecule Type`[indeces]

# Get list of unique label values
label_array <- paste0( "L", sort( unique( molecule_communities$label ) ) )
n <- length( label_array )
color_cycle <- ceiling( n / 12 )
fill_array <- ( rep( RColorBrewer::brewer.pal( 12, "Paired" ), color_cycle ) )[1:n]
paint_array <- ( mapply( rep, RColorBrewer::brewer.pal( 12, "Paired" ), color_cycle ) )[1:n]

# Make nodes df
similarity_nodes <- data.frame(
  id = paste0( "M", molecule_communities$node ),
  full_name = molecule_names,
  type = molecule_types,
  label = paste0( "L", molecule_communities$label ),
  stringsAsFactors = F )

# Cytoscape doesn't like integer IDs for some reason
similarity_edges <- filtered_similarity %>%
  transmute( source = paste0( "M", a ), target = paste0( "M", b ), similarity ) %>%
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
    RCy3::mapVisualProperty( 'node label', 'full_name', 'p' ),
    RCy3::mapVisualProperty( 'Node Border Paint', 'label', 'd', label_array, paint_array ),
    RCy3::mapVisualProperty( 'node fill color', 'label', 'd', label_array, fill_array ),
    RCy3::mapVisualProperty( 'node shape', 'type', 'd',
                             c("Protein", "Metabolite", "Lipid"),
                             c("round rectangle", "circle", "triangle") ),
    RCy3::mapVisualProperty( 'edge unselected paint', 'similarity', 'c',
                             c( -1, -0.85, 0.85, 1 ),
                             c( "#FF0000", "#F01010", "#908080", "#808090", "#1010F0", "#0000FF" ) )
  ) )

RCy3::setVisualStyle( style_name )
