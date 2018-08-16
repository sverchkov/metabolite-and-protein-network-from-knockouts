# Visualize annotated net

# Load net
annotated_net <- readRDS( "processed-data/kegg-bg-net.rds" )

# Let's just deal with the molecules in the data
nodes <- sort( unique( c( annotated_net$a, annotated_net$b ) ) )

# Load molecule ID table
molecule_id_table <- readRDS("processed-data/molecule-ids.rds" )

# Get molecule names and types
indeces <- match( nodes, molecule_id_table$id )
molecule_names <- molecule_id_table$`Molecule Name`[indeces]
molecule_types <- molecule_id_table$`Molecule Type`[indeces]

# Make nodes df
cy_nodes <- data.frame(
  id = paste0( "M", nodes ),
  full_name = molecule_names,
  type = molecule_types,
  stringsAsFactors = F )

# Cytoscape doesn't like integer IDs for some reason
cy_edges <- annotated_net %>%
  transmute( source = paste0( "M", a ), target = paste0( "M", b ), annotation ) %>%
  as.data.frame()

# Create cytoscape network with RCy3
net_id <- RCy3::createNetworkFromDataFrames(
  nodes = cy_nodes,
  edges = cy_edges,
  title = paste( "Annotated KEGG BG Network", date() ),
  collection = "H3K Networks" )

# Visual style
style_name <- "H3K Molecule Network"

RCy3::createVisualStyle(
  style.name = style_name,
  defaults = list( NODE_SHAPE="round rectangle", EDGE_COLOR = "#000000"),
  mappings = list(
    RCy3::mapVisualProperty( 'node label', 'full_name', 'p' ),
    RCy3::mapVisualProperty( 'node shape', 'type', 'd',
                             c("Protein", "Metabolite", "Lipid"),
                             c("round rectangle", "circle", "triangle") ),
    RCy3::mapVisualProperty( 'edge unselected paint', 'annotation', 'd',
                             c( "surprising", "unknown molecule" ),
                             c( "#008000", "#909090" ) )
  ) )

RCy3::setVisualStyle( style_name )
