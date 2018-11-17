# Visualize annotated net

# Load net
if ( !exists( "annotated_net") ) {
  annotated_net <- readRDS( "processed-data/annotated-filtered-cosims.rds" )
} else flog.info( "Using annotated_net object that was in environment" )

# Let's just deal with the molecules in the data
nodes <- sort( unique( c( annotated_net$a, annotated_net$b ) ) )

# Load molecule ID table
molecule_ids <- readRDS("processed-data/molecule-ids.rds" )

# Make nodes df
cy_nodes <- molecule_ids %>% filter( `Molecule ID` %in% nodes ) %>%
  select( -id ) %>%
  select( id = `Molecule ID`, full_name = `Molecule Name`, type = `Molecule Type` ) %>%
  as.data.frame()

# Make edges df
cy_edges <- annotated_net %>%
  transmute( source = a, target = b, similarity, annotation ) %>%
  as.data.frame()

# Create cytoscape network with RCy3
net_id <- RCy3::createNetworkFromDataFrames(
  nodes = cy_nodes,
  edges = cy_edges,
  title = paste( "Annotated Cosine Similarity Network", date() ),
  collection = "H3K Networks" )

# Visual style
style_name <- "H3K Annotated Molecule Network"

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
