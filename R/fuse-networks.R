## Fuse Networks

# Load things
library( "dplyr" )

annotated_molecule_edges <- readRDS( "processed-data/annotated_molecule_edges.rds" )
fgnem_ko_edges <- readRDS( "processed-data/fgnem_ko_edges.rds" )
molecule_ids <- readRDS( "processed-data/molecule-ids.rds" )
# Assume kos loaded

# Get the KO-protein map
ko_proteins <- bind_rows( lapply(
  kos,
  function( ko ) {
    molecule_ids %>%
      filter( grepl( paste0( "(^|;)", ko, "(;|$)") , `Protein Gene` ) ) %>%
      transmute( `Molecule ID`, `Protein Gene`, Knockout = ko )
  }
))

# Make ko-protein edges
ko_protein_links <-
  ko_proteins %>%
  transmute( source = Knockout, target = `Molecule ID`, interaction = "gene product" )

# Put all edges in one DF
cy_edges <- union_all( annotated_molecule_edges, fgnem_ko_edges ) %>%
  union_all( ko_protein_links ) %>%
  as.data.frame()

# Make nodes df
molecule_nodes <- union(
  distinct( annotated_molecule_edges, `Molecule ID` = source ),
  distinct( annotated_molecule_edges, `Molecule ID` = target ) )

cy_nodes <- molecule_nodes %>% left_joiin( molecule_ids ) %>%
  select( -id ) %>%
  select( id = `Molecule ID`, full_name = `Molecule Name`, type = `Molecule Type` ) %>%
  union( tibble( id = kos, full_name = id, type = `Knockout` ) ) %>%
  as.data.frame()

# Save objects
destination <- "results/fused-network"
dir.create( destination, recursive = T )

saveRDS( cy_edges, paste0( destination, "edges.rds" ) )
saveRDS( cy_nodes, paste0( destination, "nodes.rds" ) )

# Run Cytoscape
if ( askYesNo( "Make Cytoscape network?" ) ){
  # Create cytoscape network with RCy3
  net_id <- RCy3::createNetworkFromDataFrames(
    nodes = cy_nodes,
    edges = cy_edges,
    title = paste( "Annotated Cosine Similarity Network", date() ),
    collection = "H3K Networks" )
  
  # Visual style
  # TODO
}