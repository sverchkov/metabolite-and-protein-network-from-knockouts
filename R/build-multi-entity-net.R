# Make a multi-entity network based on KEGG and STRING

library("dplyr")
library( "KEGGREST" )
library( "futile.logger" )

flog.threshold(TRACE)

# Load node IDs
molecule_ids <- readRDS("processed-data/molecule-ids.rds")
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

# Load KEGG mappings
uniprot2kegg <- keggConv( "hsa", "uniprot" )

hsa2pathway <- keggLink( "pathway", "hsa" )
hsa_pathway_map <- tibble( hsa_code = names( hsa2pathway ), pathway = hsa2pathway )

hsa2enzyme <- keggLink( "enzyme", "hsa" )
hsa_enzyme_map <- tibble( hsa_code = names( hsa2enzyme ), enzyme = hsa2enzyme )

enzyme2reaction <- keggLink( "reaction", "enzyme" )
enzyme_reaction_map <- tibble( enzyme = names( enzyme2reaction ), reaction = enzyme2reaction )

compound2reaction <- keggLink( "reaction", "compound" )
compound_reaction_map <- tibble( compound = names( compound2reaction ), reaction = compound2reaction )

# Build protein and compound map
proteins <- filter( molecule_ids, .data$`Molecule Type` == "Protein" )

kegg_marked_molecules <- filter( molecule_ids, !is.na( .data$`KEGG ID` ), "" != .data$`KEGG ID` )

protein_map <- bind_rows(
  mapply(
    function ( molecule_id, uniprot_id ){
      tibble( protein_id = molecule_id, uniprot = paste0( "up:", uniprot_id ) )
    },
    proteins %>% pull( `Molecule ID` ) , strsplit( proteins %>% pull( `Molecule ID` ), "; *" ),
    SIMPLIFY = F )
) %>%
  mutate( hsa_code = uniprot2kegg[uniprot] )

compound_map <- bind_rows(
  mapply(
    function( int_id, kegg ){
      tibble( molecule_id = int_id, compound = paste0( "cpd:", kegg ) )
    },
    kegg_marked_molecules %>% pull( `Molecule ID` ), strsplit( kegg_marked_molecules %>% pull( `KEGG ID` ), "; *" ),
    SIMPLIFY = F )
)

# Pull in protein-pathway connections
protein_pathway <- distinct( protein_map, protein_id, hsa_code ) %>%
  inner_join( distinct( hsa_pathway_map, hsa_code, pathway ), by = "hsa_code" )

# Identify connector pathways
connector_pathways <- protein_pathway %>%
  distinct( pathway, protein_id ) %>%
  group_by( pathway ) %>%
  summarize( count = n() ) %>%
  ungroup() %>%
  filter( count > 1 ) %>%
  select( -count )

# Only keep connector pathways
protein_map_net <- protein_pathway %>% inner_join( connector_pathways, by = "pathway" )

## Optional: Print Protein-Pathway Map
if( askYesNo("Make a Cytoscape network for protein-pathways only?") ){
  RCy3::createNetworkFromDataFrames(
    edges = distinct( protein_map_net, source = protein_id, target = pathway) %>% as.data.frame(),
    title = paste( "Protein-Pathway KEGG Network", date() ),
    collection = "H3K Networks" )
}

warning("We need to figure out if 'enzymes' are really the correct entities to connect proteins to reactions!")
# Pull in protein-enzyme connections
protein_enzyme_edges <- inner_join(
  distinct( protein_map, protein_id, hsa_code ),
  distinct( hsa_enzyme_map, hsa_code, enzyme ),
  by = "hsa_code" ) %>%
  distinct( protein_id, enzyme ) %>%
  mutate( interaction = "protein-enzyme" )
  
# Pull in compound-reaction connections
compound_reaction_edges <- inner_join(
  distinct( compound_map, molecule_id, compound ),
  distinct( compound_reaction_map, compound, reaction ),
  by = "compound" ) %>%
  distinct( molecule_id, reaction ) %>%
  mutate( interaction = "compound-reaction" )
  
# Pull in enzyme-reaction connections
enzyme_reaction_edges <- enzyme_reaction_map %>%
  mutate( interaction = "enzyme-reaction" )

# Only keep connecting enzymes
connecting_enzymes <- union_all(
  protein_enzyme_edges %>% select( enzyme, other = protein_id, interaction ),
  distinct(
    inner_join(
      enzyme_reaction_edges %>% select( -interaction ),
      compound_reaction_edges %>% select( -interaction ),
      by = "reaction" ),
    enzyme, other = molecule_id ) ) %>%
  group_by( enzyme ) %>%
  summarize( count = n() ) %>%
  ungroup() %>%
  filter( count > 1 ) %>%
  distinct( enzyme )

# Only keep connecting reactions
connecting_reactions <- union_all(
  compound_reaction_edges %>% select( reaction, other = molecule_id, interaction ),
  enzyme_reaction_edges %>% inner_join( connecting_enzymes, by = "enzyme" ) %>%
    select( reaction, other = enzyme, interaction ) ) %>%
  group_by( reaction ) %>%
  summarize( count = n() ) %>%
  ungroup() %>%
  fiter( count > 1 ) %>%
  distinct( reaction )

## TODO: factor redundant code into functions (esp. the "keep connecting" logic)

# Pull in  STRING protein-protein connections
found_string_pp <- readRDS("temporary-data/found_string_pp.rds") %>% mutate( interaction = "protein-protein" )

#####################
# Network synthesis #
#####################

# Put together a library of entities
net_nodes <- union_all(
  tibble( id = kos, full_name = kos, type = "Knockout" ),
  molecule_ids %>% select( id = `Molecule ID`, full_name = `Molecule Name`, type = `Molecule Type` ),
  connector_pathways %>% transmute( id = pathway, full_name = pathway, type = "Pathway" ),
  connecting_reactions %>% transmute( id = reaction, full_name = reaction, type = "Reaction" ),
  connecting_enzymes %>% transmute( id = enzyme, full_name = enzyme, type = "Enzyme" ) )

# Put together all edges
net_edges <- union_all(
  transmute( ko_proteins, source = Knockout, target = `Molecule ID`, interaction = "knockout-protein" ),
  select( found_string_pp, source = a, target = b, interaction ),
  transmute( protein_pathway, source = protein_id, target = pathway, interaction = "protein-pathway" ),
  select( protein_enzyme_edges, source = protein_id, target = enzyme, interaction ),
  ## TODO other parts
)

# Save
saveRDS( net_nodes, "processed-data/augmented-bg-net-nodes.rds" )
saveRDS( net_edges, "processed-data/augmented-bg-net-edges.rds" )