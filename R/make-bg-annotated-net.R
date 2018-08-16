# Make a network of the molecules we have based on KEGG

library( "KEGGREST" )
library( "futile.logger" )

flog.threshold(TRACE)

# Load node IDs
molecule_ids <- readRDS("processed-data/molecule-ids.rds")

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
    function ( int_id, uniprot_id ){
      tibble( id = int_id, uniprot = paste0( "up:", uniprot_id ) )
    },
    proteins$id, strsplit( proteins$`Protein ID`, "; *" ),
    SIMPLIFY = F )
) %>%
  mutate( hsa_code = uniprot2kegg[uniprot] )

compound_map <- bind_rows(
  mapply(
    function( int_id, kegg ){
      tibble( id = int_id, compound = paste0( "cpd:", kegg ) )
    },
    kegg_marked_molecules$id, strsplit( kegg_marked_molecules$`KEGG ID`, "; *" ),
    SIMPLIFY = F )
)

# Build net
annotations_pp <- distinct( protein_map, a = id, a_hsa = hsa_code ) %>%
  inner_join( distinct( hsa_pathway_map, a_hsa = hsa_code, pathway ), by = "a_hsa" ) %>%
  inner_join( distinct( hsa_pathway_map, b_hsa = hsa_code, pathway ), by = "pathway" ) %>%
  inner_join( distinct( protein_map, b = id, b_hsa = hsa_code ), by = "b_hsa" ) %>%
  filter( a < b ) %>%
  group_by( a, b ) %>%
  summarize( annotation = paste( pathway, collapse = "; " ) ) %>%
  ungroup()

annotations_cc <- distinct( compound_map, a = id, a_cpd = compound ) %>%
  inner_join( distinct( compound_reaction_map, a_cpd = compound, reaction ), by = "a_cpd" ) %>%
  inner_join( distinct( compound_reaction_map, b_cpd = compound, reaction ), by = "reaction" ) %>%
  inner_join( distinct( compound_map, b = id, b_cpd = compound ), by = "b_cpd" ) %>%
  filter( a < b ) %>%
  group_by( a, b ) %>%
  summarize( annotation = paste( reaction, collapse = "; " ) ) %>%
  ungroup()

annotations_cp <- distinct( compound_map, a = id, a_cpd = compound ) %>%
  inner_join( distinct( compound_reaction_map, a_cpd = compound, reaction ), by = "a_cpd" ) %>%
  inner_join( distinct( enzyme_reaction_map, enzyme, reaction ), by = "reaction" ) %>%
  inner_join( distinct( hsa_enzyme_map, b_hsa = hsa_code, enzyme ), by = "enzyme" ) %>%
  inner_join( distinct( protein_map, b = id, b_hsa = hsa_code ), by = "b_hsa" ) %>%
  filter( a < b ) %>%
  group_by( a, b ) %>%
  summarize( annotation = paste( reaction, collapse = "; " ) ) %>%
  ungroup()

annotated_net <- bind_rows( annotations_pp, annotations_cc, annotations_cp )

saveRDS( annotated_net, "processed-data/kegg-bg-net.rds" )