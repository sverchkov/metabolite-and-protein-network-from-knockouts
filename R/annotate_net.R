#' Find surprising interactions
#' @author Yuriy Sverchkov
#' 
#' General approach:
#' * Loop through filtered similarities
#' * For each similarity, find interaction in DB to account for it
#' * If found, annotate, if not found, mark as surprising

########
# Init #
########
library( "KEGGREST" )
library( "futile.logger" )

flog.threshold(TRACE)

#############
# Load data #
#############
molecule_net <- readRDS( "processed-data/filtered-cosine-similarities.rds" )
molecule_ids <- readRDS( "processed-data/molecule-ids.rds" )

#################
# Load mappings #
#################
uniprot2kegg <- keggConv( "hsa", "uniprot" )

molecule_subset <- union( distinct( molecule_net, id = a ), distinct( molecule_net, id = b ) ) %>%
  left_join( molecule_ids, by = "id" )

proteins <- filter( molecule_subset, .data$`Molecule Type` == "Protein" )

kegg_marked_molecules <- filter( molecule_subset, !is.na( .data$`KEGG ID` ), "" != .data$`KEGG ID` )

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

hsa2pathway <- keggLink( "pathway", "hsa" )
hsa_pathway_map <- tibble( hsa_code = names( hsa2pathway ), pathway = hsa2pathway )

hsa2enzyme <- keggLink( "enzyme", "hsa" )
hsa_enzyme_map <- tibble( hsa_code = names( hsa2enzyme ), enzyme = hsa2enzyme )

enzyme2reaction <- keggLink( "reaction", "enzyme" )
enzyme_reaction_map <- tibble( enzyme = names( enzyme2reaction ), reaction = enzyme2reaction )

compound2reaction <- keggLink( "reaction", "compound" )
compound_reaction_map <- tibble( compound = names( compound2reaction ), reaction = compound2reaction )

edges_with_types <- molecule_net %>%
  left_join( distinct( molecule_subset, a = id, a_type = `Molecule Type` ), by = "a" ) %>%
  left_join( distinct( molecule_subset, b = id, b_type = `Molecule Type` ), by = "b" )

edges_pp <- filter( edges_with_types, "Protein" == a_type, "Protein" == b_type )
edges_pc <- filter( edges_with_types, "Protein" == a_type, "Protein" != b_type )
edges_cp <- filter( edges_with_types, "Protein" != a_type, "Protein" == b_type )
edges_cc <- filter( edges_with_types, "Protein" != a_type, "Protein" != b_type )

# Protein-protein edges are annotated by common pathways
found_pp <- edges_pp %>%
  inner_join( distinct( protein_map, a = id, a_hsa = hsa_code ), by = "a" ) %>%
  inner_join( distinct( protein_map, b = id, b_hsa = hsa_code ), by = "b" )

annotations_pp <- found_pp %>%
  inner_join( distinct( hsa_pathway_map, a_hsa = hsa_code, pathway ), by = "a_hsa" ) %>%
  inner_join( distinct( hsa_pathway_map, b_hsa = hsa_code, pathway ), by = c( "b_hsa", "pathway" ) ) %>%
  group_by( a, b ) %>%
  summarize( annotation = paste( pathway, collapse = "; " ) )

found_cc <- edges_cc %>%
  inner_join( distinct( compound_map, a = id, a_cpd = compound ), by = "a" ) %>%
  inner_join( distinct( compound_map, b = id, b_cpd = compound ), by = "b" )

annotations_cc <- found_cc %>%
  inner_join( distinct( compound_reaction_map, a_cpd = compound, reaction ), by = "a_cpd" ) %>%
  inner_join( distinct( compound_reaction_map, b_cpd = compound, reaction ), by = c( "b_cpd", "reaction" ) ) %>%
  group_by( a, b ) %>%
  summarize( annotation = paste( reaction, collapse = "; " ) )

found_cp <- edges_cp %>%
  inner_join( distinct( compound_map, a = id, a_cpd = compound ), by = "a" ) %>%
  inner_join( distinct( protein_map, b = id, b_hsa = hsa_code ), by = "b" )

annotations_cp <- found_cp %>%
  inner_join( distinct( compound_reaction_map, a_cpd = compound, reaction ), by = "a_cpd" ) %>%
  inner_join( distinct( hsa_enzyme_map, b_hsa = hsa_code, enzyme ), by = c( "b_hsa" ) ) %>%
  inner_join( distinct( enzyme_reaction_map, enzyme, reaction ), by = c( "enzyme", "reaction" ) ) %>%
  group_by( a, b ) %>%
  summarize( annotation = paste( reaction, collapse = "; " ) )

if ( nrow( edges_pc ) == 0 ) {
  flog.info( "Skipping protin-comound edges (none found because of molecule ordering" )
} else {
  stop( flog.fatal( "We weren't expecting protein-compound edges but we found some") )
}

annotations <- bind_rows(
  distinct( found_cc, a, b ),
  distinct( found_cp, a, b ),
  distinct( found_pp, a, b ) ) %>%
  left_join( bind_rows( annotations_cc, annotations_cp, annotations_pp ), by = c( "a", "b" ) ) %>%
  mutate( annotation = if_else( is.na( annotation ), "surprising", annotation ) )

annotated_net <- left_join( molecule_net, annotations, by = c( "a", "b" ) ) %>%
  mutate( annotation = if_else( is.na( annotation ), "unknown molecule", annotation ) )

saveRDS( annotated_net, "processed-data/annotated-filtered-cosims.rds" )
