# This file loads the relevant STRING relations
#

library(dplyr)

u2s_map <-
  read.table( "/Volumes/LaCie/yuriy/string-db/full_uniprot_2_string.04_2015.tsv",
              header = FALSE,
              sep = "\t",
              check.names = FALSE,
              strip.white = TRUE,
              comment.char = "#",
              stringsAsFactors = FALSE )

colnames( u2s_map ) <- c( "species", "uniprot_both", "string_id", "identity", "bit_score" )

u2s_map <- u2s_map %>%
  mutate(
    uniprot_acc = sapply( strsplit( uniprot_both, split = "|", fixed = TRUE ), function(x) x[1] ) )

molecule_ids <- readRDS( "processed-data/molecule-ids.rds" )

proteins <- molecule_ids %>%
  filter( `Molecule Type` == "Protein" ) %>%
  pull( `Molecule ID` )

my_uniprot <- bind_rows(
  mapply(
    function ( one, many ){
      tibble( `Molecule ID` = one, uniprot_id = many )
    },
    proteins,
    strsplit( proteins, ";", fixed = TRUE ),
    SIMPLIFY = FALSE
    )
  )

# Our uniprot IDs are accession numbers

my_string_map <-
  left_join( my_uniprot, u2s_map, by = c("uniprot_id" = "uniprot_acc" ) ) %>%
  mutate( string_long = paste0( species, ".", string_id ) )

# (Optional, for diagnostics)

# How many matches are there for each molecule

mapping_statistics <- my_string_map %>%
  group_by( `Molecule ID` ) %>%
  distinct( string_id ) %>%
  mutate( count_me = if_else( is.na( string_id ), 0L, 1L ) ) %>%
  summarise( n_matches = sum( count_me ) )

# How many matches are there for each string id

inverse_mapping_statistics <- my_string_map %>%
  group_by( string_id ) %>%
  distinct( `Molecule ID` ) %>%
  mutate( count_me = if_else( is.na( `Molecule ID` ), 0L, 1L ) ) %>%
  summarize( n_matches = sum( count_me ) )

# Load STRING interaction table
string_interactions <-
  read.table( "/Volumes/LaCie/yuriy/string-db/human.protein.actions.txt",
              header = T,
              sep = "\t",
              check.names = F,
              strip.white = T,
              stringsAsFactors = F )