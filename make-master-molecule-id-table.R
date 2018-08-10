stop("Are you sure you want to run this? It's a one-time script for making a master molecule ID table.")

library(dplyr)

spead_table <- readRDS( "processed-data/spread-table-for-similarities.rds" )
molecule_id_table <- spead_table[c("Molecule ID", "Molecule Name", "Molecule Type", "id")]
unimputed_metabolomics <- read.csv( "raw-data/20180322_H3K_Batches1and2_Metabolomics_Unimputed.csv",
                                    check.names = F, stringsAsFactors = F )
molecule_id_table_sav <- molecule_id_table # For in case we screw up
molecule_id_table <- full_join( molecule_id_table,
                                transmute( unimputed_metabolomics,
                                           Metabolites, `KEGG IDs`, `Molecule Type` = "Metabolite" ),
                                by = c( "Molecule ID" = "Metabolites", "Molecule Type" ) )
molecule_id_table_sav <- molecule_id_table # For in case we screw up
unimputed_proteomics <- read.csv( "raw-data/20180322_H3K_Batches1and2_Proteomics_Unimputed.csv",
                                  check.names = F, stringsAsFactors = F )
molecule_id_table <- full_join( molecule_id_table,
                                transmute( unimputed_proteomics,
                                           `Majority protein IDs`, `Molecule Name` = `Protein names`, `Gene names`,
                                           `Molecule ID` = `Fasta headers`, `Molecule Type` = "Protein" ),
                                by = c( "Molecule Name", "Molecule ID", "Molecule Type" ) )
molecule_id_table_sav <- molecule_id_table # For in case we screw up
lipidomics <- read.csv( "raw-data/20180322_H3K_Batches1and2_Lipidomics.csv",
                        check.names = F, stringsAsFactors = F )[1:2]
molecule_id_table <- full_join( molecule_id_table_sav,
                                transmute( lipidomics, `Molecule ID` = `Lipid`, `Molecule Type` = "Lipid", `KEGG ID` ),
                                by = c("Molecule ID", "Molecule Type" ) ) %>%
  mutate( `KEGG ID` = if_else( is.na( `KEGG ID` ) | `KEGG ID` == "", `KEGG IDs`, `KEGG ID` ) ) %>%
  select( -`KEGG IDs` ) %>% distinct()

# Remove duplicate CoQ10 entry
molecule_id_table <- filter( molecule_id_table, (`Molecule ID` != "CoQ10") | !is.na( `KEGG ID` ) )

saveRDS( molecule_id_table, "processed-data/molecule-ids.rds" )
