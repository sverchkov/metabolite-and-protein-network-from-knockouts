stop("Are you sure you want to run this? It's a one-time script for making a master molecule ID table.")

library(dplyr)

# Load all the molecule table
source("R/load-data.R")

# Load the table with my integer IDs
spead_table <- readRDS( "processed-data/spread-table-for-similarities.rds" )

# Select out the columns we care about from all the tables
molecule_id_table <- spead_table[c("Molecule ID", "Molecule Name", "Molecule Type", "id")]

lipid_ids <- lipidomics[c("Lipid", "KEGG ID")]

metabolite_ids <- bind_rows(
  metabolomics[c("Metabolites","KEGG IDs")],
  metabolomics_unimputed[c("Metabolites","KEGG IDs")]
)

protein_ids <- bind_rows(
  proteomics_unimputed[c( "Majority protein IDs", "Protein names", "Gene names", "Fasta headers" )],
  proteomics1[c( "Majority protein IDs", "Protein names", "Gene names", "Fasta headers" )],
  proteomics2[c( "Majority protein IDs", "Protein names", "Gene names", "Fasta headers" )],
  proteomics3[c( "Majority protein IDs", "Protein names", "Gene names", "Fasta headers" )],
  proteomics4[c( "Majority protein IDs", "Protein names", "Gene names", "Fasta headers" )]
)

# Remove duplicates and normalize names
molecule_table <- bind_rows(
  lipid_ids %>% distinct() %>%
    transmute( `Molecule ID` = Lipid,
               `Molecule Name` = Lipid,
               `KEGG ID`,
               `Molecule Type` = "Lipid" ),
  metabolite_ids %>% distinct() %>%
    transmute( `Molecule ID` = Metabolites,
               `Molecule Name` = Metabolites,
               `KEGG ID` = `KEGG IDs`,
               `Molecule Type` = "Metabolite"),
  protein_ids %>% distinct() %>%
    transmute( `Molecule ID` = `Fasta headers`, # Ideally we would want this in to be `Majority protein IDs`,
               `Protein ID` = `Majority protein IDs`, # Remove when above updated.
               `Molecule Name` = `Protein names`,
               `Protein Gene` = `Gene names`,
               FASTA = `Fasta headers`,
               `Molecule Type` = "Protein" )
)

# There are for some reason two CoQ10 lipids? Will want to redo everything in a workflow that keeps track of that.
# For now, will hack around to rename one of them
molecule_table <- molecule_table %>%
  mutate( `Molecule ID` = if_else( (`Molecule ID` == "CoQ10") & (`KEGG ID` == ""),
                                   "CoQ10.", `Molecule ID` ) )

# Pretty sure this one KEGG ID is wrong
molecule_table <- molecule_table %>%
  mutate( `KEGG ID` = if_else( `KEGG ID` == "C0023", "C00233", `KEGG ID` ) )

# Merge with spread_table to get integer IDs
molecule_ids <- full_join( molecule_id_table, molecule_table, by = c( "Molecule ID", "Molecule Name", "Molecule Type" ) )

saveRDS( molecule_ids, "processed-data/molecule-ids.rds" )
