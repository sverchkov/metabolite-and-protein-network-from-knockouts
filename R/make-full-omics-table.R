#' Make a full -omics table in long-form from the imputed data
#' @author Yuriy Sverchkov

library("dplyr")
library("tidyr")

# Load functions
source("R/makeLongTable.R")
source("R/tTestPValue.R")

# Load data
source("R/load-data.R")
lipidomics <- lipidomics[-123:-126] # Cleaning that maybe shold happen in load?
metabolomics <- metabolomics[-277:-278,]

# Specify structure of resultant table
standard_columns <- c(
  "Molecule ID",
  "Molecule Name",
  "Molecule Type",
  "Knockout",
  "CRISPR Replicate",
  "Biological Replicate",
  "log2FC"
)

full_omics_table <- bind_rows(
  Map( function( the_table ){
    makeLongTable( the_table, colnames( the_table )[-1:-4] ) %>%
      mutate( `Molecule ID` = `Fasta headers`, `Molecule Name` = `Protein names`, `Molecule Type` = "Protein" ) %>%
      select( standard_columns )
  }, list( proteomics1, proteomics2, proteomics3, proteomics4 ) ),
  makeLongTable( metabolomics, colnames( metabolomics )[-1:-2] ) %>%
      mutate( `Molecule ID` = Metabolites, `Molecule Name` = Metabolites, `Molecule Type` = "Metabolite" ) %>%
      select( standard_columns ),
  makeLongTable( lipidomics, colnames( lipidomics )[3:122] ) %>%
      mutate( `Molecule ID` = Lipid, `Molecule Name` = Lipid, `Molecule Type` = "Lipid" ) %>%
      select( standard_columns )
)

saveRDS( full_omics_table, "processed-data/all-omics-replicates.rds" )

tested_omics_table <- full_omics_table %>%
  group_by( `Molecule ID`, `Molecule Name`, `Molecule Type`, `Knockout` ) %>%
  summarize_at( vars( log2FC ), funs( `Mean log2FC` = mean, `SD of log2FC` = sd, `p-Value` = tTestPValue ) ) 

saveRDS( tested_omics_table, "H3KExplorer/data/full-omics-table.rds")
saveRDS( tested_omics_table, "processed-data/tested-omics-table.rds")
