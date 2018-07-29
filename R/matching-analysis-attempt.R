# What I tried to do for presentation July 30, 2018

# Load table with p-values
computed_omics_table <- readRDS( "H3KExplorer/data/full-omics-table.rds" )

# Infer calls
calls_omics <- mutate( computed_omics_table,
                       call = case_when( `p-Value` < 0.05 && abs( `Mean log2FC` > 1 ) ~ sign( `Mean log2FC` ),
                                         TRUE ~ 0 ) )

# Make table shape
calls_spread_table <- tidyr::spread( calls_omics, key = Knockout, value = call, fill = 0 )
calls_spread_table_summarized <- calls_spread_table %>%
  select( -`Mean log2FC`, -`SD of log2FC`, -`p-Value` ) %>%
  group_by( `Molecule ID`, `Molecule Name`, `Molecule Type` ) %>%
  summarize_all( sum ) %>%
  ungroup()
calls_sorted <- calls_spread_table_summarized %>% arrange( desc(`Molecule Type`), `Molecule ID` )
matrix_df <- calls_sorted %>% select( -`Molecule ID`, -`Molecule Name`, -`Molecule Type` )

ko_names <- colnames( matrix_df )
n_proteins <- ( calls_sorted %>% filter( `Molecule Type` == "Protein" ) %>% summarize( count = n() ) )$count
n_others <- ( calls_sorted %>% summarize( count = n() ) )$count - n_proteins

call_m <- t( as.matrix( matrix_df ) )

# Run matcher
matches <- matchUpstreamToDownstream(
  call.matrix = call_m,
  n.cond = length(ko_names),
  n.upstrm = n_proteins,
  n.dnstrm = n_others )

# Get table
match_table <- getMatchDataTable( matches, ko_names, calls_sorted$`Molecule ID`, types = calls_sorted$`Molecule Type` )
