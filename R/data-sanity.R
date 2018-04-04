########################################
## Data exploration and sanity checks ##
########################################

source( "R/load-data.R" )

## Proteomics
#############

## Which proteins are measured in which set?

identifying_columns <- c( "Majority protein IDs", "Protein names", "Gene names", "Fasta headers" )

protein_names <- bind_rows(
  proteomics1 %>% select( identifying_columns ) %>% mutate( Set = "1" ),
  proteomics2 %>% select( identifying_columns ) %>% mutate( Set = "2" ),
  proteomics3 %>% select( identifying_columns ) %>% mutate( Set = "3" ),
  proteomics4 %>% select( identifying_columns ) %>% mutate( Set = "4" ),
  proteomics_unimputed %>% select( identifying_columns ) %>% mutate( Set = "unimputed" ) ) %>%
  group_by_at( identifying_columns ) %>%
  summarise( Sets = paste( Set, collapse = "+" ) )

grouped_by_sets <- protein_names %>% ungroup() %>% group_by( Sets )
small_sets <- grouped_by_sets %>% summarize( count = n() ) %>% filter( count < 100 )

View( inner_join( grouped_by_sets, small_sets %>% select( Sets ) ) %>% order_by( Sets ) )

## Which KOs are in which set?

proteomics_knockouts <-
  bind_rows(
    tibble( Knockout = setdiff( colnames( proteomics1 ), identifying_columns ), `Imputed set` = 1 ),
    tibble( Knockout = setdiff( colnames( proteomics2 ), identifying_columns ), `Imputed set` = 2 ),
    tibble( Knockout = setdiff( colnames( proteomics3 ), identifying_columns ), `Imputed set` = 3 ),
    tibble( Knockout = setdiff( colnames( proteomics4 ), identifying_columns ), `Imputed set` = 4 )
  ) %>%
  inner_join(
    tibble( Knockout = setdiff( colnames( proteomics_unimputed ), identifying_columns ), `Unimputed set` = TRUE )
  )

split_KO_names <- Reduce( rbind, strsplit( proteomics_knockouts$Knockout, "-" ) )
colnames( split_KO_names ) <- c( "KO", "CRISPR replicate", "Biological replicate" )
proteoomics_knockouts <- proteomics_knockouts %>% cbind( split_KO_names )

# Are the the KOs in the unimputed the same as the KOs in 1+2+3+4 ?
if ( ! all( proteomics_knockouts$`Unimputed set` & !is.na( proteomics.knockouts$`Imputed set` ) ) ) {
  flog.warn( "Proteomics: KO set of imputed data doesn't match KO set of imputed data." )
  View( proteomics_knockouts )
}

## Metabolomics
###############

identifying_columns_m <- c( "Metabolites", "KEGG IDs" )
metabolites <- metabolomics %>% select( identifying_columns_m ) %>% mutate( `Imputed` = T ) %>%
  full_join( metabolomics_unimputed %>% select( identifying_columns_m ) %>% mutate( `Unimputed` = T ) )

View( metabolites %>% filter( is.na( Imputed ) | is.na( Unimputed ) ) )

if ( ! all( colnames( metabolomics ) == colnames( metabolomics_unimputed ) ) )
  flog.warn( "Metabolomics: KO set of imputed data doesn't match KO set of imputed data." )

metabolomics_knockouts <- setdiff( colnames( metabolomics ), identifying_columns_m )
split_KO_names <- Reduce( rbind, strsplit(metabolomics_knockouts, "-" ) )
colnames( split_KO_names ) <- c( "KO", "CRISPR replicate", "Biological replicate" )
metabolomics_knockouts <- cbind( Knockout = metabolomics_knockouts, split_KO_names )
metabolomics_knockouts <-
  data.frame( metabolomics_knockouts, stringsAsFactors = F, check.names = F ) %>%
  mutate( Metabolomics = T )

## Check which knockouts were used in which set
###############################################

master_knockout_table <- full_join(
  proteomics_knockouts %>% mutate( Proteomics = T ), metabolomics_knockouts )
