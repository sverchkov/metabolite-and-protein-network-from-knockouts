library("dplyr")

# Reader with specific settings
my.read <- function( filename ) read.csv( file = filename, stringsAsFactors = F, check.names = F )

proteomics1 <- read.csv( "raw-data/20180322_H3K_Batches1and2_Proteomics_1.csv", stringsAsFactors = F, check.names = F )
proteomics2 <- read.csv( "raw-data/20180322_H3K_Batches1and2_Proteomics_2.csv", stringsAsFactors = F, check.names = F )
proteomics3 <- read.csv( "raw-data/20180322_H3K_Batches1and2_Proteomics_3.csv", stringsAsFactors = F, check.names = F )
proteomics4 <- read.csv( "raw-data/20180322_H3K_Batches1and2_Proteomics_4.csv", stringsAsFactors = F, check.names = F )
proteomics_unimputed <- read.csv( "raw-data/20180322_H3K_Batches1and2_Proteomics_Unimputed.csv", stringsAsFactors = F, check.names = F )

metabolomics <- read.csv( "raw-data/20180322_H3K_Batches1and2_Metabolomics.csv", stringsAsFactors = F, check.names = F )
metabolomics_unimputed <- read.csv( "raw-data/20180322_H3K_Batches1and2_Metabolomics_Unimputed.csv", stringsAsFactors = F, check.names = F )

lipidomics <- my.read( "raw-data/20180322_H3K_Batches1and2_Lipidomics.csv" )

## Data exploration and sanity checks

## Which proteins are measured in which set?
protein_names <- bind_rows(
  proteomics1 %>% select( `Fasta headers` ) %>% mutate( Set = "1" ),
  proteomics2 %>% select( `Fasta headers` ) %>% mutate( Set = "2" ),
  proteomics3 %>% select( `Fasta headers` ) %>% mutate( Set = "3" ),
  proteomics4 %>% select( `Fasta headers` ) %>% mutate( Set = "4" ),
  proteomics_unimputed %>% select( `Fasta headers` ) %>% mutate( Set = "unimputed" ) ) %>%
  group_by( `Fasta headers` ) %>%
  summarise( Sets = paste( Set, collapse = "+" ) )

protein_names %>% ungroup() %>% group_by( Sets ) %>% summarize( n() )

#all_protein_names <- union( proteomics1$`Protein names`, proteomics2$`Protein names`, proteomics3$`Protein names`, proteomics4$`Protein names`, proteomics_unimputed$`Protein names` )

