library("dplyr") # Data manipulation
library("futile.logger") # Logging

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
