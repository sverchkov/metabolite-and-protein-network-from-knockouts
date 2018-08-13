# Reader with specific settings
my.read <- function( filename ) read.csv( file = filename, stringsAsFactors = F, check.names = F )

# Read
proteomics1 <- my.read( "raw-data/20180322_H3K_Batches1and2_Proteomics_1.csv" )
proteomics2 <- my.read( "raw-data/20180322_H3K_Batches1and2_Proteomics_2.csv" )
proteomics3 <- my.read( "raw-data/20180322_H3K_Batches1and2_Proteomics_3.csv" )
proteomics4 <- my.read( "raw-data/20180322_H3K_Batches1and2_Proteomics_4.csv" )
proteomics_unimputed <- my.read( "raw-data/20180322_H3K_Batches1and2_Proteomics_Unimputed.csv" )

metabolomics <- my.read( "raw-data/20180322_H3K_Batches1and2_Metabolomics.csv" )
metabolomics_unimputed <- my.read( "raw-data/20180322_H3K_Batches1and2_Metabolomics_Unimputed.csv" )

lipidomics <- my.read( "raw-data/20180322_H3K_Batches1and2_Lipidomics.csv" )
