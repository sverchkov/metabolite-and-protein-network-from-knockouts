#' Get p-value-filtered similarities
#' 
#' @author Yuriy Sverchkov

# Libraries
library("dplyr")

# Load necessary data
similarities <- readRDS( "processed-data/cosine-similarities.rds" )
spread_table <- readRDS( "processed-data/spread-table-for-similarities.rds" )

# notes: see https://stats.stackexchange.com/questions/85916/distribution-of-scalar-products-of-two-random-unit-vectors-in-d-dimensions
# for null distribution of cosine similarities.
# Having a null allows us to obtain p-values
# For network we probably need to do multiple testing correction?

# We are assuming that the non-KO columns in spread_table are Molecule ID, Molecule Name, Molecule Type, and id
n_kos <- ncol( spread_table ) - 4

# Bonferroni correction:
alpha = 0.05 / nrow( similarities )

# Cosine similarity p-Value:
cutoff <- 1 - qbeta( alpha, (n_kos-1)/2, (n_kos-1)/2 ) * 2

filtered_similarity <- similarities %>% filter( abs(similarity ) > cutoff )

saveRDS( filtered_similarity, "processed-data/filtered-cosine-similarities.rds" )