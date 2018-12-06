## Run FGNEM on RESP data
library("KnockoutNets")
library("dplyr")
library("parallel")

# Load features
feature_t_df <- readRDS("processed-data/feature_t_df.rds")

# Strip molecule names
original_matrix <- as.matrix( select( feature_t_df, -`Molecule ID` ) )

# Sample bootstrap samples
n <- 10000 # Number of molecules
b <- 100 # Number of bootstrap samples

# Parallel setup
n.cores <- detectCores()-1
#cl <- makeCluster( n.cores )

# Make data
bootstrap_matrices <- mclapply(
  1:b,
  function ( i, om, n ){
    expr <- om[ sample( nrow( om ), size = n, replace = T ), ] # should we define rownames?
    # Build object for input to FGNEM
    return(
      eg <- list( egenes = expr
                , knockdown.cols = colnames( expr )
                , lof = colnames( expr )
                , stddev = apply( expr, 2, sd, na.rm = TRUE ) ) )
  },
  original_matrix,
  n,
  mc.cores = n.cores
)

# Save
saveRDS( bootstrap_matrices, file = "temporary-data/bootstrap_matrices.rds" )

# FGNEM settings
fgnem_params = paramGen( 1.5 , 1 ) # Defaults, maybe worth changing

# Parallel execute
nem_bootstrap_results <- mclapply(
  bootstrap_matrices,
  scoreBestModelEstimate,
  params = fgnem_params,
  doTransitivity = TRUE,
  summarization = max, # or logsum
  mc.cores = n.cores
)

# Save results
saveRDS( nem_bootstrap_results, file = "results/nem_bootstrap_results.rds" )

## Unused:

# Write data to a table in the format that KnockoutNets likes
#
# Particularly:
# knockdown.cols Vector of knockdown names (one per expression matrix column)
# lof Loss-of-function genes (vector)
# stddev Vector of standard deviations (one per expression matrix column)
# expr Matrix of gene expression levels
# file Output file
# append Whether to append to the file
# write.egene.tab( knockdown.cols = colnames( expr ),
#                  lof = colnames( expr ),
#                  stddev = apply( expr, 2, sd, na.rm = TRUE ),
#                  expr,
#                  file = "results/respiration.fgnem.tab",
#                  append = FALSE)
