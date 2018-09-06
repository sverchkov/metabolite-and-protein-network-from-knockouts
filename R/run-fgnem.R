## Run FGNEM on RESP data
library("KnockoutNets")
library("dplyr")

# Load features
spread_table <- readRDS("processed-data/spread-table-for-similarities.rds")

expr = as.matrix( select( spread_table, -"Molecule ID", -"Molecule Name", -"Molecule Type", -"id" ) )
rownames( expr ) = spread_table$`Molecule ID`

# Limiter for experimentation
# expr = expr[1:10,1:10]

# Build object for input to FGNEM
eg = list( egenes = expr
         , knockdown.cols = colnames( expr )
         , lof = colnames( expr )
         , stddev = apply( expr, 2, sd, na.rm = TRUE ) )

# FGNEM settings
params = paramGen( 1.5 , 1 ) # Defaults, maybe worth changing

# Run FGNEM
results <- scoreBestModelEstimate( eg
                                 , params = params
                                 , doTransitivity = FALSE
                                 , summarization = max # or logsum
           )

# Save results
saveRDS( results, file = "results/fgnem-v1.rds" )

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
