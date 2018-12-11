## Fuse Networks

# Load things
library( "dplyr" )

annotated_molecule_edges <- readRDS("processed-data/annotated_molecule_edges.rds")
fgnem_ko_edges <- readRDS("processed-data/fgnem_ko_edges.rds")

# Assume we have the fgnem and the molecule similarity net loaded

# Get the KO-protein map

# Merge Networks