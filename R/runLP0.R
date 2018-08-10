# Libraries
library( "dplyr" )
library( "futile.logger" )
library( "CommunityInference" )

# Init
flog.threshold(TRACE)

# Funcitons
#source("R/community-inference.R")

# Load similarity
if ( !exists( "similarities" ) ) similarities <- readRDS( "processed-data/cosine-similarities.rds" )

n_kos = 29
flog.warn( "Hard coded number of KOs to %s.", n_kos )

molecule_edges <- similarities %>%
  mutate( weight = abs( qnorm( pbeta( ( similarity + 1 ) / 2, (n_kos-1)/2, (n_kos-1)/2 ) ) ) )

flog.debug( "Applying label propagation to molecule similarities" )

molecules_lp <- inferCommunitiesASLP( molecule_edges )

saveRDS( molecules_lp, file = "processed-data/molecule-communities-aslp.rds")

flog.debug( "Done!" )