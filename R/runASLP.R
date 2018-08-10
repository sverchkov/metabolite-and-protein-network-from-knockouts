# Libraries
library( "dplyr" )
library( "futile.logger" )
library( "CommunityInference" )

# Constants
molecule_adj_file <- "processed-data/molecule_adjacency.rds"

# Init
flog.threshold(TRACE)

# Funcitons
#source("R/community-inference.R")

# Load similarity
if ( !exists( "similarities" ) ) similarities <- readRDS( "processed-data/cosine-similarities.rds" )

if ( !exists( "molecule_adj" ) ){
  if ( !file.exists( molecule_adj_file ) ){

    n_kos = 29
    flog.warn( "Hard coded number of KOs to %s.", n_kos )

    molecule_edges <- similarities %>%
      mutate( weight = abs( qnorm( pbeta( ( similarity + 1 ) / 2, (n_kos-1)/2, (n_kos-1)/2 ) ) ) )
    
    par_cl <- parallel::makeCluster( parallel::detectCores()-1 )
    
    flog.debug( "Making molecule adjacency matrix..." )
    molecule_adj <- adjacencyMatFromDF( molecule_edges, cluster = par_cl )
    saveRDS( molecule_adj, molecule_adj_file )
    
    parallel::stopCluster( par_cl )
  } else
    molecule_adj <- readRDS( molecule_adj_file )
}

flog.debug( "Applying label propagation to molecule similarities" )

molecules_lp <- inferCommunitiesASLP( adj_mat = molecule_adj )

saveRDS( molecules_lp, file = "processed-data/molecule-communities-aslp.rds")

flog.debug( "Done!" )