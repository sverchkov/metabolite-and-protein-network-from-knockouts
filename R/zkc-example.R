# Karate club test

library( "dplyr" )
library( "futile.logger" )

flog.threshold(TRACE)

# Constants
collection_name <- "Zachary's Karate Club Networks"

# Load functions
source("R/community-inference.R")

# Load data
zkc <- read.csv( "toy-data/zachary-karate-club.csv", header = T )

# LP communities test
zkc_lp <- inferCommunitiesLP( zkc )

# Make Cy network
cy_nodes <- zkc_lp %>% transmute( id = paste0( "N", node ), community = paste0( "C", label ) ) %>% as.data.frame()
cy_edges <- zkc %>% transmute( source = paste0( "N", a ), target = paste0( "N", b ) ) %>% as.data.frame()

RCy3::createNetworkFromDataFrames( cy_nodes, cy_edges,
                                   title = paste( "Karate Club LP Communities", date() ),
                                   collection = collection_name )

# Get communities
zkc_members <- inferCommunitiesGreedily( zkc, full_trajectory = T )

# Prepare for network
cy_nodes <- zkc_members %>% transmute( id = paste0( "N", id ), community = paste0( "N", `Peak 1` ) ) %>% as.data.frame()
cy_edges <- zkc %>% transmute( source = paste0( "N", a ), target = paste0( "N", b ) ) %>% as.data.frame()

# Create cytoscape network with RCy3
net_id <- RCy3::createNetworkFromDataFrames(
  nodes = cy_nodes,
  edges = cy_edges,
  title = paste( "Karate Club Communities", date() ),
  collection = collection_name )
