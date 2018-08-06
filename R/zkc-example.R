# Karate club test

library( "dplyr" )
library( "futile.logger" )

flog.threshold(TRACE)

# Load functions
source("R/community-inference.R")

# Load data
zkc <- read.csv( "toy-data/zachary-karate-club.csv", header = T )

# Get communities
zkc_c <- inferCommunitiesGreedily( zkc, full_trajectory = T )

# Prepare for network
cy_nodes <- zkc_c$membership %>% transmute( id = paste0( "N", id ), community = paste0( "N", `Peak 1` ) ) %>% as.data.frame()
cy_edges <- zkc %>% transmute( source = paste0( "N", a ), target = paste0( "N", b ) ) %>% as.data.frame()

# Create cytoscape network with RCy3
net_id <- RCy3::createNetworkFromDataFrames(
  nodes = cy_nodes,
  edges = cy_edges,
  title = paste( "Karate Club Communities", date() ),
  collection = "Zachary's Karate Club Networks" )
