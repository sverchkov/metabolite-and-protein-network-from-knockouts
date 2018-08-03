#' Community detection
#' 
#' @author Yuriy Sverchkov

#' Community detection
#' 
#' A naive implementation of community detection in a fully-connected weighted graph based on
#' Newman ME, Girvan M. Finding and evaluating community structure in networks.
#' Physical review E. 2004 Feb 26;69(2):026113.
#' 
#' @param similarities a data frame with columns a, b, similarity representing the similarities between nodes.
#' We assume undirected graph, and therefore b < a.
#' @return a data frame with two columns, n - node id (taken from a, b of input) and cluster - unique cluster ID.
inferCommunities <- function( similarities ){

  edges %>% filter( community = c ) %>% summarize( sum() )
  
  # Modularity = sum_c (e_cc - a_c^2)
  # where c - communities
  #       e_xy - proportion of edges from community x to community y
  #       a_c - proportion of edges going into community c total
  
  # When two communities (x and y) are merged, the change in sum is:
  # (e_xx + e_xy + e_yy - (a_x + a_y)^2) - (e_xx - a_x^2) - (e_yy - a_y^2) =
  # e_xy - 2 a_x a_y
  # The prospective cost of merging z to the newly merged group is
  # e_zx + e_zy - 2 a_z ( a_x + a_y )
  # which is the sum of the costs of merging to a_x and a_y
  
  # Proposed algorithm:
  # compute merge costs for initial nodes
  # merge best
  # update merge costs table
  # repeat
  # either stop at first peak or report all peaks
  
  return( stop() )
}

computeModularity <- function( communities ){
  
}