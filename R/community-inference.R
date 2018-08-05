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
#' @import dplyr
#' @import futile.logger
inferCommunities <- function( similarities, save.file = NULL ){

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
  
  # IMPLEMENTATION
  
  # Compute merge costs for initial nodes, prepare community indeces
  nodes <- distinct( similarities, a )
  
  # Compute total similarity
  total <- summarize( similarities, total = sum( similarity ) )$total
  
  flog.trace( "Total similarity of graph is %s.", total )
  
  # Figure out weighted degree of nodes
  flog.trace( "Computing weighted degree of nodes..." )
  degrees <- union_all(
    similarities %>% group_by( a ) %>% summarize( deg = sum( similarity ) ) %>% ungroup() %>% rename( id = a ),
    similarities %>% group_by( b ) %>% summarize( deg = sum( similarity ) ) %>% ungroup() %>% rename( id = b ) ) %>%
    group_by( id ) %>% summarize( degree = sum( deg ) ) %>% ungroup()
  flog.trace( "Weighted degrees computed.")

  # Convert similarities to community sum change
  flog.trace( "Converting similarities to community score change..." )
  communities <- similarities %>%
    left_join( degrees, by = c("a" = "id") ) %>%
    rename( deg_a = degree ) %>%
    left_join( degrees, by = c("b" = "id") ) %>%
    rename( deg_b = degree ) %>%
    transmute( a = as.character( a ), b = as.character( b ),
               cost = similarity - 2 * deg_a * deg_b / total )
  flog.trace( "Community score change computed." )
  
  # The greedy iteration part
  repeat {
    merge_row <- filter( communities, cost == max( cost ) )
    
    if ( 1 != ( n_merge <- nrow( merge_row ) ) )
      flog.error( "The merge candidate is %s rows :/", merge_row )
    
    if ( 0 > ( m_cost <- merge_row$cost ) ){
      flog.trace("Max merge cost is %s, we're done.", m_cost )
      break;
    }
    
    m_ids <- merge_row[ c( "a", "b" ) ]
    flog.trace( "Merging %s with %s...", merge_row$a, merge_row$b )
    new_id <- paste( m_ids, collapse = " " )
    
    # Find all rows where selected communities appear
    find_communities <- mutate( communities, target = ( a %in% m_ids ) || ( b %in% m_ids ) )
    # Will remove them from the community table
    
    # Compute costs of merging with new community
    new_community_costs <-
      union_all( communities %>% filter( a %in% m_ids, !( b %in% m_ids ) ) %>% rename( src = b, tgt = a ),
                 communities %>% filter( b %in% m_ids, !( a %in% m_ids ) ) %>% rename( src = a, tgt = b ) ) %>%
      group_by( src ) %>%
      summarize( new_cost = sum( cost ) ) %>%
      transmute( a = src, b = new_id, cost = new_cost )
    
    # Add to community table
    communities <- union_all( filter( communities, !( a %in% m_ids ), !( b %in% m_ids ) ), new_community_costs )
    
    flog.trace("...merged!")
    
    if( !is.null( save.file ) )
      saveRDS( communities, save.file )
  }
  
  return ( communities )
}

computeModularity <- function( communities ){
  
}