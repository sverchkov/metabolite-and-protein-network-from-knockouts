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
inferCommunities <- function( similarities, save.file = NULL, full_trajectory = F ){

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
  membership <- union( distinct( similarities, member = a ), distinct( similarities, member = b ) ) %>%
    arrange( member )
  
  if ( convert_ids <- is.integer( membership$member ) )
    membership <- mutate( membership, id = member )
  else
    membership <- mutate( membership, id = row_number() )
  
  membership <- mutate( membership, community = id )
  next_id <- max( membership$id ) + 1
  
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
    transmute( a, b, cost = similarity - 2 * deg_a * deg_b / total )
  
  if ( convert_ids )
    communities <- communities %>%
      left_join( select( membership, member, a_id = id ), by = c("a" = "member") ) %>%
      left_join( select( membership, member, b_id = id ), by = c("b" = "member") ) %>%
      select( a = a_id, b = b_id, cost )
  
  flog.trace( "Community score change computed." )
  
  # For full trajectory
  peak_number <- 1
  prev_m_cost <- 0
  
  # The greedy iteration part
  while ( nrow( communities ) > 0 ) {
    merge_row <- filter( communities, cost == max( cost ) )
    
    if ( 1 != ( n_merge <- nrow( merge_row ) ) ){
      flog.error( "The merge candidate is %s rows :/", merge_row )
      if ( n_merge < 1 ) stop()
      if ( n_merge > 1 ) merge_row = sample_n( merge_row, 1 )
    }
    
    if ( 0 > ( m_cost <- merge_row$cost ) && 0 <= prev_m_cost ){
      if ( full_trajectory ) {
        membership[[paste("Peak", peak_number)]] <- membership$community
        peak_number <- 1 + peak_number
      } else {
        flog.trace("Max merge cost is %s, we're done.", m_cost )
        break;
      }
    }
    prev_m_cost <- m_cost
    
    m_ids <- merge_row[ c( "a", "b" ) ]
    flog.trace( "Merging %s with %s...", merge_row$a, merge_row$b )

    # Find all rows where selected communities appear
    find_communities <- mutate( communities, target = ( a %in% m_ids ) || ( b %in% m_ids ) )
    # Will remove them from the community table
    
    # Compute costs of merging with new community
    new_community_costs <-
      union_all( communities %>% filter( a %in% m_ids, !( b %in% m_ids ) ) %>% rename( src = b, tgt = a ),
                 communities %>% filter( b %in% m_ids, !( a %in% m_ids ) ) %>% rename( src = a, tgt = b ) ) %>%
      group_by( src ) %>%
      summarize( new_cost = sum( cost ) ) %>%
      transmute( a = src, b = next_id, cost = new_cost )
    
    # Add to community table
    communities <- union_all( filter( communities, !( a %in% m_ids ), !( b %in% m_ids ) ), new_community_costs )

    flog.trace("...merged!")

    # Update membership table
    membership <- mutate( membership, community = if_else( community %in% m_ids, as.integer( next_id ), community ) )
    
    # Update next id
    next_id <- 1 + next_id
    
    if( !is.null( save.file ) )
      saveRDS( communities, save.file )
  }
  
  return ( list( communities = communities, membership = membership ) )
}