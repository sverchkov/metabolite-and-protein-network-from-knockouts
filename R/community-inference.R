#' Community detection, greedy algorithm
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
#' @author Yuriy Sverchkov
inferCommunitiesGreedily <- function( similarities, save.file = NULL, full_trajectory = F, simplify = T, loop_limit = Inf ){

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
  
  my <- initCommunitySearch( similarities )
  next_id <- max( my$membership$id ) + 1
  
  # For full trajectory
  peak_number <- 1
  prev_m_cost <- 0
  
  # The greedy iteration part
  while ( nrow( my$communities ) > 0 ) {
    merge_row <- filter( my$communities, cost == max( cost ) )
    
    if ( 1 != ( n_merge <- nrow( merge_row ) ) ){
      flog.error( "The merge candidate is %s rows :/", merge_row )
      if ( n_merge < 1 ) stop()
      if ( n_merge > 1 ) merge_row = sample_n( merge_row, 1 )
    }
    
    if ( 0 > ( m_cost <- merge_row$cost ) && 0 <= prev_m_cost ){
      if ( full_trajectory ) {
        my$membership[[paste("Peak", peak_number)]] <- my$membership$community
        peak_number <- 1 + peak_number
      } else {
        flog.trace("Max merge cost is %s, we're done.", m_cost )
        break;
      }
    }
    prev_m_cost <- m_cost
    
    # We're going to call the two communities we're merging yin and yang as a shorthand
    yin <- merge_row$a
    yang <- merge_row$b
    
    flog.trace( "Merging %s with %s...", yin, yang )

    ## Compute costs of merging with new community
    right_communities <- my$communities %>%
      filter( (a == yin) || (a == yang),
              b != yin, b != yang ) %>%
      rename( src = b, tgt = a )
    left_communities <- my$communities %>%
      filter( (b == yin) || (b == yang),
              a != yin, a != yang ) %>%
      rename( src = a, tgt = b )
    
    new_community_costs <- union_all( right_communities, left_communities ) %>%
      group_by( src ) %>%
      summarize( new_cost = sum( cost ) ) %>%
      transmute( a = src, b = next_id, cost = new_cost )
    
    # Add to community table
    my$communities <- union_all( 
      filter( my$communities, a != yin, a != yang, b != yin, b != yang ),
      new_community_costs )

    flog.trace("...merged!")

    # Update membership table
    my$membership <- mutate( my$membership,
                             community = if_else(
                               (community == yin) | (community == yang),
                               as.integer( next_id ), community ) )
    
    # Update next id
    next_id <- 1 + next_id
    
    if( !is.null( save.file ) )
      saveRDS( my$communities, save.file )
    
    iteration <-
      if ( !exists( "iteration" ) ) 1
      else 1 + iteration
    if ( iteration >= loop_limit ) break;
  }
  
  if ( simplify )
    return ( my$membership )
  else
    return ( my )
}

#' Initialize community search objects
#' 
#' Initialize the communities and membership tables
initCommunitySearch <- function( similarities ) {
  # Compute merge costs for initial nodes, prepare community indeces
  membership <- union( distinct( similarities, member = a ), distinct( similarities, member = b ) ) %>%
    arrange( member )
  
  if ( convert_ids <- is.integer( membership$member ) )
    membership <- mutate( membership, id = member )
  else
    membership <- mutate( membership, id = row_number() )
  
  membership <- mutate( membership, community = id )
  
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
  
  return ( list( membership = membership, communities = communities ) )
}

voteForLabelPropagation <- function( edges, nodes ) {
  left_join( edges, nodes, by = c( "src" = "node" ) ) %>%
    group_by( tgt, label, add = F ) %>%
    summarize( vote = sum( weight ) ) %>%
    group_by( tgt, add = F ) %>%
    filter( vote == max( vote ) ) %>%
    select( node = tgt, new_label = label )
}

#' Community detection, label propagation
#' 
#' An implementation of community detection by label propagation in an undirected weighted graph based on
#' Phys Rev E 76, 036106 (2007)
#' 
#' @param unique_edges a data frame with columns a, b, weight representing the connections between nodes.
#' We assume undirected graph, and therefore b < a.
#' @return a data frame with two columns, n - node id (taken from a, b of input) and cluster - unique cluster ID.
#' @import dplyr
#' @import futile.logger
#' @author Yuriy Sverchkov
inferCommunitiesLP <- function( unique_edges, async_prop = .5 ){
  
  # For this algorithm it's more convenient to just have all edges listed twice
  flog.trace( "Converting distinct edges to bidirectional edges..." )
  edges <- union_all( select( unique_edges, src = a, tgt = b, weight ),
                      select( unique_edges, src = b, tgt = a, weight ) )
  flog.trace( "Making sure edges are unique..." )
  edges <- edges %>% distinct( src, tgt, .keep_all = T )
  
  # Create node table and initialize label table
  flog.trace( "Making node table..." )
  nodes <- distinct( edges, node = src ) %>% mutate( label = node )
  nodes_array <- nodes$node
  
  repeat {
    flog.trace( "Label propagation: Number of communities: %s.", nrow( distinct( nodes, label ) ) )
    
    # Select first batch of nodes to update
    first_batch <- nodes %>% select( node ) %>% sample_frac(async_prop )
    
    # Propagate votes from first batch
    first_batch_votes <- edges %>%
      right_join( first_batch, by = c( "tgt" = "node" ) ) %>%
      voteForLabelPropagation( nodes ) %>%
      sample_n( 1 ) %>% ungroup()
    
    # Update nodes
    nodes <- left_join( nodes, first_batch_votes, by = "node" ) %>%
      mutate( label = if_else( is.na( new_label ), label, new_label ) ) %>%
      select( node, label )
    
    # Get votes from all
    votes <- voteForLabelPropagation( edges, nodes )
    
    # Check whether we're done
    checks <- votes %>% ungroup() %>%
      left_join( nodes, by = "node" ) %>%
      group_by( node ) %>%
      summarize( concensus = any( label == new_label ) ) %>%
      ungroup() %>%
      summarize( done = all( concensus ) )
    
    if ( checks$done ) break;
    
    # Propagate votes from all
    nodes <- votes %>%
      sample_n( 1 ) %>%
      ungroup() %>%
      select( node, label = new_label )
  }
  
  return ( nodes )
}