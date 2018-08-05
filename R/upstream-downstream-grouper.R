# Maximum score multi-level network grouping,
# specialized for a KO, protein, lipid, metabolite grouping
#' @author Yuriy Sverchkov

library("futile.logger")

#' Specialized group inference
#' 
#' Matrix is +1,0,-1 calls of change of each protein, metabolite, lipid subject to knockout.
#' 
#' @param call.matrix n by p+m+l matrix
#' @param n number of knockouts
#' @param p number of proteins
#' @param m number of metabolites
#' @param l number of lipids
#' @return Named list of groups of each kind
groupFromCallsSpecialized <- function( call.matrix, n=nrow(call.matrix), p=0, m=0, l=0 ){

  if ( any( dim( pos.call.matrix ) != c( n, p+m+l ) ) )
    stop( "Call matrix dimensions don't match given numbers of proteins, metabolites, and lipids" )
  
  # Group proteins by KOs
  protein.groups <- NULL
  
  # Group metablites by proteins
  metabolite.groups <- groupDownstreamByUpstream( call.matrix[,1:(p+m)], n, p, m )
  
  # Group lipids by proteins
  lipid.groups <- groupDownstreamByUpstream( call.matrix[,c(1:p,p+m+(1:l))], n, p, l )
  
  return ( list( proteins = protein.groups, metabolites = metabolite.groups, lipids = lipid.groups ) )
}

#' Generalized upstream-downstream group inference
#' 
#' @param call.matrix
#' @param n.cond
#' @param n.upstrm
#' @param n.dnstrm
groupDownstreamByUpstream <- function(
  call.matrix,
  n.cond=nrow(call.matrix),
  n.upstrm,
  n.dnstrm )
{

  if ( any( dim( call.matrix ) != c( n.cond, n.upstrm + n.dnstrm ) ) )
    stop( "call matrix dimensions don't match given numbers of upstream and downstream things" )

  call.matrix <- as.matrix( call.matrix )
  
  unexplained <- logical( length = n.dnstrm )
  the.groups <- list()
    
  for ( change in c( "pos", "neg" ) ) { # For each effect change direction
    calls <- if ( change == "pos" )
      ( call.matrix > 0 )[1:n.cond, n.upstrm+(1:n.dnstrm) ]
    else
      call.matrix[1:n.cond, n.upstrm+(1:n.dnstrm) ] < 0
    
    flog.trace("Calls: ", calls, capture = F )
    
    for ( effect in 1:n.dnstrm ) { # For each downstream effect
      
      conds <- calls[,effect] # Select conditions where effect happens
      
      if ( any( conds ) ){
        neg.upstrm <- columnwiseAnd( call.matrix[conds,1:n.upstrm] < 0 )
        pos.upstrm <- columnwiseAnd( call.matrix[conds,1:n.upstrm] > 0 )
        the.groups <- addToGroupStruct( the.groups, neg.upstrm, pos.upstrm, effect, change, conds )
      } else unexplained[ effect ] <- T
    }
  }
  
  return ( list( the.groups, unexplained ) )
}

#' Infer data table from match matrix
#' 
#' Infer a data table from the output of [matchUpstreamToDownstream], with columns
#' groupID, KOs (concatenated string and group indicator). Molecule, Molecule type
#' @import dplyr
getMatchDataTable <- function( matches, kos, names, types ) {
  n_molecules <- length( names )
  bind_rows( Map( function (gid) {
    bind_rows( Map( function( molid ) {
      tibble( Molecule = names[molid], `Molecule Type` = types[molid] )
    }, which( matches$matches[gid,] != 0 ) ) ) %>%
      mutate( KOs = paste( kos[ matches$conditions[gid,] ], collapse = " " ),
              groupID = gid )
    
  }, 1:nrow( matches$matches )) )
}

#' Generalized upstream-downstream pattern matching
#' 
#' Generalized upstream-downstream pattern matching
#' 
#' @param call.matrix
#' @param n.cond
#' @param n.upstrm
#' @param n.dnstrm
matchUpstreamToDownstream <- function(
  call.matrix,
  n.cond=nrow(call.matrix),
  n.upstrm,
  n.dnstrm )
{
  
  if ( any( dim( call.matrix ) != c( n.cond, n.upstrm + n.dnstrm ) ) )
    stop( "call matrix dimensions don't match given numbers of upstream and downstream things" )
  
  call.matrix <- as.matrix( call.matrix )
  
  skip <- logical( length = n.dnstrm )
  
  match_matrix <- NULL
  condition_matrix <- NULL
  
  for ( effect in 1:n.dnstrm ) if ( !skip[effect] ) {
    ind <- n.upstrm + effect
    # Find matches
    matching <- matchCols( call.matrix[, ind], call.matrix )
    antiMatching <- matchCols( -call.matrix[, ind], call.matrix )
    # Update skips
    skip <- skip | (matching|antiMatching)[ n.upstrm + (1:n.dnstrm) ]
    # Update match matrix
    match_matrix <- rbind( match_matrix, as.integer( matching ) - as.integer( antiMatching ) )
    # Update condition matrix
    condition_matrix <- rbind( condition_matrix, call.matrix[, ind] != 0 )
  }
  
  return ( list( matches = match_matrix, conditions = condition_matrix ) )
}

#' Find rows matching a vector
#' 
#' Finds the rows in matrix m that are exactly equal to vector v, returns a logical vector
#' @param v
#' @param m
#' @return logical vector
matchCols <- function( v, m ) {
  apply( m, 2, function(x){ all( x == v ) } )
}

#' Columnwise and that handles input vectors sanely
columnwiseAnd <- function( matrix ){
  if ( is.null( dim( matrix ) ) ) return( matrix )
  else return( apply( matrix, 2, all ) )
}

#' Add an element to the group struct
#' 
#' @param the.groups group struct
#' @param neg binary vector of negative change calls
#' @param pos binary vector of positive change calls
#' @param effect integer indicating effect number
#' @param change string indicating effect change direction (+=pos)
#' @param conds a binary vector of the conditions that cause the change
addToGroupStruct <- function( the.groups, neg, pos, effect, change, conds ){
  
  string <- paste( actions[conds], collapse = " " )

  num <-
    if ( change == "pos" ) effect else -effect
  
  if ( ! ( string %in% names( the.groups ) ) )
    the.groups[[string]] <-
      list(
        pos.upstrm = pos,
        neg.upstrm = neg,
        effects = character() )
  
  the.groups[[string]]$effects <- union( the.groups[[string]]$effects, num )

  return ( the.groups )
}
