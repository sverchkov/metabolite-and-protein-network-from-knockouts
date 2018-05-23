# Maximum score multi-level network grouping,
# specialized for a KO, protein, lipid, metabolite grouping
#' @author Yuriy Sverchkov


#' Specialized group inference
#' 
#' Matrices are boolean calls of pos or negative change of each protein, metabolite, lipid subject to knockout.
#' 
#' @param pos.call.matrix n by p+m+l matrix
#' @param neg.call.matrix n by p+m+l matrix
#' @param n number of knockouts
#' @param p number of proteins
#' @param m number of metabolites
#' @param l number of lipids
#' @return Named list of groups of each kind
groupFromCallsSpecialized <- function( pos.call.matrix, neg.call.matrix, n=nrow(pos.call.matrix), p=0, m=0, l=0 ){

  if ( any( dim( pos.call.matrix ) != dim( neg.call.matrix ) ) )
    stop( "positive and negative score matrices have different dimensions!" )
  if ( any( dim( pos.call.matrix ) != c( n, p+m+l ) ) )
    stop( "score matrix dimensions don't match given numbers of proteins, metabolites, and lipids" )
  
  # Group proteins by KOs
  protein.groups <- NULL
  
  # Group metablites by proteins
  metabolite.groups <- groupDownstreamByUpstream(
    pos.call.matrix[,1:(p+m)], neg.call.matrix[,1:(p+m)],
    n, p, m )
  
  # Group lipids by proteins
  lipid.groups <- groupDownstreamByUpstream(
    pos.call.matrix[,c(1:p,p+m+(1:l))], neg.call.matrix[,c(1:p,p+m+(1:l))],
    n, p, l )
  
  return ( list( proteins = protein.groups, metabolites = metabolite.groups, lipids = lipid.groups ) )
}

#' Generalized upstream-downstream group inference
#' 
#' @param pos.call.matrix
#' @param neg.call.matrix
#' @param n.cond
#' @param n.upstrm
#' @param n.dnstrm
groupDownstreamByUpstream <- function(
  pos.call.matrix,
  neg.call.matrix,
  n.cond=nrow(pos.call.matrix),
  n.upstrm,
  n.dnstrm )
{

  if ( any( dim( pos.call.matrix ) != dim( neg.call.matrix ) ) )
    stop( "positive and negative call matrices have different dimensions!" )
  if ( any( dim( pos.call.matrix ) != c( n.cond, n.upstrm + n.dnstrm ) ) )
    stop( "call matrix dimensions don't match given numbers of upstream and downstream things" )

  unexplained <- logical( length = n.dnstrm )
  the.groups <- list()
    
  for ( change in c( "pos", "neg" ) ) { # For each effect change direction
    calls <- if ( change == "pos" )
      pos.call.matrix[1:n.cond, n.upstrm+(1:n.dnstrm) ]
    else
      neg.call.matrix[1:n.cond, n.upstrm+(1:n.dnstrm) ]
    
    for ( effect in 1:n.dnstrm ) { # For each downstream effect
      
      conds <- calls[,effect] # Select conditions where effect happens
      
      neg.upstrm <- apply( neg.call.matrix[calls,1:n.upstrm], 2, all )
      pos.upstrm <- apply( pos.call.matrix[calls,1:n.upstrm], 2, all )
      
      if ( !any( c(neg.upstrm, pos.upstrm) ) ) unexplained[ effect ] <- T
      
      the.groups <- addToGroupStruct( the.groups, neg.upstrm, pos.upstrm, effect, change )
    }
  }
  
  return ( list( the.groups, unexplained ) )
}

#' Add an element to the group struct
#' 
#' @param the.groups group struct
#' @param neg binary vector of negative change calls
#' @param pos binary vector of positive change calls
#' @param effect integer indicating effect number
#' @param change string indicating effect change direction (+=pos)
addToGroupStruct <- function( the.groups, neg, pos, effect, change ){
  
  string <-
    paste(
      ifelse( neg&pos, "*",
              ifelse( neg, "-",
                      ifelse( pos, "+", "." ) ) ), collapse = "" )
  
  num <-
    if ( change == pos ) effect else -effect
  
  if ( ! ( string %in% names( the.groups ) ) )
    the.groups[[string]] <-
      list(
        pos.upstrm = pos,
        neg.upstrm = neg,
        effects = character() )
  
  the.groups[[string]]$effects <- union( the.groups[[string]]$effects, num )
  
  return ( the.groups )
}
