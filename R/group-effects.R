# Intervention-set effect grouping
#' @author Yuriy Sverchkov

#' Helper function that converts a binary vector to a string of . (=FALSE) and | (=TRUE)
#' @param x vector of true/false values
#' @return string of . and |
binaryVectorAsString <- function ( x ) paste( ifelse( x, "|", "." ), collapse = "" )

#' Add a group-effect pair to the group-effect list
addGroupEffectPair <- function ( theList, group, effect ){
  key <- binaryVectorAsString( group )
  if ( !all( key %in% names( theList ) ) )
    theList[[key]] <- list( group=group, effects=c() )
  theList[[key]]$effects <- union( theList[[key]]$effects, effect )
  
  return ( theList )
}

#' Get intervention sets and effect groups for a deterministic [-1,0,1] matrix
#' @param m the effect matrix (rows are interventions, columns are effects)
getSetEffectGroupsD <- function( m ){
  groups <- list()
  for ( col in 1:ncol(m) ) {
    # Add groups in which we saw +ve effect
    groups <- addGroupEffectPair( groups, m[,col] > 0, col )
    # Add groups in which we saw -ve effect
    groups <- addGroupEffectPair( groups, m[,col] < 0, -col )
  }
  
  return( groups )
}
