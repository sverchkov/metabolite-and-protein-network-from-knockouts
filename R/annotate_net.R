#' Find surprising interactions
#' @author Yuriy Sverchkov
#' 
#' General approach:
#' * Loop through filtered similarities
#' * For each similarity, find interaction in DB to account for it
#' * If found, annotate, if not found, mark as surprising

########
# Init #
########
library( "futile.logger" )

flog.threshold(TRACE)

################
# Load bk data #
################
molecule_ids <- readRDS("processed-data/molecule-ids.rds")

########################
# Function definitions #
########################

#' Look up molecule
#' 
#' Look up database ids that correspond to this molecule integer id
#' @param id integer id
lookupMolecule <- function( id ){
  id_row <- filter( molecule_ids, id == id )
  
}

#' Look up annotation
#' 
#' Find some evidence in the background dbs that a relates to b
lookupAnnotation <- function( a, b ){
  return ( NULL )
}

#' Annotate an edge
#' 
#' Annotate an edge (as either surprising or provide evidence)
#' 
#' @param a integer ID of molecule
#' @param b integer ID of molecule
#' @return a string representing an annotation of how these molecules interact,
#' the string "surprising" if they do not interact, and the string
#' "unknown molecule" if either a or b don't have a match in the background DBs
#' @import futile.logger
annotateEdge <- function ( a, b ) {
  if ( length( a ) > 1 || length( b ) > 1 ) stop( flog.fatal( "annotateEdge can't deal with vectors" ) )
  if ( is.null( molecule_a <- lookupMolecule( a ) ) ||
       is.null( molecule_b <- lookupMolecule( b ) ) ) return ( "unknown molecule" )
  
  annotation <- lookupAnnotation( molecule_a, molecule_b )
  
  if ( is.null( annotation ) ) return ( "surprising" )
  
  return ( annotation )
}

#################
# Script to run #
#################
molecule_net <- readRDS( "processed-data/filtered-cosine-similarities.rds" )

annotated_net <- molecule_net %>% rowwise() %>% mutate( annotation = annotateEdge( a, b ) )

saveRDS( annotated_net, "processed-data/annotated-filtered-cosims.rds" )