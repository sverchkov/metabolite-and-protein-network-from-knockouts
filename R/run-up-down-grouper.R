#' Run up-down grouper
#' @author Yuriy Sverchkov

source("R/small-make-call-matrix.R")
source("R/upstream-downstream-grouper.R")

# Hacky fix
metabolite_call_matrix[ which( is.na( metabolite_call_matrix ) ) ] = 0

ud.groups <- groupDownstreamByUpstream(
  call.matrix = cbind( call_matrix, metabolite_call_matrix ),
  n.cond = length(actions),
  n.upstrm = length(effects),
  n.dnstrm = length(metabolite_effects) )

for ( group in ud.groups[[1]] ) {
  print( "up-proteins" )
  sum(group$pos.upstrm)
  #print( effects[group$pos.upstrm] )
  print( "down-proteins" )
  sum(group$neg.upstrm)
  #print( effects[group$neg.upstrm] )
  
  up.m = group$effects[ group$effects > 0 ]
  if( length( up.m > 0 ) ){
    print( "up-metabolites" )
    print( metabolite_effects[up.m] )
  }
  
  down.m = 0-group$effects[ group$effects < 0 ]
  if( length( down.m > 0 ) ){
    print( "down-metabolites" )
    print( metabolite_effects[down.m] )
  }
}