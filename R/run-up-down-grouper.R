#' Run up-down grouper
#' @author Yuriy Sverchkov

source("R/small-make-call-matrix.R")
source("R/upstream-downstream-grouper.R")

# Hacky fix
metabolite_call_matrix[ which( is.na( metabolite_call_matrix ) ) ] = 0

ud.groups.1 <- groupDownstreamByUpstream(
  call.matrix = cbind( call_matrix, metabolite_call_matrix ),
  n.cond = length(actions),
  n.upstrm = length(effects),
  n.dnstrm = length(metabolite_effects) )

for ( group in ud.groups[[1]] ) {
  print( "up-proteins" )
  print( sum(group$pos.upstrm) )
  #print( effects[group$pos.upstrm] )
  print( "down-proteins" )
  print( sum(group$neg.upstrm) )
  #print( effects[group$neg.upstrm] )
  
  mets <- as.numeric( group$effects )
  
  up.m = mets[ mets > 0 ]
  if( length( up.m > 0 ) ){
    print( "up-metabolites" )
    print( metabolite_effects[up.m] )
  }
  
  down.m = 0-mets[ mets < 0 ]
  if( length( down.m > 0 ) ){
    print( "down-metabolites" )
    print( metabolite_effects[down.m] )
  }
}

# Make a table summarizing the results
# Table columns: Group, Molecule, Molecule type, Change direction
ud_group_table <- bind_rows( Map( function( group ){
  mets <- as.numeric( group$effects )
  up.m = mets[ mets > 0 ]
  down.m = 0-mets[ mets < 0 ]
  bind_rows(
    tibble( Molecule = hit_counts$`Fasta headers`[ group$pos.upstrm ],
            `Molecule type` = "Protein",
            `Change direction` = "up" ),
    tibble( Molecule = hit_counts$`Fasta headers`[ group$neg.upstrm ],
            `Molecule type` = "Protein",
            `Change direction` = "down" ),
    tibble( Molecule = metabolite_effects[up.m],
            `Molecule type` = "Metabolite",
            `Change direction` = "up" ),
    tibble( Molecule = metabolite_effects[down.m],
            `Molecule type` = "Metabolite",
            `Change direction` = "down" )
  )
}, ud.groups.1[[1]] ), .id = "Group" )