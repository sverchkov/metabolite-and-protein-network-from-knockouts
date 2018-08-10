library("profvis")
library("dplyr")
library("futile.logger")

source("R/community-inference.R")

flog.threshold( TRACE )

abs_sim <- similarities %>% mutate( similarity = abs( similarity ) )

profvis({
  membership <- inferCommunitiesGreedily( similarities, simplify = F, loop_limit = 3 )
})