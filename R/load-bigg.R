##################################
## Loading BiGG Metabolic Model ##
##################################

library("jsonlite")
library("dplyr")

umodel <- fromJSON( "raw-data/universal_model.json" )

# umodel$reactions$metabolites is a flux matrix