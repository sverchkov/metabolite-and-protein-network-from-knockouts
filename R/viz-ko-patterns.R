library("dplyr")
library("futile.logger")

# Load tables if needed
if( !exists( tested_omics_table ) )
  tested_omics_table <- readRDS( "processed-data/tested-omics-table.rds" )

if( !exists( spread_table ) )
  spread_table <- readRDS( "processed-data/spread-table-for-similarities.rds" )

# Need a way to map cytoscape IDs to molecules
id_map <- spread_table %>% transmute( `Molecule ID`, `Molecule Name`, `Molecule Type`, CyID = paste0( "M", id ) )

flog.warn("Using hard-coded cytoscape network UID")
net_id <- 62510 # RCy3::getNetworkSuid()

node_table <- RCy3::getTableColumns( table = "node", network = net_id )

cluster_map <- left_join( id_map,
                          node_table %>% select( id, Cluster = `__glayCluster` ),
                          by = c("CyID"="id") )

# Get calls
calls_df <- tested_omics_table %>% mutate( Call = ( abs( `Mean log2FC` ) > 1 ) && ( `p-Value` < 0.05 ) )

# Summarize call stats per cluster
cluster_calls <- left_join( calls_df, cluster_map, by = "Molecule ID" ) %>%
  group_by( Cluster, Knockout ) %>%
  summarize( `Call Fraction` = mean( Call ), `Mean log2FC` = mean( `Mean log2FC`) ) %>%
  ungroup()

# Plot the profiles
library( ggplot2 )
ggplot( data = cluster_calls %>% filter( Cluster <=3 ), mapping = aes( x = Knockout, y = `Mean log2FC` ) ) +
  facet_grid( rows = vars( Cluster ), labeller = "label_both" ) +
  geom_col()

stop()
grouped_spread_table <- spread_table %>% mutate( CyID = paste0( "M", id ) )
grouped_spread_table <- left_join( grouped_spread_table, node_table, by = c("CyID"="id") )
grouped_spread_table <- grouped_spread_table %>% group_by( `__glayCluster` )
grouped_spread_table
ko_patterns <- grouped_spread_table %>% select( -`Molecule ID`, -`Molecule Name`, -`Molecule Type`, -id, -CyID, -SUID, -`shared name`, -Type, -name, -selected )
ko_patterns
ko_pattern_sum <- ko_patterns %>% summarize_all( function(x) sum(abs(x)) )
ko_pattern_sum
ko_pattern_sum <- ko_patterns %>% summarize_all( function(x) mean(abs(x)) )
ko_pattern_sum
image( ko_pattern_sum %>% select( -`__glayCluster` ) )
image( as.matrix( ko_pattern_sum %>% select( -`__glayCluster` ) ) )
hist( as.matrix( ko_pattern_sum %>% select( -`__glayCluster` ) )[] )
View( ko_pattern_sum %>% mutate_all( > 4 ) )
View( ko_pattern_sum %>% mutate_all( function(x) x > 4 ) )
