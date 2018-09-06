library( futile.logger )
library( dplyr )

stopifnot( exists( "long_table" ) )

source( "R/tTestPValue.R" )

flog.info("Aggregating across Crispr and Bio replicates")

tested_table <- long_table %>%
  filter( `Molecule ID` != "" ) %>%
  group_by( `Molecule ID`, Knockout ) %>%
  summarize( `Mean log2FC` = mean( log2FC ),
             `StDev log2FC` = sd( log2FC ),
             count = n(),
             `p value` = tTestPValue( log2FC ) ) %>%
  ungroup() %>%
  mutate( `t statistic` = `Mean log2FC` / ( `StDev log2FC` / sqrt( count ) ) )

flog.info("Saving tested table")
saveRDS( tested_table, file = "processed-data/tested_table.rds" )