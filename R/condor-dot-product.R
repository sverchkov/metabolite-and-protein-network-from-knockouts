# Compute dot product
a <- readRDS("a.rds")
b <- readRDS("b.rds")
saveRDS( sum( a*b ), "result.rds" )