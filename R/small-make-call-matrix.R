#' Data arranging helpers
#' @author Yuriy Sverchkov

library("dplyr")
library("Matrix")

# Matrix from calls table

# Assume that we have a call table 'call_table'
binary_calls <- call_table %>% filter( p.value < 0.05 ) %>% mutate( sign = sign( `Mean log2fc` ) )

actions <- unique( binary_calls$`Gene Name` )
effects <- unique( binary_calls$`Fasta headers` )

call_matrix <- sparseMatrix(
  i = match( binary_calls$`Gene Name`, actions ),
  j = match( binary_calls$`Fasta headers`, effects ),
  x = binary_calls$sign )