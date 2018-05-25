#' Get a t-test p-value and fail gracefully
#' 
#' @param x vector of samples to test
#' @author Yuriy Sverchkov
tTestPValue <- function( x ){
  tryCatch( t.test(x)$p.value, error = function( anything ) NA )
}