
#' @export
ldbudboptim <- function (z,yy1,yy0, F1,F2,F3=NULL,F4=NULL,delta=0.1,sopY=NULL){
  if (is.null(F3)) F3 <- F1
  if (is.null(F4)) F4 <- F2
  #Y <- union(yy1,yy0)
  #if (!is.null(delta)) {
  #  x.limit <- range(yy1)
  #  y.limit <- range(yy0)
  #  u <- seq(  min( x.limit , y.limit + z ) - 3*delta , max( x.limit , y.limit + z ) + 3*delta , delta  )
  #  u <- sort( union(Y,u) )
  #  yy1 <- u
  #  yy0 <- u
  #} else {
  #  yy1 <- Y
  #  yy0 <- Y
  #}
  yy1 <- range(yy1)
  yy0 <- range(yy0)
  ldbudboptimCpp(z, F1, F2, F3, F4, yy1, yy0,delta=delta, sopY=sopY)
}
