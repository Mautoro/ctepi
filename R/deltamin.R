#' @export
deltamin <- function(yy1,yy0){
  aux <- outer(yy1,yy0, function(x, y) abs(x - y))
  aux[aux==0] <- Inf
  min(aux)
}
