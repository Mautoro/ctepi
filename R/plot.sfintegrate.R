#'
#' @export
plot.sfintegrate <- function (sfi, ..., from = NULL, to = NULL, n = 501,
                              ylab = "Integrate", col.01line = "gray70") {
  x_vals <- get("x_vals", envir = environment(sfi) )
  y_vals <- get("y_vals", envir = environment(sfi) )
  if ( is.null(from) ) { from <- min(x_vals) }
  if ( is.null(to) ) { to <- max(x_vals) }
  plot.function(sfi, from = from, to = to, n = n, ..., ylab = ylab)
  abline(h = c(0), col = col.01line, lty = 2)
} 
