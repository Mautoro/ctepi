
#' Calculate ECDF with Partial Identification Bounds (considering NAs)
#'
#' This function calculates the empirical cumulative distribution function (ECDF) 
#' along with its partial identification bounds under no assumptions.
#' It is a modification of the base \code{ecdf()} function and uses the class \code{ecdfPI} 
#' to allow a specific plot method for these objects (\code{plot.ecdfPI}).
#'
#' @param x A numeric vector.
#' @param weight A numeric vector of weights. The sum of weights must be 1.
#' 
#' @return An object of class \code{ecdfPI}, which is a list with:
#' \item{ecdf}{An ECDF function of class \code{ecdf}.}
#' \item{F1}{Lower bound of the ECDF with class \code{ecdfPIbound}.}
#' \item{F2}{Upper bound of the ECDF with class \code{ecdfPIbound}.}
#'
#' @examples
#' x <- c(1, 2, 3, NA, 4, 5)
#' result <- ecdfPI(x)
#' plot(result)
#' plot(result$ecdf)
#' lines(result$F1, vertical=F, pch=19, col = "blue")
#' lines(result$F2, vertical=F, pch=19, col = "red")
#'
#' @export
ecdfPI <- function (x , weight=NULL) {
  #<-># Verifico que haya un peso por elemento de x 
  if ( !is.null(weight) ) {
    if ( length(x) != length(weight) ) stop("Error: Lengths of x and weight are differents.")
    if ( sum(weight) != 1 ) warning( paste("Sum of weight is",sum(weight)) )
  }
  
  #<-># Vector de NAs
  x.na <- x[is.na(x)]
  p.na <- length(x.na) / length(x)
  # ecdf
  x <- x[!is.na(x)]
  n <- length(x)
  if (n < 1) 
    stop("'x' must have 1 or more non-missing values")
  vals <- sort(unique(x))
  if ( is.null(weight) ) {
    Fvals <- cumsum(tabulate(match(x, vals)))/n
  } else {
    #<-># Pesos de cada valor único en x. Uso mi función aggregate_cpp
    #<-># porque es mucho más rápida y no quita valores únicos.
    weight.cum <- aggregate_cpp( weight, x, vals )
    Fvals <- cumsum( weight.cum[,2] )
  }
  ecdf2 <- approxfun(vals, Fvals, 
                     method = "constant", yleft = 0, yright = 1, f = 0, ties = "ordered")
  #<-># Cotas sin supuestos
  ecdf2x <- ecdf2(vals)
  F1vals <- ecdf2x * (1-p.na)
  F2vals <- ecdf2x * (1-p.na) + p.na
  F1 <- approxfun(vals, F1vals, method = "constant", 
                  yleft = 0, yright = 1-p.na, f = 0, ties = "ordered")
  F2 <- approxfun(vals, F2vals, method = "constant", 
                  yleft = p.na, yright = 1, f = 0, ties = "ordered")
  #<-># 
  class(ecdf2) <- c("ecdf", "stepfun", class(ecdf2))
  attr(ecdf2, "call") <- sys.call()
  class(F1) <- c("ecdfPIbound", "stepfun", class(F1))
  class(F2) <- c("ecdfPIbound", "stepfun", class(F2))
  #<-># Salida
  rval <- list( ecdf=ecdf2, F1=F1, F2=F2)
  #<-># Asigno varias variables al entorno de cada función en 'rval'
  for (fn in rval) {
    assign("nobs", n, envir = environment(fn))
    assign("vals", vals, envir = environment(fn))
    assign("Fvals", Fvals, envir = environment(fn))
    assign("F1vals", F1vals, envir = environment(fn))
    assign("F2vals", F2vals, envir = environment(fn))
    #<-># Propiedades de objeto para mi función ldboptim
    assign("yy1", vals, envir = environment(fn))
    assign("yy2", vals, envir = environment(fn))
    assign("F1yy1", F1vals, envir = environment(fn))
    assign("F2yy2", F2vals, envir = environment(fn))
  }
  #<-># Si no le asigno clase "ecdfPI" a rval, R se marea con plot.ecdfPI() porque el objeto no tiene la clase "ecdfPI" descrita después del punto
  class(rval) <- c("ecdfPI",class(rval)) 
  rval
}


