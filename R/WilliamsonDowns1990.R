#' Williamson and Downs (1990) bounds for the distribution of the difference of counterfactual outcomes 
#'
#' Implementation of the bounds for \eqn{U = Y_{1} - Y_{0}} from Williamson and Downs (1990) using partial identification bounds for the counterfactual distributions. The function returns step functions for the lower and upper bounds of the treatment effect distribution, \eqn{P(Y_{1} - Y_{0} \leq y)}, under different assumptions.
#'
#' The bounds are computed over a grid of values determined by the `delta` parameter. When `includesupport` is TRUE, the grid also incorporates all possible values of \eqn{U} based on the observed data. The `delta.ldbudboptim` parameter specifies the grid on which the minimum-convolution is evaluated.
#'
#' @param Y A numeric vector of observed outcomes.
#' @param Z A dummy numeric vector indicating the treatment assignment (1 for treatment, 0 for control).
#' @param y1.limit,y0.limit Numeric vectors of length 2, specifying the lower and upper bounds for the potential outcomes under treatment 1 and 0. If `NULL`, the bounds are automatically set based on the observed treatment outcomes.
#' @param boundaries A character vector specifying the boundary methods to use. Options include "Ignorability", "NoAssumptions", and "FourEpsilon".
#' @param delta.u A scalar defining the spacing of the sequence for \eqn{U = Y_{1} - Y_{0}} values. By default, `delta.u` is determined using `gridsize` and the empirical range of \eqn{U = Y_{1} - Y_{0}}.
#' @param delta.ldbudboptim A scalar defining the spacing for the sequence used in the minimum-convolution steps. By default, `delta.ldbudboptim` is determined using `gridsize` and the empirical range of \eqn{U = Y_{1} - Y_{0}}.
#' @param gridsize Number of elements in the sequences for bounds evaluation (values of \eqn{U = Y_{1} - Y_{0}}) and the grid for minimum-convolution. Default is 2000. `gridsize` is ignored for the respective grid if `delta.u` or `delta.ldbudboptim` is specified. 
#' @param includesupport Logical value indicating whether to include the set \eqn{\{y_1 - y_0 : y_1 \in Y[Z == 1] \wedge y_0 \in Y[Z == 0]\}} in the sequence of values for \eqn{U = Y_{1} - Y_{0}} constructed with `delta.u`. Default is `FALSE`.
#' @param addtoU Optional numeric vector containing specific values to be added to the grid for evaluating the bounds. Default is `NULL`.
#' @param parallel A logical indicating whether parallel computation should be used (if `TRUE`, the `parallel` library will be employed). Default is `FALSE`.
#' @param freeCores The number of CPU cores to reserve for other tasks when `parallel = TRUE`. Default is 1.
#' @param eps11,eps12,eps01,eps02 Scalars defining the epsilon parameters for the "FourEpsilon" boundary method. Default values are 0.
#'
#' @details
#' One way to relax ignorability assumptions is by introducing tolerance for discrepancies between identified and non-identified distributions. The `FourEpsilon` boundaries are based on the following identification assumptions: .
#'
#' \deqn{
#' -\varepsilon_{11} 
#' \leq 
#' P(Y_{z_1} \leq y | Z=z_1) - P(Y_{z_1} \leq y | Z=z_0) 
#' \leq 
#' \varepsilon_{12}
#' \,,
#' }
#' \deqn{
#' -\varepsilon_{01} 
#' \leq 
#' P(Y_{z_0} \leq y | Z=z_1) - P(Y_{z_0} \leq y | Z=z_0) 
#' \leq 
#' \varepsilon_{02}
#' \,,
#' }
#' 
#' where \eqn{\varepsilon_{11}, \varepsilon_{12}, \varepsilon_{01}, \varepsilon_{02}} correspond to \code{eps11}, \code{eps12}, \code{eps01}, and \code{eps02},, respectively.
#' 
#' Warning: The minimum-convolution strongly depends on the delta and delta.ldbudboptim parameters. If delta and delta.ldbudboptim are too large, the estimation will be poor. A reasonable value for delta is given by \code{deltamin(Y[Z==1], Y[Z==0])}. Be cautious: as delta (or delta.ldbudboptim) decreases, the computation time increases significantly.
#' 
#' 
#' 
#' @import parallel
#' @import Rcpp
#' 
#' @export
WilliamsonDowns1990 <- function(Y, Z,  y1.limit=NULL, y0.limit=NULL, 
                                boundaries = c("Ignorability","NoAssumptions","FourEpsilon") ,
                                gridsize = 2000,
                                delta.u = diff( range( Y[Z==1] ) - rev(range( Y[Z==0] )) )/gridsize , 
                                delta.ldbudboptim = diff( range( Y[Z==1] ) - rev(range( Y[Z==0] )) )/gridsize, 
                                includesupport=FALSE, addtoU=NULL,
                                parallel=FALSE, freeCores=1,
                                eps11=0,eps12=0,eps01=0,eps02=0){
  if ( is.null(y1.limit) ) {
    y1.limit <- c( min(Y[Z==1],na.rm = T) , max(Y[Z==1], na.rm = T) )
  }
  if ( is.null(y0.limit) ) {
    y0.limit <- c( min(Y[Z==0],na.rm = T) , max(Y[Z==0], na.rm = T) )
  }
  
  #<-># Lista para guardar las cotas ajustadas
  rval <- list()
  i <- 0
  
  #<-># cat(paste0("delta.u=", delta.u, "\ndelta.ldbudboptim=",delta.ldbudboptim,"\n"))
  
  #<-># Funciones y soporte de variables
  func <- ctefunctions(Y,Z,y1.limit = y1.limit, y0.limit = y0.limit,
                       eps11=eps11, eps12=eps12, eps01=eps01, eps02=eps02)
  sop <- supportU(Y,Z,y1.limit = y1.limit, y0.limit = y0.limit)
  seqX <- function(X,d) seq( X[1],X[2] , d )
  
  #<-># Valores de U = Y(1)-Y(0) sobre el que evalúo las cotas
  #<-># Está mal usar este uu porque el rango de U cambia con y1.limit y y0.limit
  #uu <- seqX( range(sop$sopU) , delta.u) 
  uu <- seqX( y1.limit - rev(y0.limit) , delta.u)
  if (includesupport){
    uu <- sort( union(uu, sop$sopU) )
  }
  #<->#uu <- sop$sopU
  
  yy1 <- sop$sopY1
  yy0 <- sop$sopY0
  
  #<-># 1. Cotas para Y(1)-Y(0) bajo el supuesto de ignorabilidad
  if ( is.element("Ignorability",boundaries) ){
    i <- i+1
    ###
    if (parallel) {
      cat(paste0(i,'. Computing boundaries with ignorability assumptions. Parallelized computation with "parallel" library.\n'))
      ldbudboptimP <- function(x) {
        ldbudboptim(x, range(yy1), range(yy0),
                    F1 = func$ecdfYZ1,
                    F2 = func$ecdfYZ0,
                    delta = delta.ldbudboptim,
                    sopY = addtoU) 
      }
      
      Ignboundsp <- parallel::mclapply( uu , ldbudboptimP , mc.cores = parallel::detectCores() - freeCores )
      
      Ignbounds <- data.frame(z = unlist(Ignboundsp)[3*c(1:length(Ignboundsp))-2],
                              ldb = unlist(Ignboundsp)[3*c(1:length(Ignboundsp))-1],
                              udb = unlist(Ignboundsp)[3*c(1:length(Ignboundsp))]
      )
    } else {
      cat(paste0(i,'. Computing boundaries with ignorability assumption.\n'))
      Ignbounds <- ldbudboptim(uu, range(yy1), range(yy0),
                               F1 = func$ecdfYZ1,
                               F2 = func$ecdfYZ0,
                               delta = delta.ldbudboptim,
                               sopY = addtoU)
    }
    
    ldbIgn <- approxfun(Ignbounds$z, Ignbounds$ldb, method = "constant", 
                        yleft = min(Ignbounds$ldb), yright = max(Ignbounds$ldb), 
                        f = 0, ties = "ordered")
    udbIgn <- approxfun(Ignbounds$z, Ignbounds$udb, method = "constant", 
                        yleft = min(Ignbounds$udb), yright = max(Ignbounds$udb), 
                        f = 0, ties = "ordered")
    
    class(ldbIgn) <- c("stepfun",class(ldbIgn))
    class(udbIgn) <- c("stepfun",class(udbIgn))
    
    rval[[i]] <- list( ldb=ldbIgn, udb=udbIgn)
    names(rval)[i] <- "Ignorability"
  }
  
  #<-># 2, Cotas para Y(1)-Y(0) sin supuestos
  if ( is.element("NoAssumptions",boundaries) ){
    i <- i+1
    if (parallel) {
      cat(paste0(i,'. Computing boundaries with no assumptions. Parallelized computation with "parallel" library.\n'))
      ldbudboptimP <- function(x) {
        ldbudboptim(x, range(yy1), range(yy0),
                    F1 = func$ecdfPIbounds$FlY1, #F1 F_Y(1)  F_X
                    F2 = func$ecdfPIbounds$FuY0, #F2 F^Y(0)  F^Y
                    F3 = func$ecdfPIbounds$FuY1, #F3 F^Y(1)  F^X
                    F4 = func$ecdfPIbounds$FlY0, #F4 F_Y(0)  F_Y
                    delta = delta.ldbudboptim,
                    sopY = addtoU) 
      }
      
      NABoundsp <- parallel::mclapply( uu , ldbudboptimP , mc.cores = parallel::detectCores() - freeCores )
      
      NABounds <- data.frame(z = unlist(NABoundsp)[3*c(1:length(NABoundsp))-2],
                             ldb = unlist(NABoundsp)[3*c(1:length(NABoundsp))-1],
                             udb = unlist(NABoundsp)[3*c(1:length(NABoundsp))]
      )
    } else {
      cat(paste0(i,'. Computing boundaries with no assumptions\n'))
      NABounds <- ldbudboptim(uu, range(yy1), range(yy0),
                              F1 = func$ecdfPIbounds$FlY1, 
                              F2 = func$ecdfPIbounds$FuY0,
                              F3 = func$ecdfPIbounds$FuY1, 
                              F4 = func$ecdfPIbounds$FlY0, 
                              delta = delta.ldbudboptim,
                              sopY = addtoU)
    }
    
    cotal <- approxfun(NABounds$z, NABounds$ldb, method = "constant", 
                       yleft = 0, yright = 1, 
                       f = 0, ties = "ordered")
    cotau <- approxfun(NABounds$z, NABounds$udb, method = "constant", 
                       yleft = 0, yright = 1, 
                       f = 0, ties = "ordered")
    
    class(cotal) <- c("stepfun",class(cotal))
    class(cotau) <- c("stepfun",class(cotau))
        
    rval[[i]] <- list( ldb=cotal, udb=cotau)
    names(rval)[i] <- "NoAssumptions"
  }
  
  
  #<-># 3. Cotas para Y(1)-Y(0) con los supuestos epsilon
  if ( is.element("FourEpsilon",boundaries) ){
    i <- i+1
    
    if (parallel){
      cat(paste0(i,'. Computing boundaries with 4-epsilon assumption. Parallelized computation with "parallel" library.\n'))
      ldbudboptimP <- function(x) {
        ldbudboptim(x, range(yy1), range(yy0),
                    F1 = func$ecdfPIbounds.eps$FlY1eps, #F1 F_Y(1)  F_X
                    F2 = func$ecdfPIbounds.eps$FuY0eps, #F2 F^Y(0)  F^Y
                    F3 = func$ecdfPIbounds.eps$FuY1eps, #F3 F^Y(1)  F^X
                    F4 = func$ecdfPIbounds.eps$FlY0eps, #F4 F_Y(0)  F_Y
                    delta = delta.ldbudboptim,
                    sopY = addtoU) 
      }
      
      e4Boundsp <- parallel::mclapply( uu , ldbudboptimP , mc.cores = parallel::detectCores() - freeCores )
      
      e4Bounds <- data.frame(z = unlist(e4Boundsp)[3*c(1:length(e4Boundsp))-2],
                             ldb = unlist(e4Boundsp)[3*c(1:length(e4Boundsp))-1],
                             udb = unlist(e4Boundsp)[3*c(1:length(e4Boundsp))]
      )
    } else {
      cat(paste0(i,'. Computing boundaries with 4-epsilon assumption\n'))
      e4Bounds <- ldbudboptim(uu, range(yy1), range(yy0),
                              F1 = func$ecdfPIbounds.eps$FlY1eps, 
                              F2 = func$ecdfPIbounds.eps$FuY0eps,
                              F3 = func$ecdfPIbounds.eps$FuY1eps, 
                              F4 = func$ecdfPIbounds.eps$FlY0eps, 
                              delta = delta.ldbudboptim,
                              sopY = addtoU)
    }
    
    cotaleps <- stats::approxfun(e4Bounds$z, e4Bounds$ldb, method = "constant", 
                                yleft = 0, yright = 1, 
                                f = 0, ties = "ordered")
    cotaueps <- stats::approxfun(e4Bounds$z, e4Bounds$udb, method = "constant", 
                                yleft = 0, yright = 1, 
                                f = 0, ties = "ordered")
    
    class(cotaleps) <- c("stepfun",class(cotaleps))
    class(cotaueps) <- c("stepfun",class(cotaueps))
        
    rval[[i]] <- list( ldb=cotaleps, udb=cotaueps)
    names(rval)[i] <- "FourEpsilon"
  }
  
  if (i==0) warning( paste0('There is not valid boundaries methods: "', paste0(boundaries, collapse = '", "'), '".\n') )
  #<-># Puedo añadir un warning que indique si hay métodos que no están implementados. Para después.
  
  #<-># debug
  #<-># i <- i+1
  #<-># rval[[i]] <- list( uu=uu )
  #<->#names(rval)[i] <- "debug"
  rval
}

