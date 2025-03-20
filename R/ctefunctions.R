#' Partial Identification bounds for counterfactual distributions
#' 
#' This function computes partial identification bounds for the counterfactual outcomes,
#' both with no assumptions, under ignorability and under the 4-epsilon assumption. 
#' It returns step functions bounds and the CTE under ignorability assumption.
#' 
#' @param Y A numeric vector of observed outcomes.
#' @param Z A binary treatment indicator (0 or 1).
#' @param y1.limit A numeric vector of length 2 specifying the lower and upper bounds of the 
#'        support for \eqn{Y_1}. Default is NULL, in which case the bounds are 
#'        based on the observed data.
#' @param y0.limit A numeric vector of length 2 specifying the lower and upper bounds of the 
#'        support for \eqn{Y_0}. Default is NULL, in which case the bounds are 
#'        based on the observed data.
#' @param eps11,eps12,eps01,eps02 Numeric values specifying the bound for 4-epsilon assumption. Default is 1.
#' @param eta11,eta12,eta01,eta02 Numeric values specifying the bound for 4-eta assumption. Default is 1.
#' @param na.rm A logical value indicating whether to remove missing values from the dataset. 
#'        Default is FALSE.
#' 
#' @details
#' One way to relax ignorability assumptions is by introducing tolerance for discrepancies between identified and non-identified distributions. The - **FourEpsilon** boundaries are based on the following identification assumptions: .
#'
#' \deqn{
#' -\varepsilon_{11} 
#' \leq 
#' P(Y_{1} \leq y | Z=1) - P(Y_{1} \leq y | Z=0) 
#' \leq 
#' \varepsilon_{12}
#' \,,
#' }
#' \deqn{
#' -\varepsilon_{01} 
#' \leq 
#' P(Y_{0} \leq y | Z=1) - P(Y_{0} \leq y | Z=0) 
#' \leq 
#' \varepsilon_{02}
#' \,,
#' }
#' 
#' where \eqn{\varepsilon_{11}, \varepsilon_{12}, \varepsilon_{01}, \varepsilon_{02}} correspond to \code{eps11}, \code{eps12}, \code{eps01}, and \code{eps02}, respectively. 
#' 
#' If the observed outcome has missing data, assumption 4-eta is considered:
#' 
#' \deqn{
#' -\eta_{11} 
#' \leq 
#' P(Y_{1} \leq y | Z=1, W=1) - P(Y_{1} \leq y | Z=1, W=0) 
#' \leq 
#' \eta_{12}
#' }
#' 
#' \deqn{
#' -\eta_{01} 
#' \leq 
#' P(Y_{0} \leq y | Z=0, W=1) - P(Y_{0} \leq y | Z=0, W=0) 
#' \leq 
#' \eta_{02}
#' }
#' 
#' where \eqn{\eta_{11}, \eta_{12}, \eta_{01}, \eta_{02}} correspond to \code{eta11}, \code{eta12}, \code{eta01}, and \code{eta02}, respectively. 
#' 
#' 
#' @return A list containing:
#' \itemize{
#'   \item \code{ecdfYZ1} A step function representing the empirical cumulative distribution 
#'         function of \eqn{Y_1} given \eqn{Z=1}.
#'   \item \code{ecdfYZ0} A step function representing the empirical cumulative distribution 
#'         function of \eqn{Y_0} given \eqn{Z=0}.
#'   \item \code{pZ1} The probability \eqn{P(Z = 1)}.
#'   \item \code{pZ0} The probability \eqn{P(Z = 0)}.
#'   \item \code{CTE} The conditional treatment effect (CTE) under ignorability assumption.
#'   \item \code{CTEplus} The positive part of the conditional treatment effect (CTE) under ignorability assumption.
#'   \item \code{CTEminus} The negative part of the conditional treatment effect (CTE) under ignorability assumption.
#'   \item \code{yz1} The sorted values of \eqn{Y} for units with \eqn{Z = 1}.
#'   \item \code{yz0} The sorted values of \eqn{Y} for units with \eqn{Z = 0}.
#'   \item \code{ecdfPIbounds} A list of partial identification bounds for the counterfactual distributions
#'         with no assumptions.
#'   \item \code{ecdfPIbounds.eps} A list of partial identification bounds for the counterfactual distributions
#'         under 4-epsilon and 4-eta assumptions.
#'   \item \code{dataset} The dataset used for the analysis.
#' }
#' 
#'
#' 
#' @examples
#' # Example usage
#' set.seed(123)
#' Y <- rnorm(100)
#' Z <- rbinom(100, 1, 0.5)
#' result <- ctefunctions(Y, Z, y1.limit=c(-1, 2), y0.limit=c(-2, 3))
#' 
#' 
#' @export
ctefunctions <- function(Y, Z, y1.limit=NULL, y0.limit=NULL, 
                      eps11=1,eps12=1,eps01=1,eps02=1,
                      eta11=1,eta12=1,eta01=1,eta02=1, na.rm = F) {
  YZ1min <- min(Y[Z==1],na.rm = T)
  YZ1max <- max(Y[Z==1],na.rm = T)
  YZ0min <- min(Y[Z==0],na.rm = T)
  YZ0max <- max(Y[Z==0],na.rm = T)
  # Verifico que el supuesto de soporte acotado no contradiga a lo observado.
  if ( is.null(y1.limit) ) {
    y1.limit <- c( YZ1min , YZ1max )
  } else if (length(y1.limit) != 2) { 
    stop(paste0("Error: y1.limit is not valid. It should be a vector of length 2.")) 
  } else if (YZ1min < y1.limit[1]) { 
    stop(paste0('Error: y1.limit[1] should be lower than min(Y[Z==1]) = ',YZ1min,'.')) 
  } else if (YZ1max > y1.limit[2]) { 
    stop(paste0('Error: y1.limit[2] should be greater than min(Y[Z==1]) = ',YZ1max,'.')) 
  }
  if ( is.null(y0.limit) ) {
    y0.limit <- c( YZ0min , YZ0max )
  } else if (length(y0.limit) != 2) { 
    stop(paste0("Error: y0.limit is not valid. It should be a vector of length 2.")) 
  } else if (YZ0min < y0.limit[1]) { 
    stop(paste0('Error: y0.limit[1] should be lower than min(Y[Z==0]) = ',YZ0min,'.')) 
  } else if (YZ0max > y0.limit[2]) { 
    stop(paste0('Error: y0.limit[2] should be greater than min(Y[Z==0]) = ',YZ0max,'.')) 
  }
  
  
  dataset <- data.frame( Y , Z )
  if (na.rm) { dataset <- dataset[ !is.na(Y) , ] }
  yz1 <- sort(Y[Z==1])
  yz0 <- sort(Y[Z==0])
  yvals <- sort(Y)
  
  #nz1 <- length(yz1)
  #nz0 <- length(yz0)
  nz1w1 <- length(yz1)
  nz0w1 <- length(yz0)
  nz1w0 <- dataset$Y[dataset$Z==1] |> is.na() |> sum()
  nz0w0 <- dataset$Y[dataset$Z==0] |> is.na() |> sum()
  pW1 <- (nz1w1+nz0w1) / nrow(dataset)
  pW0Z1 <- nz1w0 / ( nz1w0 + nz1w1 )
  pW0Z0 <- nz0w0 / ( nz0w0 + nz0w1 )
  
  pZ1 <- sum(dataset$Z==1) / nrow(dataset)
  pZ0 <- sum(dataset$Z==0) / nrow(dataset)
  
  
  vals1 <- unique(yz1)
  vals0 <- unique(yz0)
  cumsumpY1Z1W1 <- cumsum(tabulate(match(yz1, vals1)))/nz1w1
  cumsumpY0Z0W1 <- cumsum(tabulate(match(yz0, vals0)))/nz0w1
  
  ecdfYZ1 <- approxfun(vals1, cumsumpY1Z1W1 , method = "constant", 
                       yleft = 0, yright = 1, f = 0, ties = "ordered")
  ecdfYZ0 <- approxfun(vals0, cumsumpY0Z0W1, method = "constant", 
                       yleft = 0, yright = 1, f = 0, ties = "ordered")
  class(ecdfYZ1) <- c("stepfun", class(ecdfYZ1))
  class(ecdfYZ0) <- c("stepfun", class(ecdfYZ0))
  
  CTEvals <- ecdfYZ0(yvals) - ecdfYZ1(yvals)
  CTE <- approxfun(yvals, CTEvals , 
                   method = "constant", yleft = 0, yright = 0, f = 0, ties = "ordered")
  i <- CTEvals>=0
  CTEplusvals <- CTEminusvals <- CTEvals
  CTEplusvals[!i] <- 0
  CTEminusvals[i] <- 0
  CTEplus <- approxfun( yvals , CTEplusvals , 
                        method = "constant", yleft = 0, yright = 0, f = 0, ties = "ordered")
  CTEminus <- approxfun(yvals, -CTEminusvals , 
                        method = "constant", yleft = 0, yright = 0, f = 0, ties = "ordered")
  class(CTE) <- c("stepfun", class(CTE))
  class(CTEplus) <- c("stepfun", class(CTEplus))
  class(CTEminus) <- c("stepfun", class(CTEminus))
  
  y_valsF1 <- pZ1 * cumsumpY1Z1W1
  y_valsF2 <- pZ0 + y_valsF1
  y_valsF3 <- pZ0 * cumsumpY0Z0W1
  y_valsF4 <- pZ1 + y_valsF3
  F1 <- approxfun(vals1, y_valsF1 , 
                  method = "constant", yleft = 0, yright = pZ1, f = 0, ties = "ordered")  
  F2 <- approxfun(vals1, y_valsF2, 
                  method = "constant", yleft = pZ0, yright = 1, f = 0, ties = "ordered")
  F3 <- approxfun(vals0, y_valsF3, 
                  method = "constant", yleft = 0, yright = pZ0, f = 0, ties = "ordered")  
  F4 <- approxfun(vals0, y_valsF4, 
                  method = "constant", yleft = pZ1, yright = 1, f = 0, ties = "ordered")
  
  y_valsF1[vals1>=y1.limit[2]] <- 1
  FlY1 <- approxfun( c(y1.limit[1],vals1,y1.limit[2]) , c( pZ1 * cumsumpY1Z1W1[1], pZ1 * cumsumpY1Z1W1 , 1) , 
                     method = "constant", yleft = 0, yright = 1, f = 0, ties = "ordered")
  FuY1 <- approxfun( c(y1.limit[1],vals1,y1.limit[2]) , c(pZ0,pZ0 + pZ1 * cumsumpY1Z1W1,1) , 
                     method = "constant", yleft = 0, yright = 1, f = 0, ties = "ordered")
  FlY0 <- approxfun( c(y0.limit[1],vals0,y0.limit[2]) , c(pZ0 * cumsumpY0Z0W1[1], pZ0 * cumsumpY0Z0W1,1), 
                     method = "constant", yleft = 0, yright = 1, f = 0, ties = "ordered")
  FuY0 <- approxfun( c(y0.limit[1],vals0,y0.limit[2]) , c(pZ1,pZ1 + pZ0 * cumsumpY0Z0W1,1), 
                     method = "constant", yleft = 0, yright = 1, f = 0, ties = "ordered")
  
    
  y_valsf1eps <- cumsumpY1Z1W1 * (1-pW0Z1) * pZ1  +  pmax( 0 , cumsumpY1Z1W1-eta12) * pW0Z1 * pZ1  +
    pmax( 0 , -eps12 + cumsumpY1Z1W1 * (1-pW0Z1) + pW0Z1 * pmax( 0 , cumsumpY1Z1W1-eta12) ) * pZ0
  
  y_valsf2eps <- cumsumpY1Z1W1 * (1-pW0Z1) * pZ1  +  pmin( 1 , cumsumpY1Z1W1+eta11) * pW0Z1 * pZ1  +
    pmin( 1 , eps11 + cumsumpY1Z1W1 * (1-pW0Z1) + pW0Z1 * pmin( 1 , cumsumpY1Z1W1+eta11) ) * pZ0
  
  y_valsf3eps <- cumsumpY0Z0W1 * (1-pW0Z0) * pZ0  +  pmax( 0 , cumsumpY0Z0W1-eta02) * pW0Z0 * pZ0 +
    pmax( 0 , -eps01 + cumsumpY0Z0W1 * (1-pW0Z0) + pW0Z0 * pmax( 0 , cumsumpY0Z0W1-eta02) ) * pZ1
  
  y_valsf4eps <- cumsumpY0Z0W1 * (1-pW0Z0) * pZ0  +  pmin( 1 , cumsumpY0Z0W1+eta01) * pW0Z0 * pZ0 +
    pmin( 1 , eps02 + cumsumpY0Z0W1 * (1-pW0Z0) + pW0Z0 * pmin( 1 , cumsumpY0Z0W1+eta01) ) * pZ1
  
  F1epsilon <- approxfun( vals1 , y_valsf1eps, method = "constant", 
                          yleft = min(y_valsf1eps), yright = max(y_valsf1eps), f = 0, ties = "ordered")
  F2epsilon <- approxfun( vals1 , y_valsf2eps, method = "constant", 
                          yleft = min(y_valsf2eps), yright = max(y_valsf2eps), f = 0, ties = "ordered")
  F3epsilon <- approxfun(vals0, y_valsf3eps, 
                         method = "constant", yleft = min(y_valsf3eps), yright = max(y_valsf3eps), f = 0, ties = "ordered")
  F4epsilon <- approxfun(vals0, y_valsf4eps, 
                         method = "constant", yleft = min(y_valsf4eps), yright = max(y_valsf4eps), f = 0, ties = "ordered")
  
  FlY1eps <- approxfun( c(y1.limit[1],vals1,y1.limit[2]) , c( y_valsf1eps[1], y_valsf1eps , 1) , 
                        method = "constant", yleft = 0, yright = 1, f = 0, ties = "ordered")
  FuY1eps <- approxfun( c(y1.limit[1],vals1, y1.limit[2]) , c( min(y_valsf2eps) , y_valsf2eps , 1 ) , 
                        method = "constant", yleft = 0, yright = 1, f = 0, ties = "ordered")
  FlY0eps <- approxfun( c(y0.limit[1],vals0,y0.limit[2]) , c( y_valsf3eps[1], y_valsf3eps ,1), 
                        method = "constant", yleft = 0, yright = 1, f = 0, ties = "ordered")
  FuY0eps <- approxfun( c(y0.limit[1],vals0, y0.limit[2]) , c( min(y_valsf4eps) , y_valsf4eps , 1), 
                        method = "constant", yleft = 0, yright = 1, f = 0, ties = "ordered")
  
  class(F1) <- c("ecdf", "stepfun", class(F1))
  class(F2) <- c("ecdf", "stepfun", class(F2))
  class(F3) <- c("ecdf", "stepfun", class(F3))
  class(F4) <- c("ecdf", "stepfun", class(F4))
  class(FlY1) <- c("ecdf", "stepfun", class(FlY1))
  class(FuY1) <- c("ecdf", "stepfun", class(FuY1))
  class(FlY0) <- c("ecdf", "stepfun", class(FlY0))
  class(FuY0) <- c("ecdf", "stepfun", class(FuY0))
  class(F1epsilon) <- c("ecdf", "stepfun", class(F1epsilon))
  class(F2epsilon) <- c("ecdf", "stepfun", class(F2epsilon))
  class(F3epsilon) <- c("ecdf", "stepfun", class(F3epsilon))
  class(F4epsilon) <- c("ecdf", "stepfun", class(F4epsilon))
  class(FlY1eps) <- c("ecdf", "stepfun", class(FlY1eps))
  class(FuY1eps) <- c("ecdf", "stepfun", class(FuY1eps))
  class(FlY0eps) <- c("ecdf", "stepfun", class(FlY0eps))
  class(FuY0eps) <- c("ecdf", "stepfun", class(FuY0eps))
  
  list( ecdfYZ1 = ecdfYZ1 , 
        ecdfYZ0 = ecdfYZ0 ,
        pZ1=pZ1,
        pZ0=pZ0,
        CTE = CTE ,
        CTEplus = CTEplus ,
        CTEminus = CTEminus ,
        yz1 = yz1,
        yz0 = yz0,
        ecdfPIbounds = list(F1=F1,F2=F2,F3=F3,F4=F4,
                       FlY1=FlY1, FuY1=FuY1, FlY0=FlY0, FuY0=FuY0),
        ecdfPIbounds.eps = list(F1epsilon=F1epsilon, F2epsilon=F2epsilon,
                           F3epsilon=F3epsilon, F4epsilon=F4epsilon,
                           FlY1eps=FlY1eps, FuY1eps=FuY1eps,
                           FlY0eps=FlY0eps, FuY0eps=FuY0eps),
        dataset = dataset
  )
}
