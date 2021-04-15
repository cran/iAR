#' Simulate from an IAR(1) Model
#'
#' Simulates an IAR(1) Time Series Model.
#'
#' @param phi A coefficient of IAR(1) model. A value between 0 and 1
#' @param n Length of the output time series. A strictly positive integer.
#' @param sT Array with observational times.
#'
#' @return A list with the following components:
#' \itemize{
#' \item{times}{ Array with observation times.}
#' \item{series}{ Array with simulated IAR(1) data.}
#' }
#' @export
#' @references
#' \insertRef{Eyheramendy_2018}{iAR}
#'
#' @seealso
#'
#' \code{\link{gentime}}
#'
#' @examples
#'
#' set.seed(6714)
#' st<-gentime(n=100)
#' y<-IAR.sample(phi=0.99,n=100,st)
#' y<-y$series
#' plot(st,y,type='l')
#'
IAR.sample=function(phi,n=100,sT){
  Sigma=matrix(0,ncol=n,nrow=n)
  for( i in 1:n){
    d<-sT[i]-sT[i:n]
    Sigma[i,i:n]=phi**abs(d)
    Sigma[i:n,i]=Sigma[i,i:n]
  }
  b=eigen(Sigma)
  V <- b$vectors
  A=V %*% diag(sqrt(b$values)) %*% t(V)
  e=rnorm(n)
  y=as.vector(A%*%e)
  out=list(series=y,times=sT)
  return(out)
}
