#' Generating Irregularly spaced times
#'
#' Function to generate irregularly spaced times from a mixture of exponential distributions.
#'
#' @param n A positive integer. Length of observations times.
#' @param lambda1 Mean (1/rate) of the first exponential distribution.
#' @param lambda2 Mean (1/rate) of the second exponential distribution.
#' @param p1 Weight of the first exponential distribution.
#' @param p2 Weight of the second exponential distribution.
#'
#' @return Array with irregularly spaced observations times
#' @export
#' @references
#' \insertRef{Eyheramendy_2018}{iAR}
#'
#' @seealso
#'
#' \code{\link{IAR.sample}}
#'
#' @examples
#' st<-gentime(n=100)
gentime<-function(n,lambda1=130,lambda2=6.5,p1=0.15,p2=0.85)
{
  dT<-rep(0,n)
  a<-sample(c(lambda1,lambda2),size=n,prob=c(p1,p2),replace=TRUE)
  dT[which(a==lambda1)]=rexp(length(which(a==lambda1)),rate=1/lambda1)
  dT[which(a==lambda2)]=rexp(length(which(a==lambda2)),rate=1/lambda2)
  sT=cumsum(dT)
  return(sT)
}
