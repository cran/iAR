#' Generating Irregularly spaced times
#'
#' Function to generate irregularly spaced times from a mixture of exponential distributions.
#'
#' @param n A positive integer. Length of observation times.
#' @param distribution Distribution of the observation times that will be generated. Default value is "expmixture" for a mixture of exponential distributions. Alternative distributions are "uniform", "exponential" and "gamma".
#' @param lambda1 Mean (1/rate) of the exponential distribution or the first exponential distribution in a mixture of exponential distributions.
#' @param lambda2 Mean (1/rate) of the second exponential distribution in a mixture of exponential distributions.
#' @param p1 Weight of the first exponential distribution in a mixture of exponential distributions.
#' @param p2 Weight of the second exponential distribution in a mixture of exponential distributions.
#' @param a Shape parameter of a gamma distribution or lower limit of the uniform distribution.
#' @param b Scale parameter of a gamma distribution or upper limit of the uniform distribution.
#'
#' @return Array with irregularly spaced observations times
#' @export
#' @references
#' \insertRef{Eyheramendy_2018}{iAR}
#'
#' @seealso
#'
#' \code{\link{IARsample}}
#'
#' @examples
#' st<-gentime(n=100)
#' st<-gentime(n=100,distribution="uniform")
#' st<-gentime(n=100,distribution="gamma",a=1,b=1)
#' st<-gentime(n=100,distribution="exponential",lambda1=1)
gentime<-function (n, distribution="expmixture",lambda1 = 130, lambda2 = 6.5, p1 = 0.15, p2 = 0.85,a=0,b=1)
{
  if(distribution=="expmixture")
  {
    dT <- rep(0, n)
    a <- sample(c(lambda1, lambda2), size = n, prob = c(p1, p2),
                replace = TRUE)
    dT[which(a == lambda1)] = rexp(length(which(a == lambda1)),
                                   rate = 1/lambda1)
    dT[which(a == lambda2)] = rexp(length(which(a == lambda2)),
                                   rate = 1/lambda2)
    sT = cumsum(dT)
  }
  if(distribution=="uniform")
  {
    sT=cumsum(runif(n,a,b))
  }
  if(distribution=="exponential")
  {
    sT=cumsum(rexp(n,rate=1/lambda1))
  }
  if(distribution=="gamma")
  {
    sT=cumsum(rgamma(n,shape=a,scale=b))
  }
  return(sT)
}
