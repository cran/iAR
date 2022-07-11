#' Interpolation predictor of the iAR Model
#'
#' One step interpolation predictor of the iAR Model
#'
#' @param x A given phi coefficient of the IAR model.
#' @param st Array with the irregular observational times.
#' @param y Array with the time series observations.
#' @param grid Array with the times in which the interpolation predictor will be computed.
#' @param delta Array with the measurements error standard deviations.
#' @param zero.mean logical; if true, the array y has zero mean; if false, y has a mean different from zero.
#' @param standarized logical; if true, the array y is standarized; if false, y contains the raw time series.
#'
#' @return A list with the following components:
#' \itemize{
#' \item{yhat}{One step interpolation predictor of the iAR model for each given time.}
#' \item{grid}{Times in which the interpolation predictor was computed.}
#' }
#' @export
#' @references
#' \insertRef{Elorrieta_2021}{iAR}
#'
#' @seealso
#'
#' \code{\link{gentime}}, \code{\link{IARsample}}, \code{\link{IARloglik}}
#'
#' @examples
#' n=300
#' set.seed(6713)
#' st<-gentime(n)
#' y<-IARsample(phi=0.9,st=st,n=n)
#' y<-y$series
#' phi=IARkalman(y=y,st=st)$phi
#' napos=10
#' y0=y
#' yest=IARsmoothing(x=phi,st=st[-napos],y=y[-napos],grid=st[napos])
#' print(yest)
IARsmoothing = function(x,st,y,grid,delta=0,zero.mean='FALSE', standarized = "TRUE")
{
    sigma = 1
    mu = 0
    if (sum(delta) == 0) {
        delta = rep(0, length(y))
    }
    if (standarized == "FALSE")
        sigma = var(y)
    if (zero.mean == "FALSE")
        mu = mean(y)
    grid=grid[which(grid>min(st) & grid<max(st))]
    yhat=rep(0,length(grid))
    phi=x
    for(i in 1:length(grid))
    {
        pos=which(grid[i]<st)[1]
        times=c(st[pos-1],grid[i],st[pos])
        delta=diff(times)
        alpha=phi**(delta[1])*(1-phi**(2*delta[2]))/(1-phi**(2*(delta[2]+delta[1])))
        beta=phi**(delta[2])-alpha*phi**(delta[2]+delta[1])
        yhat[i]=alpha*y[pos]+beta*y[pos+1]
    }
    if (standarized == "FALSE")
        yhat=yhat*sigma
    if (zero.mean == "FALSE")
        yhat=yhat+mu
    return(list(yhat=yhat,grid=grid))
}
