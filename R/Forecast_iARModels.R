#' Forecast from iAR package model's
#'
#' Forecast with any of the models available in the iAR package
#'
#' @param phi Autocorrelation coefficient estimated by the method specified.
#' @param y Array with the time series observations.
#' @param st Array with the observational times.
#' @param tAhead The time ahead for which the forecast is required.
#' @param model model to be used for the forecast. The default is to use the iAR model. Other models available are "iAR-T", "iAR-Gamma", "CiAR" and "BiAR".
#' @param mu Level parameter of the IAR-Gamma process. A positive value.
#' @param phiI Imaginary parameter of CIAR model or Cross-correlation parameter of BIAR model.
#' @param nu degrees of freedom parameter of iAR-T model.
#' @param level significance level for the confidence interval. The default value is 95.
#'
#' @return A dataframe with the following columns:
#' \itemize{
#' \item{tAhead}{ The time ahead used for the forecast.}
#' \item{forecast}{ Point forecast in the time ahead required.}
#' \item{stderror}{ Standard error of the forecast.}
#' \item{lowerCI}{ Lower limit of the confidence interval.}
#' \item{upperCI}{ Upper limit of the confidence interval.}
#' }
#' @export
#' @references
#' \insertRef{Eyheramendy_2018}{iAR}
#'
#' @seealso
#'
#' \code{\link{gentime}}, \code{\link{IARforecast}}, \code{\link{IARgforecast}}, \code{\link{IARforecast}}, \code{\link{BIARforecast}}
#'
#' @examples
#' st <- gentime(n=200,lambda1=15,lambda2=2)
#' y  <- IARsample(phi=0.9,n=200,st=st)
#' model<-IARloglik(y=y$series,st=st)
#' phi=model$phi
#' forIAR<-IARforecast(phi=phi,y$series,st=st,tAhead=c(1.3),standardized=FALSE,zero.mean=FALSE)
#' forIAR
#' forIAR<-Forecast_iARModels(phi=phi,y=y$series,st=st,tAhead=c(1.3,2.6))
#' forIAR
Forecast_iARModels<-function(phi,y,st,tAhead,model="iAR",mu=NULL,phiI=NULL,nu=NULL,level=95)
{
  if(model == "iAR")
  {
    y=as.numeric(y)
    sigma = sqrt(var(y))
    mu = mean(y)
    y1=(y-mu)/sigma
    fore=IARforecast(phi=phi,y=y1,st=st,tAhead=tAhead)
    fore=fore*sigma+mu
    error=sigma*(1-phi**(2*tAhead))
  }
  if(model == "iAR-T")
  {
    fore=IARforecast(phi=phi,y=y,st=st,tAhead=tAhead,standardized=FALSE,zero.mean=FALSE)
    error=((nu-2)/nu)*sigma*(1-phi**(2*tAhead))
  }
  if(model == "iAR-Gamma")
  {
    y1=y
    if(min(y)<0)
      y1=y-min(y)
    fore=IARgforecast(phi=phi,mu=mu,y=y1,st=st,tAhead=tAhead)
    if(min(y)<0)
      fore=fore+min(y)
    error=sigma*(1-phi**(2*tAhead))
  }
  if(model == "CiAR")
  {
    if(is.null(phiI))
      phiI=0
    sigma = sqrt(var(as.numeric(y)))
    mu = mean(as.numeric(y))
    y1=y
    fore=CIARforecast(phiR=phi,phiI=phiI,y1=y1,st=st,tAhead=tAhead)$forecast
    Mod=sqrt(phi^2+phiI^2)
    error=sigma*(1-Mod**(2*tAhead))
  }
  if(model == "BiAR")
  {
    mu = apply(y,1,mean)
    y1=y
    y2=y1[2,]
    y1=y1[1,]
    fore=BIARforecast(phiR=phi,phiI=phiI,y1=y1,y2=y2,t=st,tAhead=tAhead)$forecast
    sigma=matrix(apply(y,1,sd),nrow=2)
    Mod=sqrt(phi^2+phiI^2)
    error=sigma*(1-Mod**(2*tAhead))
  }
  upper=fore+qnorm(level/100)*error
  lower=fore-qnorm(level/100)*error
  fin=data.frame(tAhead=tAhead,forecast=fore,stderror=error,lowerCI=lower,upperCI=upper)
  return(fin)
}
