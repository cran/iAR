#' Pairing two irregularly observed time series
#'
#' Pairing the observational times of two irregularly observed time series
#'
#' @param lc1 data frame with three columns corresponding to the first irregularly observed time series. The columns must be ordered as follow: First the observational times, second the measures of each time, and third the measurement errors.
#' @param lc2 data frame with three columns corresponding to the second irregularly observed time series. The columns must be ordered as follow: First the observational times, second the measures of each time, and third the measurement errors.
#' @param tol tolerance parameter. Minimum time gap to consider that two observations have measured at different times.
#'
#' @return A list with the following components:
#' \itemize{
#' \item{n}{ Number of observations paired by their observational times.}
#' \item{par}{Data Frame with the paired datasets.}
#' }
#' @export
#' @references
#' \insertRef{Elorrieta_2021}{iAR}
#'
#' @seealso
#'
#' \code{\link{cvnovag}}, \code{\link{cvnovar}}, \code{\link{BIARkalman}}
#'
#'
#' @examples
#' data(cvnovag)
#' data(cvnovar)
#' pargr=pairingits(cvnovag,cvnovar,tol=0.1)
pairingits<-function(lc1,lc2,tol=0.1)
{
  t1=lc1[,1]
  t2=lc2[,1]
  A=cbind(c(t1,t2),c(rep(1,length(t1)),rep(2,length(t2))),c(1:length(t1),1:length(t2)))
  A=A[order(A[,1]),]
  fin=NULL
  i=2
  while(i<=dim(A)[1])
  {
    if(A[i-1,2]!=A[i,2])
    {
      dt=diff(c(A[i,1],A[i-1,1]))
      if(abs(dt)<tol)
      {
        if(A[i,2]>A[i-1,2])
          par=c(lc1[A[i-1,3],],lc2[A[i,3],])
        if(A[i,2]<A[i-1,2])
          par=c(lc1[A[i,3],],lc2[A[i-1,3],])
        i=i+1
      }
      else
      {
        if(A[i-1,2]==1)
          par=c(lc1[A[i-1,3],],rep(NA,3))
        if(A[i-1,2]==2)
          par=c(rep(NA,3),lc2[A[i-1,3],])
      }
    }
    else
    {
      if(A[i-1,2]==1 & A[i,2]==1)
        par=c(lc1[A[i-1,3],],rep(NA,3))
      if(A[i-1,2]==2)
        par=c(rep(NA,3),lc2[A[i-1,3],])
    }
    fin=rbind(fin,c(unlist(par)))
    i=i+1
  }
  i=i-1
  if(i==dim(A)[1])
  {
    if(A[i-1,2]==A[i,2] | (A[i-1,2]!=A[i,2] & abs(dt)>=tol))
    {
      if(A[i,2]==1)
        par=c(lc1[A[i,3],],rep(NA,3))
      if(A[i,2]==2)
        par=c(rep(NA,3),lc2[A[i,3],])
      fin=rbind(fin,c(unlist(par)))
    }
  }
  prop=dim(na.omit(fin))[1]
  length(na.omit(unique(fin[,1])))
  length(na.omit(unique(fin[,4])))
  return(list(n=prop,par=fin))
}
