#' combination function resulting in (non-negative) test statistics
#'
#' @param p vector or matrix of p-values
#' @param method either "Stouffer" or "Fisher" combination test
#' @param func userdefined function. Is set to zero if no user defined function is given
#' @returns List of test statistics for each dimension
#' @export
#'
fk=function(p,method,func=NULL)
{
  if (!is.null(func))
  {
    if(is.matrix(p)){}
    return(func(p))
  }

  else
  {
    if (method=="fisher") #product test (Fisher) #results in 1-p-value
    {
      if(is.matrix(p))
        return(pchisq(-2 * rowSums(log(p)), df = 2 * ncol(p))) else
          return(pchisq(-2 * sum(log(p)), df = 2 * length(p)))
    }
    if (method=="stouffer")   #inverse normal (Stouffer) #results in 1-p-value
    {
      if(is.matrix(p))
        return(pnorm(rowSums(qnorm(1-p))/sqrt(ncol(p)))) else
          return(pnorm(sum(qnorm(1-p))/sqrt(length(p))))
    }
  }
}
