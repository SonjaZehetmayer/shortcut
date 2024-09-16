#rm(list = ls())


#Example 1
runs2=10^6 #number of simulation runs
h=5 #upper bound for dimensions in null distribution
alp=0.05 #significance level
cores<-detectCores() - 1 #number of cores for parallelization
GFI<-create.null("fisher",alp,k)

cons.test(c(0.01,0.05),method="fisher",GFI,alp)
# 0.00456 0.05000
cons.test(c(0.01,0.1),method="fisher",GFI,alp)
#0.00847 0.10000
cons.test(c(0.01,0.05,0.1,0.2,0.3),method="fisher",GFI,alp)
#0.00449 1.00000 1.00000 1.00000 1.00000




#Example 2
runs2=10^6 #number of simulation runs
h=5 #upper bound for dimensions in null distribution
alp=0.05 #significance level
cores<-detectCores() - 1 #number of cores for parallelization


means<-function(p)
{
  if(is.matrix(p))
    return(1-rowSums(p)/ncol(p))
  return(1-mean(p))
}


GM<-create.null(alp=0.05,h=5,func=means)


cons.test(p=c(0.01,0.05),Gnull=GM,alp=0.05,func=means)
#0.0022 0.0500
cons.test(p=c(0.01,0.1),Gnull=GM,alp=0.05,func=means)
# 0.00647 0.10000
cons.test(p=c(0.01,0.05,0.1,0.2,0.3),Gnull=GM,alp=alp,func=means)
#0.00087 0.00529 1.00000 1.00000 1.00000
