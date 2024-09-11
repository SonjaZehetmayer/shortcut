#fk....combination function resulting in (non-negative) test statistics
#p vector or matrix of p-values
#method combination test method; "fisher" or "stouffer".


library(matrixStats)

#rm(list = ls())



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

means<-function(p)
{
  if(is.matrix(p)) 
    return(1-rowSums(p)/ncol(p))#return(rep(1,length=ncol(p)))#
  return(1-mean(p))
}

fk(c(0.01,0.2),method="stouffer")
fk(c(0.01,0.2),func=means)
r<-matrix(c(0.1,0.2,0.5,0.6,0.1,0.1),3,byrow=T)
fk(r,"stouffer",func=means)
#Example
#fk(c(0.01,0.2),"fisher")
#fk(c(0.01,0.2),"stouffer")



#original from library(survcomp)
#cp <- pchisq(-2 * sum(log(p)), df = 2 * length(p), lower.tail = FALSE)

#original from library(survcomp)
#z <- qnorm(p, lower.tail = FALSE)
#cp <- pnorm(sum( z)/sqrt(length(p)), lower.tail = FALSE)
#1-cp


#simes (hommel)
# fk=function(p) if(is.matrix(p)) 
# {
#   m=ncol(p)
#   r=nrow(p)
#   s=(m/1:m)
#   pranks=rowRanks(p,ties.method="first")+matrix(rep(1:r-1,m),nrow=r)*m
#   return(1-colMins(t(matrix(as.vector(t(p))[order(as.vector(t(pranks)))], ncol=m,byrow=T))*s))
# } else
#     {m=length(p)
#       return(1-min(sort(p)*m/(1:m)))
#     }




create.null<-function(method="stouffer",alp,h,runs=10^5,eps=10^(-10),func=NULL) ## ad Null distribution
{
  cores<-detectCores() - 1 
  
  create_function <- function(f,eps) 
  {
    force(f)  # Ensure that f is evaluated when the factory is called
    function(x) 1 - f(x-eps) # eps is needed due to rounding errors
  }
  
  ## Compute the G[[k]] for k=1,...,h 
  #compute the null distribution for dimension h (??)
  #initialize G
  G=list(function (x) 1-x)
  
  pint<-function(p,method,alp,h)  #Algorithm 2 in Paper (see below, there it is repeated)
  { 
    p_sort=sort(p)
    l=length(p)
    L=1:l
    h=l
    K=L[1]
    pK=list(p_sort[1])
    while(h>=2)
    {
      K=c(K,L[h])
      k=length(K)
      ps=p_sort[K]
      TK=fk(ps,method,func) # use here the combination function. #TK is test statistic
      if(pK[[k-1]]>alp) TK=0
      pK[[k]]=G[[k]](TK)
      h=h-1
    }
    return(unlist(pK[l]))
  }


#null<-function(h,alp,runs,eps,method,cores)
#{
for(k in 2:h)
{
  p=matrix(runif(k*runs),nrow=runs) #matrix of random p-values under null
  TK=fk(p,method,func)#combination function - test stat=1-p (original)
  pranks=rowRanks(p,ties.method="first")+matrix(rep(1:runs-1,k),nrow=runs)*k
  #ranks per row + absolutes ranking... ergibt Ranking innerhalb Zeile aber Rankingzahl insgesamt.
  psort=matrix(as.vector(t(p))[order(as.vector(t(pranks)))], ncol=k,byrow=T) #nun sind die Zahlen pro Zeile sortiert
  pwis=matrix(psort[,-2],nrow=runs) #lass in der Matrix 2. groessten weg;
  rows_list <- split(pwis, row(pwis)) #erzeugt eine Liste
  result_list=mclapply(rows_list, pint, method=method,alp=alp,h=h,mc.cores = cores) #wendet nun combination function an. WIE KANN HIER G ANGEWENDET WERDEN,
  #result_list=mclapply(rows_list, sum,mc.cores = cores) #wendet nun combination function an. WIE KANN HIER G ANGEWENDET WERDEN, 
  #WENN ES ERST ERZEUGT WIRD?  --> wird oben definiert: G=list(function (x) 1-x)
  pKwis <- do.call(rbind, result_list) #jetzt wird wieder matrix draus.
  Tc=ifelse(pKwis <= alp,TK,0) #modified test stat
  f=ecdf(Tc) #kann mir jeweils Tc plotten, dann sehe ich alles.
  G[[k]]= create_function(f,eps)  #hier wird funktion erzeugt, die dann in pint angewendet wird (closure) - in jeden Schritt wird Liste G erweitert.
}
 return(G)
}


#create.null("stouffer",alp=0.05,h=5)
#null(5,0.05,10000,10^(-10),"stouffer",cores=39)
#null()



## Algorithm 2 in Zehetmayer, Koenig, Posch (2024)
## Computes a consonant test (p-value) for a specific intersection hypothesis iteratively
#shortcut
#
#
#p vector or matrix of p-values
#method combination test method; "fisher" or "stouffer".

pint<-function(p,method="stouffer",Gnull,alp,func=NULL)
{ 
  p_sort=sort(p)
  l=length(p)
  L=1:l
  h=l
  K=L[1]
  pK=list(p_sort[1])
  while(h>=2)
  {
    K=c(K,L[h])
    k=length(K)
    ps=p_sort[K]
    TK=fk(ps,method,func) # use here the combination function. #TK is test statistic
    if(pK[[k-1]]>alp) TK=0
    pK[[k]]=Gnull[[k]](TK)
    h=h-1
  }
  return(unlist(pK[l]))
}

#Example
#pint(c(0.1,0.01),"stouffer",GST,alp)
#pint(c(0.1,0.01),"fisher",GFI,alp)



### Algorithm 3 in Zehetmayer, Koenig, Posch (2024)
# Shortcut for testing m intersection hyp.
#eigentliche Funktion!  #calculates adjusted consonant p-values
cons.test = function(p,method="stouffer",Gnull,alp,func=NULL)
{
  m = length(p)
  p_adj = numeric(m)
  order_raw=order(p)
  for (i in 1:m) {
    sorted_p <- p[order_raw]
    ith_smallest <- sorted_p[i]
    # Get indices of all elements larger or equal to the i-th smallest element
    indices <- which(p >= ith_smallest)
    p_adj[i] = pint(p[indices],method,Gnull,alp,func) #brauch hier #G; Sollte ich vorher berechnen.
    if (i > 1) {
      p_adj[i] = max(p_adj[i], p_adj[i - 1])
    }
  }
  return(p_adj[order(order_raw)])
}


#Example
#p=c(.01,.1,.1)
#GST<-create.null("stouffer",0.05,5)
#cons.test(p,method="stouffer",GST,0.05)
#cons.test(c(0.1,0.01,0.1,0.1,0.1),"stouffer",GST,0.05)
#GFI<-create.null("fisher",0.05,5)
#cons.test(p,method="fisher",GFI,0.05)
#cons.test(c(0.1,0.01,0.1,0.1,0.1),"stouffer",GFI,0.05)

#GFI[[5]](1)
#hist(Tc)
#x<-environment(GFI[[2]])
#ls(x)
#get(f,x)
