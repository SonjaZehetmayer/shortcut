#fk....combination function resulting in (non-negative) test statistics
#p vector or matrix of p-values
#method combination test method; "fisher" or "stouffer".


library(matrixStats)





fk=function(p,method)
{
  if (method=="fisher") #product test (Fisher) #results in 1-p-value
  {
    if(is.matrix(p))
      return(pchisq(-2 * rowSums(log(p)), df = 2 * length(p))) else
        return(pchisq(-2 * sum(log(p)), df = 2 * length(p))) 
  }
  if (method=="stouffer")   #inverse normal (Stouffer) #results in 1-p-value
  {
    if(is.matrix(p))
      return(pnorm(rowSums(qnorm(1-p))/sqrt(ncol(p)))) else
        return(pnorm(sum(qnorm(1-p))/sqrt(length(p))))
  }
}

#Example
fk(c(0.01,0.2),"fisher")
fk(c(0.01,0.2),"stouffer")



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



## Algorithm 2 in Zehetmayer, Koenig, Posch (2024)
## Computes a consonant test (p-value) for a specific intersection hypothesis iteratively
#shortcut
#
#
#p vector or matrix of p-values
#method combination test method; "fisher" or "stouffer".


pint=function(p,method)
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
    TK=fk(ps,method) # use here the combination function. #TK is test statistic
    if(pK[[k-1]]>alp) TK=0
    pK[[k]]=G[[k]](TK)
    h=h-1
  }
  return(unlist(pK[l]))
}

#Example
#pint(c(0.1,0.01),"stouffer")



# ad Null distribution

cores<-detectCores() - 1 

create_function <- function(f,eps) {
  force(f)  # Ensure that f is evaluated when the factory is called
  function(x) 1 - f(x-eps) # eps is needed due to rounding errors
}


## Compute the G[[k]] for k=1,...,h 
#compute the null distribution for dimension h (??)


###############################WAS IST G?????????????


#initialize G
G=list(function (x) 1-x)

#null<-function(h,alp,runs,eps,method,cores)
#{
for(k in 2:h)
{
  p=matrix(runif(k*runs),nrow=runs) #matrix of random p-values under null
  TK=fk(p,method)#combination function - test stat=1-p (original)
  pranks=rowRanks(p,ties.method="first")+matrix(rep(1:runs-1,k),nrow=runs)*k
  #ranks per row + absolutes ranking... ergibt Ranking innerhalb Zeile aber Rankingzahl insgesamt.
  psort=matrix(as.vector(t(p))[order(as.vector(t(pranks)))], ncol=k,byrow=T) #nun sind die Zahlen pro Zeile sortiert - ist schneller als Schleife, sortiert Matrix zeilenweise
  pwis=matrix(psort[,-2],nrow=runs) #lass in der Matrix 2. groessten weg;
  rows_list <- split(pwis, row(pwis)) #erzeugt eine Liste wegen mcapply
  result_list=mclapply(rows_list, pint, method=method,mc.cores = cores) #wendet nun combination function an. WIE KANN HIER G ANGEWENDET WERDEN,
  #result_list=mclapply(rows_list, sum,mc.cores = cores) #wendet nun combination function an. WIE KANN HIER G ANGEWENDET WERDEN, 
  #WENN ES ERST ERZEUGT WIRD?  --> wird oben definiert: G=list(function (x) 1-x)
  pKwis <- do.call(rbind, result_list) #jetzt wird wieder matrix draus.
  Tc=ifelse(pKwis <= alp,TK,0) #modified test stat
  f=ecdf(Tc)
  G[[k]]= create_function(f,eps)  #hier wird funktion erzeugt, die dann in pint angewendet wird (closure) - in jeden Schritt wird Liste G erweitert.
}
# return(G)
#}


#null(5,0.05,10000,10^(-10),"stouffer",cores=39)
#null()



### Algorithm 3 in Zehetmayer, Koenig, Posch (2024)
# Shortcut for testing m intersection hyp.
#eigentliche Funktion!  #calculated adjusted consonant p-values
cons_test = function(p_raw,method)
{
  m = length(p_raw)
  p_adj = numeric(m)
  order_raw=order(p_raw)
  for (i in 1:m) {
    sorted_p_raw <- p_raw[order_raw]
    ith_smallest <- sorted_p_raw[i]
    # Get indices of all elements larger or equal to the i-th smallest element
    indices <- which(p_raw >= ith_smallest)
    p_adj[i] = pint(p_raw[indices],method) #brauch hier #G; Sollte ich vorher berechnen.
    if (i > 1) {
      p_adj[i] = max(p_adj[i], p_adj[i - 1])
    }
  }
  return(p_adj[order(order_raw)])
}


#Example
p_raw=c(.01,.1,.1)
cons_test(p_raw,method="stouffer")
cons_test(p_raw,method="fisher")


