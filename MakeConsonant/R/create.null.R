#' Generation of null distribution of combination function to compute empirical p-values
#'
#' @param method either "Stouffer" or "Fisher" combination test
#' @param func userdefined function. Is set to zero if no user defined function is given
#' @param alp significance level
#' @param h number of dimensions of null distribution
#' @param runs number of random variables under null hypothesis
#' @param eps auxiliary parameter
#' @returns list of null distributions
#' @export
#'


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


