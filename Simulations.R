#rm(list = ls())





##############################################################################
###
### Simulation 1: Type 1 error rate of the intersection hypothesis test 
### under the null hypothesis

#source("Shortcut.R")

runs2=10^6# simulation runs
k=5 #number of hypotheses in the intersection
alp=0.05
cores<-detectCores() - 1


#Create null distribution
#ps=matrix(runif(k*runs2),nrow=runs2)
GST<-create.null("stouffer",alp,k)
#rows_list <- split(ps, row(ps))
#result_list=mclapply(rows_list, pint,method="stouffer",GST, alp,mc.cores = cores) #combination test
#pad <- do.call(rbind, result_list)
#mean(pad<alp) #adjusted p-value

GFI<-create.null("fisher",alp,k)
#result_list=mclapply(rows_list, pint,method="fisher",GFI, alp,mc.cores = cores) #combination test
#pad <- do.call(rbind, result_list)
#mean(pad<alp) #adjusted p-value

p<-c(0.01,0.05)
cons.test(c(0.01,0.05),method="fisher",GFI,alp)
cons.test(c(0.01,0.1),method="fisher",GFI,alp)

#1-pnorm(qnorm(1-alp/3)-2.6525)
#1-pnorm(qnorm(1-alp/5)-2.851)




### Simulation 2: consonant multiple test with shortcut
# vector of effect sizes, maximum length 
# (= number of hypotheses) is h
#mu=c(0,0,2.851,2.851,2.851)#
mu=c(2.8,0,0)#,0)#,0,0,0)
set.seed(2)
#ps=matrix(1-pnorm(rnorm(length(mu)*runs2,mean=mu)),nrow=runs2,byrow=T)

ps=matrix(1-pnorm(rnorm(length(mu)*runs2,mean=mu)),nrow=runs2,byrow=T)
rows_list <- split(ps, row(ps))

#Fisher
result_list=mclapply(rows_list,cons.test , method="fisher",GFI,alp,mc.cores = cores) #adjusted p-value
pad <- do.call(rbind, result_list)

#Fisher.closed
#eigentliches closed program
#provides all combinations of hypotheses (1:n)
all_comb <- function(n){
  do.call(rbind,
          lapply(X = 1:n, FUN = function(i){
            t(combn(x = 1:n, m = i, FUN = function(x) c(x, rep(0,n-i))))
          })
  )
}

#schnelleres Programm
all_comb_fast <- function(n) {
  # Anzahl der Zeilen vorab berechnen
  total_rows <- sum(choose(n, 1:n))
  
  # Leere Matrix mit der richtigen Größe vorab erstellen
  result <- matrix(0, nrow = total_rows, ncol = n)
  
  # Zähler für die Zeilenposition
  row_index <- 1
  
  # Für jede Kombination von Länge i (1 bis n)
  for (i in 1:n) {
    combs <- combn(1:n, i)
    num_combs <- ncol(combs)
    
    # Fülle jede Kombination auf die restlichen Positionen mit Nullen auf
    result[row_index:(row_index + num_combs - 1), 1:i] <- t(combs)
    
    # Aktualisiere den Zeilenindex
    row_index <- row_index + num_combs
  }
  
  return(result)
}

#aus library survcomp
combine.test<-function (p, weight, method = c("fisher", "z.transform", "logit"), hetero = FALSE, na.rm = FALSE) 
{
  if (hetero) {
    stop("function to deal with heterogeneity is not implemented yet!")
  }
  method <- match.arg(method)
  na.ix <- is.na(p)
  if (any(na.ix) && !na.rm) {
    stop("missing values are present!")
  }
  if (all(na.ix)) {
    return(NA)
  }
  p <- p[!na.ix]
  k <- length(p)
  if (k == 1) {
    return(p)
  }
  if (missing(weight)) {
    weight <- rep(1, k)
  }
  switch(method, fisher = {
    cp <- pchisq(-2 * sum(log(p)), df = 2 * k, lower.tail = FALSE)
  }, z.transform = {
    z <- qnorm(p, lower.tail = FALSE)
    cp <- pnorm(sum(weight * z)/sqrt(sum(weight^2)), lower.tail = FALSE)
  }, logit = {
    tt <- (-sum(log(p/(1 - p))))/sqrt(k * pi^2 * (5 * k + 
                                                    2)/(3 * (5 * k + 4)))
    cp <- pt(tt, df = 5 * k + 4, lower.tail = FALSE)
  })
  return(cp)
}


fisher.closed<-function(p,alpha)
{
  m<-length(p)
  #allcombs<-lapply(1:m,function(x) {combinations(m,x)})
  #bincombinations<-bincomb(m)[-1,]
  allcombs<-all_comb_fast(m)
  bincombinations<-matrix(0,ncol=m,nrow=nrow(allcombs))
  
  for (i in 1:nrow(bincombinations))
    bincombinations[i,c(allcombs[i,])]<-TRUE
  #bincombinations[cbind(1:nrow(bincombinations),t(allcombs))]<-TRUE
  
  #Vektor fuer abgelehnte Hyp and p-values
  rejected<-rep(0,(2^m - 1))
  pall<-vector(length=(2^m - 1))  #rep(NA,(2^m - 1))#
  pall[1:m]<-p
  rejected[1:m]<-p<=alpha
  
  if (sum(rejected[1:m])>0) #nur wenn mind 1 Hyp<=alpha geht es weiter
  {
    for (j in (m+1):(2^m - 1))  
      pall[j]<-combine.test(p[allcombs[j,]],method="fisher")#combine.test(selection,method="fisher")
  }
  
  rres<-(bincombinations*pall)<=alpha
  result<-apply(rres,2,sum)!=(2^m - 1)  #FALSE..H1, TRUE..H0
  result
}

fisher.closed.optimized <- function(p, alpha) {
  m <- length(p)
  allcombs <- all_comb_fast(m)  # Verwende die optimierte all_comb Funktion
  bincombinations <- matrix(0, ncol = m, nrow = nrow(allcombs))
  
  # Direkte Zuordnung ohne Schleife
  bincombinations[cbind(1:nrow(bincombinations), t(allcombs))] <- TRUE
  
  # Vektor für abgelehnte Hypothesen und p-Werte
  rejected <- rep(0, (2^m - 1))
  pall <- rep(NA, (2^m - 1))
  
  # Setze die ersten m-Werte von p direkt
  pall[1:m] <- p
  rejected[1:m] <- p <= alpha
  
  # Überprüfe, ob mindestens eine Hypothese abgelehnt wurde
  if (sum(rejected[1:m]) > 0) {
    # Wende die combine.test Funktion auf alle Kombinationen an, die noch nicht berechnet wurden
    selected_combinations <- (m+1):(2^m - 1)
    pall[selected_combinations] <- apply(allcombs[selected_combinations, , drop = FALSE], 1, function(comb) {
      combine.test(p[comb > 0], method = "fisher")
    })
  }
  
  # Überprüfe, welche Kombinationen abgelehnt wurden
  rres <- (bincombinations * pall) <= alpha
  result <- apply(rres, 2, sum) != (2^m - 1)  # FALSE..H1, TRUE..H0
  
  return(result)
}


fisher.closed.optimized(c(.1,.2),0.05)

pfcl<-pad

start.time <- Sys.time()

for (i in  1:dim(ps)[1])
pfcl[i,]<-fisher.closed.optimized(ps[i,],0.05)

end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken



pow=colMeans(pad<alp)
pow1<-colMeans(1-pfcl)
pow
pow1


any_rej=mean(1-rowProds((pad>alp))) # any rejection (including nulls)
FWER=ifelse(length(which(mu==0))>1,mean(1-rowProds((pad[,which(mu==0)]>alp))),
            mean(pad[,which(mu==0)]<alp))
pow
any_rej
FWER



result_listS=mclapply(rows_list,cons.test , method="stouffer",GST,alp,mc.cores = cores) #adjusted p-value
padS <- do.call(rbind, result_listS)
powS=colMeans(padS<alp)
any_rejS=mean(1-rowProds((padS>alp))) # any rejection (including nulls)
FWERS=ifelse(length(which(mu==0))>1,mean(1-rowProds((padS[,which(mu==0)]>alp))),
            mean(padS[,which(mu==0)]<alp))
powS
any_rejS
FWERS




######
p1<-seq(0,1,0.01)
p2<-p1
ps<-expand.grid(p1,p2)
ps<-as.matrix(ps)
rows_list <- split(ps, row(ps))

result_listS=mclapply(rows_list,cons.test , method="stouffer",GST,alp,mc.cores = cores) #adjusted p-value
padS <- do.call(rbind, result_listS)
plot(padS[,1],padS[,2])

powS=colMeans(padS<alp)
any_rejS=mean(1-rowProds((padS>alp))) # any rejection (including nulls)
FWERS=ifelse(length(which(mu==0))>1,mean(1-rowProds((padS[,which(mu==0)]>alp))),
             mean(padS[,which(mu==0)]<alp))
powS
any_rejS
FWERS

#mu=c(0,0,0,0,0)


###############user defined Example 3
GM<-create.null(alp=0.05,h=5,func=means)
GM[[1]](c(.1,.2))
fk(c(0.01,0.05),func=means)
cons.test(p=c(0.01,0.05),Gnull=GM,alp=0.05,func=means)
cons.test(p=c(0.01,0.1),Gnull=GM,alp=0.05,func=means)
GM[[2]](c(.1,.2))
#[1] 0.09611 0.09611



#Example 3
#eigentliches closed program

