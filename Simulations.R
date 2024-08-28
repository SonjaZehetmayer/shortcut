#rm(list = ls())

#initialization
#source("fk.R")
#source("pint.R")

#source(null.R)
#source(constest.R)

h=5   #gibt an, wieviele Dimensionen fuer Nullverteilung
alp=0.05
runs=10^5
eps=10^(-10)
method="fisher"
cores=detectCores() - 1 #nr of cores to use

source("Shortcut.R")



#example 1
#pint(c(0.1,0.01),"stouffer")

#example 2
#pint(c(0.1,0.01),"fisher")



##############################################################################
###
### Simulation 1: Type 1 error rate of the intersection hypothesis test 
### under the null hypothesis

runs2=10^6# simulation runs
k=2 #number of hypotheses in the intersection

ps=matrix(runif(k*runs2),nrow=runs2)
rows_list <- split(ps, row(ps))
result_list=mclapply(rows_list, pint,method=method, mc.cores = cores) #combination test
pad <- do.call(rbind, result_list)
mean(pad<alp) #adjusted p-value


#1-pnorm(qnorm(1-alp/3)-2.6525)
#1-pnorm(qnorm(1-alp/5)-2.851)




### Simulation 2: consonant multiple test with shortcut
runs2=10^6 # simulation runs
# vector of effect sizes, maximum length 
# (= number of hypotheses) is h
#mu=c(0,2.5,2.5,2.5)#c(0,2.48998,2.48998)#2.851 
mu=c(0,0,2.851,2.851,2.851)#
#mu=c(0,0)#,0,0,0)
k=length(mu) 
ps=matrix(1-pnorm(rnorm(k*runs2,mean=mu)),nrow=runs2,byrow=T)
rows_list <- split(ps, row(ps))
result_list=mclapply(rows_list,cons_test , method=method,mc.cores = cores) #adjusted p-value
pad <- do.call(rbind, result_list)
pow=colMeans(pad<alp)
any_rej=mean(1-rowProds((pad>alp))) # any rejection (including nulls)
FWER=ifelse(length(which(mu==0))>1,mean(1-rowProds((pad[,which(mu==0)]>alp))),
            mean(pad[,which(mu==0)]<alp))
pow
any_rej
FWER
