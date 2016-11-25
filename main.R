rm(list=ls())
source('/Users/mbenigni/Public/repositories/IVCC/IVCC_Features.R')
source('simFunctions.R')
load('CJTC_M.Rdata')

alpha_seq=seq(0,1,by=.025)
replicates=100
dens=trans=dia=dyad=triad=as.data.frame(array(NA,dim=c(length(alpha_seq),replicates)),
                                        row.names=as.character(alpha_seq))

min.clust.par=20
                                        

s=proc.time()
obj=simIteration(nodeList=a$userID,g=m,min.clust=min.clust.par,alpha=.25,replace=TRUE,add=TRUE)
d=(proc.time()-s)

d*replicates*length(alpha_seq)/3600

for(i in 1:length(alpha)){
  for(j in 1:replicates){
    obj=simIteration(nodeList=a$userID,g=m,min.clust=min.clust.par,alpha=alpha_seq[i],replace=TRUE,add=TRUE)
    dens[i,j]=obj$dense
    trans[i,j]=obj$trans
    dia[i,j]=obj$dia
    dyad[i,j]=obj$dyad
    triad[i,j]=obj$triad
  }
}

save(alpha_seq,min.clust.par,dens,trans,dia,dyad,triad,file='VirtualExperiment1_replacement.RData')