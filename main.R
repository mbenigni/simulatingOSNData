rm(list=ls())
source('~/Public/repositories/IVCC/IVCC_Features.R')
source('simFunctions_working.R')
load('CJTC_M.RData')

M=graph.empty(n=length(a$userID),directed=FALSE)
V(M)$name=as.character(a$userID)
M=add_edges(M,edges=rbind(as.character(m[,1]),as.character(m[,2])))
membershipM=as.vector(membership(cluster_louvain(M)))

cDF=as.data.frame(ftable(membershipM))
groupsM=sum(cDF$Freq>1)

alpha_seq=seq(.025,.5,by=.025)
replicates=2
dens=trans=dia=dyad=triad=groups=users=as.data.frame(array(NA,dim=c(length(alpha_seq),replicates)),
                                        row.names=as.character(alpha_seq))

min.clust.par=100
for(i in 1:length(alpha_seq)){
  for(j in 1:replicates){
    obj=simIteration(nodeList=a$userID,g=m,min.clust=min.clust.par,alpha=alpha_seq[i],replace=FALSE,add=TRUE,membershipM)
    dens[i,j]=obj$dens
    trans[i,j]=obj$trans
    dia[i,j]=obj$dia
    dyad[i,j]=obj$dyad
    triad[i,j]=obj$triad
    users[i,j]=obj$users
    groups[i,j]=obj$groups
  }
}

save(alpha_seq,min.clust.par,dens,trans,dia,dyad,triad,file='VirtualExperiment2_add.RData')

alpha_seq=seq(.025,.5,by=.025)
replicates=100
dens=trans=dia=dyad=triad=as.data.frame(array(NA,dim=c(length(alpha_seq),replicates)),
                                        row.names=as.character(alpha_seq))

for(i in 1:length(alpha_seq)){
  for(j in 1:replicates){
    obj=simIteration(nodeList=a$userID,g=m,min.clust=min.clust.par,alpha=alpha_seq[i],replace=FALSE,add=FALSE)
    dens[i,j]=obj$dens
    trans[i,j]=obj$trans
    dia[i,j]=obj$dia
    dyad[i,j]=obj$dyad
    triad[i,j]=obj$triad
    users[i,j]=obj$users
    groups[i,j]=obj$groups
  }
}

save(alpha_seq,min.clust.par,dens,trans,dia,dyad,triad,file='VirtualExperiment2_subtract.RData')


alpha_seq=seq(.025,1,by=.025)
replicates=100
dens=trans=dia=dyad=triad=as.data.frame(array(NA,dim=c(length(alpha_seq),replicates)),
                                        row.names=as.character(alpha_seq))

min.clust.par=20
for(i in 1:length(alpha_seq)){
  for(j in 1:replicates){
    obj=simIteration(nodeList=a$userID,g=m,min.clust=min.clust.par,alpha=alpha_seq[i],replace=TRUE,add=TRUE)
    dens[i,j]=obj$dens
    trans[i,j]=obj$trans
    dia[i,j]=obj$dia
    dyad[i,j]=obj$dyad
    triad[i,j]=obj$triad
    users[i,j]=obj$users
    groups[i,j]=obj$groups
  }
}

save(alpha_seq,min.clust.par,dens,trans,dia,dyad,triad,file='VirtualExperiment1_replacement.RData')
