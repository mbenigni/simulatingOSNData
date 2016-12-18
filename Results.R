rm(list=ls())
load("CJTC_M.RData")
library(igraph)



M=graph.empty(n=length(a$userID),directed=FALSE)
V(M)$name=as.character(a$userID)
M=add_edges(M,edges=rbind(as.character(m[,1]),as.character(m[,2])))

gdens=graph.density(M)
gtrans=transitivity(M)
gdia=diameter(M)
gdyad=count_motifs(M,2)
gtriad=count_motifs(M,3)

membershipM=as.vector(membership(cluster_louvain(M)))

cDF=as.data.frame(ftable(membershipM))
groupsM=sum(cDF$Freq>1)
usersM=sum(cDF$Freq[cDF$Freq>1])
rm(M,a,m)

#Virtual Experiment 2: subtract
load('VirtualExperiment2_subtract.RData')

png('figures/VE2S_users_percent.png',
    width=1000,height=800,res=100)
plot(alpha_seq,apply(users,1,mean),type='l',
     ylab='co-clustered users',xlab='Percent of Removed Edges',
     main='Virtual Experiment 2B (Randomly Added Edges):\n percent edges removed vs co-clustered users')
     for(i in 1:dim(users)[2]){
       points(alpha_seq,users[,i],pch='.',col=alpha('gray',alpha=.5))
     }

dev.off()

png('figures/VE2S_4panel.png',
    width=1400,height=1000,res=100)
par(mfcol=c(2,2))
plot(apply(trans,1,mean),apply(users,1,mean),type='l',
     ylab='co-clustered users',xlab='Transitivity',
     main='Transitivity')
for(i in 1:dim(users)[2]){
  points(trans[,i],users[,i],pch='.',col=alpha('gray',alpha=.75))
}

lines(apply(trans,1,mean),apply(users,1,mean),lwd=2)
points(apply(trans,1,mean)[1],apply(users,1,mean)[1],pch=21,cex=1.5,bg='grey')

plot(apply(dia,1,mean),apply(users,1,mean),type='l',
     ylab='co-clustered users',xlab='Diameter',
     main='Diameter')
for(i in 1:dim(users)[2]){
  points(dia[,i],users[,i],pch='.',col=alpha('gray',alpha=.75))
}

lines(apply(dia,1,mean),apply(users,1,mean),lwd=2)
points(apply(dia,1,mean)[1],apply(users,1,mean)[1],pch=21,cex=1.5,bg='grey')

plot(apply(triad,1,mean),apply(users,1,mean),type='l',
     ylab='co-clustered users',xlab='Triad Count',
     main='Triad Count')
for(i in 1:dim(users)[2]){
  points(triad[,i],users[,i],pch='.',col=alpha('gray',alpha=.75))
}

lines(apply(triad,1,mean),apply(users,1,mean),lwd=2)
points(apply(triad,1,mean)[1],apply(users,1,mean)[1],pch=21,cex=1.5,bg='grey')

plot(apply(dyad,1,mean),apply(users,1,mean),type='l',
     ylab='co-clustered users',xlab='Dyad Count',
     main='Dyad Count')
for(i in 1:dim(users)[2]){
  points(dyad[,i],users[,i],pch='.',col=alpha('gray',alpha=.75))
}

lines(apply(dyad,1,mean),apply(users,1,mean),lwd=2)
points(apply(dyad,1,mean)[1],apply(users,1,mean)[1],pch=21,cex=1.5,bg='grey')
dev.off()




#Virtual Experiment 2: add
load('VirtualExperiment2_add.RData')

png('figures/VE2A_users_percent.png',
    width=1000,height=800,res=100)
plot(alpha_seq,apply(users,1,mean),type='l',
     ylab='co-clustered users',xlab='Percent Added Edges',
     main='Virtual Experiment 2B (Randomly Added Edges):\n percent edges added vs co-clustered users')
for(i in 1:dim(users)[2]){
  points(alpha_seq,users[,i],pch='.',col=alpha('gray',alpha=.5))
}

dev.off()

png('figures/VE2A_4panel.png',
    width=1400,height=1000,res=100)
par(mfcol=c(2,2))
plot(apply(trans,1,mean),apply(users,1,mean),type='l',
     ylab='co-clustered users',xlab='Transitivity',
     main='Transitivity')
for(i in 1:dim(users)[2]){
  points(trans[,i],users[,i],pch='.',col=alpha('gray',alpha=.75))
}

lines(apply(trans,1,mean),apply(users,1,mean),lwd=2)
points(apply(trans,1,mean)[1],apply(users,1,mean)[1],pch=21,cex=1.5,bg='grey')

plot(apply(dia,1,mean),apply(users,1,mean),type='l',
     ylab='co-clustered users',xlab='Diameter',
     main='Diameter')
for(i in 1:dim(users)[2]){
  points(dia[,i],users[,i],pch='.',col=alpha('gray',alpha=.75))
}

lines(apply(dia,1,mean),apply(users,1,mean),lwd=2)
points(apply(dia,1,mean)[1],apply(users,1,mean)[1],pch=21,cex=1.5,bg='grey')

plot(apply(triad,1,mean),apply(users,1,mean),type='l',
     ylab='co-clustered users',xlab='Triad Count',
     main='Triad Count')
for(i in 1:dim(users)[2]){
  points(triad[,i],users[,i],pch='.',col=alpha('gray',alpha=.75))
}

lines(apply(triad,1,mean),apply(users,1,mean),lwd=2)
points(apply(triad,1,mean)[1],apply(users,1,mean)[1],pch=21,cex=1.5,bg='grey')

plot(apply(dyad,1,mean),apply(users,1,mean),type='l',
     ylab='co-clustered users',xlab='Dyad Count',
     main='Dyad Count')
for(i in 1:dim(users)[2]){
  points(dyad[,i],users[,i],pch='.',col=alpha('gray',alpha=.75))
}

lines(apply(dyad,1,mean),apply(users,1,mean),lwd=2)
points(apply(dyad,1,mean)[1],apply(users,1,mean)[1],pch=21,cex=1.5,bg='grey')
dev.off()




