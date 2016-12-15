# What if we build an iterative model? One that continues to sub-group louvain groups until we meet 
# some minimal group size threshold.  Synthetic Nested Louvain Block Models.
# keep M as an edge list, just continue to work the diagonal in nested fashion
clusterList=function(vector,min.clust){
  cDF=as.data.frame(ftable(vector))
  names(cDF)=c('cluster','count')
  if(any(cDF$count>=min.clust)){
    cDF=cDF[order(cDF$count,decreasing=TRUE),]
    cDF$cluster=as.numeric(as.character(cDF$cluster))
    smallClusters=cDF$cluster[cDF$count<min.clust]
    cDF=cDF[cDF$count>=min.clust,]
    return(list(cDF=cDF,smallClusters=smallClusters)) 
  }else{
    return(list(cDF=NULL,smallClusters=cDF$cluster))
  }
  
}

getEdges=function(nodesA,nodesB,g,alpha,replace=TRUE,add=TRUE){
  indecies=(g$Source %in% nodesA & g$Target %in% nodesB)|(g$Source %in% nodesB & g$Target %in% nodesA)
  edges=g[indecies,]
  ecount=sum(indecies)
  n=ceiling(alpha*ecount)
  random.edges=as.data.frame(cbind(sample(nodesA,n,replace=TRUE),sample(nodesB,n,replace=TRUE)),stringsAsFactors=FALSE)
  
  if(replace){
    edges[sample(1:ecount,n),]=random.edges
  }else{
    if(add){names(random.edges)=c('Source','Target')
            edges=rbind(edges,random.edges)
            }else{edges=edges[sample(1:ecount,ecount-n),]}
            }
  return(list(edges=edges,ind=indecies))
}



get_clusters=function(min.clust,g,nodeList,alpha,replace,add){
    
  g=g[g$Source %in% nodeList & g$Target %in% nodeList,]
  edge.count=dim(g)[1]
  G=graph.edgelist(cbind(as.character(g$Source),as.character(g$Target)),directed=FALSE)
  
  #GET CLUSTERS
  nDF=data.frame(node=V(G)$name,cluster=as.vector(membership(cluster_louvain(G))),stringsAsFactors=FALSE)
  cL=clusterList(vector=nDF$cluster,min.clust=min.clust)

  smallNodes=nDF$node[nDF$cluster %in% cL[[2]]]
  
  clusters=lapply(cL$cDF$cluster,function(x) nDF$node[nDF$cluster==x])
  temp_small_edges=sum(g$Source %in% smallNodes | g$Target %in% smallNodes)
  # NOTE NEED TO ACCOUNT FOR WHEN THERE IS ONLY ONE CLUSTER
  #update edges
  
  # edges from small clusters
  eobj=getEdges(nodesA=smallNodes,nodesB=smallNodes,g,alpha=alpha,replace=replace,add=add)
  
  if(length(clusters)<=1){
        return(list(edges_out=eobj$edges,clusters='STOP')) 
  }else{
        if(replace){g[eobj$ind,]=eobj$edges}else{
          g=g[!eobj$ind,]
          g=rbind(g,eobj$edges)
        }
        for(i in 1:length(clusters)){
          eobj=getEdges(nodesA=smallNodes,nodesB=clusters[[i]],g,alpha=alpha,replace=replace,add=add)
          if(replace){g[eobj$ind,]=eobj$edges}else{
            g=g[!eobj$ind,]
            g=rbind(g,eobj$edges)
          }
        }
        
        
        # off diagonal edges
        for(i in length(clusters):2){
          for(j in (i-1):1){
            eobj=getEdges(nodesA=clusters[[j]],nodesB=clusters[[i]],g,alpha=alpha,replace=replace,add=add)
            if(replace){g[eobj$ind,]=eobj$edges}else{
              g=g[!eobj$ind,]
              g=rbind(g,eobj$edges)
            }
          }
        }
 
        for(i in 1:length(clusters)){
          indecies=g$Source %in% clusters[[i]] & g$Target %in% clusters[[i]]
          g=g[!indecies,]
        } 
        
        return(list(edges_out=g,clusters=clusters))}
}

# there is a problem within the while loop in the obj$clusters

simIteration=function(nodeList,g,min.clust,alpha,replace,add,membershipM){
  
  clustOBJ=get_clusters(min.clust,g,nodeList,alpha,replace,add)
  
  queue=clustOBJ$clusters
  edges=clustOBJ$edges
  
  # should be able to parrelelize this by queu length
  
  while(length(queue)>0){
    nodes=queue[[1]]
    queue[[1]]=NULL
    if(length(nodes)>0){
      obj=get_clusters(min.clust,g,nodeList=nodes,alpha,replace,add)
                      
      edges=rbind(edges,obj$edges)
      
      if(length(obj$clusters)>0){
        n=length(queue)
        
        for(i in 1:length(obj$clusters)){
          if(length(obj$clustList)>0){
            queue[[n+i]]=obj$clustList[[i]]}else{
              queue[[n+i]]=NULL  
            }
        }
      }} 
    
  }
  
  
  S=graph.empty(n=length(a$userID),directed=FALSE)
  V(S)$name=as.character(a$userID)
  S=add_edges(S,edges=rbind(as.character(edges[,1]),as.character(edges[,2])))
  tl=testLouvain(S,membership=membershipM)
  
  
  
  
  
  return(list(dens=graph.density(S),
              trans=transitivity(S),
              dia=diameter(S),
              dyad=count_motifs(S,2),
              triad=count_motifs(S,3),
              users=tl$users,
              groups=tl$groups))}


testLouvain=function(G,membership){
  t=as.data.frame(ftable(paste(membership,as.vector(membership(cluster_louvain(G))),sep='_')))
  indecies=t$Freq>1
  users=sum(t$Freq[indecies])
  groups=sum(indecies)
  return(list(users=users,groups=groups))
}