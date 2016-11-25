# What if we build an iterative model? One that continues to sub-group louvain groups until we meet 
# some minimal group size threshold.  Synthetic Nested Louvain Block Models.
# keep M as an edge list, just continue to work the diagonal in nested fashion


get_clusters=function(nodeList,g,min.clust=5,alpha,replace,add){
  
  getEdges=function(nodesA,nodesB,g,alpha,replace=TRUE,add=TRUE){
    original.edges=g[(g[,1] %in% nodesA & g[,2] %in% nodesB)|(g[,1] %in% nodesB & g[,2] %in% nodesA),]
    ecount=dim(original.edges)[1]
    n=ceiling(alpha*ecount)
    random.edges=edges=as.data.frame(cbind(sample(nodesA,n,replace=TRUE),sample(nodesB,n,replace=TRUE)))
    names(original.edges)=names(random.edges)=c('Source','Target')
    if(replace){
      edges=rbind(original.edges[-1*sample(1:ecount,n),],random.edges)
    }else{
      if(add){
        edges=rbind(original.edges,random.edges)
      }else{
        edges=original.edges[sample(1:ecount,ecount-n),]
      }
    }
    return(edges)
  }
  g=g[(g$Source %in% nodeList & g$Target %in% nodeList),]
  G=graph.edgelist(cbind(as.character(g$Source),as.character(g$Target)),directed=FALSE)
  nDF=data.frame(node=V(G)$name,cluster=as.vector(membership(cluster_louvain(G))),stringsAsFactors=FALSE)
  cDF=as.data.frame(ftable(nDF$cluster))
  names(cDF)=c('cluster','count')
  cDF$cluster=as.numeric(as.character(cDF$cluster))
  newClust=max(cDF$cluster)+1
  smallClusters=cDF$cluster[cDF$count<min.clust]
  smallNodes=nDF$node[nDF$cluster %in% smallClusters]
  nDF=nDF[!(nDF$node %in% smallNodes),]
  clist=unique(nDF$cluster)
  
  if(length(clist)>1){
    clustList=lapply(unique(nDF$cluster),function(x) nDF$node[nDF$cluster==x])}else{
    clustList=list()
    clustList[[1]]=nDF$node
    }
  if(length(clustList)==0){
    clustList=NULL
    edges=getEdges(nodesA=nodeList,nodesB=nodeList,g=g,alpha,replace,add)
  }else{
   
  clustList[[length(clustList)+1]]=smallNodes
  
  edges=getEdges(nodesA=clustList[[length(clustList)]],nodesB=clustList[[length(clustList)]],g=g,alpha,replace,add)
  cdim=length(clustList)
  for(i in 1:(cdim-1)){
    for(j in (i+1):cdim){
      edges=rbind(edges,
                  getEdges(nodesA=clustList[[i]],
                           nodesB=clustList[[j]],g=g,alpha,replace,add))
    }
  }
  clustList[[cdim]]=NULL
  }
  obj=list(clustList=clustList,edgeList=edges)
  return(obj)
}

simIteration=function(nodeList,g,min.clust,alpha,replace,add){
  clustOBJ=get_clusters(nodeList,g,min.clust,alpha,replace,add)
  queue=clustOBJ$clustList
  edges=clustOBJ$edgeList
  
  # should be able to parrelelize this by queu length
  
  while(length(queue)>0){
    nodes=queue[[1]]
    queue[[1]]=NULL
    if(length(nodes)>0){
      obj=get_clusters(nodeList=nodes,g=m,min.clust,alpha,replace,add)
      edges=rbind(edges,obj$edgeList)
      
      if(length(obj$clustList)>0){
        n=length(queue)
        
        for(i in 1:length(obj$clustList)){
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
  
  
  
  
  
  
  return(list(dens=graph.density(S),
              trans=transitivity(S),
              dia=diameter(S),
              dyad=count_motifs(S,2),
              triad=count_motifs(S,3)))}

