#!/usr/bin/Rscript




#Helper function to find number of genes in sub-dendrogram
childsize<-function(dend){
  sapply(dend,function(subdend) length(labels(subdend)))
}


#Function to traverse dendorgram and compute the ratio score for each gene
#params: dend - a dendrogram object
#        genelist - list of the genes in the pathway/set of interest
#        minClusterSize - minimal size of cluster to evalute. The
#        recursion will stop when encountering a cluster of that size or less
#        depth - depth of the recursion - defaults to 0 and updated internally 
#        with each recursive call
#        childindex - string identifiy the cluster currenly being evaluated - updated internally
#        genePoolsize - number of genes in the dendrogram being evaluated
#        subtreeFilterFunction - function used to evalute if the recursion halting condition is satisfied
ratioScoreDendrogramFilter<-function(dend,genelist,depth=0,
childindex='',genePoolSize,subtreeFilterFunction){
  childFilter<-sapply(dend,subtreeFilterFunction,genelist=genelist)	
  if(all(childFilter==FALSE)){
    return(ratioScoreCluster(dend,depth,genelist,childindex,
    genePoolSize=genePoolSize))
  }
  else{
    child_index<-which(childFilter==TRUE)
    child_tables<-lapply(child_index,function(i)
    ratioScoreDendrogramFilter(dend=dend[[i]],genelist = genelist,
    depth=depth+1,childindex = paste(childindex,i,sep = '.'),
    genePoolSize=genePoolSize,subtreeFilterFunction=subtreeFilterFunction))
    all_child_table<-do.call('rbind',child_tables)
    self_table<-ratioScoreCluster(dend,depth,genelist,childindex,
    genePoolSize=genePoolSize)
    return(rbind(self_table,all_child_table))
  }
}

#subtreeFilterFunction - stops the recursion when reaching a cluster
#with less then minClusterSize and minCount gold standard genes
subtreeFilter<-function(subtree,genelist,minClusterSize=3,minCount=2){
  clusterSize<-length(labels(subtree))
  gs_count<-length(intersect(labels(subtree),genelist))
  return(clusterSize>=minClusterSize & gs_count>=minCount)
}

#Compute the ratio score and hypergeometic p-value for each cluster
ratioScoreCluster<-function(dend,depth,genelist=gs,childindex,genePoolSize){
  members<-labels(dend)
  cluster_size<-length(members)
  gs_in_clust<-intersect(members,genelist)
  cluster_score<-length(gs_in_clust)
  cluster_ratio<-cluster_score/cluster_size
  clusterRS<-data.frame(Gene=members,
  Depth=rep(depth,cluster_size),
  Childindex=rep(childindex,cluster_size),
  Cluster_size=rep(cluster_size,cluster_size),
  n_gs_genes_in_clust=cluster_score,
  RS=rep(cluster_ratio,cluster_size),
  gold.standard.genes.in.cluster=paste(gs_in_clust,collapse = ','),
  stringsAsFactors = FALSE)
  return(clusterRS)
}
