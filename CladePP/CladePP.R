#!/usr/bin/Rscript

#This script implements the main methods of the cladePP algorithm.
#It takes as the following inputs:
#clusters - a rds file with a named list of hclust object 
#representing hierarchial clustering of the NPP matrix, one object per clade.
#The names of the list elements should correspond to the clades.
# genelist - A txt file with the gene symbols of the genes of interest.
# outfile - Path to file to which the output shall be written.
# The output is a table indicating for each gene with which gold standard genes it is clustered,
#in which clade and its MRS score.

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
library(optparse,quietly = TRUE)
library(plyr,quietly = TRUE)
library(dplyr,quietly = TRUE)
})
options<-list(
  make_option('--clusters',action='store',type='character'),
  make_option('--genelist',action='store',type='character'),
  make_option('--outfile',action='store',type='character')
  )

parser<-OptionParser(option_list=options)
args<-parse_args(parser)

path_to_clusters<-args$clusters
path_to_genelist<-args$genelist
outfile<-args$outfile




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
  gold.standard.genes.in.cluster=paste(gs_in_clust,collapse = ','))
  return(clusterRS)
}

cluster_list<-readRDS(path_to_clusters)
gs<-readLines(path_to_genelist)



MRSlist<-lapply(names(cluster_list),function(clade){
  clustobj<-cluster_list[[clade]]
  dendobj<-as.dendrogram(clustobj)
  RS<-ratioScoreDendrogramFilter(dend=dendobj,
  genelist = gs, genePoolSize = length(labels(dendobj)),
  subtreeFilterFunction=subtreeFilter)
  MRS<-RS %>% group_by(Gene) %>% filter(RS==max(RS))
  MRS$Clade<-clade
  return(MRS)
})

MRStable<-do.call('rbind',MRSlist)

MRStable.allclades<-MRStable[!duplicated(paste(MRStable$Gene,MRStable$Clade,sep='.')),]


opt_clust<-MRStable %>% group_by(Gene) %>% filter(RS==max(RS))
opt_clust_dedup<- opt_clust %>% group_by(Gene) %>% filter(Depth==min(Depth))
opt_clust_dedup<-opt_clust_dedup[!duplicated(opt_clust_dedup$Gene),]

opt_clust_dedup$Rank<-rank(opt_clust_dedup$RS,ties.method = 'random')
opt_clust_dedup$in.gold.standard<-sapply(opt_clust_dedup$Gene, function (g) ifelse(g%in%gs,'Yes','No'))
colnames(opt_clust_dedup)<-sub('RS','Score',colnames(opt_clust_dedup))
final<-opt_clust_dedup[,c('Gene','Rank','Score','in.gold.standard','gold.standard.genes.in.cluster','Clade')]

write.table(file=outfile,final,row.names = FALSE,quote = FALSE,sep='\t')