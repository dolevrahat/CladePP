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

options(stringsAsFactors = FALSE,warn=1)
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


where_am_i<-function(){
	location<-sub('^--file=','',grep('^--file=',commandArgs(trailingOnly = FALSE),value=TRUE))
	return(dirname(location))
}

source(sprintf('%s/utils.R',where_am_i()))

parser<-OptionParser(option_list=options)
args<-parse_args(parser)

path_to_clusters<-args$clusters
path_to_genelist<-args$genelist
outfile<-args$outfile



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
  print(str(MRS))
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
final<-opt_clust_dedup[,c('Gene','Score','in.gold.standard','gold.standard.genes.in.cluster','Clade')]

write.table(file=outfile,final,row.names = FALSE,quote = FALSE,sep='\t')
