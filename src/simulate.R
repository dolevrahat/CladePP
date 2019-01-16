#!/usr/bin/Rscript
options(warn=1)  

#This script generates MRS scoe for nSim simulated  gene sets.
# It takes as parameters the number of gene sets to be simulated (n), the size of each gene set (group_size),
#and the path to file containing the hclust objects.
#The file is assumed to contain a list of a hclust objects. Each element of the list contains a anoter list with two  hierarhcial clustering
#of a subset of the NPP matrix corresponding to a different clade, one hclust object in each list element was generated
#using complete linkage, the other using average linkage. Linkage type is indicated by the nest list names.
#The type of linkaeg to be used is determined by the linkage argyment.
# The output is a rds file with a list of tables, where the ith table indicates for each gene with which genes in the ith simulated gene set it is clustered,
#in which clade and its MRS score. 


suppressPackageStartupMessages({
library('dplyr')
library('optparse')
})

options<-list(make_option('--n',type='integer',action='store',help='Number of simulated gene sets to be generated'),
              make_option('--outfile',type='character',action='store',help='Path to file to which the output shall be written'),
              make_option('--genelist',type='character',action='store',help='txt file with list of genes of interest'),
              make_option('--clusters',type='character',action='store',
              help='path to file storing the hclust objects of the NPP matrix')
)
		
where_am_i<-function(){
  location<-sub('^--file=','',grep('^--file=',commandArgs(trailingOnly = FALSE),value=TRUE))
  return(dirname(location))
}

source(sprintf('%s/utils.R',where_am_i()))
       
parser<-OptionParser(option_list=options)
args<-parse_args(parser)	
print(args)
nSim<-args$n
outdir<-args$outdir


clusters<-readRDS(args$clusters)
genepool<-clusters[[1]]$labels

gs<-readLines(args$genelist)
group_size<-length(gs)

perm_list<-lapply(1:nSim,function(i){
  randgs<-sample(genepool,group_size)
  MRSlist<-lapply(names(clusters),function(clade){
  clustobj<-clusters[[clade]]
  dendobj<-as.dendrogram(clustobj)
  tryCatch(
  RS<-ratioScoreDendrogramFilter(dend=dendobj,
	genelist = randgs, genePoolSize = length(labels(dendobj)),
	subtreeFilterFunction=subtreeFilter)	
	)
  if(exists('RS')){
  MRS<-RS %>% group_by(Gene) %>% filter(RS==max(RS))
  MRS$Clade<-clade
  }
  else{
	MRS=NA
  } 
  return(MRS)
})
  MRStable<-do.call('rbind',MRSlist)
  opt_clust<-MRStable %>% group_by(Gene) %>% filter(RS==max(RS)) %>% filter %>% filter(row_number()==1) 
  return(list(MRStable=opt_clust,random_gs=randgs))

})


saveRDS(file=args$outfile,perm_list)


