#!/usr/bin/Rscript

# This script is meant to be run after generating a MRS table using the MRS_master.R script
# and a series of files containing MRS tables for simulated gene sets using the MRS_permutation_tests.R script.
# This scrips acceptes in the permdir argument the directory holding the results of the simulated genes sets
# generates using MRS_permutation_tests.R.
# The MSRtable argument contains the path to the MRS table generated under the "real" gold standard gene set.
# The script then adds to MRStable two columns containing for each gene the empirircal p-value of the MRS score
# estimated from the simulated gene sets and the q-value adjusting the p-value for multiple hypothesis testing.
# The new table is saved to a file whose name is specficied in the outfile argument.


suppressPackageStartupMessages({
library('qvalue')
library('dplyr')
library('optparse')
})

options<-list(
  make_option('--permdir',action='store',type='character'),
  make_option('--MRStable',action='store',type='character'),
  make_option('--outfile',action='store',type='character')
)
parser<-OptionParser(option_list=options)
args<-parse_args(parser)

MRSdf<-read.table(args$MRStable,stringsAsFactors=FALSE,header=TRUE,sep='\t')


permfiles<-list.files(pattern='.rds',path=args$permdir,full.names=TRUE)


perm_max_MRS_list<-lapply(permfiles,function(file){
  try({
  perm<-readRDS(file)
  rsv<-unlist(lapply(perm,function(x) x[['MRStable']]$RS))
  rm(perm)
  gc()
  return(rsv)
})})

perm_MRS_vec<-unlist(perm_max_MRS_list)
n_obs<-length(perm_MRS_vec)

t<-as.data.frame(table(perm_MRS_vec))


t$val<-as.numeric(as.character(t$perm_MRS_vec))
t$cum<-cumsum(t$Freq)
t$cum<-t$cum-t$Freq
t$F<-t$cum/n_obs
t$pval<-1-t$F
scoreCol<-ifelse(any(grepl('RS',colnames(MRSdf)))==TRUE,'RS','Score')


MRSdf$p.int<-cut(MRSdf[,scoreCol],t$val,labels=FALSE)+1
MRSdf$p.int[which(MRSdf$RS<=min(t$val))]<-1
MRSdf$p.int[which(MRSdf$RS>max(t$val))]<-nrow(t)+1
p_val_vec<-c(1,t$pval[-1],0)
MRSdf$p.value<-p_val_vec[MRSdf$p.int]
MRSdf$FDR<-p.adjust(MRSdf$p.value,method='BH')
try({
MRSdf$q.value<-qvalue(MRSdf$p.value)[['qvalues']]
})
MRSdf<-MRSdf[,which(colnames(MRSdf)!='p.int')]
print(colnames(MRSdf))
write.table(file=args$outfile,
     MRSdf,row.names = FALSE,quote=FALSE,sep='\t')
save(file='mergedpertest.rdata',MRSdf,t,n_obs)



