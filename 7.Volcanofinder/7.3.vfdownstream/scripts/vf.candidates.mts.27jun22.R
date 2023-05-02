# Mon 27 Jun 2022 11:20:50 CEST
# follows from vf.overlap.skov.20jun22.R
# determine if/which mutations fall in the intersect of (vf >0.95 outliers with strict skov outliers) 
   # ie if any protein coding mts within the candidate genes identified 

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
# module load R/4.0.1 
require(data.table)
options(stringsAsFactors=F)
library(tidyverse)
library(GenomicRanges)
library(GenomicFeatures) 

# read in the candidates (intersection of vf 0.95 outliers - strict skov outliers - genes of gtf)
load(file="/scratch/devel/hpawar/admix/volcanofinder/output/merged/em.el.vf.skov.intersect.genes", verbose=T)

# A) assess variant types
find_mts_fun<-function(in_genes,i) {
x<-(data.frame(in_genes[[i]])[1])
chrom<- as.numeric(gsub('chr', '', x))
vep<-read.table(paste("/scratch/devel/shan/gorillas/vep/chr",chrom,".txt",sep=""))
vep$chr<-paste('chr',chrom,sep='')
test<-vep[,1]
test<-do.call(rbind,strsplit(test,split="_"))
vep_ranges<-GRanges(seqnames=vep$chr,ranges=IRanges(start=as.numeric(test[,2]),end=as.numeric(test[,2])+1,names=vep$chr),strand=rep("*",length(vep[,1])),
   variant=(vep[,7]),impact=(vep[,14]))

# 2) overlap the vep data with the candidate genic region (intersect of vf-skov-gtf)
y<- findOverlaps(in_genes[[i]],vep_ranges)
y1<-unique(as.data.frame(y)[,2])

# extract mutations from vep data within range of overlap with the candidate genic region
out_genes<-list()
for (ind in (1:length(y1))) {
out_genes[[ind]]<-data.frame(vep_ranges[y1[[ind]]])
}

mts<-do.call(rbind,out_genes)

return(mts)
}



# B) other aspects of metadata

process_meta_fun<-function(input) {
tes<-as.data.frame(input[,7])
colnames(tes)<-c("impact")
Split <- strsplit(tes$impact, ";", fixed = TRUE)
Ncol <- vapply(Split, length, 1L)
M <- matrix(NA_character_, nrow = nrow(tes),
            ncol = max(Ncol), 
            dimnames = list(NULL, paste0("V", sequence(max(Ncol)))))
M[cbind(rep(1:nrow(tes), Ncol), 
        sequence(Ncol))] <- unlist(Split, use.names = FALSE)
M1<-as.data.frame(M)

# counts frequencies of entries for each column
counts_M1<-apply(M1, 2, table) 
return(counts_M1)
}


mts_per_cand_fun<-function(scen, i){
mts<-find_mts_fun(scen,i)
# calculate frequency of variant types of the mutations found
count_types <- apply(mts[6], 2, table)
mis_mt<-mts[ which(mts$variant=='missense_variant'),]
imp_mt<-process_meta_fun(mts) [[1]]
biot_mt<-process_meta_fun(mts) [[7]]
gene_mt<-process_meta_fun(mts) [[10]]
return(list(mts,count_types,mis_mt,imp_mt,biot_mt,gene_mt))
}


process_mts_cand<-function(scen){
out_mts<-list()
for (ind in (1:length(scen))) {
out_mts[[ind]]<-mts_per_cand_fun(scen, ind)
}
return(out_mts)
}


em.mts<-process_mts_cand(e_cand[[1]]) 
el.mts<-process_mts_cand(e_cand[[2]]) 

e.mts<-list(em.mts,el.mts)

save(e.mts, file="/scratch/devel/hpawar/admix/volcanofinder/output/merged/em.el.vf.skov.intersect.genes.mts")
