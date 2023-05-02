# Thu 21 Jul 2022 16:51:47 CEST
# determine if/which mutations fall in the intersect of (vf >0.95 outliers with strict skov outliers) 
   # ie if any protein coding mts within the candidate genes identified 

#-----------------------------------------------------------------------------------------------------------------------

# module load R/4.0.1 
require(data.table)
options(stringsAsFactors=F)
library(tidyverse)
library(GenomicRanges)
library(GenomicFeatures) 
# read in the candidates (intersection of vf 0.95 outliers - strict skov outliers - genes of gtf)
load(file="/scratch/devel/hpawar/admix/volcanofinder/output/merged/em.el.vf.s*.skov.intersect.genes.21jul22", verbose=T)

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


#-----------------------------------------------------------------------------------------------------------------------

em_mts<-list()
for (ind in (1:length(e_cand_fin))) {
em_mts[[ind]]<-find_mts_fun(e_cand_fin, ind)
}


save(em_mts, file="/scratch/devel/hpawar/admix/volcanofinder/output/merged/em.el.vf.s*.skov.intersect.genes.mts.21jul22")


#-----------------------------------------------------------------------------------------------------------------------

# calculate frequency of variant types of the mutations found

assess_types_fun<-function(input_df) {
mts<-input_df
count_types <- apply(mts[6], 2, table)
mis_mt<-mts[ which(mts$variant=='missense_variant'),]
return(list(count_types,mis_mt))
}

typ<-list()
for (ind in (1:length(e_cand_fin))) {
typ[[ind]]<-assess_types_fun(em_mts[[ind]])
}

save(typ, file="/scratch/devel/hpawar/admix/volcanofinder/output/merged/em.el.vf.s*.skov.intersect.genes.mts.meta.21jul22")


q()


#-----------------------------------------------------------------------------------------------------------------------
