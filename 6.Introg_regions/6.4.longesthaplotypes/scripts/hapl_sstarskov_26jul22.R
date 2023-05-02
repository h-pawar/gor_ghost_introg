# Tue 26 Jul 2022 12:23:10 CEST
# generate haplotype networks for intersect regions of length >=100kb
  # first identify how many such regions exist

#-----------------------------------------------------------------------------------------------------------------------
## module load gcc/6.3.0 openssl/1.0.2q  R/4.0.1  BCFTOOLS/1.12
require(data.table)
options(stringsAsFactors=F)
library(tidyverse)
options(scipen=100)
library(adegenet)
library(mgcv)
library(GenomicRanges)
library(tidyr)
library(ggplot2)

#-----------------------------------------------------------------------------------------------------------------------
load(file=paste("/scratch/devel/hpawar/admix/sstar/simul_postabc/stepwisesegs/final_8apr22/check_22apr22/glm",sep=""),verbose=T)

proc_proportion<-function(nput, ci) {

chroms=1:23 #there is no chrom 23.star file -> this is effectively 1:22 (reading in the autosomal data) 

cn1<-list(c("GBG","GBB","GGG","GGD"))

starout<-list()
for (chrom in (chroms)) {
starout[[chrom]]<-read.table(paste("/scratch/devel/hpawar/admix/sstar/out.subsp.comp/",cn1[[1]][[nput]],"_",chrom,".star",sep=""),sep="\t",header=T)[,c(9,4,10,52,12,13,7,1,2)]
}

starout[[9]]<-read.table(paste("/scratch/devel/hpawar/admix/sstar/out.subsp.comp/",cn1[[1]][[nput]],"_newT_9.star",sep=""),sep="\t",header=T)[,c(9,4,10,52,12,13,7,1,2)] 

staroutgbb<-do.call(rbind,starout)

staroutgbb$ind_id[which(staroutgbb$ind_id=="Gbg-Tumani")]<-"Gorilla_beringei_graueri-Tumani"  

staroutgbb[,8]<-paste(staroutgbb[,8],staroutgbb[,9])
indiv.gbb<-unique(staroutgbb[,7])


starperind.gbb<-list()
for (ind in (1:length(indiv.gbb))) {
  starperind.gbb[[ind]]<-list()
  allstars<- staroutgbb[which(staroutgbb[,7]==indiv.gbb[ind] & staroutgbb[,4]>=33333 &staroutgbb[,2]>=30),c(1,8,2)]
  allstars[,3]<-as.numeric(jitter(allstars[,3]))
  newdatA=data.frame(sS=allstars[,3])
  out.pred<-predict(mods[[nput]][[ci]],newdata=newdatA,type="response")
  indval<-cbind(allstars,out.pred) 
  starperind.gbb[[ind]]<-indval[which(indval[,1]>indval[,4]),,drop=F]
  }
iva<-c()
for (ind in 1:length(indiv.gbb)) { iva[ind]<-nrow(starperind.gbb[[ind]]) }

converttoranges_function<-function(nput) {
test<-starperind.gbb[[nput]][,c(1,2)]
test1<-separate(test, chrom, into = c("chr", "winstart"), sep = " (?=[^ ]+$)")
test1$winstart<-as.numeric(test1$winstart)
test1$winend<-test1$winstart + 40000
fG<-GRanges(seqnames=test1[,2],ranges=IRanges(start=as.numeric(test1[,3]),end=as.numeric(test1[,4]),names=test1[,2]),sstar=test1[,1],strand=rep("*",length(test1[,1])))
return(fG)
}


allids<-list()
for (ind in (1:length(indiv.gbb))) {
allids[[ind]]<-converttoranges_function(ind)
}

reduceallids<-list()
for (ind in (1:length(indiv.gbb))) {
reduceallids[[ind]]<-reduce(allids[[ind]])
}

return(reduceallids)
}

#-----------------------------------------------------------------------------------------------------------------------

# read in sstar outlier data for easterns, 99% CI
# GBG, CI 2 (99%)
sstar_GBG_99<-proc_proportion(1,2)

# GBB, CI 2 (99%)
sstar_GBB_99<-proc_proportion(2,2)

#-----------------------------------------------------------------------------------------------------------------------

# skov HMM bed file with the putative introgressed fragments ( ~2% with strict cutoff)
sk<-read.table("/scratch/devel/mkuhlwilm/arch/skov/gor_fragments_strict.bed")

sk_ids_gbb<-unique(sk$V4)[1:12]
sk_ids_gbg<-unique(sk$V4)[13:21]

sk<-subset(sk, V1 != "chrX")

#-----------------------------------------------------------------------------------------------------------------------

skov_proc_proportion<-function(lids) {

sk_per_id<-list()
for (ind in (1:length(lids))) {
sk_per_id[[ind]]<-subset(sk, V4 == lids[ind])
}

skov_converttoranges_function<-function(nput) {
sk<-sk_per_id[[nput]]  
skranges<-GRanges(seqnames=sk[,1],ranges=IRanges(start=as.numeric(sk[,2]),end=as.numeric(sk[,3]),names=sk[,1]),strand=rep("*",length(sk[,1])))
return(skranges)
}

sk_allids<-list()
for (ind in (1:length(lids))) {
sk_allids[[ind]]<-skov_converttoranges_function(ind)
}

reduce_sk_allids<-list()
for (ind in (1:length(lids))) {
reduce_sk_allids[[ind]]<-reduce(sk_allids[[ind]])
}
return(reduce_sk_allids)
}

#-----------------------------------------------------------------------------------------------------------------------

sk_regions_gbb<-skov_proc_proportion(sk_ids_gbb)
sk_regions_gbg<-skov_proc_proportion(sk_ids_gbg)

sk$V5<-sk$V3-sk$V2
sk<-sk[sk$V5 >= 40000, ]
sk_40kbregions_gbb<-skov_proc_proportion(sk_ids_gbb)
sk_40kbregions_gbg<-skov_proc_proportion(sk_ids_gbg)
#-----------------------------------------------------------------------------------------------------------------------

rm_xchr_skov<-function(fG){
aut_sk_regions_gbb<-list()
for (ind in (1:length(fG))) {
test_sk<-fG[[ind]]  
aut_sk_regions_gbb[[ind]]<-test_sk[seqnames(test_sk) != "chrX"] }
return(aut_sk_regions_gbb)
}


autosomes_sk_gbb<-rm_xchr_skov(sk_regions_gbb)
autosomes_sk_gbg<-rm_xchr_skov(sk_regions_gbg)


autosomes_40sk_gbb<-rm_xchr_skov(sk_40kbregions_gbb)
autosomes_40sk_gbg<-rm_xchr_skov(sk_40kbregions_gbg)

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

# overlap skov outliers with S* outliers


calc_prop_fun<-function(fG,fB){
fG<-reduce(fG);fB<-reduce(fB)
pos.overlaps<-reduce(subsetByOverlaps(fG,fB))
prop.bp.x<-sum(reduce(subsetByOverlaps(fG,fB))@ranges@width)/sum(fG@ranges@width) # proportion of overlapping bp
prop.frag.x<-length(reduce(subsetByOverlaps(fG,fB)))/length(fG) # proportion of overlapping fragments
prop.bp.y<-sum(reduce(subsetByOverlaps(fG,fB))@ranges@width)/sum(fB@ranges@width) # proportion of overlapping bp
prop.frag.y<-length(reduce(subsetByOverlaps(fG,fB)))/length(fB) # proportion of overlapping fragments
out.x<-(cbind(prop.bp.x, prop.frag.x))
out.y<-(cbind(prop.bp.y, prop.frag.y))
return(rbind(out.x,out.y))
}
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
# process output of overlaps : regions overlapping, length dist, proportion (overlapping bp & fragments) - order of query matters
process_tooverlap_fun<-function(sstar_subs,skov_subs){
ranges_2<-reduce(subsetByOverlaps(skov_subs,sstar_subs))
lengths_2<-summary(reduce(subsetByOverlaps(skov_subs,sstar_subs))@ranges@width)/1000
prop_2<-calc_prop_fun(skov_subs,sstar_subs)
return(list(ranges_2,lengths_2,prop_2))
}

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

intersect_function<-function(sstar_subs,skov_subs) {
out_gbb_95<-list()
for (ind in (1:length(sstar_subs))) {
out_gbb_95[[ind]]<-process_tooverlap_fun(sstar_subs[[ind]],skov_subs[[ind]])[[1]]
}
return(out_gbb_95)
}

ov_gbb_99<-intersect_function(sstar_GBB_99,autosomes_sk_gbb)
ov_gbg_99<-intersect_function(sstar_GBG_99,autosomes_sk_gbg)
ov_gbb40_99<-intersect_function(sstar_GBB_99,autosomes_40sk_gbb)
ov_gbg40_99<-intersect_function(sstar_GBG_99,autosomes_40sk_gbg)

#-----------------------------------------------------------------------------------------------------------------------

save(ov_gbb_99,file="/scratch/devel/hpawar/admix/overlap.s*.skov/overlap_objects/ov_gbb_99")
save(ov_gbg_99,file="/scratch/devel/hpawar/admix/overlap.s*.skov/overlap_objects/ov_gbg_99")
save(ov_gbb40_99,file="/scratch/devel/hpawar/admix/overlap.s*.skov/overlap_objects/ov_gbb40_99")
save(ov_gbg40_99,file="/scratch/devel/hpawar/admix/overlap.s*.skov/overlap_objects/ov_gbg40_99")

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------


extract_longwin_fun<-function(scen,x){
longestwindowsperid<-list()
for (ind in (1:length(scen))) {
longestwindowsperid[[ind]]<-scen[[ind]][which(width(scen[[ind]])>x),]
}
return(longestwindowsperid)
}


l_gbb_100k<-extract_longwin_fun(ov_gbb_99,100000)
l_gbb_250k<-extract_longwin_fun(ov_gbb_99,250000)
l_gbg_100k<-extract_longwin_fun(ov_gbg_99,100000)
l_gbg_250k<-extract_longwin_fun(ov_gbg_99,250000)


longwindows<-list(l_gbb_100k,l_gbb_250k,l_gbg_100k,l_gbg_250k)

save(longwindows,file=paste("/scratch/devel/hpawar/admix/overlap.s*.skov/hapl/longwindows",sep="")) 

