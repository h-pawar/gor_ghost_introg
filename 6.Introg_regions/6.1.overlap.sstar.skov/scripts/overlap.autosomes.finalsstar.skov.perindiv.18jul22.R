# Mon 18 Jul 2022 09:54:24 CEST
# overlap S* windows of 99% CI with strict skov hmm autosomal outliers - calculate per eastern individual

#-----------------------------------------------------------------------------------------------------------------------
## module load gcc/6.3.0 openssl/1.0.2q R/4.0.1
require(data.table)
options(stringsAsFactors=F)
options(scipen=100)
library(mgcv)
library(GenomicRanges)
library(tidyr)
library(ggbio)
library(pheatmap)
library(ggplot2)
library(ggpubr)

#-----------------------------------------------------------------------------------------------------------------------
# glm generated for the S* CIs calc from simulated data (null simulations - gor demog, no ghost)
load(file=paste("/scratch/devel/hpawar/admix/sstar/simul_postabc/stepwisesegs/final_8apr22/check_22apr22/glm",sep=""),verbose=T)

#-----------------------------------------------------------------------------------------------------------------------
# function to read in sstar outlier data per scenario & per CI & split this into s* windows per target individual

proc_proportion<-function(nput, ci) {

# read in the empirical data
chroms=1:23 #there is no chrom 23.star file -> this is effectively 1:22 (reading in the autosomal data) 

  # update to per subspecies comparisons
# 1) WL outgroup, EL ingroup (GGG outg, GBG ing) = "GBG"
# 2) WL outgroup, EM ingroup (GGG outg, GBB ing) = "GBB"
# 3) EL outgroup, WL ingroup (GBG outg, GGG ing) = "GGG"
# 4) EL outgroup, WC ingroup (GBG outg, GGD ing) = "GGD"


cn1<-list(c("GBG","GBB","GGG","GGD"))

starout<-list()
for (chrom in (chroms)) {
starout[[chrom]]<-read.table(paste("/scratch/devel/hpawar/admix/sstar/out.subsp.comp/",cn1[[1]][[nput]],"_",chrom,".star",sep=""),sep="\t",header=T)[,c(9,4,10,52,12,13,7,1,2)]
}
# read in data for s* applied to empirical data
starout[[9]]<-read.table(paste("/scratch/devel/hpawar/admix/sstar/out.subsp.comp/",cn1[[1]][[nput]],"_newT_9.star",sep=""),sep="\t",header=T)[,c(9,4,10,52,12,13,7,1,2)] 

staroutgbb<-do.call(rbind,starout)

staroutgbb$ind_id[which(staroutgbb$ind_id=="Gbg-Tumani")]<-"Gorilla_beringei_graueri-Tumani"  

#get values per individual - for each comparison
staroutgbb[,8]<-paste(staroutgbb[,8],staroutgbb[,9])
indiv.gbb<-unique(staroutgbb[,7]) # col 7 = ind_id


starperind.gbb<-list()
for (ind in (1:length(indiv.gbb))) {
  starperind.gbb[[ind]]<-list()
  # only use windows where there is enough data (33,333 bp out of 40,000, with at least 30 segregating sites)
  allstars<- staroutgbb[which(staroutgbb[,7]==indiv.gbb[ind] & staroutgbb[,4]>=33333 &staroutgbb[,2]>=30),c(1,8,2)]
  allstars[,3]<-as.numeric(jitter(allstars[,3]))
  # create a data frame of the same segregating sites:
  newdatA=data.frame(sS=allstars[,3])
  # predict S* vals (given the segregating sites) at given ci
  out.pred<-predict(mods[[nput]][[ci]],newdata=newdatA,type="response")
  indval<-cbind(allstars,out.pred)
  # which windows lie outside the expectation for the ci (3) 
  starperind.gbb[[ind]]<-indval[which(indval[,1]>indval[,4]),,drop=F]
  }

## you can see that each individual has a different number of significant regions:
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
# 1) read in sstar outlier data, per scenario, per CI (splits per individual & converts to reduced granges objects)

# read in sstar outlier data for easterns, 99% CI
# GBG, CI 1 (99%)
sstar_GBG_99<-proc_proportion(1,2)

# GBB, CI 1 (99%)
sstar_GBB_99<-proc_proportion(2,2)

#-----------------------------------------------------------------------------------------------------------------------

# 2) read in skov strict outliers & using function below split per inidivdual & generate reduced granges objects
# skov HMM bed file with the putative introgressed fragments ( ~2% with strict cutoff)
sk<-read.table("/scratch/devel/mkuhlwilm/arch/skov/gor_fragments_strict.bed")

# individuals for GBB & GBG
sk_ids_gbb<-unique(sk$V4)[1:12]
sk_ids_gbg<-unique(sk$V4)[13:21]

#-----------------------------------------------------------------------------------------------------------------------

# convert skov outliers to granges objects per individual
skov_proc_proportion<-function(lids) {

# split skov data -> skov fragments per individual
sk_per_id<-list()
for (ind in (1:length(lids))) {
sk_per_id[[ind]]<-subset(sk, V4 == lids[ind])
}

# convert skov per id also to granges 
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

# 2.1) filter strict skov regions by length - 40kb cutoff
sk$V5<-sk$V3-sk$V2
sk<-sk[sk$V5 >= 40000, ]
sk_40kbregions_gbb<-skov_proc_proportion(sk_ids_gbb)
sk_40kbregions_gbg<-skov_proc_proportion(sk_ids_gbg)

#-----------------------------------------------------------------------------------------------------------------------
# filter skov to only autosomes

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
#-----------------------------------------------------------------------------------------------------------------------

# calculate proportion of bp overlapping in intersect

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
# pairwise comparison of outlier windows : diff number of ranges obtained depending on which is query & which labelled as source
# ie query s* outliers by skov bed files & vice versa
ranges_2<-reduce(subsetByOverlaps(skov_subs,sstar_subs))
# length distribution of overlapping regions
lengths_2<-summary(reduce(subsetByOverlaps(skov_subs,sstar_subs))@ranges@width)/1000
# proportion overlapping
prop_2<-calc_prop_fun(skov_subs,sstar_subs)
return(list(ranges_2,lengths_2,prop_2))
}

#-----------------------------------------------------------------------------------------------------------------------

# perform intersect b/n s* & skov outliers for all ids
intersect_function<-function(sstar_subs,skov_subs) {
# 1) compare GBB - 95% CI for S*
#sstar_subs #sk_regions_gb
out_gbb_95<-list()
for (ind in (1:length(sstar_subs))) {
out_gbb_95[[ind]]<-process_tooverlap_fun(sstar_subs[[ind]],skov_subs[[ind]])
}


prop_overlap_gbb_95<-list()
for (ind in (1:length(sstar_subs))) {
prop_overlap_gbb_95[[ind]]<-out_gbb_95[[ind]][[3]][1,]
}

prop_df<-do.call(rbind, prop_overlap_gbb_95)

# length dists in kb
length_overlap_gbb_95<-list()
for (ind in (1:length(sstar_subs))) {
length_overlap_gbb_95[[ind]]<-out_gbb_95[[ind]][[2]]
}

length_df<- do.call(rbind,length_overlap_gbb_95)

return(list(out_gbb_95,prop_df,length_df))
}
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------


intersect_function(sstar_GBB_99,autosomes_sk_gbb)[[2]]
intersect_function(sstar_GBG_99,autosomes_sk_gbg)[[2]]

intersect_function(sstar_GBB_99,autosomes_40sk_gbb)[[2]]
intersect_function(sstar_GBG_99,autosomes_40sk_gbg)[[2]]

#-----------------------------------------------------------------------------------------------------------------------


as.data.frame(intersect_function(sstar_GBB_99,autosomes_sk_gbb)[[2]][,1])
as.data.frame(intersect_function(sstar_GBG_99,autosomes_sk_gbg)[[2]][,1])
as.data.frame(intersect_function(sstar_GBB_99,autosomes_40sk_gbb)[[2]][,1])
as.data.frame(intersect_function(sstar_GBG_99,autosomes_40sk_gbg)[[2]][,1])


#-----------------------------------------------------------------------------------------------------------------------

# GBG, CI 1 (95%)
sstar_GBG_95<-proc_proportion(1,1)
# GBB, CI 1 (95%)
sstar_GBB_95<-proc_proportion(2,1)

# GBG, CI 2 (99.5%)
sstar_GBG_99.5<-proc_proportion(1,3)
# GBB, CI 2 (99%)
sstar_GBB_99.5<-proc_proportion(2,3)


scen_GBB<-list(sstar_GBB_95,sstar_GBB_99,sstar_GBB_99.5)
scen_GBG<-list(sstar_GBG_95,sstar_GBG_99,sstar_GBG_99.5)


re_calc_overlap<-function(scen,skov_version){
overlap40kb_gbb<-list()
for (ind in (1:length(scen))) {
overlap40kb_gbb[[ind]]<-intersect_function(scen[[ind]],skov_version)[[2]]}
return(overlap40kb_gbb)
}


oprop_gbb<-re_calc_overlap(scen_GBB,autosomes_sk_gbb)
oprop_gbg<-re_calc_overlap(scen_GBG,autosomes_sk_gbg)

oprop_40gbb<-re_calc_overlap(scen_GBB,autosomes_40sk_gbb)
oprop_40gbg<-re_calc_overlap(scen_GBG,autosomes_40sk_gbg)


#-----------------------------------------------------------------------------------------------------------------------

