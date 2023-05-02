# Thu 21 Jul 2022 12:22:29 CEST
#-----------------------------------------------------------------------------------------------------------------------

# 6. Output file format - http://degiorgiogroup.fau.edu/Manual_VolcanoFinder_v1.0.pdf

#Each row represents the calculation of a log likelihood ratio test statistic for adaptive introgression at a given grid position. 
# col 1 = location = test site position
# col 2 = LR = log likelihood ratio test stat for adaptive intro
# col 3 = alpha
# col 4 = D = genetic dist b.n smpled & donor - set by the user

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
# module load gcc/6.3.0 openssl/1.0.2q  R/4.0.1 
require(data.table)
options(stringsAsFactors=F)
library(tidyverse)
options(scipen=100)
library(mgcv)
library(GenomicRanges)
library(GenomicFeatures)
library(tidyr)
library(ggbio)
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

# 1.1) read in volcanofinder output & extract outliers
vfout<-list()
for (i in 1:22) {
test<-read.table(paste("/scratch/devel/hpawar/admix/volcanofinder/output/merged/4_",i,"_test",sep=""),sep="\t",header=T)
test$chr<-paste('chr',i,sep='')
vfout[[i]]<-test
}
vfoutall<-do.call(rbind,vfout)
# extract top subsets of volcanofinder scores
#  > 95% outliers
vf_95<-vfoutall[which(vfoutall[,2]>=quantile(vfoutall$LR, probs = 0.95)),]


#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

vf_95_ranges<-GRanges(seqnames=vf_95[,5],ranges=IRanges(start=as.numeric(vf_95[,1]),end=as.numeric(vf_95[,1]+1),names=vf_95[,5]),strand=rep("*",length(vf_95[,1])),LR=as.numeric(vf_95[,2]))

#-----------------------------------------------------------------------------------------------------------------------

# overlap s* 99% outliers with strict skov outliers

load(file=paste("/scratch/devel/hpawar/admix/sstar/simul_postabc/stepwisesegs/final_8apr22/check_22apr22/glm",sep=""),verbose=T)

#-----------------------------------------------------------------------------------------------------------------------

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
# GBG, CI 1 (99%)
sstar_GBG_99<-proc_proportion(1,2)

# GBB, CI 1 (99%)
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


calc_prop_fun<-function(fG,fB){
fG<-reduce(fG);fB<-reduce(fB) 
pos.overlaps<-reduce(subsetByOverlaps(fG,fB)) 
prop.bp.x<-sum(reduce(subsetByOverlaps(fG,fB))@ranges@width)/sum(fG@ranges@width) 
prop.frag.x<-length(reduce(subsetByOverlaps(fG,fB)))/length(fG) 
prop.bp.y<-sum(reduce(subsetByOverlaps(fG,fB))@ranges@width)/sum(fB@ranges@width)
prop.frag.y<-length(reduce(subsetByOverlaps(fG,fB)))/length(fB)
out.x<-(cbind(prop.bp.x, prop.frag.x))
out.y<-(cbind(prop.bp.y, prop.frag.y))
return(rbind(out.x,out.y))
}
#-----------------------------------------------------------------------------------------------------------------------

process_tooverlap_fun<-function(sstar_subs,skov_subs){
ranges_2<-reduce(subsetByOverlaps(skov_subs,sstar_subs))
lengths_2<-summary(reduce(subsetByOverlaps(skov_subs,sstar_subs))@ranges@width)/1000
prop_2<-calc_prop_fun(skov_subs,sstar_subs)
return(list(ranges_2,lengths_2,prop_2))
}


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

gtf="/scratch/project/devel/mkuhlwilm/chimpdiv/index/Homo_sapiens.GRCh37.75_c.gtf"

test<-makeTxDbFromGFF(gtf, format='gtf')

human.genes = genes(test)

autosome.hum.genes<-human.genes[seqnames(human.genes) == c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22")]

#-----------------------------------------------------------------------------------------------------------------------


overlap_vfoutl_skov2<-function(vf_obj,sk_obj) {
out_overlaps_95<-list()
out_1<-list()
for (i in (1:length(sk_obj))) {
out_overlaps_95[[i]]<-subsetByOverlaps(vf_obj,sk_obj[[i]]) 
out_1[[i]]<-intersect(autosome.hum.genes, out_overlaps_95[[i]], ignore.strand=TRUE)
}
return(out_1)
}

proc_genes2<-function(overlap_obj) {
out_1<-overlap_obj  
x<- findOverlaps(out_1,autosome.hum.genes) 
x1<-unique(as.data.frame(x)[,2])
out_genes<-list()
for (ind in (1:length(x1))) {
out_genes[[ind]]<-autosome.hum.genes[x1[[ind]]]
}

y<-unlist(out_genes)

return(y)
}

output_genes<-function(skov_obj){

ov2<-overlap_vfoutl_skov2(vf_95_ranges,skov_obj)

check_overlaps<-list()
for (a in (1:length(ov2))) {
check_overlaps[[a]]<-length(ov2[[a]])
}

hold_genes2<-list()
for (a in (1:length(ov2))) {
hold_genes2[[a]]<-proc_genes2(ov2[[a]])
}

y<-unlist(hold_genes2)

return(y)
}


o_gbg<-output_genes(ov_gbg_99)



output_genes_gbb<-function(skov_obj){

ov2<-overlap_vfoutl_skov2(vf_95_ranges,skov_obj)

check_overlaps<-list()
for (a in (1:length(ov2))) {
check_overlaps[[a]]<-length(ov2[[a]])
}

ov2<-ov2[-which(check_overlaps==0)]

hold_genes2<-list()
for (a in (1:length(ov2))) {
hold_genes2[[a]]<-proc_genes2(ov2[[a]])
}

y<-unlist(hold_genes2)

return(y)
}

o_gbb<-output_genes_gbb(ov_gbb_99)

#-----------------------------------------------------------------------------------------------------------------------

unique(o_gbb)
unique(o_gbg)

#-----------------------------------------------------------------------------------------------------------------------

l1<-list(unique(o_gbb), unique(o_gbg))
e_cand_fin<-unique(unlist(l1))
save(e_cand_fin, file="/scratch/devel/hpawar/admix/volcanofinder/output/merged/em.el.vf.s*.skov.intersect.genes.21jul22")


