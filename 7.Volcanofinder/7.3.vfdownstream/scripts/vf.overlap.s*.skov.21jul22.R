# Thu 21 Jul 2022 12:22:29 CEST
# RUN INTERACTIVELY **

# amended from vf.overlap.skov.20jun22.R

# Mon 20 Jun 2022 08:58:34 CEST
# 0) extract top subsets of volcanofinder scores
# 1) overlap with skov fragments

#-----------------------------------------------------------------------------------------------------------------------

# asked MK re QC 

# Thu, 16 Jun, 17:58 
# MK
#I think the very first steps are to look at the top 10 peaks, since they
#are very nicely looking peaks, and ask:
#- do they overlap with S*/Skov fragments?
#- do they contain interesting genes?
#- do these genes contain interesting mutations?

#-----------------------------------------------------------------------------------------------------------------------

# follows from 
# /scratch/devel/hpawar/admix/volcanofinder/scripts/vf.test.autosomes.blocks.arr # to run volcanofinder per fragment 
#/Volumes/"Ultra USB 3.0"/IBE/further.analysis.feb.2020/gorillas/abc/simul_postabc/abc+ghost/modelcomp/volcanofinder/plot.vf.25may22.R # merge the completed reps, & generate diagnostic plots of histograms etc
#-----------------------------------------------------------------------------------------------------------------------
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
#BiocManager::install(c("GenomicRanges", "plyranges", "HelloRangesData"))
#library(plyranges) # may not be necessary here
library(GenomicFeatures) # if need to read in the gtf
library(tidyr)
library(ggbio)
#library(valr)
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

# 1.1) read in volcanofinder output & extract outliers
vfout<-list()
for (i in 1:22) {
test<-read.table(paste("/scratch/devel/hpawar/admix/volcanofinder/output/merged/4_",i,"_test",sep=""),sep="\t",header=T)
# chromosomes need to be specified in the same way in both data objects (the skov & the volcanofinder)
test$chr<-paste('chr',i,sep='')
vfout[[i]]<-test
}
vfoutall<-do.call(rbind,vfout)
# extract top subsets of volcanofinder scores
#  > 95% outliers
vf_95<-vfoutall[which(vfoutall[,2]>=quantile(vfoutall$LR, probs = 0.95)),]

#> head(vf_95)
#     location       LR        alpha        D  chr
#1573  1634282 24.36087 4.814070e-05 0.019332 chr1
#1574  1635282 25.34131 5.027065e-05 0.019332 chr1
# marker every 10kb

# Q - which cutoff to take?
#> nrow(vfoutall[which(vfoutall[,2]>=quantile(vfoutall$LR, probs = 0.95)),])
#[1] 131633
#> nrow(vfoutall[which(vfoutall[,2]>=quantile(vfoutall$LR, probs = 0.99)),])
#[1] 26327
#> nrow(vfoutall[which(vfoutall[,2]>=quantile(vfoutall$LR, probs = 0.995)),])
#[1] 13164

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

# 1.2) convert vf outliers to genomic ranges format 

vf_95_ranges<-GRanges(seqnames=vf_95[,5],ranges=IRanges(start=as.numeric(vf_95[,1]),end=as.numeric(vf_95[,1]+1),names=vf_95[,5]),strand=rep("*",length(vf_95[,1])),LR=as.numeric(vf_95[,2]))

#> vf_95_ranges
#GRanges object with 131633 ranges and 1 metadata column:
#        seqnames            ranges strand |        LR
#           <Rle>         <IRanges>  <Rle> | <numeric>
#   chr1     chr1   1634282-1634283      * |   24.3609
#   chr1     chr1   1635282-1635283      * |   25.3413
#   chr1     chr1   1636282-1636283      * |   25.9850
#   chr1     chr1   1637282-1637283      * |   26.4573
#   chr1     chr1   1638282-1638283      * |   26.8564
#    ...      ...               ...    ... .       ...
#  chr22    chr22 49709960-49709961      * |   46.6849
#  chr22    chr22 49710960-49710961      * |   40.8915
#  chr22    chr22 49713960-49713961      * |   32.6973
#  chr22    chr22 49714960-49714961      * |   32.5112
#  chr22    chr22 49715960-49715961      * |   28.0831
#  -------
#  seqinfo: 22 sequences from an unspecified genome; no seqlengths

  # coudl also add the alpha & D metadata columns? - or this info is lost anyway when intersecting..

#tail(vf_95[order(vf_95$LR),], n=10)
# ie this is the same peak - want to extract 10 distinct peaks



#> tail(vf_95[order(vf_95$LR),])
#        location       LR        alpha        D chr
#1728079 55128436 1144.024 3.797502e-05 0.077329  11
#1728084 55133436 1146.978 3.804284e-05 0.077329  11

# orders of magnitude higher than max log likelihood ratio of table d1? of setter et al?
#https://storage.googleapis.com/plos-corpus-prod/10.1371/journal.pgen.1008867/2/pgen.1008867.s004.pdf?X-Goog-Algorithm=GOOG4-RSA-SHA256&X-Goog-Credential=wombat-sa%40plos-prod.iam.gserviceaccount.com%2F20220620%2Fauto%2Fstorage%2Fgoog4_request&X-Goog-Date=20220620T072904Z&X-Goog-Expires=86400&X-Goog-SignedHeaders=host&X-Goog-Signature=7f0fd55d3e62fd00c5e3160502a140108bb4c52447f82cf6c4bc32b1c18e3a1f01ab1546e0a5d2abfb290b21c4d97dff2a916fe4e7bf40b69d5daea8eda5d8a9a03a6d2445469f4e16fff7e7833573fc3a3b1669f63c05d1fc6fd813dfeab52583dfd8196da1f0264328d999248ca00b25398b4a5a3da69df164d646acf81e5abb524c3582ce86e3ab5b6ce1327fc6e8b1f5243e236c4b87dbaa79eb21fb667cce6b34fa090312797bbd2db761e6eb03d97899d13822f743e48c245513c79a454aa73bf0783345d3cce2e799b80f3fa570e10d2ef52c1aecf772725bf1524678eae1581c1080c2d306dc13050d6de7ef39156e4c0a65789cdef8d68df5a7fb83


#-----------------------------------------------------------------------------------------------------------------------

# overlap s* 99% outliers with strict skov outliers


# glm generated for the S* CIs calc from simulated data (null simulations - gor demog, no ghost)
load(file=paste("/scratch/devel/hpawar/admix/sstar/simul_postabc/stepwisesegs/final_8apr22/check_22apr22/glm",sep=""),verbose=T)

#-----------------------------------------------------------------------------------------------------------------------
# function to read in sstar outlier data per scenario & per CI & split this into s* windows per target individual

proc_proportion<-function(nput, ci) {
  # for nput=1 # GBG

# read in the empirical data
chroms=1:23 #there is no chrom 23.star file -> this is effectively 1:22 (reading in the autosomal data) 
# will need to specify either GB or GG
#cn1<-list(c("GB","GG"))

  # update to per subspecies comparisons
# 1) WL outgroup, EL ingroup (GGG outg, GBG ing) = "GBG"
# 2) WL outgroup, EM ingroup (GGG outg, GBB ing) = "GBB"
# 3) EL outgroup, WL ingroup (GBG outg, GGG ing) = "GGG"
# 4) EL outgroup, WC ingroup (GBG outg, GGD ing) = "GGD"


cn1<-list(c("GBG","GBB","GGG","GGD"))


# 1)  GBG, CI 3, 99.5%

# for nput=1 # GBG

#nput=1

starout<-list()
for (chrom in (chroms)) {
starout[[chrom]]<-read.table(paste("/scratch/devel/hpawar/admix/sstar/out.subsp.comp/",cn1[[1]][[nput]],"_",chrom,".star",sep=""),sep="\t",header=T)[,c(9,4,10,52,12,13,7,1,2)]
}
# read in data for s* applied to chr 9 for target individuals plus newly processed tumani (whose chr 9 had sparse data issues)
starout[[9]]<-read.table(paste("/scratch/devel/hpawar/admix/sstar/out.subsp.comp/",cn1[[1]][[nput]],"_newT_9.star",sep=""),sep="\t",header=T)[,c(9,4,10,52,12,13,7,1,2)] 

staroutgbb<-do.call(rbind,starout)
# ensure all tumani fragments are named consistently
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

#iva
#[1] 900 842 898 876 883 800 864 818 876
#> length(iva)
#[1] 9

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

# (1) read in skov strict outliers & using function below split per inidivdual & generate reduced granges objects

# skov HMM bed file with the putative introgressed fragments ( ~2% with strict cutoff)
sk<-read.table("/scratch/devel/mkuhlwilm/arch/skov/gor_fragments_strict.bed")
#head(sk)
#    V1       V2       V3                                V4
#1 chr1  5199000  5219000 Gorilla_beringei_beringei-Bwiruka
#2 chr1  5402000  5460000 Gorilla_beringei_beringei-Bwiruka

#> nrow(sk)
#[1] 15691

# convert startpos & endpos from int to numeric values
#sk[,c(2:3)]<-as.data.frame(lapply(sk[,c(2:3)], as.numeric)) # not necessary - b/c doing this step in the convert to granges function below

#head(sk)
#    V1       V2       V3                                V4
#1 chr1  5199000  5219000 Gorilla_beringei_beringei-Bwiruka
#2 chr1  5402000  5460000 Gorilla_beringei_beringei-Bwiruka

# individuals for GBB & GBG
sk_ids_gbb<-unique(sk$V4)[1:12]
sk_ids_gbg<-unique(sk$V4)[13:21]

sk<-subset(sk, V1 != "chrX")

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
# perhaps already reduced? in the sk_allids step?
return(reduce_sk_allids)
}

#-----------------------------------------------------------------------------------------------------------------------

sk_regions_gbb<-skov_proc_proportion(sk_ids_gbb) # works
sk_regions_gbg<-skov_proc_proportion(sk_ids_gbg)

# 2.1) filter strict skov regions by length - 40kb cutoff
sk$V5<-sk$V3-sk$V2
sk<-sk[sk$V5 >= 40000, ]
sk_40kbregions_gbb<-skov_proc_proportion(sk_ids_gbb) # works
sk_40kbregions_gbg<-skov_proc_proportion(sk_ids_gbg)
#-----------------------------------------------------------------------------------------------------------------------
# there are X chr in the skov outliers but not in the s* -> not a fair comparison

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

# Wed 20 Jul 2022 17:35:51 CEST
calc_prop_fun<-function(fG,fB){
fG<-reduce(fG);fB<-reduce(fB) # amendment - Wed 23 Feb 2022 17:28:27 CET
pos.overlaps<-reduce(subsetByOverlaps(fG,fB)) ## this may be more useful output? (gives ranges which overlap ie the exact positions) (yes - these are the unique overlapping regions)
## shoudl calc proportion overlapping both ways (ie of A how much overlaps with B + vice versa)
## & calc propotion of only unique overlapping regions
prop.bp.x<-sum(reduce(subsetByOverlaps(fG,fB))@ranges@width)/sum(fG@ranges@width) # proportion of overlapping bp
prop.frag.x<-length(reduce(subsetByOverlaps(fG,fB)))/length(fG) # proportion of overlapping fragments
prop.bp.y<-sum(reduce(subsetByOverlaps(fG,fB))@ranges@width)/sum(fB@ranges@width) # proportion of overlapping bp
prop.frag.y<-length(reduce(subsetByOverlaps(fG,fB)))/length(fB) # proportion of overlapping fragments
## shoudl combine the following into 1 output? b/c also need to output for vice versa B in A
out.x<-(cbind(prop.bp.x, prop.frag.x))
out.y<-(cbind(prop.bp.y, prop.frag.y))
#return(pos.overlaps)
return(rbind(out.x,out.y))
}
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
# process output of overlaps : regions overlapping, length dist, proportion (overlapping bp & fragments) - order of query matters
process_tooverlap_fun<-function(sstar_subs,skov_subs){
# pairwise comparison of outlier windows : diff number of ranges obtained depending on which is query & which labelled as source
# ie query s* outliers by skov bed files & vice versa
#ranges_1<-reduce(subsetByOverlaps(sstar_subs,skov_subs)) 
ranges_2<-reduce(subsetByOverlaps(skov_subs,sstar_subs))
# length distribution of overlapping regions
#lengths_1<-summary(reduce(subsetByOverlaps(sstar_subs,skov_subs))@ranges@width)/1000
lengths_2<-summary(reduce(subsetByOverlaps(skov_subs,sstar_subs))@ranges@width)/1000
# proportion overlapping
#col1 = proportion of overlapping bp for unique overlapping regions
#col2 = proportion of overlapping fragments for unique overlapping regions
#prop_1<-calc_prop_fun(sstar_subs,skov_subs)
prop_2<-calc_prop_fun(skov_subs,sstar_subs)
#return(list(ranges_1,ranges_2,lengths_1,lengths_2,prop_1,prop_2))
# only output ranges here? - RUN THROUGH INTERACTIVELY TMRW *******
return(list(ranges_2,lengths_2,prop_2))
}

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

# SEE IF THIS WORKS & GIVES EQUIVALENT RESULTS TO PREVIOUS FN ******** yes, but only outputs the ranges (not the proportions or the lengths)
# or only output ranges_2
intersect_function<-function(sstar_subs,skov_subs) {
# 1) compare GBB - 95% CI for S*
#sstar_subs #sk_regions_gb
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
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

# 3) read in the gtf (gene info)
# use makeTxDbFromGFF to create the TxDB object 
gtf="/scratch/project/devel/mkuhlwilm/chimpdiv/index/Homo_sapiens.GRCh37.75_c.gtf"

test<-makeTxDbFromGFF(gtf, format='gtf')
#Import genomic features from the file as a GRanges object ... 
#OK
#Prepare the 'metadata' data frame ... OK
#Make the TxDb object ... 
#OK
#Warning message:
#In .get_cds_IDX(mcols0$type, mcols0$phase) :
#  The "phase" metadata column contains non-NA values for features of type
#  stop_codon. This information was ignored.

#str(test)
#Reference class 'TxDb' [package "GenomicFeatures"] with 5 fields

human.genes = genes(test)

# only retain autosomal chr
autosome.hum.genes<-human.genes[seqnames(human.genes) == c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22")]

#> autosome.hum.genes
#GRanges object with 2519 ranges and 1 metadata column:
#                  seqnames              ranges strand |         gene_id
#                     <Rle>           <IRanges>  <Rle> |     <character>
#  ENSG00000002822     chr7     1855429-2272878      - | ENSG00000002822
#  ENSG00000004700    chr12   21621845-21654603      - | ENSG00000004700
#  ENSG00000004779    chr16   23592323-23607677      - | ENSG00000004779
#  ENSG00000005156    chr17   33307513-33332083      + | ENSG00000005156
#  ENSG00000006114    chr17   35874900-35969544      - | ENSG00000006114
#              ...      ...                 ...    ... .             ...
#  ENSG00000273221     chr1 111727037-111727683      - | ENSG00000273221
#  ENSG00000273293     chr7 149578448-149578669      - | ENSG00000273293
#  ENSG00000273348    chr18   19927809-19928215      - | ENSG00000273348
#  ENSG00000273398     chr2   68358370-68488362      - | ENSG00000273398
#  ENSG00000273450    chr10   96074664-96075084      - | ENSG00000273450
#  -------
#  seqinfo: 265 sequences from an unspecified genome; no seqlengths

#-----------------------------------------------------------------------------------------------------------------------


# 4) perform the overlap b/n skov & vf outliers


# to intersect volcanofinder quantile with fragments from the s*-skov overlap of each individual
overlap_vfoutl_skov2<-function(vf_obj,sk_obj) {
out_overlaps_95<-list()
out_1<-list()
for (i in (1:length(sk_obj))) {
# overlap the volcanofinder quantile with the skov fragments per individual
out_overlaps_95[[i]]<-subsetByOverlaps(vf_obj,sk_obj[[i]]) 
# intersect this overlap with human gene gtf annotation
out_1[[i]]<-intersect(autosome.hum.genes, out_overlaps_95[[i]], ignore.strand=TRUE)
}
return(out_1)
}

proc_genes2<-function(overlap_obj) {
out_1<-overlap_obj  
#out_1<-intersect(autosome.hum.genes, overlap_obj, ignore.strand=TRUE)
x<- findOverlaps(out_1,autosome.hum.genes) 
x1<-unique(as.data.frame(x)[,2])
# genes in the intersect of the 95% volcanofinder outliers with the strict skov outliers for individual 1
out_genes<-list()
for (ind in (1:length(x1))) {
out_genes[[ind]]<-autosome.hum.genes[x1[[ind]]]
}

y<-unlist(out_genes)

return(y)
}

output_genes<-function(skov_obj){

ov2<-overlap_vfoutl_skov2(vf_95_ranges,skov_obj)

# check lengths of the object
check_overlaps<-list()
for (a in (1:length(ov2))) {
check_overlaps[[a]]<-length(ov2[[a]])
}

# 0s have no overlap in ranges b/n volcanofinder & skov
# remove this from ov2 object 
#ov2<-ov2[-which(check_overlaps==0)]

# loop for all individuals which have an overlap b/n vf & skov
hold_genes2<-list()
for (a in (1:length(ov2))) {
hold_genes2[[a]]<-proc_genes2(ov2[[a]])
}

# now works
y<-unlist(hold_genes2)

return(y)
}

#o_gbb<-output_genes(ov_gbb_99)
o_gbg<-output_genes(ov_gbg_99)
#o_gbb40<-output_genes(ov_gbb40_99)
#o_gbg40<-output_genes(ov_gbg40_99)

#o_gbb<-output_genes(ov_gbb_99)
#Error in h(simpleError(msg, call)) : 
#  error in evaluating the argument 'i' in selecting a method for function '[': subscript out of bounds


o_gbg # gives output
[[1]]
GRanges object with 1 range and 1 metadata column:
                  seqnames            ranges strand |         gene_id
                     <Rle>         <IRanges>  <Rle> |     <character>
  ENSG00000212127    chr12 11090005-11324172      - | ENSG00000212127
  -------
  seqinfo: 265 sequences from an unspecified genome; no seqlengths

[[2]]
GRanges object with 1 range and 1 metadata column:
                  seqnames          ranges strand |         gene_id
                     <Rle>       <IRanges>  <Rle> |     <character>
  ENSG00000112902     chr5 9035138-9546187      - | ENSG00000112902
  -------
  seqinfo: 265 sequences from an unspecified genome; no seqlengths

[[3]]
GRanges object with 1 range and 1 metadata column:
                  seqnames            ranges strand |         gene_id
                     <Rle>         <IRanges>  <Rle> |     <character>
  ENSG00000212127    chr12 11090005-11324172      - | ENSG00000212127
  -------
  seqinfo: 265 sequences from an unspecified genome; no seqlengths

[[4]]
GRanges object with 1 range and 1 metadata column:
                  seqnames              ranges strand |         gene_id
                     <Rle>           <IRanges>  <Rle> |     <character>
  ENSG00000228215     chr1 191844625-191980390      + | ENSG00000228215
  -------
  seqinfo: 265 sequences from an unspecified genome; no seqlengths

[[5]]
GRanges object with 1 range and 1 metadata column:
                  seqnames            ranges strand |         gene_id
                     <Rle>         <IRanges>  <Rle> |     <character>
  ENSG00000212127    chr12 11090005-11324172      - | ENSG00000212127
  -------
  seqinfo: 265 sequences from an unspecified genome; no seqlengths

[[6]]
GRanges object with 1 range and 1 metadata column:
                  seqnames              ranges strand |         gene_id
                     <Rle>           <IRanges>  <Rle> |     <character>
  ENSG00000160145     chr3 123798870-124445172      + | ENSG00000160145
  -------
  seqinfo: 265 sequences from an unspecified genome; no seqlengths

[[7]]
GRanges object with 1 range and 1 metadata column:
                  seqnames              ranges strand |         gene_id
                     <Rle>           <IRanges>  <Rle> |     <character>
  ENSG00000119927    chr10 113909624-113975135      - | ENSG00000119927
  -------
  seqinfo: 265 sequences from an unspecified genome; no seqlengths

[[8]]
GRanges object with 1 range and 1 metadata column:
                  seqnames            ranges strand |         gene_id
                     <Rle>         <IRanges>  <Rle> |     <character>
  ENSG00000212127    chr12 11090005-11324172      - | ENSG00000212127
  -------
  seqinfo: 265 sequences from an unspecified genome; no seqlengths

[[9]]
GRanges object with 1 range and 1 metadata column:
                  seqnames              ranges strand |         gene_id
                     <Rle>           <IRanges>  <Rle> |     <character>
  ENSG00000119927    chr10 113909624-113975135      - | ENSG00000119927
  -------
  seqinfo: 265 sequences from an unspecified genome; no seqlengths

[[10]]
GRanges object with 1 range and 1 metadata column:
                  seqnames            ranges strand |         gene_id
                     <Rle>         <IRanges>  <Rle> |     <character>
  ENSG00000212127    chr12 11090005-11324172      - | ENSG00000212127
  -------
  seqinfo: 265 sequences from an unspecified genome; no seqlengths

[[11]]
GRanges object with 1 range and 1 metadata column:
                  seqnames              ranges strand |         gene_id
                     <Rle>           <IRanges>  <Rle> |     <character>
  ENSG00000160145     chr3 123798870-124445172      + | ENSG00000160145
  -------
  seqinfo: 265 sequences from an unspecified genome; no seqlengths

[[12]]
GRanges object with 1 range and 1 metadata column:
                  seqnames            ranges strand |         gene_id
                     <Rle>         <IRanges>  <Rle> |     <character>
  ENSG00000212127    chr12 11090005-11324172      - | ENSG00000212127
  -------
  seqinfo: 265 sequences from an unspecified genome; no seqlengths

[[13]]
GRanges object with 1 range and 1 metadata column:
                  seqnames              ranges strand |         gene_id
                     <Rle>           <IRanges>  <Rle> |     <character>
  ENSG00000119927    chr10 113909624-113975135      - | ENSG00000119927
  -------
  seqinfo: 265 sequences from an unspecified genome; no seqlengths

[[14]]
GRanges object with 1 range and 1 metadata column:
                  seqnames            ranges strand |         gene_id
                     <Rle>         <IRanges>  <Rle> |     <character>
  ENSG00000212127    chr12 11090005-11324172      - | ENSG00000212127
  -------
  seqinfo: 265 sequences from an unspecified genome; no seqlengths

[[15]]
GRanges object with 1 range and 1 metadata column:
                  seqnames            ranges strand |         gene_id
                     <Rle>         <IRanges>  <Rle> |     <character>
  ENSG00000212127    chr12 11090005-11324172      - | ENSG00000212127
  -------
  seqinfo: 265 sequences from an unspecified genome; no seqlengths

[[16]]
GRanges object with 1 range and 1 metadata column:
                  seqnames            ranges strand |         gene_id
                     <Rle>         <IRanges>  <Rle> |     <character>
  ENSG00000105997     chr7 27145803-27192200      - | ENSG00000105997
  -------
  seqinfo: 265 sequences from an unspecified genome; no seqlengths

[[17]]
GRanges object with 1 range and 1 metadata column:
                  seqnames            ranges strand |         gene_id
                     <Rle>         <IRanges>  <Rle> |     <character>
  ENSG00000212127    chr12 11090005-11324172      - | ENSG00000212127
  -------
  seqinfo: 265 sequences from an unspecified genome; no seqlengths

# add to spreadsheet *

# for the gbb:
output_genes_gbb<-function(skov_obj){

ov2<-overlap_vfoutl_skov2(vf_95_ranges,skov_obj)

# check lengths of the object
check_overlaps<-list()
for (a in (1:length(ov2))) {
check_overlaps[[a]]<-length(ov2[[a]])
}

# 0s have no overlap in ranges b/n volcanofinder & skov
# remove this from ov2 object 
ov2<-ov2[-which(check_overlaps==0)]

# loop for all individuals which have an overlap b/n vf & skov
hold_genes2<-list()
for (a in (1:length(ov2))) {
hold_genes2[[a]]<-proc_genes2(ov2[[a]])
}

# now works
y<-unlist(hold_genes2)

return(y)
}

# MG now fine
o_gbb<-output_genes_gbb(ov_gbb_99)


> o_gbb
[[1]]
GRanges object with 1 range and 1 metadata column:
                  seqnames          ranges strand |         gene_id
                     <Rle>       <IRanges>  <Rle> |     <character>
  ENSG00000112902     chr5 9035138-9546187      - | ENSG00000112902
  -------
  seqinfo: 265 sequences from an unspecified genome; no seqlengths

[[2]]
GRanges object with 1 range and 1 metadata column:
                  seqnames            ranges strand |         gene_id
                     <Rle>         <IRanges>  <Rle> |     <character>
  ENSG00000212127    chr12 11090005-11324172      - | ENSG00000212127
  -------
  seqinfo: 265 sequences from an unspecified genome; no seqlengths

[[3]]
GRanges object with 1 range and 1 metadata column:
                  seqnames              ranges strand |         gene_id
                     <Rle>           <IRanges>  <Rle> |     <character>
  ENSG00000160145     chr3 123798870-124445172      + | ENSG00000160145
  -------
  seqinfo: 265 sequences from an unspecified genome; no seqlengths

[[4]]
GRanges object with 1 range and 1 metadata column:
                  seqnames            ranges strand |         gene_id
                     <Rle>         <IRanges>  <Rle> |     <character>
  ENSG00000087128     chr4 69313167-69363322      + | ENSG00000087128
  -------
  seqinfo: 265 sequences from an unspecified genome; no seqlengths

[[5]]
GRanges object with 1 range and 1 metadata column:
                  seqnames            ranges strand |         gene_id
                     <Rle>         <IRanges>  <Rle> |     <character>
  ENSG00000087128     chr4 69313167-69363322      + | ENSG00000087128
  -------
  seqinfo: 265 sequences from an unspecified genome; no seqlengths

[[6]]
GRanges object with 1 range and 1 metadata column:
                  seqnames            ranges strand |         gene_id
                     <Rle>         <IRanges>  <Rle> |     <character>
  ENSG00000212127    chr12 11090005-11324172      - | ENSG00000212127
  -------
  seqinfo: 265 sequences from an unspecified genome; no seqlengths

[[7]]
GRanges object with 1 range and 1 metadata column:
                  seqnames            ranges strand |         gene_id
                     <Rle>         <IRanges>  <Rle> |     <character>
  ENSG00000212127    chr12 11090005-11324172      - | ENSG00000212127
  -------
  seqinfo: 265 sequences from an unspecified genome; no seqlengths

[[8]]
GRanges object with 1 range and 1 metadata column:
                  seqnames            ranges strand |         gene_id
                     <Rle>         <IRanges>  <Rle> |     <character>
  ENSG00000087128     chr4 69313167-69363322      + | ENSG00000087128
  -------
  seqinfo: 265 sequences from an unspecified genome; no seqlengths

[[9]]
GRanges object with 1 range and 1 metadata column:
                  seqnames            ranges strand |         gene_id
                     <Rle>         <IRanges>  <Rle> |     <character>
  ENSG00000212127    chr12 11090005-11324172      - | ENSG00000212127
  -------
  seqinfo: 265 sequences from an unspecified genome; no seqlengths



#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------



> length(o_gbb)
[1] 9
> length(unique(o_gbb)) # add the unique ones to the spreadsheet
[1] 4


length(o_gbg)
length(unique(o_gbg))
> length(o_gbg)
[1] 17
> length(unique(o_gbg))
[1] 6

#-----------------------------------------------------------------------------------------------------------------------

unique(o_gbb)
unique(o_gbg)


# have added these to the spreadsheet - https://docs.google.com/spreadsheets/d/1yLq5JE5OQC-yZ6E6aBDyz_4LZZwLf6PRlO0IrLcqCZo/edit#gid=1264452875

> unique(o_gbb)
[[1]]
GRanges object with 1 range and 1 metadata column:
                  seqnames          ranges strand |         gene_id
                     <Rle>       <IRanges>  <Rle> |     <character>
  ENSG00000112902     chr5 9035138-9546187      - | ENSG00000112902
  -------
  seqinfo: 265 sequences from an unspecified genome; no seqlengths

[[2]]
GRanges object with 1 range and 1 metadata column:
                  seqnames            ranges strand |         gene_id
                     <Rle>         <IRanges>  <Rle> |     <character>
  ENSG00000212127    chr12 11090005-11324172      - | ENSG00000212127
  -------
  seqinfo: 265 sequences from an unspecified genome; no seqlengths

[[3]]
GRanges object with 1 range and 1 metadata column:
                  seqnames              ranges strand |         gene_id
                     <Rle>           <IRanges>  <Rle> |     <character>
  ENSG00000160145     chr3 123798870-124445172      + | ENSG00000160145
  -------
  seqinfo: 265 sequences from an unspecified genome; no seqlengths

[[4]]
GRanges object with 1 range and 1 metadata column:
                  seqnames            ranges strand |         gene_id
                     <Rle>         <IRanges>  <Rle> |     <character>
  ENSG00000087128     chr4 69313167-69363322      + | ENSG00000087128
  -------
  seqinfo: 265 sequences from an unspecified genome; no seqlengths

> unique(o_gbg)
[[1]]
GRanges object with 1 range and 1 metadata column:
                  seqnames            ranges strand |         gene_id
                     <Rle>         <IRanges>  <Rle> |     <character>
  ENSG00000212127    chr12 11090005-11324172      - | ENSG00000212127
  -------
  seqinfo: 265 sequences from an unspecified genome; no seqlengths

[[2]]
GRanges object with 1 range and 1 metadata column:
                  seqnames          ranges strand |         gene_id
                     <Rle>       <IRanges>  <Rle> |     <character>
  ENSG00000112902     chr5 9035138-9546187      - | ENSG00000112902
  -------
  seqinfo: 265 sequences from an unspecified genome; no seqlengths

[[3]]
GRanges object with 1 range and 1 metadata column:
                  seqnames              ranges strand |         gene_id
                     <Rle>           <IRanges>  <Rle> |     <character>
  ENSG00000228215     chr1 191844625-191980390      + | ENSG00000228215
  -------
  seqinfo: 265 sequences from an unspecified genome; no seqlengths

[[4]]
GRanges object with 1 range and 1 metadata column:
                  seqnames              ranges strand |         gene_id
                     <Rle>           <IRanges>  <Rle> |     <character>
  ENSG00000160145     chr3 123798870-124445172      + | ENSG00000160145
  -------
  seqinfo: 265 sequences from an unspecified genome; no seqlengths

[[5]]
GRanges object with 1 range and 1 metadata column:
                  seqnames              ranges strand |         gene_id
                     <Rle>           <IRanges>  <Rle> |     <character>
  ENSG00000119927    chr10 113909624-113975135      - | ENSG00000119927
  -------
  seqinfo: 265 sequences from an unspecified genome; no seqlengths

[[6]]
GRanges object with 1 range and 1 metadata column:
                  seqnames            ranges strand |         gene_id
                     <Rle>         <IRanges>  <Rle> |     <character>
  ENSG00000105997     chr7 27145803-27192200      - | ENSG00000105997
  -------
  seqinfo: 265 sequences from an unspecified genome; no seqlengths

#-----------------------------------------------------------------------------------------------------------------------
#e_cand_fin<-list( unique(o_gbb), unique(o_gbg))

l1<-list(unique(o_gbb), unique(o_gbg))

unique(unlist(l1))

e_cand_fin<-unique(unlist(l1))

# save these R objects, so don't have to execute this step when intersecting with mts

# already executed:
save(e_cand_fin, file="/scratch/devel/hpawar/admix/volcanofinder/output/merged/em.el.vf.s*.skov.intersect.genes.21jul22")

# perhaps rather than this, it would be better to write out only the unique genes?
# ie reduce to 7 unique genes*
