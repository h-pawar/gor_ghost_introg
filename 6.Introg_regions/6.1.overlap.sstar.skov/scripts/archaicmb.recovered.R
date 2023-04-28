# Archaic Mb of the high-confidence putative introgressed regions (overlap of S* -Skov HMM outliers)

#module load gcc/6.3.0 openssl/1.0.2q  R/4.0.1  BCFTOOLS/1.14
require(data.table)
options(stringsAsFactors=F)
library(tidyverse)
options(scipen=100)
library(adegenet)
library(GenomicRanges)
load(file="/scratch/devel/hpawar/admix/overlap.s*.skov/overlap_objects/ov_gbb_99",verbose=T)
load(file="/scratch/devel/hpawar/admix/overlap.s*.skov/overlap_objects/ov_gbg_99",verbose=T)

gg<-c(ov_gbg_99[[1]],ov_gbg_99[[2]],ov_gbg_99[[3]],ov_gbg_99[[4]],ov_gbg_99[[5]],ov_gbg_99[[6]],ov_gbg_99[[7]],ov_gbg_99[[8]],ov_gbg_99[[9]])
gb<-c(ov_gbb_99[[1]],ov_gbb_99[[2]],ov_gbb_99[[3]],ov_gbb_99[[4]],ov_gbb_99[[5]],ov_gbb_99[[6]],ov_gbb_99[[7]],ov_gbb_99[[8]],ov_gbb_99[[9]],ov_gbb_99[[10]],ov_gbb_99[[11]],ov_gbb_99[[12]])

sum(reduce(gg)@ranges@width)/1000000
#[1] 115.1618
# ~115 Mbp in GBG
sum(reduce(gb)@ranges@width)/1000000
#[1] 152.5841
# ~153 Mbp in GBB
sum(reduce(c(gb,gg))@ranges@width)/1000000
#[1] 213.9295
# ~214 Mbp


sum(reduce(c(gb,gg))@ranges@width)/2881033286*100
#[1] 7.425443

# for action of the genome (minus regions without data) 
# 1) read in the informative regions (ie with callable sites)

informative<-read.table("/scratch/devel/hpawar/admix/overlap.s*.skov/pairwisediff/callableprop.txt")

# convert to granges
informranges<-GRanges(seqnames=informative[,1],ranges=IRanges(start=as.numeric(informative[,2]),end=as.numeric(informative[,2]+40000),names=informative[,1]),strand=rep("*",length(informative[,1])),prop=as.numeric(informative[,3]))
informative[,1]<-sapply("chr",paste, informative[,1], sep="")
informranges<-GRanges(seqnames=informative[,1],ranges=IRanges(start=as.numeric(informative[,2]),end=as.numeric(informative[,2]+40000),names=informative[,1]),strand=rep("*",length(informative[,1])),prop=as.numeric(informative[,3]))

sum(reduce(informranges)@ranges@width)/1000000



sum(reduce(gg)@ranges@width)/sum(reduce(informranges)@ranges@width)
sum(reduce(gb)@ranges@width)/sum(reduce(informranges)@ranges@width)
sum(reduce(c(gb,gg))@ranges@width)/sum(reduce(informranges)@ranges@width)
