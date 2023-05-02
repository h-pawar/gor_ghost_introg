# generate random regions of equal length which have sufficient callable sites

#-----------------------------------------------------------------------------------------------------------------------
#module load gcc/6.3.0 openssl/1.0.2q  R/4.0.1  BCFTOOLS/1.14
require(data.table)
options(stringsAsFactors=F)
library(tidyverse)
options(scipen=100)
library(adegenet)
library(GenomicRanges)
library(tidyr)
library(phangorn) 
library(ape) 
library('pegas')
library(valr)


load(file="/scratch/devel/hpawar/admix/overlap.s*.skov/overlap_objects/ov_gbb_99",verbose=T)
load(file="/scratch/devel/hpawar/admix/overlap.s*.skov/overlap_objects/ov_gbg_99",verbose=T)


#-----------------------------------------------------------------------------------------------------------------------

# 1) read in the informative regions (ie with callable sites)

informative<-read.table("/scratch/devel/hpawar/admix/overlap.s*.skov/pairwisediff/callableprop.txt")
=
informranges<-GRanges(seqnames=informative[,1],ranges=IRanges(start=as.numeric(informative[,2]),end=as.numeric(informative[,2]+40000),names=informative[,1]),strand=rep("*",length(informative[,1])),prop=as.numeric(informative[,3]))
informative[,1]<-sapply("chr",paste, informative[,1], sep="")
informranges<-GRanges(seqnames=informative[,1],ranges=IRanges(start=as.numeric(informative[,2]),end=as.numeric(informative[,2]+40000),names=informative[,1]),strand=rep("*",length(informative[,1])),prop=as.numeric(informative[,3]))


#-----------------------------------------------------------------------------------------------------------------------

hg19<-fread('/home/devel/marcmont/scratch/Harvi/geneflow/test/test.reconst.ids/ora_1_1/proportion/hg19.autosomes.chrom.sizes.txt')
hg19.lengths<-apply(hg19[,2],2,function(x) as.numeric(as.character(x)))

hg19<-hg19[c(1:18,20,19,22,21),]
colnames(hg19)<-c("chrom","size") 

#-----------------------------------------------------------------------------------------------------------------------



test_random<-list()

replace_random_function<-function(scen,i){

repeat{

replace_random_1<-bed_random(hg19, length = (scen[i]@ranges@width-1), n=1)

x<-data.frame(matrix(unlist(replace_random_1), ncol = length(replace_random_1), byrow = TRUE))

test_s<-GRanges(seqnames=x[,1],ranges=IRanges(start=as.numeric(x[,2]),end=as.numeric(x[,3]),names=x[,1]),strand=rep("*",length(x[,1])))


if(length(subsetByOverlaps(test_s,informranges, invert = TRUE))==0) {
test_random[[i]]<-test_s
break
}
}
return(test_s)
}



#-----------------------------------------------------------------------------------------------------------------------
hold_new_ran<-list()
for (i in (1:length(ov_gbb_99[[1]]))) {
hold_new_ran[[i]]<-replace_random_function(ov_gbb_99[[1]],i)}

rand_one<-unlist(hold_new_ran)



#random_18oct22<-list() # once
random_18oct22[[1]]<-rand_one

save(random_18oct22,file=paste("/scratch/devel/hpawar/admix/overlap.s*.skov/pairwisediff/random_18oct22",sep=""))


#-----------------------------------------------------------------------------------------------------------------------


# 1) generate random regions of sufficient callable sites, with length distribution equivalent to one ind
random_oneid<-function(scen,id){
hold_new_ran<-list()
for (i in (1:length(scen[[id]]))) {
hold_new_ran[[i]]<-replace_random_function(scen[[id]],i)}
rand_out<-unlist(GRangesList(hold_new_ran))
return(rand_out)}


rem<-seq(2:12)+1

# 2) for all gbb ids (MG)
for (i in 1:length(rem)){
val<-rem[i]
random_18oct22[[val]]<-random_oneid(ov_gbb_99,val)
}


save(random_18oct22,file=paste("/scratch/devel/hpawar/admix/overlap.s*.skov/pairwisediff/random_18oct22",sep=""))
#-----------------------------------------------------------------------------------------------------------------------


# 3) for all gbg ids (EL)
for (i in 1:length(ov_gbg_99)){
val<-12+i
random_18oct22[[val]]<-random_oneid(ov_gbg_99,i)
}


save(random_18oct22,file=paste("/scratch/devel/hpawar/admix/overlap.s*.skov/pairwisediff/random_18oct22",sep=""))



q()

