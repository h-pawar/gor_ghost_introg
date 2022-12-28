# regenerate random regions of equal length which do have sufficient callable sites

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
library(valr) # need to generate random regions


load(file="/scratch/devel/hpawar/admix/overlap.s*.skov/overlap_objects/ov_gbb_99",verbose=T)
load(file="/scratch/devel/hpawar/admix/overlap.s*.skov/overlap_objects/ov_gbg_99",verbose=T)


#-----------------------------------------------------------------------------------------------------------------------

# 1) read in the informative regions (ie with callable sites)

informative<-read.table("/scratch/devel/hpawar/admix/overlap.s*.skov/pairwisediff/callableprop.txt")

# convert to granges
informranges<-GRanges(seqnames=informative[,1],ranges=IRanges(start=as.numeric(informative[,2]),end=as.numeric(informative[,2]+40000),names=informative[,1]),strand=rep("*",length(informative[,1])),prop=as.numeric(informative[,3]))
informative[,1]<-sapply("chr",paste, informative[,1], sep="")
informranges<-GRanges(seqnames=informative[,1],ranges=IRanges(start=as.numeric(informative[,2]),end=as.numeric(informative[,2]+40000),names=informative[,1]),strand=rep("*",length(informative[,1])),prop=as.numeric(informative[,3]))

#GRanges object with 32561 ranges and 1 metadata column:
#        seqnames              ranges strand |      prop
#           <Rle>           <IRanges>  <Rle> | <numeric>
#  chr10    chr10       160000-200000      * |  0.542675
#  chr10    chr10       240000-280000      * |  0.688575

#-----------------------------------------------------------------------------------------------------------------------
> ov_gbb_99[[1]]
GRanges object with 232 ranges and 0 metadata columns:
        seqnames            ranges strand
           <Rle>         <IRanges>  <Rle>
    [1]     chr1   5402000-5460000      *
    [2]     chr1   9786000-9884000      *
    [3]     chr1 10650000-10698000      *
    [4]     chr1 34780000-34960000      *
    [5]     chr1 37295000-37456000      *
    ...      ...               ...    ...
  [228]    chr21 40388000-40519000      *
  [229]    chr21 44500000-44619000      *
  [230]    chr22 18040000-18140000      *
  [231]    chr22 29378000-29461000      *
  [232]    chr22 36864000-36978000      *
  -------
  seqinfo: 22 sequences from an unspecified genome; no seqlengths
#-----------------------------------------------------------------------------------------------------------------------
# lengths of hg19
hg19<-fread('/home/devel/marcmont/scratch/Harvi/geneflow/test/test.reconst.ids/ora_1_1/proportion/hg19.autosomes.chrom.sizes.txt')

hg19.lengths<-apply(hg19[,2],2,function(x) as.numeric(as.character(x)))

# hg19 chr in the correct order
hg19<-hg19[c(1:18,20,19,22,21),]
colnames(hg19)<-c("chrom","size") 

#-----------------------------------------------------------------------------------------------------------------------

 ov_gbb_99[[1]][1]

# here scen=ov_gbb_99[[i]]


test_random<-list()

replace_random_function<-function(scen,i){

repeat{

replace_random_1<-bed_random(hg19, length = (scen[i]@ranges@width-1), n=1)

x<-data.frame(matrix(unlist(replace_random_1), ncol = length(replace_random_1), byrow = TRUE))

# convert to granges - then this will be what will be being intersected
test_s<-GRanges(seqnames=x[,1],ranges=IRanges(start=as.numeric(x[,2]),end=as.numeric(x[,3]),names=x[,1]),strand=rep("*",length(x[,1])))


if(length(subsetByOverlaps(test_s,informranges, invert = TRUE))==0) {
test_random[[i]]<-test_s
break
}
}
return(test_s)
}


#out<-replace_random_function(ov_gbb_99[[1]],1)
# out@ranges
#IRanges object with 1 range and 0 metadata columns:
#            start       end     width
#        <integer> <integer> <integer>
#  chr11  97707637  97765637     58001

# generate like this to generate equivalent regions for all the introg regions for the first ind

#-----------------------------------------------------------------------------------------------------------------------
hold_new_ran<-list()
for (i in (1:length(ov_gbb_99[[1]]))) {
hold_new_ran[[i]]<-replace_random_function(ov_gbb_99[[1]],i)}

rand_one<-unlist(hold_new_ran)


#> width(ov_gbb_99[[1]][232])
#[1] 114001
#> width(rand_one[232])
#[1] 114001


#random_18oct22<-list() # once
random_18oct22[[1]]<-rand_one

save(random_18oct22,file=paste("/scratch/devel/hpawar/admix/overlap.s*.skov/pairwisediff/random_18oct22",sep=""))


#-----------------------------------------------------------------------------------------------------------------------

# for rest of gbb ids

# 1) generate random regions of sufficient callable sites, with length distribution equivalent to one ind
random_oneid<-function(scen,id){
hold_new_ran<-list()
for (i in (1:length(scen[[id]]))) {
hold_new_ran[[i]]<-replace_random_function(scen[[id]],i)}
rand_out<-unlist(GRangesList(hold_new_ran))
return(rand_out)}


#x<-random_oneid(ov_gbb_99,2)


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

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
