# Mon 20 Mar 2023 15:50:03 CET
#get random regions, like before, and ask how many of the sEs in these regions are genic etc., and do this, like before, 100 times
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

# random regions have amended following relevant sections of  /scratch/devel/hpawar/admix/sstar/scripts/simul_postabc/downstream/gene.density.s*.skov.overlap.21jul22.R
#-----------------------------------------------------------------------------------------------------------------------

## module load gcc/6.3.0 openssl/1.0.2q R/4.0.1
#library(GenomicScores)
require(data.table)
options(stringsAsFactors=F)
options(scipen=100)
library(mgcv)
library(GenomicRanges)
library(tidyr)
library(ggbio)
#library(pheatmap)
#library(ggplot2)
#library(ggpubr)
library(valr)
library(GenomicFeatures)

#-----------------------------------------------------------------------------------------------------------------------


# read in introgressed regions
load(file="/scratch/devel/hpawar/admix/overlap.s*.skov/overlap_objects/ov_gbb_99",verbose=T)
load(file="/scratch/devel/hpawar/admix/overlap.s*.skov/overlap_objects/ov_gbg_99",verbose=T)

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

# (2) read in gtf of gene coordinates from hg19 reference (protein coding bp)
    # initially just use the human annotation (rather than filtering by orthologs)

# following - https://www.biostars.org/p/169171/ to convert gtf file -> granges object

# need 2 objects : 1 with coordinates of genes, 2nd object of outlier regions

# generate object- with coordinates of genes

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
#Warning message:","chr17","chr18","chr19","chr20","chr21","chr22")]
#In e1 == Rle(e2) :
#  longer object length is not a multiple of shorter object length

# data structure - contains gene coor
#str(human.genes)
#-----------------------------------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------------------------------
# lengths of hg19
hg19<-fread('/home/devel/marcmont/scratch/Harvi/geneflow/test/test.reconst.ids/ora_1_1/proportion/hg19.autosomes.chrom.sizes.txt')

hg19.lengths<-apply(hg19[,2],2,function(x) as.numeric(as.character(x)))

# hg19 chr in the correct order
hg19<-hg19[c(1:18,20,19,22,21),]
colnames(hg19)<-c("chrom","size") 
# or whether need to read in the genome itself?
  # path to hg19 ref

#-----------------------------------------------------------------------------------------------------------------------


# random regions of sufficient callability already generated -> read these in & assess protein coding bp in these regions *
    # will need to generate 100 * such sets for this analysis **

# random regions: 1-12 MG, 13-21 EL
#load(file="/scratch/devel/hpawar/admix/overlap.s*.skov/pairwisediff/random_18oct22",verbose=T)

# but now need to generate 100 iterations for each random rep **
    #  yes will need to amend generate.randomreg.sufficallable.18oct22.R # to generate 100 times *
        # & query this with the genes *
# run through the original script interactively to see where to amend **

# WILL NEED TO RUN THROUGH INTERACTIVELY<, TO CHECK IS WORKING *

# 1) read in the informative regions (ie with callable sites)

informative<-read.table("/scratch/devel/hpawar/admix/overlap.s*.skov/pairwisediff/callableprop.txt")

# convert to granges
informranges<-GRanges(seqnames=informative[,1],ranges=IRanges(start=as.numeric(informative[,2]),end=as.numeric(informative[,2]+40000),names=informative[,1]),strand=rep("*",length(informative[,1])),prop=as.numeric(informative[,3]))
informative[,1]<-sapply("chr",paste, informative[,1], sep="")
informranges<-GRanges(seqnames=informative[,1],ranges=IRanges(start=as.numeric(informative[,2]),end=as.numeric(informative[,2]+40000),names=informative[,1]),strand=rep("*",length(informative[,1])),prop=as.numeric(informative[,3]))

replace_random_function<-function(scen,i){
test_random<-list()

repeat{

replace_random_1<-bed_random(hg19, length = (scen[i]@ranges@width-1), n=1)

x<-data.frame(matrix(unlist(replace_random_1), ncol = length(replace_random_1), byrow = TRUE))

# convert to granges - then this will be what will be being intersected
test_s<-GRanges(seqnames=x[,1],ranges=IRanges(start=as.numeric(x[,2]),end=as.numeric(x[,3]),names=x[,1]),strand=rep("*",length(x[,1])))

# ie random regions need to overlap with informative ranges (ie obtain random regions of sufficient callable sites)
if(length(subsetByOverlaps(test_s,informranges, invert = TRUE))==0) {
test_random[[i]]<-test_s
break
}
}
return(test_s)
}
#-----------------------------------------------------------------------------------------------------------------------

# 1) generate random regions of sufficient callable sites, with length distribution equivalent to one ind
random_oneid<-function(scen,id){
hold_new_ran<-list()
for (i in (1:length(scen[[id]]))) {
hold_new_ran[[i]]<-replace_random_function(scen[[id]],i)}
rand_out<-unlist(GRangesList(unlist(hold_new_ran)))
return(rand_out)}


#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

gor_reg=read.table("/scratch/devel/pesteller/gor_reg_4Harvi/Gorilla.re_annotation.with_Non-RE.hg19.bed")
gor_reg<-as.data.frame(gor_reg)

# sE regions
gor_sereg<-gor_reg[(which(gor_reg$V7=="sE")),]

# shoudl first filter gor_sereg by autosomes only ***
#gor_allreg<-gor_reg[-(which(gor_reg$V7=="P/Non-re" | gor_reg$V7=="E/Non-re"| gor_reg$V7=="Non-re")),]

#nrow(gor_sereg)
#[1] 16311

filt<- unique(gor_sereg$V1)[c(5,19:25,28,30,32,33)]

gor_sereg<-gor_sereg[-(which(gor_sereg$V1==filt[[1]] | gor_sereg$V1==filt[[2]] | gor_sereg$V1==filt[[3]]  | gor_sereg$V1==filt[[4]]   | gor_sereg$V1==filt[[5]]    | gor_sereg$V1==filt[[6]]     | gor_sereg$V1==filt[[7]]      | gor_sereg$V1==filt[[8]]       | gor_sereg$V1==filt[[9]]        | gor_sereg$V1==filt[[10]]         | gor_sereg$V1==filt[[11]]  | gor_sereg$V1==filt[[12]])),]

# nrow(gor_sereg)
#[1] 16295

# add column of genic/promoter interacting etc
gor_seregranges<-GRanges(seqnames=gor_sereg[,1],ranges=IRanges(start=as.numeric(gor_sereg[,2]),end=as.numeric(gor_sereg[,3]),names=gor_sereg[,1]),strand=rep("*",length(gor_sereg[,1])),re=(gor_sereg[,7]),pos=(gor_sereg[,8]))



#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------


# CHECK IF THESE FUNCTIONS WORK **


# where here scen would be the random rep 1 of multiple ids
proc_perid_function<-function(scen,id) {
rtes<-findOverlaps(gor_seregranges,scen[[id]])
y1<-unique(as.data.frame(rtes)[,1])
test<-data.frame(gor_seregranges[y1])
out<-table(test$pos)
return(out)}

proc_p<-function(scen){
pos_perid<-list()
for (i in 1:length(scen)){
print(i)
pos_perid[[i]]<-proc_perid_function(scen,i)}
out<-do.call(rbind,pos_perid)
return(out)}

#-----------------------------------------------------------------------------------------------------------------------


# run throguh interactively - think need this format - for this iteration
process_onerandomrep<-function(scen){
random_rep1_MG_allids<-list()
for (val in 1:length(scen)){
random_rep1_MG_allids[[val]]<-random_oneid(scen,val)
}

#intersect_withgenes_function(random_rep1_MG_allids) # think do not need this step...
tes<-proc_p(random_rep1_MG_allids)
return(tes)
}

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

#process_onerandomrep(ov_gbb_99) # see if runs fine...


#      EiE  gE  P PiE prE
# [1,]  60   8 10  14  60
# [2,]  12 131 12  12  14
# [3,]   3 100 15  32   8
# [4,]   1  88 10  13  13
# [5,]   3 121 10  21  16
# [6,]   8 106 11  24  10
# [7,]   4 111 11  19  14
# [8,]   4 103 10  21  10
# [9,]   4 110  9  28  17
#[10,]   2  49  3  14   7
#[11,]   3  87 13  19  10
#[12,]   7 117  8  28   8

# yes output is generated


# see if this is useful.. -> coudl be good as an initial sanity check
# coudl be good to output also the random regions
#process_onerandomrep<-function(scen){
#random_rep1_MG_allids<-list()
#for (val in 1:length(scen)){
#random_rep1_MG_allids[[val]]<-random_oneid(scen,val)
#}

#intersect_withgenes_function(random_rep1_MG_allids) # think do not need this step...
#tes<-proc_p(random_rep1_MG_allids)
#tes1<-list(random_rep1_MG_allids,tes)
#return(tes1)
#}

#test<-process_onerandomrep(ov_gbb_99)

# then check if the functions - give the expected output *
#[[2]]
#      EiE  gE  P PiE prE
# [1,]   7  91  3  19  12
# [2,]  11 128 13  13  14
# [3,]   5 117  9  24  17
# [4,]   2  91  9  15  11
# [5,]   3 116 11  20  11
# [6,]   5 110 20  16  16
# [7,]  11 119 16  14  14
# [8,]   5 106 13  22  18
# [9,]   3 129 13  23   9
#[10,]   4  99  8  17   7
#[11,]   7 102  6  20  11
#[12,]   8 104  9  22  15


#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------


# but would only need random regions for MG *

ran_setype_gbb_20mar23<-list()
for (i in 1:100){
val<-i
print(val)
ran_setype_gbb_20mar23[[val]]<-process_onerandomrep(ov_gbb_99)}

# is taking a while to run - may be better to send as a job,..
save(ran_setype_gbb_20mar23,file=paste("/scratch/devel/hpawar/admix/overlap.s*.skov/revisions_8feb23/ran_setype_gbb_20mar23"))

#-----------------------------------------------------------------------------------------------------------------------
