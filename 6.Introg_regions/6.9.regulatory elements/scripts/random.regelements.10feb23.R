# Fri 10 Feb 2023 16:43:23 CET
# generate random regions of equivalent length & callability
# & assess for regulatory regions

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

# have amended following relevant sections of  /scratch/devel/hpawar/admix/sstar/scripts/simul_postabc/downstream/gene.density.s*.skov.overlap.21jul22.R
#-----------------------------------------------------------------------------------------------------------------------

## module load gcc/6.3.0 openssl/1.0.2q R/4.0.1
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
# empirical overlap b/n sstar-skov
load(file="/scratch/devel/hpawar/admix/overlap.s*.skov/overlap_objects/ov_gbb_99",verbose=T)
load(file="/scratch/devel/hpawar/admix/overlap.s*.skov/overlap_objects/ov_gbg_99",verbose=T)


#-----------------------------------------------------------------------------------------------------------------------

# gorilla annotated regulatory regions (generated in garcia-perez et al, & lifted over to hg19 by PEC)

gor_reg=read.table("/scratch/devel/pesteller/gor_reg_4Harvi/Gorilla.re_annotation.with_Non-RE.hg19.bed")

gor_reg<-as.data.frame(gor_reg)

gor_reg=read.table("/scratch/devel/pesteller/gor_reg_4Harvi/Gorilla.re_annotation.with_Non-RE.hg19.bed")

gor_reg<-as.data.frame(gor_reg)


# PEC: Only focus of column 7, Column 7 is the consensus RE of the SPECIES --> The one column you want

#nrow(gor_reg)
#[1] 68441
#-----------------------------------------------------------------------------------------------------------------------

# retain only those with annotated regulatory function in both replicates *
gor_allreg<-gor_reg[-(which(gor_reg$V7=="P/Non-re" | gor_reg$V7=="E/Non-re"| gor_reg$V7=="Non-re")),]

gor_allregranges<-GRanges(seqnames=gor_allreg[,1],ranges=IRanges(start=as.numeric(gor_allreg[,2]),end=as.numeric(gor_allreg[,3]),names=gor_allreg[,1]),strand=rep("*",length(gor_allreg[,1])),re=(gor_allreg[,7]))



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

#hold_new_ran<-list()
#for (i in (1:length(ov_gbb_99[[1]]))) {
#hold_new_ran[[i]]<-replace_random_function(ov_gbb_99[[1]],i)}

#rand_one<-unlist(hold_new_ran)


#GRangesList(unlist(hold_new_ran))

#unlist(GRangesList(unlist(hold_new_ran)))
#GRanges object with 232 ranges and 0 metadata columns:
#        seqnames              ranges strand
#           <Rle>           <IRanges>  <Rle>
#   chr3     chr3   66777810-66835810      *
#  chr12    chr12 132627514-132725514      *
#  chr10    chr10   80698431-80746431      *
# this is the desired format *

# one random rep : would need to loop over these functions 100 times
#-----------------------------------------------------------------------------------------------------------------------


# 1) generate random regions of sufficient callable sites, with length distribution equivalent to one ind
random_oneid<-function(scen,id){
hold_new_ran<-list()
for (i in (1:length(scen[[id]]))) {
hold_new_ran[[i]]<-replace_random_function(scen[[id]],i)}
rand_out<-unlist(GRangesList(unlist(hold_new_ran)))
return(rand_out)}


# testing



#x<-random_oneid(ov_gbb_99,1)

proc_ran_1id<-function(scen,id){
x<-random_oneid(scen,id)
rant<-intersect(gor_allregranges, x, ignore.strand=TRUE)
# output rancounts: how many bp in gorilla annotated regulatory regions 
rancounts<-cbind(sum(reduce(rant)@ranges@width), sum(x@ranges@width))

rtes<-findOverlaps(gor_allregranges,x)
y1<-unique(as.data.frame(rtes)[,1])
# extract mutations from vep data within range of overlap with the candidate genic region
out_re<-list()
for (ind in (1:length(y1))) {
out_re[[ind]]<-data.frame(gor_allregranges[y1[[ind]]])
}
reg_info_1id<-do.call(rbind,out_re)
paths = by(reg_info_1id, reg_info_1id[,"re"], function(x) x)
counts_perel<-list()
for (ind in (1:length(paths))) {
counts_perel[[ind]]<-cbind(sum(paths[[ind]][,4]-1),  unique(paths[[ind]][[6]]))}
p<-as.data.frame(do.call(rbind, counts_perel))
p[,1]<-as.numeric(p[,1])
p[,3]<-p[,1]/sum(x@ranges@width)
# ratio of bp per regulatory element (# bp per reg element / total length of random region)
p[,c(2:3)] # output this
return(list(x,rancounts,p[,c(2:3)]))
}



process_onerandomrep<-function(scen){
random_rep1_allids<-list()
for (val in 1:length(scen)){
random_rep1_allids[[val]]<-proc_ran_1id(scen,val)
}
return(random_rep1_allids)
}

#-----------------------------------------------------------------------------------------------------------------------
# yes works for 1 rep - this section was run successfully, now commenting out to run for all ids & all random reps 

#t_o<-proc_ran_1id(ov_gbb_99,1)
# check if works, may not be necessary to output also the random regions.. (the x object: that may be making it too heavy?)


# str( test_all1pop)
#List of 12
# $ :List of 3
#  ..$ :Formal class 'GRanges' [package "GenomicRanges"] with 7 slots


#out_re<-list()
#for (ind in (1:length(y1))) {
#out_re[[ind]]<-data.frame(gor_allregranges[y1[[ind]]])
#}


#process_onerandomrep<-function(scen){
#random_rep1_allids<-list()
#for (val in 1:length(scen)){
#random_rep1_allids[[val]]<-proc_ran_1id(scen,val)
#}
#return(random_rep1_allids)
#}


#test_all1pop<-process_onerandomrep(ov_gbb_99)
# check if this runs 
    # think an echo here would be useful.. to check is going to the next step

# may need to think of a better output format, but in principle this shoudl be fine...



#per_pro<-list()
#for (val in 1:length(test_all1pop)){
#per_pro[[val]]<-test_all1pop[[val]][[2]]
#}
#do.call(rbind,per_pro)

#       [,1]     [,2]
# [1,] 1415274 28826232
# [2,] 1802994 34622287
# [3,] 1741899 34215269
# [4,] 1327560 29988231
# [5,] 1903651 32824258
# [6,] 1561659 32871263
# [7,] 1174557 31186259
# [8,] 1712152 33132257
# [9,] 1793768 34663269
#[10,]  942029 19783171
#[11,] 1053681 26814224
#[12,] 1687588 32461259


#output this 



#vper_pro<-list()
#for (val in 1:length(test_all1pop)){
#vper_pro[[val]]<-test_all1pop[[val]][[3]]
#}
# directly output as list
# but this may be more manageable


#list(do.call(rbind,per_pro),vper_pro)
# add this as list element 1 & proceed ...

# could also output the random regions...

#do.call(cbind, vper_pro)[,c(1,seq(2,length(input)*2,2))]
#Error in data.frame(..., check.names = FALSE) : 
#  arguments imply differing number of rows: 9, 8
# & deal with this further downstream

#rper_pro<-list()
#for (val in 1:length(test_all1pop)){
#rper_pro[[val]]<-test_all1pop[[val]][[1]]
#}

#list(rper_pro,do.call(rbind,per_pro),vper_pro)

# i think this output is more manageable *


#/scratch/devel/hpawar/admix/abc/simul/test/ghost/revisions_8feb23/
#mkdir -p /scratch/devel/hpawar/admix/overlap.s*.skov/revisions_8feb23

# just once 
##ran_reg_gbb_10oct23<-list()
##ran_reg_gbb_10oct23[[1]]<-list(rper_pro,do.call(rbind,per_pro),vper_pro)

##save(ran_reg_gbb_10oct23,file=paste("/scratch/devel/hpawar/admix/overlap.s*.skov/revisions_8feb23/ran_reg_gbb_10oct23"))

# write these last steps as a function **
# & write to run for rest of the 100 reps for MG

# then for 1-100 for EL


#q()



#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
# could potentially send as a job?



out_onerandomrep_fun<-function(scen){

test_all1pop<-process_onerandomrep(scen)


per_pro<-list()
for (val in 1:length(test_all1pop)){
per_pro[[val]]<-test_all1pop[[val]][[2]]
}

vper_pro<-list()
for (val in 1:length(test_all1pop)){
vper_pro[[val]]<-test_all1pop[[val]][[3]]
}

rper_pro<-list()
for (val in 1:length(test_all1pop)){
rper_pro[[val]]<-test_all1pop[[val]][[1]]
}

return(list(rper_pro,do.call(rbind,per_pro),vper_pro))
}


# tes the function
#tf<-out_onerandomrep_fun(ov_gbb_99)
# in progress **


# to run **
#ran_reg_gbb_10oct23[[2]]<-tf # also worked
#save(ran_reg_gbb_10oct23,file=paste("/scratch/devel/hpawar/admix/overlap.s*.skov/revisions_8feb23/ran_reg_gbb_10oct23"))


# For MG
# in rest of jobs load this R object & run *
load(file=paste("/scratch/devel/hpawar/admix/overlap.s*.skov/revisions_8feb23/ran_reg_gbb_10oct23"),verbose=T)

# ran rep 3 interactively - fixed the issue *

rem<-seq(4:100)+2


# RUN IN PROGRESS (TAKING A WHILE..)  - send as a job

for (i in 1:length(rem)){
val<-rem[i]
print(val)
ran_reg_gbb_10oct23[[val]]<-out_onerandomrep_fun(ov_gbb_99)}

# is taking a while to run - may be beetter to send as a job,..
save(ran_reg_gbb_10oct23,file=paste("/scratch/devel/hpawar/admix/overlap.s*.skov/revisions_8feb23/ran_reg_gbb_10oct23"))

#-----------------------------------------------------------------------------------------------------------------------




# TO RUN ***

# for EL

ran_reg_gbg_10oct23<-list()
# only once 
save(ran_reg_gbg_10oct23,file=paste("/scratch/devel/hpawar/admix/overlap.s*.skov/revisions_8feb23/ran_reg_gbg_10oct23"))


for (i in 1:100){
val<-i
print(val)
ran_reg_gbg_10oct23[[val]]<-out_onerandomrep_fun(ov_gbg_99)}

save(ran_reg_gbg_10oct23,file=paste("/scratch/devel/hpawar/admix/overlap.s*.skov/revisions_8feb23/ran_reg_gbg_10oct23"))

# i think this shoudl work 

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------


q()

