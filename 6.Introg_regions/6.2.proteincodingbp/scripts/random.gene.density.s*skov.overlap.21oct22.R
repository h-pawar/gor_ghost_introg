# Thu 20 Oct 2022 16:35:59 CEST
# generate random regions of sufficient callable sites & recalc proportion of protein coding bp in these random regions

#-----------------------------------------------------------------------------------------------------------------------

# NEEDS AMENDING/REWRITING ** done

# me

#Analyses using random regions: 
#- bp overlap in introgressed vs random regions
#- gene density in introgressed vs random regions
#- pcas, nj trees of introgressed vs random regions (for all introgressed regions, for private regions)
#- pairwise differences in introgressed vs random regions.


#Fri 14 Oct 2022 11:06:55 CEST
# MK
#ok, then that was all. 
#Gene density and pairwise differences are most relevant. bp overlap of regions is quite neutral, 
#but once you generated random regions easy to test again. As I said, PCAs and NJ trees patterns should not be affected, 
#only branch lengths or PC values, but that is not of priority for now.

#-----------------------------------------------------------------------------------------------------------------------


# 6.2) calculate % protein coding bp in overlap (sstar-skov)
# scripts
# calc for empirical & random in the same script *
# /scratch/devel/hpawar/admix/sstar/scripts/simul_postabc/downstream/gene.density.s*.skov.overlap.21jul22.R
# /scratch/devel/hpawar/admix/sstar/scripts/simul_postabc/downstream/gene.density.s*.skov.overlap.21jul22.arr
# p value calculated in: plot.gene.density.skov.20jul22.R 

# ran interactively - to generate violin plots of means of random regions
#/Volumes/"Ultra USB 3.0"/IBE/further.analysis.feb.2020/gorillas/abc/simul_postabc/plot.gene.density.skov.20jul22.R
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

# ie rewrite section for random reps **

# see whether need to send as a job, or if can run interactively () - think may be too slow to run interactively for 100 iterations


# amend relevant sections of  /scratch/devel/hpawar/admix/sstar/scripts/simul_postabc/downstream/gene.density.s*.skov.overlap.21jul22.R

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


#x<-random_oneid(ov_gbb_99,2)

#x
#GRanges object with 287 ranges and 0 metadata columns:
#        seqnames              ranges strand
#           <Rle>           <IRanges>  <Rle>
#  chr15    chr15   46441889-46571889      *
#  chr13    chr13   93013420-93128420      *
 #  chr5     chr5 147491893-147604893      *
 #  chr2     chr2 234783442-235038442      *
 #  chr6     chr6   90957804-91061804      *
 #   ...      ...                 ...    ...
 #  chr1     chr1 186780857-186969857      *
 #  chr5     chr5 103948190-104006190      *
 #  chr1     chr1     1799531-1912531      *
 #  chr3     chr3   25815457-26035457      *
 #  chr2     chr2 175985443-176344443      *
 # -------
 # seqinfo: 22 sequences from an unspecified genome; no seqlengths

#rem<-seq(2:12)+1

# 2) for all gbb ids (MG)
#for (val in 1:12){
#random_18oct22[[val]]<-random_oneid(ov_gbb_99,val)
#}

#-----------------------------------------------------------------------------------------------------------------------

#random_rep1_MG_allids<-list()
#for (val in 1:12){
#random_rep1_MG_allids[[val]]<-random_oneid(ov_gbb_99,val)
#}

#tes<-intersect_withgenes_function(random_rep1_MG_allids)
# in principle this should work.. - has now worked

#x1<-do.call(rbind,tes[[2]])
#y1<-x1[,1]/x1[,2]

#> y1
# [1] 0.02050771 0.02808278 0.03016045 0.01881241 0.02990550 0.04963271
# [7] 0.03148582 0.01726580 0.03832830 0.01103468 0.01891198 0.03306175
#> length(y1)
#[1] 12

# ie this is the same output as previously

#-----------------------------------------------------------------------------------------------------------------------
intersect_withgenes_function<-function(empirical_id) {
intersect_id<-list()
counts_id<-list()
for (ind in (1:length(empirical_id))) {
intersect_id[[ind]]<-intersect(autosome.hum.genes, empirical_id[[ind]], ignore.strand=TRUE)
counts_id[[ind]]<-cbind(sum(reduce(intersect_id[[ind]])@ranges@width), sum(empirical_id[[ind]]@ranges@width))
}
return(list(intersect_id,counts_id))
}
#-----------------------------------------------------------------------------------------------------------------------

# write as a function to run 100 times, & shoudl cat out which rep & id is in the iteration *


process_onerandomrep<-function(scen){
random_rep1_MG_allids<-list()
for (val in 1:length(scen)){
random_rep1_MG_allids[[val]]<-random_oneid(scen,val)
}
tes<-intersect_withgenes_function(random_rep1_MG_allids)
x1<-do.call(rbind,tes[[2]])
y1<-x1[,1]/x1[,2]
return(y1)
}

#-----------------------------------------------------------------------------------------------------------------------
#- has run fine
# perhpas see what looks like for 1:2
#reps_10<-list()
#for (i in (1:2)) {
#print(i)   
#reps_10[[i]]<-process_onerandomrep(ov_gbb_99)}



#> reps_10
#[[1]]
# [1] 0.02704863 0.02953176 0.02692377 0.02793623 0.02214579 0.02543060
# [7] 0.03628739 0.04218756 0.02990315 0.02407920 0.03192421 0.02582999

#[[2]]
# [1] 0.025573027 0.023078458 0.029442791 0.036790433 0.023455184 0.024085658
# [7] 0.017842987 0.026205881 0.034167955 0.003779222 0.058166442 0.036460200


# looks fine

# if fine scale to 100 reps & to EL also ***

# once
#ran_gbb_22oct22<-list()
#save(ran_gbb_22oct22,file=paste("/scratch/devel/hpawar/admix/overlap.s*.skov/ran_gbb_22oct22"))


#ran_gbb_22oct22[[1]]<-reps_10[[1]]
#ran_gbb_22oct22[[2]]<-reps_10[[2]]
#-----------------------------------------------------------------------------------------------------------------------


# For MG
# in rest of jobs load this R object & run *
load(file=paste("/scratch/devel/hpawar/admix/overlap.s*.skov/ran_gbb_22oct22"),verbose=T)

rem<-seq(3:100)+2

for (i in 1:length(rem)){
val<-rem[i]
print(val)
ran_gbb_22oct22[[val]]<-process_onerandomrep(ov_gbb_99)}
save(ran_gbb_22oct22,file=paste("/scratch/devel/hpawar/admix/overlap.s*.skov/ran_gbb_22oct22"))

#-----------------------------------------------------------------------------------------------------------------------

# For EL
# run once
ran_gbg_22oct22<-list()
save(ran_gbg_22oct22,file=paste("/scratch/devel/hpawar/admix/overlap.s*.skov/ran_gbg_22oct22"))

for (i in 1:100){
print(i)
ran_gbg_22oct22[[i]]<-process_onerandomrep(ov_gbg_99)}

save(ran_gbg_22oct22,file=paste("/scratch/devel/hpawar/admix/overlap.s*.skov/ran_gbg_22oct22"))


# may need to send as a job *
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
