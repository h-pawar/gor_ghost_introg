
# Fri 24 Feb 2023 16:33:26 CET
# AMEND BELOW re polyphen


# Thu 23 Feb 2023 10:07:58 CET

# AMEND BELOW ***
# generate random regions of equivalent length & callability
# & assess for high/low gerp mutations


# have extracted gerp muts of interest which are polymorphic in E & intersect this with the random regions


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
# empirical overlap b/n sstar-skov - introgressed regions
load(file="/scratch/devel/hpawar/admix/overlap.s*.skov/overlap_objects/ov_gbb_99",verbose=T)
load(file="/scratch/devel/hpawar/admix/overlap.s*.skov/overlap_objects/ov_gbg_99",verbose=T)


# read in polyphen scores
load(file=paste("/scratch/devel/hpawar/admix/overlap.s*.skov/revisions_8feb23/polyphen/damage_polyphen"),verbose=T)
load(file=paste("/scratch/devel/hpawar/admix/overlap.s*.skov/revisions_8feb23/polyphen/benign_polyphen"),verbose=T)


damage_polyphen_scoresList<-GRangesList(damage_polyphen) 
all_damage_polyphen<-(unlist(damage_polyphen_scoresList))

benign_polyphen_scoresList<-GRangesList(benign_polyphen) 
all_benign_polyphen<-(unlist(benign_polyphen_scoresList))

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



# AMEND THIS FUNCTION *

#x<-random_oneid(ov_gbb_99,1)

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------



tesratiogerp_perind<-function(scen,gerp_rang1,gerp_rang2){
# take hits per ind / width of introg regions per ind : ie high gerp bp in introgressed regions 
        # how many high impact gerp mutations which are polymorphic in E gorillas
hold<-list()
for (ind in (1:length(scen))) {
hold[[ind]]<-length(findOverlaps(gerp_rang1,scen[[ind]]))/sum(width(scen[[ind]])-1)
}

# how many high impact gerp mutations 
hold1<-list()
for (ind in (1:length(scen))) {
hold1[[ind]]<-length(findOverlaps(gerp_rang1,scen[[ind]]))
}

# how many low impact gerp mutations
hold2<-list()
for (ind in (1:length(scen))) {
hold2[[ind]]<-length(findOverlaps(gerp_rang2,scen[[ind]]))
}

# of number of gerp sites identified - what proportion are high gerp scores
p_counts<-list()
for (ind in (1:length(scen))) {
p_counts[[ind]]<-hold1[[ind]]/(hold1[[ind]]+hold2[[ind]])}

pout<-list(unlist(hold),unlist(p_counts))

return(pout)

}


process_onerandomrep<-function(scen){
random_rep1_MG_allids<-list()
for (val in 1:length(scen)){
random_rep1_MG_allids[[val]]<-random_oneid(scen,val)
}
tes<-tesratiogerp_perind(random_rep1_MG_allids,all_damage_polyphen,all_benign_polyphen)
x1<- cbind(mean(tes[[1]]), mean(tes[[2]]))
return(x1)
}



#test_all1pop<-process_onerandomrep(ov_gbb_99)
    # if runs fine, then amend rest of script & send as a job *


#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------


# Mon 27 Feb 2023 10:30:33 CET - this section has already been run


# only once
#ran_polyphen_gbb_24feb23<-list()
#save(ran_polyphen_gbb_24feb23,file=paste("/scratch/devel/hpawar/admix/overlap.s*.skov/revisions_8feb23/polyphen/ran_polyphen_gbb_24feb23"))


# RUN IN PROGRESS (TAKING A WHILE..) 

#for (i in 1:100){
#val<-i
#print(val)
#ran_polyphen_gbb_24feb23[[val]]<-process_onerandomrep(ov_gbb_99)}

# is taking a while to run - may be beetter to send as a job,..
#save(ran_polyphen_gbb_24feb23,file=paste("/scratch/devel/hpawar/admix/overlap.s*.skov/revisions_8feb23/polyphen/ran_polyphen_gbb_24feb23"))

#-----------------------------------------------------------------------------------------------------------------------




# TO RUN ***

# for EL

ran_polyphen_gbg_24feb23<-list()
# only once 
#save(ran_polyphen_gbg_24feb23,file=paste("/scratch/devel/hpawar/admix/overlap.s*.skov/revisions_8feb23/sift/ran_polyphen_gbg_24feb23"))


for (i in 1:100){
val<-i
print(val)
ran_polyphen_gbg_24feb23[[val]]<-process_onerandomrep(ov_gbg_99)}

save(ran_polyphen_gbg_24feb23,file=paste("/scratch/devel/hpawar/admix/overlap.s*.skov/revisions_8feb23/polyphen/ran_polyphen_gbg_24feb23"))

# i think this shoudl work 

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
