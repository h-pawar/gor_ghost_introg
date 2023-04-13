#Mon 20 Mar 2023 09:25:35 CET
# generate random regions of equivalent length & callability
# assess linsight scores - mutational tolerance in noncoding regions 
    # follows from /Volumes/"Ultra USB 3.0"/IBE/further.analysis.feb.2020/gorillas/abc/revisions_feb23/linsight.14mar23.R



#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

# random regions have amended following relevant sections of  /scratch/devel/hpawar/admix/sstar/scripts/simul_postabc/downstream/gene.density.s*.skov.overlap.21jul22.R
#-----------------------------------------------------------------------------------------------------------------------

## module load gcc/6.3.0 openssl/1.0.2q R/4.0.1
library(GenomicScores)
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

# read in linsight gscore object
load(file=paste("/scratch/devel/hpawar/admix/overlap.s*.skov/revisions_8feb23/linsight/linsight_hg19"),verbose=T)


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

lin_perind<-function(scen,id){
gr1sco <- gscores(linsight, scen[[id]])
# retain regions where there is a linsight score
    # some introg regions annotated with NA instead of a linsight score
gr1sco<- gr1sco[!is.na((mcols(gr1sco))),]
# subset this object by high (>0.8), low (<0.0448) scores, & mean val across all introg regions for whcih there is a score 

# output list of 
# 1 & 2) counts of bp with high/low linsight scores,  : raw counts
# 3) mean linsihgt score - across regions where a score is available, 
# 4) proportion of high linsight bp across total introg bp for this ind, : proportions
# 5) proportion of low linsight bp across total introg bp for this ind : proportions
out<-list( sum(width(gr1sco[ gr1sco$default  > 0.8 ])-1),
   sum(width(gr1sco[ gr1sco$default < 0.0448 ])-1),
   mean((mcols(gr1sco))[,1]),
 (sum(width(gr1sco[ gr1sco$default > 0.8 ])-1))/ (sum(width(scen[[id]])-1)),
 (sum(width(gr1sco[ gr1sco$default < 0.0448 ])-1))/ (sum(width(scen[[id]])-1))
   )

return(out)}



proc_lin<-function(scen){
hold<-list()
for (ind in (1:length(scen))) {
print(ind)
hold[[ind]]<-lin_perind(scen,ind)
}
return(hold)}

#-----------------------------------------------------------------------------------------------------------------------

# changing output format
#process_onerandomrep<-function(scen){
#random_rep1_MG_allids<-list()
#for (val in 1:length(scen)){
#random_rep1_MG_allids[[val]]<-random_oneid(scen,val)
#}
#tes<-proc_lin(random_rep1_MG_allids)
#return(tes)
#}

#-----------------------------------------------------------------------------------------------------------------------

# check this - if fine, send as job - yes is fine
process_onerandomrep<-function(scen){
random_rep1_MG_allids<-list()
for (val in 1:length(scen)){
random_rep1_MG_allids[[val]]<-random_oneid(scen,val)
}
tes<-proc_lin(random_rep1_MG_allids)
tes1<-do.call(rbind,tes)
return(tes1)
}

#-----------------------------------------------------------------------------------------------------------------------

#process_onerandomrep(ov_gbb_99)
    # not sure why is hanging.. runs fine but takes time

# [,1] [,2]   [,3]       [,4] [,5]       
# [1,] 0    236000 0.06987195 0    0.008187053
# [2,] 0    82000  0.07171473 0    0.002368436
# [3,] 0    425000 0.07067952 0    0.01242145 

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------



ran_linsight_gbb_20mar23<-list()
for (i in 1:100){
val<-i
print(val)
ran_linsight_gbb_20mar23[[val]]<-process_onerandomrep(ov_gbb_99)}

# is taking a while to run - may be better to send as a job,..
save(ran_linsight_gbb_20mar23,file=paste("/scratch/devel/hpawar/admix/overlap.s*.skov/revisions_8feb23/linsight/ran_linsight_gbb_20mar23"))

#-----------------------------------------------------------------------------------------------------------------------

# for EL

ran_linsight_gbg_20mar23<-list()
for (i in 1:100){
val<-i
print(val)
ran_linsight_gbg_20mar23[[val]]<-process_onerandomrep(ov_gbg_99)}

save(ran_linsight_gbg_20mar23,file=paste("/scratch/devel/hpawar/admix/overlap.s*.skov/revisions_8feb23/linsight/ran_linsight_gbg_20mar23"))
#-----------------------------------------------------------------------------------------------------------------------
