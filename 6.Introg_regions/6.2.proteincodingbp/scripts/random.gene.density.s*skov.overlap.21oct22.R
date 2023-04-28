# Thu 20 Oct 2022 16:35:59 CEST
# generate random regions of sufficient callable sites & recalc proportion of protein coding bp in these random regions

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

## module load gcc/6.3.0 openssl/1.0.2q R/4.0.1
require(data.table)
options(stringsAsFactors=F)
options(scipen=100)
library(mgcv)
library(GenomicRanges)
library(tidyr)
library(ggbio)
library(valr)
library(GenomicFeatures)

#-----------------------------------------------------------------------------------------------------------------------
# empirical overlap b/n sstar-skov
load(file="/scratch/devel/hpawar/admix/overlap.s*.skov/overlap_objects/ov_gbb_99",verbose=T)
load(file="/scratch/devel/hpawar/admix/overlap.s*.skov/overlap_objects/ov_gbg_99",verbose=T)


#-----------------------------------------------------------------------------------------------------------------------

# (2) read in gtf of gene coordinates from hg19 reference (protein coding bp)
    # initially just use the human annotation (rather than filtering by orthologs)

# use makeTxDbFromGFF to create the TxDB object 
gtf="/scratch/project/devel/mkuhlwilm/chimpdiv/index/Homo_sapiens.GRCh37.75_c.gtf"

test<-makeTxDbFromGFF(gtf, format='gtf')

human.genes = genes(test)

# only retain autosomal chr
autosome.hum.genes<-human.genes[seqnames(human.genes) == c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22")]

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
# lengths of hg19
hg19<-fread('/home/devel/marcmont/scratch/Harvi/geneflow/test/test.reconst.ids/ora_1_1/proportion/hg19.autosomes.chrom.sizes.txt')

hg19.lengths<-apply(hg19[,2],2,function(x) as.numeric(as.character(x)))

# hg19 chr in the correct order
hg19<-hg19[c(1:18,20,19,22,21),]
colnames(hg19)<-c("chrom","size") 

#-----------------------------------------------------------------------------------------------------------------------

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

# ran these lines interactively
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


