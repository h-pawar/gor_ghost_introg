#Mon 20 Mar 2023 09:25:35 CET
# generate random regions of equivalent length & callability
# assess linsight scores - mutational tolerance in noncoding regions 

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
gtf="/scratch/project/devel/mkuhlwilm/chimpdiv/index/Homo_sapiens.GRCh37.75_c.gtf"

test<-makeTxDbFromGFF(gtf, format='gtf')

human.genes = genes(test)

# only retain autosomal chr
autosome.hum.genes<-human.genes[seqnames(human.genes) == c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22")]

#-----------------------------------------------------------------------------------------------------------------------
# lengths of hg19
hg19<-fread('/home/devel/marcmont/scratch/Harvi/geneflow/test/test.reconst.ids/ora_1_1/proportion/hg19.autosomes.chrom.sizes.txt')

hg19.lengths<-apply(hg19[,2],2,function(x) as.numeric(as.character(x)))


hg19<-hg19[c(1:18,20,19,22,21),]
colnames(hg19)<-c("chrom","size") 

#-----------------------------------------------------------------------------------------------------------------------

# 1) read in the informative regions (ie with callable sites)

informative<-read.table("/scratch/devel/hpawar/admix/overlap.s*.skov/pairwisediff/callableprop.txt")


informranges<-GRanges(seqnames=informative[,1],ranges=IRanges(start=as.numeric(informative[,2]),end=as.numeric(informative[,2]+40000),names=informative[,1]),strand=rep("*",length(informative[,1])),prop=as.numeric(informative[,3]))
informative[,1]<-sapply("chr",paste, informative[,1], sep="")
informranges<-GRanges(seqnames=informative[,1],ranges=IRanges(start=as.numeric(informative[,2]),end=as.numeric(informative[,2]+40000),names=informative[,1]),strand=rep("*",length(informative[,1])),prop=as.numeric(informative[,3]))

replace_random_function<-function(scen,i){
test_random<-list()

repeat{

replace_random_1<-bed_random(hg19, length = (scen[i]@ranges@width-1), n=1)

x<-data.frame(matrix(unlist(replace_random_1), ncol = length(replace_random_1), byrow = TRUE))

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
    # some regions annotated with NA instead of a linsight score
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
#-----------------------------------------------------------------------------------------------------------------------

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


ran_linsight_gbb_20mar23<-list()
for (i in 1:100){
val<-i
print(val)
ran_linsight_gbb_20mar23[[val]]<-process_onerandomrep(ov_gbb_99)}

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
