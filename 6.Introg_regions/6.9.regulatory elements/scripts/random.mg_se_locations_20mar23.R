# Mon 20 Mar 2023 15:50:03 CET

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

gor_reg=read.table("/scratch/devel/pesteller/gor_reg_4Harvi/Gorilla.re_annotation.with_Non-RE.hg19.bed")
gor_reg<-as.data.frame(gor_reg)

# sE regions
gor_sereg<-gor_reg[(which(gor_reg$V7=="sE")),]

filt<- unique(gor_sereg$V1)[c(5,19:25,28,30,32,33)]

gor_sereg<-gor_sereg[-(which(gor_sereg$V1==filt[[1]] | gor_sereg$V1==filt[[2]] | gor_sereg$V1==filt[[3]]  | gor_sereg$V1==filt[[4]]   | gor_sereg$V1==filt[[5]]    | gor_sereg$V1==filt[[6]]     | gor_sereg$V1==filt[[7]]      | gor_sereg$V1==filt[[8]]       | gor_sereg$V1==filt[[9]]        | gor_sereg$V1==filt[[10]]         | gor_sereg$V1==filt[[11]]  | gor_sereg$V1==filt[[12]])),]

gor_seregranges<-GRanges(seqnames=gor_sereg[,1],ranges=IRanges(start=as.numeric(gor_sereg[,2]),end=as.numeric(gor_sereg[,3]),names=gor_sereg[,1]),strand=rep("*",length(gor_sereg[,1])),re=(gor_sereg[,7]),pos=(gor_sereg[,8]))



#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------


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


process_onerandomrep<-function(scen){
random_rep1_MG_allids<-list()
for (val in 1:length(scen)){
random_rep1_MG_allids[[val]]<-random_oneid(scen,val)
}

tes<-proc_p(random_rep1_MG_allids)
return(tes)
}

#-----------------------------------------------------------------------------------------------------------------------

# only need random regions for MG here

ran_setype_gbb_20mar23<-list()
for (i in 1:100){
val<-i
print(val)
ran_setype_gbb_20mar23[[val]]<-process_onerandomrep(ov_gbb_99)}

save(ran_setype_gbb_20mar23,file=paste("/scratch/devel/hpawar/admix/overlap.s*.skov/revisions_8feb23/ran_setype_gbb_20mar23"))

#-----------------------------------------------------------------------------------------------------------------------
