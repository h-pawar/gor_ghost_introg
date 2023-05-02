# Fri 10 Feb 2023 16:43:23 CET
# generate random regions of equivalent length & callability
# & assess for regulatory regions


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

# gorilla annotated regulatory regions (generated in garcia-perez et al, & lifted over to hg19 by PEC)

gor_reg=read.table("/scratch/devel/pesteller/gor_reg_4Harvi/Gorilla.re_annotation.with_Non-RE.hg19.bed")

gor_reg<-as.data.frame(gor_reg)

gor_reg=read.table("/scratch/devel/pesteller/gor_reg_4Harvi/Gorilla.re_annotation.with_Non-RE.hg19.bed")

gor_reg<-as.data.frame(gor_reg)


# Column 7 is the consensus RE of the SPECIES 

#-----------------------------------------------------------------------------------------------------------------------

# retain only those with annotated regulatory function in both replicates *
gor_allreg<-gor_reg[-(which(gor_reg$V7=="P/Non-re" | gor_reg$V7=="E/Non-re"| gor_reg$V7=="Non-re")),]

gor_allregranges<-GRanges(seqnames=gor_allreg[,1],ranges=IRanges(start=as.numeric(gor_allreg[,2]),end=as.numeric(gor_allreg[,3]),names=gor_allreg[,1]),strand=rep("*",length(gor_allreg[,1])),re=(gor_allreg[,7]))



#-----------------------------------------------------------------------------------------------------------------------


gtf="/scratch/project/devel/mkuhlwilm/chimpdiv/index/Homo_sapiens.GRCh37.75_c.gtf"

test<-makeTxDbFromGFF(gtf, format='gtf')
human.genes = genes(test)

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

proc_ran_1id<-function(scen,id){
x<-random_oneid(scen,id)
rant<-intersect(gor_allregranges, x, ignore.strand=TRUE)
# output rancounts: how many bp in gorilla annotated regulatory regions 
rancounts<-cbind(sum(reduce(rant)@ranges@width), sum(x@ranges@width))

rtes<-findOverlaps(gor_allregranges,x)
y1<-unique(as.data.frame(rtes)[,1])
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


#tf<-out_onerandomrep_fun(ov_gbb_99)
#ran_reg_gbb_10oct23[1]]<-tf
           
           
# For MG # ran reps 1-3 interactively
# in rest of jobs load this R object & run *
load(file=paste("/scratch/devel/hpawar/admix/overlap.s*.skov/revisions_8feb23/ran_reg_gbb_10oct23"),verbose=T)



rem<-seq(4:100)+2

for (i in 1:length(rem)){
val<-rem[i]
print(val)
ran_reg_gbb_10oct23[[val]]<-out_onerandomrep_fun(ov_gbb_99)}

save(ran_reg_gbb_10oct23,file=paste("/scratch/devel/hpawar/admix/overlap.s*.skov/revisions_8feb23/ran_reg_gbb_10oct23"))

#-----------------------------------------------------------------------------------------------------------------------

# for EL

ran_reg_gbg_10oct23<-list()

save(ran_reg_gbg_10oct23,file=paste("/scratch/devel/hpawar/admix/overlap.s*.skov/revisions_8feb23/ran_reg_gbg_10oct23"))


for (i in 1:100){
val<-i
print(val)
ran_reg_gbg_10oct23[[val]]<-out_onerandomrep_fun(ov_gbg_99)}

save(ran_reg_gbg_10oct23,file=paste("/scratch/devel/hpawar/admix/overlap.s*.skov/revisions_8feb23/ran_reg_gbg_10oct23"))


#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------


q()

