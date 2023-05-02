# Mon 20 Mar 2023 10:15:33 CET
# focus on the sEs - & check if they are genic enhancers/promoter-interacting/enhancer-interacting 

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
#-----------------------------------------------------------------------------------------------------------------------

# read in annotated gorilla regulatory elements (generated in garcia-perez et al, & lifted over to hg19 by PEC)

gor_reg=read.table("/scratch/devel/pesteller/gor_reg_4Harvi/Gorilla.re_annotation.with_Non-RE.hg19.bed")
gor_reg<-as.data.frame(gor_reg)

#-----------------------------------------------------------------------------------------------------------------------

# sE regions
gor_sereg<-gor_reg[(which(gor_reg$V7=="sE")),]

# add column of genic/promoter interacting etc
gor_seregranges<-GRanges(seqnames=gor_sereg[,1],ranges=IRanges(start=as.numeric(gor_sereg[,2]),end=as.numeric(gor_sereg[,3]),names=gor_sereg[,1]),strand=rep("*",length(gor_sereg[,1])),re=(gor_sereg[,7]),pos=(gor_sereg[,8]))


#-----------------------------------------------------------------------------------------------------------------------
# only output the intersect regions
intersect_withgenes_function<-function(empirical_id) {
intersect_id<-list()
counts_id<-list()
for (ind in (1:length(empirical_id))) {
intersect_id[[ind]]<-intersect(gor_seregranges, empirical_id[[ind]], ignore.strand=TRUE)
}
return(intersect_id)
}

mg_se<-intersect_withgenes_function(ov_gbb_99)
mg_se_all<-reduce(unlist(GRangesList(mg_se)))

#-----------------------------------------------------------------------------------------------------------------------

# 1) looking at the population level

rtes<-findOverlaps(gor_seregranges,mg_se_all)

y1<-unique(as.data.frame(rtes)[,1])
test<-data.frame(gor_seregranges[y1])


table(test$pos)

nrow(test)

sum(table(test$pos))

sum(table(test$pos))/ nrow(test)

#-----------------------------------------------------------------------------------------------------------------------

# 2) position of strong enhancers per individual

proc_perid_function<-function(id) {
rtes<-findOverlaps(gor_seregranges,mg_se[[id]])
y1<-unique(as.data.frame(rtes)[,1])
test<-data.frame(gor_seregranges[y1])
out<-table(test$pos)
return(out)}


pos_perid<-list()
for (i in 1:length(mg_se)){
print(i)
pos_perid[[i]]<-proc_perid_function(i)}

x<-do.call(rbind,pos_perid)
colMeans(x)

