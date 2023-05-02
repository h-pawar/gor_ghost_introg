# Wed  8 Feb 2023 17:25:01 CET
# check the putative introgressed regions for gorilla regulatory elements which PEC & RG annotated in the LCLs project. 
    # published in https://www.nature.com/articles/s41467-021-23397-1  # relevant file is  Supplementary Data 1

#-----------------------------------------------------------------------------------------------------------------------
# PEC converted excel to bed file & performed liftover to hg19
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

# read in annotated gorilla regulatory elements (generated in garcia-perez et al, & lifted over to hg19 by PEC)
gor_reg=read.table("/scratch/devel/pesteller/gor_reg_4Harvi/Gorilla.re_annotation.with_Non-RE.hg19.bed")
gor_reg<-as.data.frame(gor_reg)
#Column 7 is the consensus RE of the SPECIES
#-----------------------------------------------------------------------------------------------------------------------
# retain only those with annotated regulatory function in both replicates *
gor_allreg<-gor_reg[-(which(gor_reg$V7=="P/Non-re" | gor_reg$V7=="E/Non-re"| gor_reg$V7=="Non-re")),]

gor_allregranges<-GRanges(seqnames=gor_allreg[,1],ranges=IRanges(start=as.numeric(gor_allreg[,2]),end=as.numeric(gor_allreg[,3]),names=gor_allreg[,1]),strand=rep("*",length(gor_allreg[,1])),re=(gor_allreg[,7]))

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

intersect_withgenes_function<-function(empirical_id) {
intersect_id<-list()
counts_id<-list()
for (ind in (1:length(empirical_id))) {
intersect_id[[ind]]<-intersect(gor_allregranges, empirical_id[[ind]], ignore.strand=TRUE)
# output - total number of introgressed bp in a regulatory element;  total introgressed bp in this individual
counts_id[[ind]]<-cbind(sum(reduce(intersect_id[[ind]])@ranges@width), sum(empirical_id[[ind]]@ranges@width))
}
return(list(intersect_id,do.call(rbind,counts_id)))
}



# for MG
test<-intersect_withgenes_function(ov_gbb_99)

# pop level estimate for MG: introg bp in regulatory element - normalised by total introg bp in this individual
mean(test[[2]][,1]/test[[2]][,2])
sd(test[[2]][,1]/test[[2]][,2])

#-----------------------------------------------------------------------------------------------------------------------


eltest<-intersect_withgenes_function(ov_gbg_99)


mean(eltest[[2]][,1]/eltest[[2]][,2])
sd(eltest[[2]][,1]/eltest[[2]][,2])

#-----------------------------------------------------------------------------------------------------------------------
mg_r<-cbind(mean(test[[2]][,1]/test[[2]][,2]),sd(test[[2]][,1]/test[[2]][,2]))
el_r<-cbind(mean(eltest[[2]][,1]/eltest[[2]][,2]),sd(eltest[[2]][,1]/eltest[[2]][,2]))


#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
# introg bp per reg element
# need to annotate overlap with which reg elements identified

# input = test (MG)
# input = eltest (EL) 
# id -> for each id of each scenario

bp_perreg_fun<-function(input,id){

rtes<-findOverlaps(gor_allregranges,input[[1]][[id]])

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

p[,3]<-p[,1]/input[[2]][id,2]

return(p[,c(2:3)])
}


# run for all individuals per scenario
pop_perreg_fun<-function(input){
o_per<-list()
for (ind in (1:nrow(input[[2]]))) {
o_per[[ind]]<-bp_perreg_fun(input,ind)}
return(o_per)
}


p_mg<-pop_perreg_fun(test)
p_el<-pop_perreg_fun(eltest)

# then take pop-wise means of introg bp per reg element 

process_fun<-function(input){

p1_o<- do.call(cbind, input)[,c(1,seq(2,length(input)*2,2))]
p2_o<-as.data.frame(lapply((p1_o)[,-1],as.numeric))

q_per<-list()
for (ind in (1:nrow(p2_o))) {
q_per[[ind]]<- mean(as.numeric(p2_o[ind,]))}

return(as.data.frame(cbind(p1_o[,1],do.call(rbind,q_per))))
}

process_fun(p_mg)
process_fun(p_el)


#-----------------------------------------------------------------------------------------------------------------------

mg_p <- process_fun(p_mg)
el_p <- process_fun(p_el)

emp_reg_mg<-list(mg_r,mg_p)
emp_reg_el<-list(el_r,el_p)

save(emp_reg_mg,file=paste("/scratch/devel/hpawar/admix/overlap.s*.skov/revisions_8feb23/emp_reg_gbb_8feb23"))
save(emp_reg_el,file=paste("/scratch/devel/hpawar/admix/overlap.s*.skov/revisions_8feb23/emp_reg_gbg_8feb23"))

#-----------------------------------------------------------------------------------------------------------------------


q()

