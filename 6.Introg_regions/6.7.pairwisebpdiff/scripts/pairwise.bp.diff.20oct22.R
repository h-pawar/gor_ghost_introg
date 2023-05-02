# Thu 20 Oct 2022 12:06:49 CEST
# calculate pairwise bp diff per individual
#  E-W, E-E, W-W
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#module load gcc/6.3.0 openssl/1.0.2q  R/4.0.1  BCFTOOLS/1.14
require(data.table)
options(stringsAsFactors=F)
library(tidyverse)
options(scipen=100)
library(adegenet)
library(GenomicRanges)
library(tidyr)
library(phangorn) 
library(ape) 
library('pegas')
ptype=(commandArgs(TRUE))
ptype=as.numeric(as.character(ptype))

load(file="/scratch/devel/hpawar/admix/overlap.s*.skov/overlap_objects/ov_gbb_99",verbose=T)
load(file="/scratch/devel/hpawar/admix/overlap.s*.skov/overlap_objects/ov_gbg_99",verbose=T)


#-----------------------------------------------------------------------------------------------------------------------

samples.in.vcf<-system(paste("bcftools query -l /scratch/devel/mkuhlwilm/arch/subsets/3_gorilla_22.vcf.gz",sep=""),intern=T) # -l, --list-samples: list sample names and exit

identifiers<-cbind(samples.in.vcf, c(rep("GBB",12), rep("GBG",9), rep("GGD",1),rep("GGG",27)))
colnames(identifiers)<-c("ids","pop")
#-----------------------------------------------------------------------------------------------------------------------
e<-seq(1:21)
w<-seq(22:49)
#-----------------------------------------------------------------------------------------------------------------------

# per region
count_diff_introg_regions<-function(inp,chrom,ind){

testsnps<-system(paste("bcftools view -m2 -M2 -v snps -r ",inp," /scratch/devel/mkuhlwilm/arch/subsets/3_gorilla_",chrom,".vcf.gz | bcftools query -f '%CHROM %POS %REF %ALT [%GT ]\n' ",sep=""),intern=T) 
testsnps<-do.call(rbind,strsplit(testsnps,split=" ")) 
# 2) randomly sample heterozygotes for the GTs per site - equal probability 0 or 1 - 0.5
aut_oneid<-testsnps

for(id in (5:ncol(aut_oneid))) {
aut_oneid[,id][which(aut_oneid[,id]=="0|1")]<-sample(c(0,1),size=length(which(aut_oneid[,id]=="0|1")),replace=T,prob=c(0.5,0.5))
}

aut_oneid[which(aut_oneid=="0|0")]<-0
aut_oneid[which(aut_oneid=="1|1")]<-1

aut_oneid<-as.data.frame(aut_oneid)
diff.e<-aut_oneid[,-c(1:4)]

diff.e[,]<-sapply(diff.e[,],as.numeric)

#-----------------------------------------------------------------------------------------------------------------------

# across one introg region (all rows) count the number of sites where the GT of id x differs from id y
one_site_diff<-function(x,y){
out<-sum(diff.e[,x]!=diff.e[,y])
return(out)
}

first_eid_all_w_diff<-function(z){
wtes<-list()
for(i in 1:length(w)) { 
wtes[[i]]<-one_site_diff(z, (w+21)[i])}
return(unlist(wtes))
}

ewtes<-list()
for(i in 1:length(e)) { 
ewtes[[i]]<-first_eid_all_w_diff(i)}


#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

# matrix of non-duplicate & non-self pairwise comparisons
count_diff<-function(colstart,colend){
reps=list()
for(i in colstart:colend) {
hold=list()
for(j in colstart:colend) { 
hold[[j]]<-sum(diff.e[,i]!=diff.e[,j])
}
reps[[i]]<-do.call(rbind,hold)
}

a<- data.frame(matrix(unlist(reps), ncol = max(lengths(reps)), byrow = TRUE))

a1<-as.matrix(a)

upper.tri(a1, diag = FALSE)

a1[lower.tri(a1,diag=TRUE)] <- NA # to remove the self comparisons & duplicate comparisons -> left only with unique comparisons
out_count<- unlist(as.list(a1)[!is.na(as.list(a1))]) 
return(out_count)
}

# all possible eastern comparisons
tese<-count_diff(1,21)
# all possible western comparisons
tesw<-count_diff(22,49)
#-----------------------------------------------------------------------------------------------------------------------

pout<-list( ewtes, tese, tesw)
return(pout)
}


#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
chr9_count_diff_introg_regions<-function(inp,chrom,ind){

testsnps<-system(paste("bcftools view -m2 -M2 -v snps -r ",inp," /scratch/devel/hpawar/admix/sstar/vcf_24jun22/3_gor.chr9.vcf.gz | bcftools query -f '%CHROM %POS %REF %ALT [%GT ]\n' ",sep=""),intern=T) 
testsnps<-do.call(rbind,strsplit(testsnps,split=" "))  
# 2) randomly sample heterozygotes for the GTs per site - equal probability 0 or 1 - 0.5
aut_oneid<-testsnps

for(id in (5:ncol(aut_oneid))) {
aut_oneid[,id][which(aut_oneid[,id]=="0|1")]<-sample(c(0,1),size=length(which(aut_oneid[,id]=="0|1")),replace=T,prob=c(0.5,0.5))
}

aut_oneid[which(aut_oneid=="0|0")]<-0
aut_oneid[which(aut_oneid=="1|1")]<-1


aut_oneid<-as.data.frame(aut_oneid)
diff.e<-aut_oneid[,-c(1:4)]

diff.e[,]<-sapply(diff.e[,],as.numeric)

#-----------------------------------------------------------------------------------------------------------------------

# across one introg region (all rows) count the number of sites where the GT of id x differs from id y
one_site_diff<-function(x,y){
out<-sum(diff.e[,x]!=diff.e[,y])
return(out)
}

first_eid_all_w_diff<-function(z){
wtes<-list()
for(i in 1:length(w)) { 
wtes[[i]]<-one_site_diff(z, (w+21)[i])}
return(unlist(wtes))
}

ewtes<-list()
for(i in 1:length(e)) { 
ewtes[[i]]<-first_eid_all_w_diff(i)}


#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

# matrix of non-duplicate & non-self pairwise comparisons
count_diff<-function(colstart,colend){
reps=list()
for(i in colstart:colend) {
hold=list()
for(j in colstart:colend) { 
hold[[j]]<-sum(diff.e[,i]!=diff.e[,j])
}
reps[[i]]<-do.call(rbind,hold)
}

a<- data.frame(matrix(unlist(reps), ncol = max(lengths(reps)), byrow = TRUE))

a1<-as.matrix(a)

upper.tri(a1, diag = FALSE)

a1[lower.tri(a1,diag=TRUE)] <- NA # to remove the self comparisons & duplicate comparisons -> left only with unique comparisons
out_count<- unlist(as.list(a1)[!is.na(as.list(a1))]) 
return(out_count)
}

# all possible eastern comparisons
tese<-count_diff(1,21)
# all possible western comparisons
tesw<-count_diff(22,49)
#-----------------------------------------------------------------------------------------------------------------------

pout<-list( ewtes, tese, tesw)

return(pout)
}

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------


# 1) calculate for individual 1 ** (for multiple individuals go to the next section)

counts_id1<-list()
for(i in 1:length(ov_gbb_99[[1]])) { 
  tryCatch({
    print(i)
testreg<-paste(as.character(ov_gbb_99[[1]][i]@seqnames),":",ov_gbb_99[[1]][i]@ranges@start,"-", as.data.frame(ov_gbb_99[[1]][i]@ranges)[,2],sep="")
testreg_chr<- as.numeric(strsplit( as.character(ov_gbb_99[[1]][i]@seqnames), split = "chr")[[1]][[2]])
counts_id1[[i]]<-count_diff_introg_regions(testreg,testreg_chr,1)
}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}



which(sapply(counts_id1, is.null))

chr9_regions<-which(sapply(counts_id1, is.null)) # b/c all null reps correspond to chr 9 regions 

for (i in 1:length(chr9_regions)) {
ptype<-as.numeric(chr9_regions[[i]])
  tryCatch({
    print(i)

testreg<-paste(as.character(ov_gbb_99[[1]][ptype]@seqnames),":",ov_gbb_99[[1]][ptype]@ranges@start,"-", as.data.frame(ov_gbb_99[[1]][ptype]@ranges)[,2],sep="")

testreg_chr<-9
counts_id1[[ptype]]<-chr9_count_diff_introg_regions(testreg,testreg_chr,1)
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}



#-----------------------------------------------------------------------------------------------------------------------
# process e-w
process_ew<-function(ind){
first_ew<-list()
for(i in 1:length(counts_id1)) { 
first_ew[[i]]<-counts_id1[[i]][[1]][[ind]]}
return(colSums(do.call(rbind,first_ew)))
}

all_ew<-list()
for(i in 1:length(counts_id1[[1]][[1]])) { 
all_ew[[i]]<-process_ew(i)}

#-----------------------------------------------------------------------------------------------------------------------
# process e-e, w-w

process_same<-function(ind){
same<-list()
for(i in 1:length(counts_id1)) { 
same[[i]]<-counts_id1[[i]][[ind]]}
return(colSums(do.call(rbind,same)))
}

e_e_id1<-process_same(2)
w_w_id1<-process_same(3)


proc<-list(all_ew,e_e_id1,w_w_id1)

gor_counts<-list()
save(gor_counts,file=paste("/scratch/devel/hpawar/admix/overlap.s*.skov/pairwisediff/gor_counts",sep=""))

load(file=paste("/scratch/devel/hpawar/admix/overlap.s*.skov/pairwisediff/gor_counts",sep=""))

ptype=1
gor_counts[[ptype]]<-proc
save(gor_counts,file=paste("/scratch/devel/hpawar/admix/overlap.s*.skov/pairwisediff/gor_counts",sep=""))


#-----------------------------------------------------------------------------------------------------------------------

process_oneid<-function(scen,id){

# 1) calculate for individual 1
counts_id1<-list()
for(i in 1:length(scen[[id]])) { 
  tryCatch({
    print(i)
testreg<-paste(as.character(scen[[id]][i]@seqnames),":",scen[[id]][i]@ranges@start,"-", as.data.frame(scen[[id]][i]@ranges)[,2],sep="")
testreg_chr<- as.numeric(strsplit( as.character(scen[[id]][i]@seqnames), split = "chr")[[1]][[2]])
counts_id1[[i]]<-count_diff_introg_regions(testreg,testreg_chr,id)
}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}


# query chrom 9 for these reps
chr9_regions<-which(sapply(counts_id1, is.null)) # b/c all null reps correspond to chr 9 regions 

for (i in 1:length(chr9_regions)) {
ptype<-as.numeric(chr9_regions[[i]])
  tryCatch({
    print(i)

testreg<-paste(as.character(scen[[id]][ptype]@seqnames),":",scen[[id]][ptype]@ranges@start,"-", as.data.frame(scen[[id]][ptype]@ranges)[,2],sep="")

testreg_chr<-9
counts_id1[[ptype]]<-chr9_count_diff_introg_regions(testreg,testreg_chr,1)
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

#-----------------------------------------------------------------------------------------------------------------------
# process e-w
process_ew<-function(ind){
first_ew<-list()
for(i in 1:length(counts_id1)) { 
first_ew[[i]]<-counts_id1[[i]][[1]][[ind]]}
return(colSums(do.call(rbind,first_ew)))
}

all_ew<-list()
for(i in 1:length(counts_id1[[1]][[1]])) { 
all_ew[[i]]<-process_ew(i)}

#-----------------------------------------------------------------------------------------------------------------------
# process e-e, w-w

process_same<-function(ind){
same<-list()
for(i in 1:length(counts_id1)) { 
same[[i]]<-counts_id1[[i]][[ind]]}
return(colSums(do.call(rbind,same)))
}

e_e_id1<-process_same(2)
w_w_id1<-process_same(3)

proc<-list(all_ew,e_e_id1,w_w_id1)

return(proc)
}
#-----------------------------------------------------------------------------------------------------------------------



o<-process_oneid(ov_gbb_99,2)

gor_counts[[2]]<-o


rem<-seq(3:12)+2

for (i in 1:length(rem)){
val<-rem[i]
gor_counts[[val]]<-process_oneid(ov_gbb_99,val)
}

save(gor_counts,file=paste("/scratch/devel/hpawar/admix/overlap.s*.skov/pairwisediff/gor_counts",sep=""))



# then for EL
for (i in 1:8){
val<-12+i
gor_counts[[val]]<-process_oneid(ov_gbg_99,i)
}

save(gor_counts,file=paste("/scratch/devel/hpawar/admix/overlap.s*.skov/pairwisediff/gor_counts",sep=""))


#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
# for tumani run interactively
scen=ov_gbg_99
id=9

counts_id1<-list()
for(i in 1:length(scen[[id]])) { 
  tryCatch({
    print(i)
testreg<-paste(as.character(scen[[id]][i]@seqnames),":",scen[[id]][i]@ranges@start,"-", as.data.frame(scen[[id]][i]@ranges)[,2],sep="")
testreg_chr<- as.numeric(strsplit( as.character(scen[[id]][i]@seqnames), split = "chr")[[1]][[2]])
counts_id1[[i]]<-count_diff_introg_regions(testreg,testreg_chr,id)
}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}



chr9_regions<-which(sapply(counts_id1, is.null)) 

#-----------------------------------------------------------------------------------------------------------------------
# process e-w
process_ew<-function(ind){
first_ew<-list()
for(i in 1:length(counts_id1)) { 
first_ew[[i]]<-counts_id1[[i]][[1]][[ind]]}
return(colSums(do.call(rbind,first_ew)))
}

all_ew<-list()
for(i in 1:length(counts_id1[[1]][[1]])) { 
all_ew[[i]]<-process_ew(i)}

#-----------------------------------------------------------------------------------------------------------------------
# process e-e, w-w

process_same<-function(ind){
same<-list()
for(i in 1:length(counts_id1)) { 
same[[i]]<-counts_id1[[i]][[ind]]}
return(colSums(do.call(rbind,same)))
}

e_e_id1<-process_same(2)
w_w_id1<-process_same(3)

proc<-list(all_ew,e_e_id1,w_w_id1)

gor_counts[[21]]<-proc
save(gor_counts,file=paste("/scratch/devel/hpawar/admix/overlap.s*.skov/pairwisediff/gor_counts",sep=""))


