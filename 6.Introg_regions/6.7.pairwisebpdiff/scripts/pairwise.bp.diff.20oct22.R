# Thu 20 Oct 2022 12:06:49 CEST
# REWRITE THE BELOW - done


#-----------------------------------------------------------------------------------------------------------------------
#Tue 18 Oct 2022 14:22:04 CEST

# rewrite calculation for pairwise bp diff - go from the individual level
#  E-W, E-E, W-W
#-----------------------------------------------------------------------------------------------------------------------

# Tue 18 Oct 2022 12:48:46 CEST
#Should I instead be calculating per individual, ie eastern individual 1 vs all eastern individuals vs all western individuals etc?
# MK
#yes, that is what I mean. 
#The problem is that you can't compare such short regions with random regions, 
#or of course you can but it's very difficult to make sense of it, 
#especially if the random regions are not selected for ~matching local mutation rates.
# If you calculate the sum for each individual, you reduce the complexity of the data drastically, 
# and with the larger number of regions local effects get smoothened.
#Then, for plotting you should use the same scale for X and Y axis, 
#to see deviations from the diagonal, and draw the diagonal as well. 
#The most important information is in E-W (exclusively), while E-E and W-W are there for comparison.

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
#ptype=1 # for troubleshooting purposes
ptype=as.numeric(as.character(ptype))

load(file="/scratch/devel/hpawar/admix/overlap.s*.skov/overlap_objects/ov_gbb_99",verbose=T)
load(file="/scratch/devel/hpawar/admix/overlap.s*.skov/overlap_objects/ov_gbg_99",verbose=T)


#-----------------------------------------------------------------------------------------------------------------------
# samples & their populations
# for ids
samples.in.vcf<-system(paste("bcftools query -l /scratch/devel/mkuhlwilm/arch/subsets/3_gorilla_22.vcf.gz",sep=""),intern=T) # -l, --list-samples: list sample names and exit

# add population identifier
identifiers<-cbind(samples.in.vcf, c(rep("GBB",12), rep("GBG",9), rep("GGD",1),rep("GGG",27)))
colnames(identifiers)<-c("ids","pop")
#-----------------------------------------------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
e<-seq(1:21)
w<-seq(22:49)

#-----------------------------------------------------------------------------------------------------------------------
# Tue 18 Oct 2022 14:53:23 CEST
# me
# - when I take the introgressed regions for eastern individual 1, I only calculate id 1 vs rest of the easterns (or westerns), 
# rather than also id 2 vs rest, id 3 vs rest (these should rather be calculated when reading in their specific introgressed regions)? 
# Is this correct?

# MK
#yes, this is correct.
#-----------------------------------------------------------------------------------------------------------------------

# per region
count_diff_introg_regions<-function(inp,chrom,ind){

testsnps<-system(paste("bcftools view -m2 -M2 -v snps -r ",inp," /scratch/devel/mkuhlwilm/arch/subsets/3_gorilla_",chrom,".vcf.gz | bcftools query -f '%CHROM %POS %REF %ALT [%GT ]\n' ",sep=""),intern=T) 
testsnps<-do.call(rbind,strsplit(testsnps,split=" "))  ## ie list merged vertically into df

# 2) randomly sample heterozygotes for the GTs per site - equal probability 0 or 1 - 0.5
aut_oneid<-testsnps

for(id in (5:ncol(aut_oneid))) {
aut_oneid[,id][which(aut_oneid[,id]=="0|1")]<-sample(c(0,1),size=length(which(aut_oneid[,id]=="0|1")),replace=T,prob=c(0.5,0.5))
}

aut_oneid[which(aut_oneid=="0|0")]<-0
aut_oneid[which(aut_oneid=="1|1")]<-1


# convert columns to numeric
aut_oneid<-as.data.frame(aut_oneid)
diff.e<-aut_oneid[,-c(1:4)]


# convert cols from character to numeric
diff.e[,]<-sapply(diff.e[,],as.numeric)

#head(diff.e)
#  V5 V6 V7 V8 V9 V10 V11 V12 V13 V14 V15 V16 V17 V18 V19 V20 V21 V22 V23 V24
#1  1  1  1  1  1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1
#2  1  1  1  1  1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1

#-----------------------------------------------------------------------------------------------------------------------

# as a function, to output raw vals for each of the eastern ids vs all western ids 
# E id 1 vs all westerns, E id 2 vs all westerns, etc

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
# calc all easterns against each other, then all westerns against each other
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
# convert to df
a<- data.frame(matrix(unlist(reps), ncol = max(lengths(reps)), byrow = TRUE))
# convert to matrix
a1<-as.matrix(a)
# extract upper triangle only
upper.tri(a1, diag = FALSE)
# set vals of lower triangle to nas
a1[lower.tri(a1,diag=TRUE)] <- NA # to remove the self comparisons & duplicate comparisons -> left only with unique comparisons
out_count<- unlist(as.list(a1)[!is.na(as.list(a1))]) # ie output this, rather than the sums
return(out_count)
}

# all possible eastern comparisons
tese<-count_diff(1,21)
# all possible western comparisons
tesw<-count_diff(22,49)
#-----------------------------------------------------------------------------------------------------------------------


# perhaps output the raw counts? yes the raw counts for each pairwise comparison
pout<-list( ewtes, tese, tesw)
return(pout)
}







#-----------------------------------------------------------------------------------------------------------------------

#i=1
#testreg=paste(as.character(ov_gbb_99[[1]][i]@seqnames),":",ov_gbb_99[[1]][i]@ranges@start,"-", as.data.frame(ov_gbb_99[[1]][i]@ranges)[,2],sep="")
#testreg_chr=as.numeric(strsplit( as.character(ov_gbb_99[[1]][i]@seqnames), split = "chr")[[1]][[2]])
#count_diff_introg_regions(testreg,testreg_chr,1)
# works as expected -> amend for rest & start running **

#-----------------------------------------------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
chr9_count_diff_introg_regions<-function(inp,chrom,ind){

testsnps<-system(paste("bcftools view -m2 -M2 -v snps -r ",inp," /scratch/devel/hpawar/admix/sstar/vcf_24jun22/3_gor.chr9.vcf.gz | bcftools query -f '%CHROM %POS %REF %ALT [%GT ]\n' ",sep=""),intern=T) 
testsnps<-do.call(rbind,strsplit(testsnps,split=" "))  ## ie list merged vertically into df

# 2) randomly sample heterozygotes for the GTs per site - equal probability 0 or 1 - 0.5
aut_oneid<-testsnps

for(id in (5:ncol(aut_oneid))) {
aut_oneid[,id][which(aut_oneid[,id]=="0|1")]<-sample(c(0,1),size=length(which(aut_oneid[,id]=="0|1")),replace=T,prob=c(0.5,0.5))
}

aut_oneid[which(aut_oneid=="0|0")]<-0
aut_oneid[which(aut_oneid=="1|1")]<-1


# convert columns to numeric
aut_oneid<-as.data.frame(aut_oneid)
diff.e<-aut_oneid[,-c(1:4)]


# convert cols from character to numeric
diff.e[,]<-sapply(diff.e[,],as.numeric)

#-----------------------------------------------------------------------------------------------------------------------

# as a function, to output raw vals for each of the eastern ids vs all western ids 
# E id 1 vs all westerns, E id 2 vs all westerns, etc

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
# calc all easterns against each other, then all westerns against each other
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
# convert to df
a<- data.frame(matrix(unlist(reps), ncol = max(lengths(reps)), byrow = TRUE))
# convert to matrix
a1<-as.matrix(a)
# extract upper triangle only
upper.tri(a1, diag = FALSE)
# set vals of lower triangle to nas
a1[lower.tri(a1,diag=TRUE)] <- NA # to remove the self comparisons & duplicate comparisons -> left only with unique comparisons
out_count<- unlist(as.list(a1)[!is.na(as.list(a1))]) # ie output this, rather than the sums
return(out_count)
}

# all possible eastern comparisons
tese<-count_diff(1,21)
# all possible western comparisons
tesw<-count_diff(22,49)
#-----------------------------------------------------------------------------------------------------------------------


# perhaps output the raw counts? yes the raw counts for each pairwise comparison
pout<-list( ewtes, tese, tesw)

return(pout)
}

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
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

# which(sapply(counts_id1, is.null))
#[1] 145 146 147

# query chrom 9 for these reps

# these are the identifiers (which rows)
# extract also which are the specific regions to query
  # & sep function for chr 9

chr9_regions<-which(sapply(counts_id1, is.null)) # b/c all null reps correspond to chr 9 regions 

for (i in 1:length(chr9_regions)) {
ptype<-as.numeric(chr9_regions[[i]])
  tryCatch({
    print(i)
# send in parallel for each of the 1494 regions
testreg<-paste(as.character(ov_gbb_99[[1]][ptype]@seqnames),":",ov_gbb_99[[1]][ptype]@ranges@start,"-", as.data.frame(ov_gbb_99[[1]][ptype]@ranges)[,2],sep="")
#testreg_chr<- as.numeric(strsplit( as.character(ov_gbb_99[[1]][ptype]@seqnames), split = "chr")[[1]][[2]])
testreg_chr<-9
counts_id1[[ptype]]<-chr9_count_diff_introg_regions(testreg,testreg_chr,1)
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

#which(sapply(counts_id1, is.null))
#integer(0)


#str(counts_id1)
#List of 232
# $ :List of 3
#  ..$ :List of 21
#  .. ..$ : int [1:28] 22 26 21 21 19 28 22 23 26 24 ...
#  .. ..$ : int [1:28] 11 21 16 22 20 21 9 20 13 13 ...
#  .. ..$ : int [1:28] 22 30 25 23 21 28 20 23 26 20 ...


# now take the sum across regions (or wait until have across all introg regions from all ids) 
  # could output like this, ie the raw counts?
    # more manageable to take the sums for the same pairwise comparisons across these 232 introg regions

# for individual 1
#counts_id1[[1]]
#[[1]][[20]]
# [1]  7 21 12 24 22 21  7 22  9 11 12  9 22 26 25 16 14 15 16 17 23 15 21 25 16
#[26] 15 11 25

#[[1]][[21]]
# [1]  8 20 13 23 21 20  8 21 10 10 13 10 23 25 26 15 15 16 15 16 22 16 22 24 15
#[26] 16 12 24


#[[2]]
#  [1] 23 18 15 23 12 17 22 25 20 17 21 12 11 14 21 23  0 15 12 25 12 21 24 15 16
# [26]  5 16 24 22 11 12 11 16 13 11 15 18 27 12 15 10 15 27  5 16 20 15 12 15 16
# [51] 11 15 11  8 12 18 13  8 13 20 13 13 15 14 14 10 18  9 22 21 30 21  9 31 18
# [76] 34 22 22 18  5 18 17 26 17  5 27 14 30 18 18  4 17  8 21 20 29 20  8 30 17
#[101] 33 21 21  3  3 21  4 17 16 25 16  4 26 13 29 17 17  7  5  6 22 25 20 17  0
#[126] 21 25  5 16 10 16 20 30 26 29 25 19  8 21 20 29 20  8 30 17 33 21 21  5  5
#[151]  2  4 29 20  3 16 15 26 15  3 25 12 28 16 16  6  4  5  5 26  5 17  8 21 20
#[176] 29 20  8 30 17 33 21 21  1  3  2  6 29  4  5 20  7 20 19 28 19  7 29 16 32
#[201] 20 20  2  2  5  7 28  7  6  3

#[[3]]
#  [1] 20 13 17 25 17 18 23 15 18  2 20 22 17 15 17  6 18 11 23 21 16 23 17 20 12
# [26] 10 25 21  8 18 13 29 27 20  8 27  8 16 13 19 17 16 10 19 10  9 25 18 28 26
# [51] 21  7 24 13 13  6 18 11 27 25 18  4 25  6  8  7 19 21 18 30 28 23 17 24 17
# [76] 17 22 15 27 19 16  6  8 15 25 16 29 21 30 27 30 26 24 19 15 17 20 24 19 28
#[101] 24 29 26 19 11 15 15 16 18 16 19 11 20 17 11 16 13 22 20 27 13 21 16 24 22
#[126] 23 13 18 13 17 16 13 22 26 25 22 14 22 15 25 23 20 10 19 16 14 15 12 15 27
#[151] 26 19 17 15 17 14 18 18 17 13 24 17 11 18 13 22 14 21 12 22 19 16 22 19 25
#[176] 23 24 12 17 18 16 17 12 15 27 26 19 13 12 19 24 14 19  5  3 18 22 11 28 20
#[201] 27 26 29  9 18 15 23 24 19 24 14 18 15 15 13 24 12 17 18 18 17 16 27 19 22
#[226] 15 15 20 19 20 12 20 20 19 13 11 22 20  9 24 18 23 22 19 15 18 19 19 20 23
#[251] 20 10 16 26 16 21  7  7 16 24 13 30 22 29 28 31  5 14 19 25 26 17 26  4 16
#[276] 12 13 21 18 24 22 19 11 20 15 11 16 11 14 26 25 16 14 11 16  7 23 21 17 25
#[301] 12 20 15 17 15 18 12 23 16 14 15 14 23 21 22 15 17 16 17 20 14 10 18 18 19
#[326]  8 18 13 27 25 20  8 25  8 10 13  6 19 27 28 15 13 14 15 18 26 14 22 28 15
#[351] 14 26 14 21  5  3 18 24 11 30 20 29 28 31 11 20 17 23 26 21 26  2 12 12  6
#[376] 25 16 26


# taking the sum across these regions may be a more manageable output?
#counts_id1[[1]][[1]][[1]]+ counts_id1[[2]][[1]][[1]]
# [1]  82 103  84  97  83  77 101  96 112  99  94 112  91 106  88 107  97  97  96
#[20]  90 104  97 107 109  87  82  94  93


#first_ew<-list()
#for(i in 1:length(counts_id1)) { 
#first_ew[[i]]<-counts_id1[[i]][[1]][[1]]}
#> ncol(do.call(rbind,first_ew))
#[1] 28
#> length(counts_id1[[1]][[1]][[1]])
#[1] 28

# ie output this *
# for each of the comparisons
#colSums(do.call(rbind,first_ew))
# [1] 25498 25706 25681 25833 25747 25306 25971 25642 25698 25781 25560 25886
#[13] 25735 25699 25780 25825 25485 26011 25920 25738 25741 25844 25670 25735
#[25] 25548 25518 25865 25806


# for region 1
#counts_id1[[1]][[1]] # is the e-w  of length 21 (21 list objects) - counts_id1[[i]][[1]][[1]] to counts_id1[[i]][[1]][[21]]
#counts_id1[[1]][[2]] # e-e length 1
#counts_id1[[1]][[3]] # e-w length 1

# process e-w
#process_ew<-function(ind){
#first_ew<-list()
#for(i in 1:length(counts_id1)) { 
#first_ew[[i]]<-counts_id1[[i]][[1]][[ind]]}
#return(colSums(do.call(rbind,first_ew)))
#}

#process_ew(1)
#> process_ew(1)
# [1] 25498 25706 25681 25833 25747 25306 25971 25642 25698 25781 25560 25886
#[13] 25735 25699 25780 25825 25485 26011 25920 25738 25741 25844 25670 25735
#[25] 25548 25518 25865 25806
## works fine

#all_ew<-list()
#for(i in 1:length(counts_id1[[1]][[1]])) { 
#all_ew[[i]]<-process_ew(i)}

# this as one of the list outputs *

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
# length(process_same(2))
#[1] 210

# length(counts_id1[[1]][[2]])
#[1] 210

#> length(counts_id1[[1]][[3]])
#[1] 378
#> length(process_same(3))
#[1] 378


proc<-list(all_ew,e_e_id1,w_w_id1)

gor_counts<-list()
save(gor_counts,file=paste("/scratch/devel/hpawar/admix/overlap.s*.skov/pairwisediff/gor_counts",sep=""))

load(file=paste("/scratch/devel/hpawar/admix/overlap.s*.skov/pairwisediff/gor_counts",sep=""))

ptype=1
gor_counts[[ptype]]<-proc
save(gor_counts,file=paste("/scratch/devel/hpawar/admix/overlap.s*.skov/pairwisediff/gor_counts",sep=""))


# continue with amending the below ***
  # will need to add these processing steps before sendign to proc ***


# AMEND BELOW 
#-----------------------------------------------------------------------------------------------------------------------

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
# send in parallel for each of the 1494 regions
testreg<-paste(as.character(scen[[id]][ptype]@seqnames),":",scen[[id]][ptype]@ranges@start,"-", as.data.frame(scen[[id]][ptype]@ranges)[,2],sep="")
#testreg_chr<- as.numeric(strsplit( as.character(scen[[1]][ptype]@seqnames), split = "chr")[[1]][[2]])
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


# function works for id 2
o<-process_oneid(ov_gbb_99,2)

gor_counts[[2]]<-o

# run for rest of gbb ids 
rem<-seq(3:12)+2

for (i in 1:length(rem)){
val<-rem[i]
gor_counts[[val]]<-process_oneid(ov_gbb_99,val)
}

save(gor_counts,file=paste("/scratch/devel/hpawar/admix/overlap.s*.skov/pairwisediff/gor_counts",sep=""))


# TESTING HERE **

# then for EL
for (i in 1:8){
val<-12+i
gor_counts[[val]]<-process_oneid(ov_gbg_99,i)
}

save(gor_counts,file=paste("/scratch/devel/hpawar/admix/overlap.s*.skov/pairwisediff/gor_counts",sep=""))


#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
# run tumani interactively * (b/c no overlap regions on chr 9 -> function breaks)
# for tumani run interactively
scen=ov_gbg_99
id=9

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
#integer(0)

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

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------


