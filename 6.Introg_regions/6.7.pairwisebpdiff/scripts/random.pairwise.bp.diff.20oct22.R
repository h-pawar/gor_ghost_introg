#Thu 20 Oct 2022 14:38:47 CEST

# AMEND THE BELOW **


#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

#Tue 18 Oct 2022 17:35:40 CEST
# for random regions generated with - /Volumes/"Ultra USB 3.0"/IBE/further.analysis.feb.2020/gorillas/abc/simul_postabc/generate.randomreg.sufficallable.18oct22.R

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

# random regions: 1-12 MG, 13-21 EL
load(file="/scratch/devel/hpawar/admix/overlap.s*.skov/pairwisediff/random_18oct22",verbose=T)


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

#head(diff.e)
#  V5 V6 V7 V8 V9 V10 V11 V12 V13 V14 V15 V16 V17 V18 V19 V20 V21 V22 V23 V24
#1  1  1  1  1  1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1
#2  1  1  1  1  1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1


#ie eastern individual 1 vs all eastern individuals vs all western individuals etc - yes
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

# only once
random_counts<-list()
save(random_counts,file=paste("/scratch/devel/hpawar/admix/overlap.s*.skov/pairwisediff/random_counts",sep=""))

#random_counts[[1]]<-proc


# in every job
load(file=paste("/scratch/devel/hpawar/admix/overlap.s*.skov/pairwisediff/random_counts",sep=""))
#random_counts[[ptype]]<-proc
#save(random_counts,file=paste("/scratch/devel/hpawar/admix/overlap.s*.skov/pairwisediff/random_counts",sep=""))



#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------------------------------

process_oneid<-function(scen,id){

# 1) calculate for individual 1
randcounts_id1<-list()
for(i in 1:length(scen[[id]])) { 
  tryCatch({
    print(i)
testreg<-paste(as.character(scen[[id]][i]@seqnames),":",scen[[id]][i]@ranges@start,"-", as.data.frame(scen[[id]][i]@ranges)[,2],sep="")
testreg_chr<- as.numeric(strsplit( as.character(scen[[id]][i]@seqnames), split = "chr")[[1]][[2]])
randcounts_id1[[i]]<-count_diff_introg_regions(testreg,testreg_chr,id)
}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}


# query chrom 9 for these reps
randchr9_regions<-which(sapply(randcounts_id1, is.null)) # b/c all null reps correspond to chr 9 regions 

for (i in 1:length(randchr9_regions)) {
ptype<-as.numeric(randchr9_regions[[i]])
  tryCatch({
    print(i)
# send in parallel for each of the 1494 regions
testreg<-paste(as.character(scen[[id]][ptype]@seqnames),":",scen[[id]][ptype]@ranges@start,"-", as.data.frame(scen[[id]][ptype]@ranges)[,2],sep="")
#testreg_chr<- as.numeric(strsplit( as.character(scen[[1]][ptype]@seqnames), split = "chr")[[1]][[2]])
testreg_chr<-9
randcounts_id1[[ptype]]<-chr9_count_diff_introg_regions(testreg,testreg_chr,1)
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

#-----------------------------------------------------------------------------------------------------------------------
# process e-w
process_ew<-function(ind){
first_ew<-list()
for(i in 1:length(randcounts_id1)) { 
first_ew[[i]]<-randcounts_id1[[i]][[1]][[ind]]}
return(colSums(do.call(rbind,first_ew)))
}

all_ew<-list()
for(i in 1:length(randcounts_id1[[1]][[1]])) { 
all_ew[[i]]<-process_ew(i)}

#-----------------------------------------------------------------------------------------------------------------------
# process e-e, w-w

process_same<-function(ind){
same<-list()
for(i in 1:length(randcounts_id1)) { 
same[[i]]<-randcounts_id1[[i]][[ind]]}
return(colSums(do.call(rbind,same)))
}

e_e_id1<-process_same(2)
w_w_id1<-process_same(3)

proc<-list(all_ew,e_e_id1,w_w_id1)
return(proc)
}
#-----------------------------------------------------------------------------------------------------------------------

# only once
#random_counts<-list()
#save(random_counts,file=paste("/scratch/devel/hpawar/admix/overlap.s*.skov/pairwisediff/random_counts",sep=""))


o<-process_oneid(random_18oct22,1)
random_counts[[1]]<-o

save(random_counts,file=paste("/scratch/devel/hpawar/admix/overlap.s*.skov/pairwisediff/random_counts",sep=""))


# AMEND BELOW **

# run for rest of gbb ids 
rem<-seq(2:12)+1

for (i in 1:length(rem)){
val<-rem[i]
random_counts[[val]]<-process_oneid(random_18oct22,val)
}

save(random_counts,file=paste("/scratch/devel/hpawar/admix/overlap.s*.skov/pairwisediff/random_counts",sep=""))
# ie 1:12 : sums for gbb ids : col 1 = id vs easterns, col 2 = id vs westerns 

# check random_counts look fine
# check str of random_18oct22

#-----------------------------------------------------------------------------------------------------------------------

# then run for gbg ids 

#for (i in 13:21){
#val<-i
#random_counts[[val]]<-process_oneid(random_18oct22,val)
#}

# error here
# need to break the function for 14 -> b/c this random rep did not have any regions on chr 9 

for (i in 13:14){
val<-i
random_counts[[val]]<-process_oneid(random_18oct22,val)
}


#check_chr9<-list()
#for (i in 15:21){
#check_chr9[[i]]<-which(unique(random_18oct22[[i]]@seqnames)=="chr9")}
# ie function should work for the remaining reps (all have regions on chr 9)


for (i in 15:21){
val<-i
random_counts[[val]]<-process_oneid(random_18oct22,val)
}
save(random_counts,file=paste("/scratch/devel/hpawar/admix/overlap.s*.skov/pairwisediff/random_counts",sep=""))


#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

scen=random_18oct22
id=14

# 1) calculate for individual 1
randcounts_id1<-list()
for(i in 1:length(scen[[id]])) { 
  tryCatch({
    print(i)
testreg<-paste(as.character(scen[[id]][i]@seqnames),":",scen[[id]][i]@ranges@start,"-", as.data.frame(scen[[id]][i]@ranges)[,2],sep="")
testreg_chr<- as.numeric(strsplit( as.character(scen[[id]][i]@seqnames), split = "chr")[[1]][[2]])
randcounts_id1[[i]]<-count_diff_introg_regions(testreg,testreg_chr,id)
}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}


# query chrom 9 for these reps
randchr9_regions<-which(sapply(randcounts_id1, is.null)) # b/c all null reps correspond to chr 9 regions 

randchr9_regions
#integer(0)
# for id =14 no regions on chr 9 -> breaking hte function when running

#-----------------------------------------------------------------------------------------------------------------------
# process e-w
process_ew<-function(ind){
first_ew<-list()
for(i in 1:length(randcounts_id1)) { 
first_ew[[i]]<-randcounts_id1[[i]][[1]][[ind]]}
return(colSums(do.call(rbind,first_ew)))
}

all_ew<-list()
for(i in 1:length(randcounts_id1[[1]][[1]])) { 
all_ew[[i]]<-process_ew(i)}

#-----------------------------------------------------------------------------------------------------------------------
# process e-e, w-w

process_same<-function(ind){
same<-list()
for(i in 1:length(randcounts_id1)) { 
same[[i]]<-randcounts_id1[[i]][[ind]]}
return(colSums(do.call(rbind,same)))
}

e_e_id1<-process_same(2)
w_w_id1<-process_same(3)

proc<-list(all_ew,e_e_id1,w_w_id1)


random_counts[[14]]<-proc
save(random_counts,file=paste("/scratch/devel/hpawar/admix/overlap.s*.skov/pairwisediff/random_counts",sep=""))

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

#Thu 20 Oct 2022 15:40:22 CEST

# plot empirical vs random pairwise bp counts

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
library('ggplot2')

#-----------------------------------------------------------------------------------------------------------------------

load(file=paste("/scratch/devel/hpawar/admix/overlap.s*.skov/pairwisediff/random_counts",sep=""),verbose=T)
load(file=paste("/scratch/devel/hpawar/admix/overlap.s*.skov/pairwisediff/gor_counts"),verbose=T)


#-----------------------------------------------------------------------------------------------------------------------

#> str(random_counts)
#List of 21
# $ :List of 3
#  ..$ :List of 21
#  .. ..$ : num [1:28] 17221 17132 17646 17417 17553 ...
#  .. ..$ : num [1:28] 16838 16843 17457 17246 17228 ...

# for the 21 easterns: 21 sets of introgressed regions

# for e_w
# random_counts[[1]][[1]]
#[[1]]
# [1] 17221 17132 17646 17417 17553 17097 17681 17648 17575 17268 17401 17679
#[13] 17259 17605 17325 17673 17249 17360 17413 17226 17447 17346 17377 17475
#[25] 17251 17157 17030 17492

# random_counts[[1]][[1]][[1]] to random_counts[[1]][[1]][[21]]

# e-e
#random_counts[[1]][[2]]

# w-w
#random_counts[[1]][[2]]

process_all<-function(scen){

all_process_ew<-function(scen,ind){
all_ew<-list()
# ie length should be number of windows **
for(i in 1:length(scen)) { 
all_ew[[i]]<-scen[[i]][[1]][[ind]]}
return(colSums(do.call(rbind,all_ew)))
}

all_ew<-list()
for(i in 1:length(scen)) { 
all_ew[[i]]<-all_process_ew(scen,i)}


#-----------------------------------------------------------------------------------------------------------------------
# process e-e, w-w

process_same<-function(scen,ind){
same<-list()
for(i in 1:length(scen)) { 
same[[i]]<-scen[[i]][[ind]]}
return(colSums(do.call(rbind,same)))
}

# ie this is fine for e-e, w-w
alle_e<-process_same(scen,2)
allw_w<-process_same(scen,3)


return(list(all_ew,alle_e,allw_w))
}
#-----------------------------------------------------------------------------------------------------------------------

random_out<-process_all(random_counts)

load(file=paste("/scratch/devel/hpawar/admix/overlap.s*.skov/pairwisediff/gor_counts"),verbose=T)
emp_out<-process_all(gor_counts)


# lists of 3 : first list is per ind for the 21 ids (E id1 vs all westerns, E id2 vs all westerns etc)
# lists 2,3: all eastern pairwise comparisons, all western pairwise comparisons 

# now convert to dfs & add identifiers of whcih comparison **

#-----------------------------------------------------------------------------------------------------------------------

# head(as.data.frame(random_out[[2]]))
#  random_out[[2]]
#1          199306
#2          197528
#3          136626



# for E-E
eastern<-as.data.frame(cbind(as.data.frame(emp_out[[2]]),as.data.frame(random_out[[2]])))
eastern$id<-"E-E"
colnames(eastern)<-c("empirical","random","comparison")

# for W-W
western<-as.data.frame(cbind(as.data.frame(emp_out[[3]]),as.data.frame(random_out[[3]])))
western$id<-"W-W"
colnames(western)<-c("empirical","random","comparison")

#> head(western)
#  empirical random comparison
#1    332965 323791        W-W
#2    329353 323275        W-W

#-----------------------------------------------------------------------------------------------------------------------

# for E-W

#> ncol(as.data.frame(random_out[[1]]))
#[1] 21
# think this shoudl rather be sequential, since sep comparisons

#as.data.frame(random_out[[1]][[1]])

proc_ew<-list()
for(i in 1:length(random_out[[1]])) { 
proc_ew[[i]]<-as.data.frame(random_out[[1]][[i]])}


e_w_df<-do.call(rbind,proc_ew)



empproc_ew<-list()
for(i in 1:length(emp_out[[1]])) { 
empproc_ew[[i]]<-as.data.frame(emp_out[[1]][[i]])}

emp_e_w_df<-do.call(rbind,empproc_ew)

comp_ew<-as.data.frame(cbind(emp_e_w_df,e_w_df))
comp_ew$id<-"E-W"
colnames(comp_ew)<-c("empirical","random","comparison")


pairwise<-rbind(eastern,western,comp_ew)

summary(pairwise)
#107640 
#507731

pdf(paste("/scratch/devel/hpawar/admix/overlap.s*.skov/pairwisediff/plot.pairwisediff.20oct22.pdf",sep="")) 
ggplot(pairwise, aes(x=empirical, y=random, color=comparison)) +
  geom_point() + 
  geom_abline(intercept = 0 , slope = 1) +
  theme_classic() + xlim(107640,507731) + ylim(107640,507731)

dev.off()
#-----------------------------------------------------------------------------------------------------------------------


scp -r hpawar@172.16.10.20:/scratch/devel/hpawar/admix/overlap.s\*.skov/pairwisediff/plot.pairwisediff.20oct22.pdf  /Users/harvi/Downloads/gorilla_postabc/overlap.sstar.skov/plots



#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

# Sun 23 Oct 2022 10:01:08 CEST
# replot (changing cols, increasing size of axes) 


# plot empirical vs random pairwise bp counts

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
library('ggplot2')

#-----------------------------------------------------------------------------------------------------------------------

load(file=paste("/scratch/devel/hpawar/admix/overlap.s*.skov/pairwisediff/random_counts",sep=""),verbose=T)
load(file=paste("/scratch/devel/hpawar/admix/overlap.s*.skov/pairwisediff/gor_counts"),verbose=T)


#-----------------------------------------------------------------------------------------------------------------------

#> str(random_counts)
#List of 21
# $ :List of 3
#  ..$ :List of 21
#  .. ..$ : num [1:28] 17221 17132 17646 17417 17553 ...
#  .. ..$ : num [1:28] 16838 16843 17457 17246 17228 ...

# for the 21 easterns: 21 sets of introgressed regions

# for e_w
# random_counts[[1]][[1]]
#[[1]]
# [1] 17221 17132 17646 17417 17553 17097 17681 17648 17575 17268 17401 17679
#[13] 17259 17605 17325 17673 17249 17360 17413 17226 17447 17346 17377 17475
#[25] 17251 17157 17030 17492

# random_counts[[1]][[1]][[1]] to random_counts[[1]][[1]][[21]]

# e-e
#random_counts[[1]][[2]]

# w-w
#random_counts[[1]][[2]]

process_all<-function(scen){

all_process_ew<-function(scen,ind){
all_ew<-list()
# ie length should be number of windows **
for(i in 1:length(scen)) { 
all_ew[[i]]<-scen[[i]][[1]][[ind]]}
return(colSums(do.call(rbind,all_ew)))
}

all_ew<-list()
for(i in 1:length(scen)) { 
all_ew[[i]]<-all_process_ew(scen,i)}


#-----------------------------------------------------------------------------------------------------------------------
# process e-e, w-w

process_same<-function(scen,ind){
same<-list()
for(i in 1:length(scen)) { 
same[[i]]<-scen[[i]][[ind]]}
return(colSums(do.call(rbind,same)))
}

# ie this is fine for e-e, w-w
alle_e<-process_same(scen,2)
allw_w<-process_same(scen,3)


return(list(all_ew,alle_e,allw_w))
}
#-----------------------------------------------------------------------------------------------------------------------

random_out<-process_all(random_counts)

load(file=paste("/scratch/devel/hpawar/admix/overlap.s*.skov/pairwisediff/gor_counts"),verbose=T)
emp_out<-process_all(gor_counts)


# lists of 3 : first list is per ind for the 21 ids (E id1 vs all westerns, E id2 vs all westerns etc)
# lists 2,3: all eastern pairwise comparisons, all western pairwise comparisons 

# now convert to dfs & add identifiers of whcih comparison **

#-----------------------------------------------------------------------------------------------------------------------

# head(as.data.frame(random_out[[2]]))
#  random_out[[2]]
#1          199306
#2          197528
#3          136626



# for E-E
eastern<-as.data.frame(cbind(as.data.frame(emp_out[[2]]),as.data.frame(random_out[[2]])))
eastern$id<-"E-E"
colnames(eastern)<-c("empirical","random","comparison")

# for W-W
western<-as.data.frame(cbind(as.data.frame(emp_out[[3]]),as.data.frame(random_out[[3]])))
western$id<-"W-W"
colnames(western)<-c("empirical","random","comparison")

#> head(western)
#  empirical random comparison
#1    332965 323791        W-W
#2    329353 323275        W-W

#-----------------------------------------------------------------------------------------------------------------------

# for E-W

#> ncol(as.data.frame(random_out[[1]]))
#[1] 21
# think this shoudl rather be sequential, since sep comparisons

#as.data.frame(random_out[[1]][[1]])

proc_ew<-list()
for(i in 1:length(random_out[[1]])) { 
proc_ew[[i]]<-as.data.frame(random_out[[1]][[i]])}


e_w_df<-do.call(rbind,proc_ew)



empproc_ew<-list()
for(i in 1:length(emp_out[[1]])) { 
empproc_ew[[i]]<-as.data.frame(emp_out[[1]][[i]])}

emp_e_w_df<-do.call(rbind,empproc_ew)

comp_ew<-as.data.frame(cbind(emp_e_w_df,e_w_df))
comp_ew$id<-"E-W"
colnames(comp_ew)<-c("empirical","random","comparison")


pairwise<-rbind(eastern,western,comp_ew)

summary(pairwise)
#107640 
#507731

pdf(paste("/scratch/devel/hpawar/admix/overlap.s*.skov/pairwisediff/plot.pairwisediff.20oct22.v1.pdf",sep="")) 
ggplot(pairwise, aes(x=empirical, y=random, color=comparison)) +
  geom_point() + 
  scale_color_manual(values=c('#619E9A','#7b778c','#F09164')) +
  geom_abline(intercept = 0 , slope = 1) + 
  theme_classic() + xlim(107640,507731) + ylim(107640,507731)  + ylab("Random")  +  xlab("Empirical") + 
  labs(color='Comparison') +
  theme(axis.title = element_text(size = 15))  +  theme(axis.title = element_text(size = 15)) + 
  theme(legend.title = element_text(size=15), legend.text = element_text(size=10))


dev.off()
#-----------------------------------------------------------------------------------------------------------------------


scp -r hpawar@172.16.10.20:/scratch/devel/hpawar/admix/overlap.s\*.skov/pairwisediff/plot.pairwisediff.20oct22.v1.pdf  /Users/harvi/Downloads/gorilla_postabc/overlap.sstar.skov/plots
