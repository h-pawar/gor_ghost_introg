#Thu 20 Oct 2022 14:38:47 CEST

# random pairwise bp diff
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

# random regions: 1-12 MG, 13-21 EL
load(file="/scratch/devel/hpawar/admix/overlap.s*.skov/pairwisediff/random_18oct22",verbose=T)


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
out_count<- unlist(as.list(a1)[!is.na(as.list(a1))])=
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

# only once
random_counts<-list()
save(random_counts,file=paste("/scratch/devel/hpawar/admix/overlap.s*.skov/pairwisediff/random_counts",sep=""))

# in every job
load(file=paste("/scratch/devel/hpawar/admix/overlap.s*.skov/pairwisediff/random_counts",sep=""))
#random_counts[[ptype]]<-proc
#save(random_counts,file=paste("/scratch/devel/hpawar/admix/overlap.s*.skov/pairwisediff/random_counts",sep=""))




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


randchr9_regions<-which(sapply(randcounts_id1, is.null)) 

for (i in 1:length(randchr9_regions)) {
ptype<-as.numeric(randchr9_regions[[i]])
  tryCatch({
    print(i)

testreg<-paste(as.character(scen[[id]][ptype]@seqnames),":",scen[[id]][ptype]@ranges@start,"-", as.data.frame(scen[[id]][ptype]@ranges)[,2],sep="")

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


# run for rest of gbb ids 
rem<-seq(2:12)+1

for (i in 1:length(rem)){
val<-rem[i]
random_counts[[val]]<-process_oneid(random_18oct22,val)
}

save(random_counts,file=paste("/scratch/devel/hpawar/admix/overlap.s*.skov/pairwisediff/random_counts",sep=""))


#-----------------------------------------------------------------------------------------------------------------------

# then run for gbg ids 


for (i in 13:14){
val<-i
random_counts[[val]]<-process_oneid(random_18oct22,val)
}


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


randchr9_regions<-which(sapply(randcounts_id1, is.null)) 

randchr9_regions

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

process_all<-function(scen){

all_process_ew<-function(scen,ind){
all_ew<-list()

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

# for E-E
eastern<-as.data.frame(cbind(as.data.frame(emp_out[[2]]),as.data.frame(random_out[[2]])))
eastern$id<-"E-E"
colnames(eastern)<-c("empirical","random","comparison")

# for W-W
western<-as.data.frame(cbind(as.data.frame(emp_out[[3]]),as.data.frame(random_out[[3]])))
western$id<-"W-W"
colnames(western)<-c("empirical","random","comparison")

# for E-W

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
