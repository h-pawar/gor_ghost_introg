# generate random regions for nj trees
#-----------------------------------------------------------------------------------------------------------------------
#module load gcc/6.3.0 openssl/1.0.2q  R/4.0.1  BCFTOOLS/1.14
require(data.table)
options(stringsAsFactors=F)
library(tidyverse)
options(scipen=100)
library(adegenet)
#library(mgcv)
library(GenomicRanges)
#BiocManager::install(c("GenomicRanges", "plyranges", "HelloRangesData"))
#library(plyranges) # may not be necessary here
#library(GenomicFeatures) # if need to read in the gtf
library(tidyr)
#library(ggbio)
library(ggplot2)
#-----------------------------------------------------------------------------------------------------------------------


load(file="/scratch/devel/hpawar/admix/overlap.s*.skov/overlap_objects/ov_gbb_99",verbose=T)
load(file="/scratch/devel/hpawar/admix/overlap.s*.skov/overlap_objects/ov_gbg_99",verbose=T)
load(file="/scratch/devel/hpawar/admix/overlap.s*.skov/overlap_objects/ov_gbb40_99",verbose=T)
load(file="/scratch/devel/hpawar/admix/overlap.s*.skov/overlap_objects/ov_gbg40_99",verbose=T)

#-----------------------------------------------------------------------------------------------------------------------
library(phangorn) 
library(ape) 
library('pegas')
library(valr) 
#-----------------------------------------------------------------------------------------------------------------------
# samples & their populations
# for ids
samples.in.vcf<-system(paste("bcftools query -l /scratch/devel/mkuhlwilm/arch/subsets/3_gorilla_22.vcf.gz",sep=""),intern=T) # -l, --list-samples: list sample names and exit

# add population identifier
identifiers<-cbind(samples.in.vcf, c(rep("GBB",12), rep("GBG",9), rep("GGD",1),rep("GGG",27)))
colnames(identifiers)<-c("ids","pop")
#-----------------------------------------------------------------------------------------------------------------------


# A) generate random regions of equal length distribution as the empirical data for each individual
# (B) intersect w vcfs -> matrix of gts -> pca

#-----------------------------------------------------------------------------------------------------------------------
# lengths of hg19
hg19<-fread('/home/devel/marcmont/scratch/Harvi/geneflow/test/test.reconst.ids/ora_1_1/proportion/hg19.autosomes.chrom.sizes.txt')

hg19.lengths<-apply(hg19[,2],2,function(x) as.numeric(as.character(x)))

# hg19 chr in the correct order
hg19<-hg19[c(1:18,20,19,22,21),]
colnames(hg19)<-c("chrom","size") 

#-----------------------------------------------------------------------------------------------------------------------

# gives lengths in bp
calc_winlengths<-function(fG){
win_lengths<-list()
for (ind in (1:length(fG))) {
win_lengths[[ind]]<-table(fG[[ind]]@ranges@width-1)}
return(win_lengths)
}
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

# generate random regions of a given length distribution (matching that of the fG input for each individual li)
# run this function per individual
test_proc_randomregions<-function(fG, li) {
# for one individual (li)
y<-as.data.frame(calc_winlengths(fG)[[li]])
# convert from factor -> character -> numeric
y$Var1<-as.numeric(as.character(y$Var1))

# write as a function - go through each of lengths & n for this individual
test_randomreg<-list()
for (ind in (1:nrow(y))) {
#test_randomreg[[ind]]<-bed_random(hg19, length = y[ind,1], n =  y[ind,2], seed = 10104) 
# run without specifying the seed
test_randomreg[[ind]]<-bed_random(hg19, length = y[ind,1], n =  y[ind,2]) 
}


convert_tib_df<-function(lin) {
x<- test_randomreg[[lin]]
x<-as.data.frame(x)
# convert to granges - then this will be what will be being intersected - later 
#test_s<-GRanges(seqnames=x[,1],ranges=IRanges(start=as.numeric(x[,2]),end=as.numeric(x[,3]),names=x[,1]),strand=rep("*",length(x[,1])))
#return(test_s)
return(x)
}

all_windowcategories<-list()
for (ind in (1:nrow(y))) {
all_windowcategories[[ind]]<-convert_tib_df(ind) }

#do.call("rbind", list(dataframes to merge))
x1<-do.call("rbind", all_windowcategories)
# convert to granges - then this will be what will be being intersected
test_s<-GRanges(seqnames=x1[,1],ranges=IRanges(start=as.numeric(x1[,2]),end=as.numeric(x1[,3]),names=x1[,1]),strand=rep("*",length(x1[,1])))

return(test_s)
}

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

generate_random<-function(fG) {
random_reg_x<-list()
for (ind in (1:length(fG))) {
random_reg_x[[ind]]<-test_proc_randomregions(fG,ind) }
return(random_reg_x)
}

#-----------------------------------------------------------------------------------------------------------------------


sk_ran_gbb<-generate_random(ov_gbb_99) # generates for 1 rep (random regions for each id)
sk_ran_gbg<-generate_random(ov_gbg_99) # generates for 1 rep (random regions for each id)
sk_ran_gbb_40<-generate_random(ov_gbb40_99) # generates for 1 rep (random regions for each id)
sk_ran_gbg_40<-generate_random(ov_gbg40_99) # generates for 1 rep (random regions for each id)
#-----------------------------------------------------------------------------------------------------------------------

# write these out to sep dir *
#/scratch/devel/hpawar/admix/overlap.s*.skov/trees/random

#[hpawar@cnb1 ~]$ mkdir -p /scratch/devel/hpawar/admix/overlap.s*.skov/trees/random/GBB/gbb
#[hpawar@cnb1 ~]$ mkdir -p /scratch/devel/hpawar/admix/overlap.s*.skov/trees/random/GBB/gbb.40
#[hpawar@cnb1 ~]$ mkdir -p /scratch/devel/hpawar/admix/overlap.s*.skov/trees/random/GBG/gbg
#[hpawar@cnb1 ~]$ mkdir -p /scratch/devel/hpawar/admix/overlap.s*.skov/trees/random/GBG/gbg.40

# functions with paths outputting to the random data directories
# granges of introg regions for 1 individual -> bed file 
  # contains introg regions across all autosomes for the specified individual
random_granges_tobed<-function(scen,ind,SPE,spe){

df <- data.frame(seqnames=seqnames(scen[[ind]]),
  start=start(scen[[ind]])-1,
  end=end(scen[[ind]]))

# convert to same structure of sk
df[,1]<-as.character(df[,1])
df[,2]<-as.integer(df[,2])

a=ind
tmp<-paste("/scratch/devel/hpawar/admix/overlap.s*.skov/trees/random/",SPE,"/",spe,".",a,".random.tmp.bed",sep="")

write.table(df,tmp,sep="\t",row.names=F,col.names=F,quote=F) # should write this out only once

}



#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

# run this before the function below?

random_allids<-function(scen,SPE,spe) {
# run for all individuals - may be better to split into separate jobs??
for (ind in (1:length(scen))) {
random_granges_tobed(scen,ind,SPE,spe)}
}

random_allids(sk_ran_gbb,"GBB","gbb")
random_allids(sk_ran_gbg,"GBG","gbg")
random_allids(sk_ran_gbb_40,"GBB","gbb.40")
random_allids(sk_ran_gbg_40,"GBG","gbg.40")

# need to regenerate for GBB - done, but need to check if they look ok
# for GBG random regions look ok

# random reps have generated - Fri 29 Jul 2022 08:48:10 CEST
# ie this section has been executed *  - so can send the below to a sep file & send as a job on the cluster


q()
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------


# 2/6/22 - random_granges_tobed - function & output looks fine - 

#>  head(read.table(tmp))
#     V1        V2        V3
#1  chr1  67179803  67180804
#2 chr20  14127100  14128101
#3  chr2 235120469 235122470

#-----------------------------------------------------------------------------------------------------------------------

#bed file -> intersect with vcf, extract biallelic snps only -> output chr, pos, gts

random_intersect_introgbed_regions<-function(ind,SPE,spe,chrom){
a=ind
tmp<-paste("/scratch/devel/hpawar/admix/overlap.s*.skov/trees/random/",SPE,"/",spe,".",a,".random.tmp.bed",sep="")

chrom=chrom
# filter by biallelic snps only, query by the putative introg regions for id 1 & output chr, pos & gts
testsnps<-system(paste("bcftools view -m2 -M2 -v snps -R ",tmp," /scratch/devel/mkuhlwilm/arch/subsets/3_gorilla_",chrom,".vcf.gz | bcftools query -f '%CHROM %POS [%GT ]\n' ",sep=""),intern=T) 
testsnps<-do.call(rbind,strsplit(testsnps,split=" "))  ## ie list merged vertically into df

return(testsnps)
}

#-----------------------------------------------------------------------------------------------------------------------

random_test_gts_topca<-function(ind,SPE,spe) {

#random_granges_tobed(scen,ind,SPE,spe)

# extract biallelic GTs for introg regions across all autosomes for the first individual
# loop over all chr - to obtain testsnps df for all chr
autosome_sk_oneid<-list()
for (chr in (1:22)) {  
autosome_sk_oneid[[chr]]<-random_intersect_introgbed_regions(ind,SPE,spe,chr)
}
aut_oneid<-do.call(rbind,autosome_sk_oneid) 

# randomly sample heterozygotes for the GTs per site - equal probability 0 or 1 - 0.5
for(id in (5:ncol(aut_oneid))) {
aut_oneid[,id][which(aut_oneid[,id]=="0|1")]<-sample(c(0,1),size=length(which(aut_oneid[,id]=="0|1")),replace=T,prob=c(0.5,0.5))
}

aut_oneid[which(aut_oneid=="0|0")]<-0
aut_oneid[which(aut_oneid=="1|1")]<-1

return(aut_oneid)
}


generate_tree_perid<-function(ind,SPE,spe) {
# run through for 1 individual first
#sbset<-test_gts_tonj(sk_regions_gbb,1,"GBB","gbb") # if is fine, then can amend to run for all individuals # yes works

sbset<-random_test_gts_topca(ind,SPE,spe)
# MK:
#2) maybe the reason why K80 didnt work was that you ended up 
#with monomorphic sites in your dataset, I suggest to remove these:
ssb<-matrix(as.numeric(sbset[,-c(1:4)]),ncol=ncol(sbset)-4)
#head(ssb)
#     [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12] [,13] [,14]
#[1,]    1    1    1    1    1    1    1    1    1     1     1     1     1     1
#[2,]    1    1    1    0    0    0    1    1    1     1     1     1     0     0
sr<-rowSums(ssb)
r2<-which(sr==0)
sbset2<-sbset[-r2,] # remove sites where all ids are homozy for the ref allele

#3) You may want to retain only a subset of the Homo-Gorilla differences, 
#as it is just for rooting and there are plenty of such sites:
r1<-which(sr[-r2]==49) # sample sites where all gor ids are fixed derived (relative to the human ref)
sbset2<-sbset2[-sample(r1,size=round(length(r1)*.9)),]
#head(sbset2)
#     [,1]   [,2]      [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12]
#[1,] "chr1" "5200698" "T"  "C"  "1"  "1"  "1"  "0"  "0"  "0"   "1"   "1"  
#[2,] "chr1" "5200750" "G"  "T"  "1"  "1"  "1"  "0"  "0"  "0"   "1"   "0"  
#[3,] "chr1" "5200756" "T"  "C"  "0"  "0"  "0"  "1"  "0"  "1"   "1"   "1"  
 
#-----------------------------------------------------------------------------------------------------------------------
# numeric GTs -> letters
hold_rows<-list()
for (i in (1:nrow(sbset2))) { 
tes<-sbset2[i,]
ref<-tes[3]
alt<-tes[4]
tes[tes == 0] <- paste0(ref)
tes[tes == 1] <- paste0(alt)
hold_rows[[i]]<-tes
}
#-----------------------------------------------------------------------------------------------------------------------

gts<-do.call(rbind,hold_rows)

#> head(gts)
#     [,1]   [,2]      [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12]
#[1,] "chr1" "5200698" "T"  "C"  "C"  "C"  "C"  "T"  "T"  "T"   "C"   "C"  
#[2,] "chr1" "5200750" "G"  "T"  "T"  "T"  "T"  "G"  "G"  "G"   "T"   "G"  

gts_2<-t(gts[,c(3,5:ncol(gts))])
#str(gts_2)
# chr [1:50, 1:383816] "A" "T" "T" "T" "T" "T" "T" "T" "T" "T" "T" "T" "T" ...

sample<-c('Homo_sapiens',samples.in.vcf)
rownames(gts_2)<-sample

db_2<-as.DNAbin(gts_2)

# MK
#4) Now, the dm works with K80:
#> dm_2 = dist.dna(db_2,pairwise.deletion=T,model="K80")

dm_2 = dist.dna(db_2,pairwise.deletion=F,model="K80")

# head(dm_2)
#[1] 0.4978397 0.4821156 0.4840672 0.4907302 0.4880527 0.4841037

outg<-c('Homo_sapiens')

treeNJ = NJ(dm_2);  fit = pml(treeNJ, data=phyDat(gts_2)); fitJC = optim.pml(fit, TRUE)
#treeNJ
#Phylogenetic tree with 50 tips and 48 internal nodes.

#-----------------------------------------------------------------------------------------------------------------------
rootedNJ <- ladderize(root(treeNJ, outgroup=outg, resolve.root=TRUE))
#-----------------------------------------------------------------------------------------------------------------------
#rootedNJ
#Phylogenetic tree with 50 tips and 49 internal nodes.
#-----------------------------------------------------------------------------------------------------------------------

# alternate way of doing the bootstrapping
#boots <- boot.phylo(phy=rootedNJ,
#                    x=db_2,FUN=function(bs_db_2) nj(dist.dna(bs_db_2, model="K80")),
#                    trees=TRUE)

# Add bootstraps to PHYLO object
#rootedNJ$node.label<-boots

save(rootedNJ, file=paste("/scratch/devel/hpawar/admix/overlap.s*.skov/trees/random/",SPE,"/",spe,".id.",ind,".random.rootedNJploutg.Robj.tree",sep=""))


#bss = bootstrap.pml(fitJC, bs=100, optNni=TRUE, control = pml.control(trace = 0)) # this step was causing recursive errors in the cluster

# drop the outgroup before doing the bootstrapping
rootedNJ <- drop.tip(rootedNJ, "Homo_sapiens")
#bss = bootstrap.pml(fitJC, bs=10, optNni=TRUE, control = pml.control(trace = 0)) # this works
bss = bootstrap.pml(fitJC, bs=100, optNni=TRUE, control = pml.control(trace = 0)) # generate 100 bs replicates instead of 10

random_out<-list(rootedNJ,bss)

save(random_out, file=paste("/scratch/devel/hpawar/admix/overlap.s*.skov/trees/random/",SPE,"/",spe,".id.",ind,".random.Robj.tree",sep=""))

#write.tree(rootedNJ, file=paste("/scratch/devel/hpawar/admix/overlap.s*.skov/trees/random/",SPE,"/",spe,".id.",ind,".random.tree",sep="")  )

}

generate_tree_perid(il,SPEC,spec)


#nj_allids<-function(scen,SPE,spe) {
# run for all individuals - may be better to split into separate jobs??
#for (ind in (1:length(scen))) {
#generate_tree_perid(scen,ind,SPE,spe)}
#}

# where need to specify il
#generate_tree_perid(ov_gbb_99,il,"GBB","gbb")
#generate_tree_perid(ov_gbb40_99,il,"GBB","gbb40")
#generate_tree_perid(ov_gbg_99,il,"GBG","gbg")
#generate_tree_perid(ov_gbg40_99,il,"GBG","gbg40")

#nj_allids(ov_gbb_99,"GBB","gbb")
#nj_allids(ov_gbb40_99,"GBB","gbb.40")
#nj_allids(ov_gbg_99,"GBG","gbg")
#nj_allids(ov_gbg40_99,"GBG","gbg.40")


