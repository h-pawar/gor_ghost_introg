# Tue  4 Oct 2022 10:01:25 CEST
# generate nj trees of private regions (intersect sstar-skov regions which are found only in one id of the pop)

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
#-----------------------------------------------------------------------------------------------------------------------

args = commandArgs(trailingOnly=TRUE)

il=args[1] 
SPEC=args[2]
spec=args[3] 

# il=1
# SPE="GBB"
# spe="gbb"

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

#bed file -> intersect with vcf, extract biallelic snps only -> output chr, pos, gts

intersect_introgbed_regions<-function(ind,SPE,spe,chrom){
a=ind
tmp<-paste("/scratch/devel/hpawar/admix/overlap.s*.skov/privateregionsperid/",SPE,"/",spe,".",a,".private.tmp.bed",sep="")

chrom=chrom
# filter by biallelic snps only, query by the putative introg regions for id 1 & output chr, pos, ref allele, alt allele & gts
#testsnps<-system(paste("bcftools view -m2 -M2 -v snps -R ",tmp," /scratch/devel/mkuhlwilm/arch/subsets/3_gorilla_",chrom,".vcf.gz | bcftools query -f '%CHROM %POS [%GT ]\n' ",sep=""),intern=T) 
testsnps<-system(paste("bcftools view -m2 -M2 -v snps -R ",tmp," /scratch/devel/mkuhlwilm/arch/subsets/3_gorilla_",chrom,".vcf.gz | bcftools query -f '%CHROM %POS %REF %ALT [%GT ]\n' ",sep=""),intern=T) 
testsnps<-do.call(rbind,strsplit(testsnps,split=" "))  ## ie list merged vertically into df

return(testsnps)
}

#-----------------------------------------------------------------------------------------------------------------------
# samples & their populations
# for ids
samples.in.vcf<-system(paste("bcftools query -l /scratch/devel/mkuhlwilm/arch/subsets/3_gorilla_22.vcf.gz",sep=""),intern=T) # -l, --list-samples: list sample names and exit

# add population identifier
identifiers<-cbind(samples.in.vcf, c(rep("GBB",12), rep("GBG",9), rep("GGD",1),rep("GGG",27)))
colnames(identifiers)<-c("ids","pop")
#-----------------------------------------------------------------------------------------------------------------------

test_gts_topca<-function(ind,SPE,spe) {

#granges_tobed(scen,ind,SPE,spe) # already generated

# extract biallelic GTs for introg regions across all autosomes for the first individual
# loop over all chr - to obtain testsnps df for all chr
autosome_sk_oneid<-list()
for (chr in (1:22)) {  
autosome_sk_oneid[[chr]]<-intersect_introgbed_regions(ind,SPE,spe,chr)
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


# run through for 1 individual first
#sbset<-test_gts_topca(sk_regions_gbb,1,"GBB","gbb") #


generate_tree_perid<-function(ind,SPE,spe) {
# run through for 1 individual first
#sbset<-test_gts_tonj(sk_regions_gbb,1,"GBB","gbb") # if is fine, then can amend to run for all individuals # yes works

sbset<-test_gts_topca(ind,SPE,spe)
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


save(rootedNJ, file=paste("/scratch/devel/hpawar/admix/overlap.s*.skov/privateregionsperid/",SPE,"/nj/",spe,".id.",ind,".rootedNJploutg.Robj.tree",sep=""))

# drop the outgroup before doing the bootstrapping
rootedNJ <- drop.tip(rootedNJ, "Homo_sapiens")
#bss = bootstrap.pml(fitJC, bs=10, optNni=TRUE, control = pml.control(trace = 0)) # this works
bss = bootstrap.pml(fitJC, bs=100, optNni=TRUE, control = pml.control(trace = 0)) # generate 100 bs replicates instead of 10


#bss = bootstrap.pml(fitJC, bs=100, optNni=TRUE, control = pml.control(trace = 0)) # this step was causing recursive errors in the cluster

# Add bootstraps to PHYLO object
#rootedNJ$node.label<-boots

#[hpawar@cnb4 privateregionsperid]$ mkdir -p /scratch/devel/hpawar/admix/overlap.s*.skov/trees/GBB
#[hpawar@cnb4 privateregionsperid]$ mkdir -p /scratch/devel/hpawar/admix/overlap.s*.skov/trees/GBG

out<-list(rootedNJ,bss)

#write.tree(rootedNJ, file=paste("/scratch/devel/hpawar/admix/overlap.s*.skov/trees/",SPE,"/",spe,".id.",ind,".tree",sep="")  )
# gives segmenation faults:
#test<-read.tree('/scratch/devel/hpawar/admix/overlap.s*.skov/trees/GBB/gbb.id.1.tree')
# *** caught segfault ***
#address 0x7ffd8126b000, cause 'memory not mapped'
#Segmentation fault

save(out, file=paste("/scratch/devel/hpawar/admix/overlap.s*.skov/privateregionsperid/",SPE,"/nj/",spe,".id.",ind,".Robj.tree",sep=""))

}


generate_tree_perid(il,SPEC,spec)