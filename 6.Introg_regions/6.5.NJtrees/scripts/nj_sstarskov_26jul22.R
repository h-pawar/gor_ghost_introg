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
#library(ggplot2)

#-----------------------------------------------------------------------------------------------------------------------

args = commandArgs(trailingOnly=TRUE)

il=args[1] 
SPEC=args[2]
spec=args[3] 

#-----------------------------------------------------------------------------------------------------------------------

# read in the intersect regions directly
# overlap s* 99% outliers with strict skov outliers


# glm generated for the S* CIs calc from simulated data (null simulations - gor demog, no ghost)
#load(file=paste("/scratch/devel/hpawar/admix/sstar/simul_postabc/stepwisesegs/final_8apr22/check_22apr22/glm",sep=""),verbose=T)
	# issue this needs to be in module load gcc/6.3.0 openssl/1.0.2q  R/4.0.1  : installed adegenet in newer r version

#-----------------------------------------------------------------------------------------------------------------------
# function to read in sstar outlier data per scenario & per CI & split this into s* windows per target individual

#proc_proportion<-function(nput, ci) {
  # for nput=1 # GBG

# read in the empirical data
#chroms=1:23 #there is no chrom 23.star file -> this is effectively 1:22 (reading in the autosomal data) 
# will need to specify either GB or GG
#cn1<-list(c("GB","GG"))

  # update to per subspecies comparisons
# 1) WL outgroup, EL ingroup (GGG outg, GBG ing) = "GBG"
# 2) WL outgroup, EM ingroup (GGG outg, GBB ing) = "GBB"
# 3) EL outgroup, WL ingroup (GBG outg, GGG ing) = "GGG"
# 4) EL outgroup, WC ingroup (GBG outg, GGD ing) = "GGD"


#cn1<-list(c("GBG","GBB","GGG","GGD"))


# 1)  GBG, CI 3, 99.5%

# for nput=1 # GBG

#nput=1

#starout<-list()
#for (chrom in (chroms)) {
#starout[[chrom]]<-read.table(paste("/scratch/devel/hpawar/admix/sstar/out.subsp.comp/",cn1[[1]][[nput]],"_",chrom,".star",sep=""),sep="\t",header=T)[,c(9,4,10,52,12,13,7,1,2)]
#}
# read in data for s* applied to chr 9 for target individuals plus newly processed tumani (whose chr 9 had sparse data issues)
#starout[[9]]<-read.table(paste("/scratch/devel/hpawar/admix/sstar/out.subsp.comp/",cn1[[1]][[nput]],"_newT_9.star",sep=""),sep="\t",header=T)[,c(9,4,10,52,12,13,7,1,2)] 

#staroutgbb<-do.call(rbind,starout)
# ensure all tumani fragments are named consistently
#staroutgbb$ind_id[which(staroutgbb$ind_id=="Gbg-Tumani")]<-"Gorilla_beringei_graueri-Tumani"  

#get values per individual - for each comparison
#staroutgbb[,8]<-paste(staroutgbb[,8],staroutgbb[,9])
#indiv.gbb<-unique(staroutgbb[,7]) # col 7 = ind_id


#starperind.gbb<-list()
#for (ind in (1:length(indiv.gbb))) {
#  starperind.gbb[[ind]]<-list()
  # only use windows where there is enough data (33,333 bp out of 40,000, with at least 30 segregating sites)
 # allstars<- staroutgbb[which(staroutgbb[,7]==indiv.gbb[ind] & staroutgbb[,4]>=33333 &staroutgbb[,2]>=30),c(1,8,2)]
 # allstars[,3]<-as.numeric(jitter(allstars[,3]))
  # create a data frame of the same segregating sites:
 # newdatA=data.frame(sS=allstars[,3])
  # predict S* vals (given the segregating sites) at given ci
 # out.pred<-predict(mods[[nput]][[ci]],newdata=newdatA,type="response")
 # indval<-cbind(allstars,out.pred)
  # which windows lie outside the expectation for the ci (3) 
#  starperind.gbb[[ind]]<-indval[which(indval[,1]>indval[,4]),,drop=F]
#  }

## you can see that each individual has a different number of significant regions:
#iva<-c()
#for (ind in 1:length(indiv.gbb)) { iva[ind]<-nrow(starperind.gbb[[ind]]) }

#iva
#[1] 900 842 898 876 883 800 864 818 876
#> length(iva)
#[1] 9

#converttoranges_function<-function(nput) {
#test<-starperind.gbb[[nput]][,c(1,2)]
#test1<-separate(test, chrom, into = c("chr", "winstart"), sep = " (?=[^ ]+$)")
#test1$winstart<-as.numeric(test1$winstart)
#test1$winend<-test1$winstart + 40000
#fG<-GRanges(seqnames=test1[,2],ranges=IRanges(start=as.numeric(test1[,3]),end=as.numeric(test1[,4]),names=test1[,2]),sstar=test1[,1],strand=rep("*",length(test1[,1])))
#return(fG)
#}


#allids<-list()
#for (ind in (1:length(indiv.gbb))) {
#allids[[ind]]<-converttoranges_function(ind)
#}

#reduceallids<-list()
#for (ind in (1:length(indiv.gbb))) {
#reduceallids[[ind]]<-reduce(allids[[ind]])
#}

#return(reduceallids)
#}

#-----------------------------------------------------------------------------------------------------------------------
# 1) read in sstar outlier data, per scenario, per CI (splits per individual & converts to reduced granges objects)

# read in sstar outlier data for easterns, 99% CI
# GBG, CI 2 (99%)
#sstar_GBG_99<-proc_proportion(1,2)

# GBB, CI 2 (99%)
#sstar_GBB_99<-proc_proportion(2,2)
#-----------------------------------------------------------------------------------------------------------------------

# (1) read in skov strict outliers & using function below split per inidivdual & generate reduced granges objects

# skov HMM bed file with the putative introgressed fragments ( ~2% with strict cutoff)
#sk<-read.table("/scratch/devel/mkuhlwilm/arch/skov/gor_fragments_strict.bed")
#head(sk)
#    V1       V2       V3                                V4
#1 chr1  5199000  5219000 Gorilla_beringei_beringei-Bwiruka
#2 chr1  5402000  5460000 Gorilla_beringei_beringei-Bwiruka

#> nrow(sk)
#[1] 15691

# convert startpos & endpos from int to numeric values
#sk[,c(2:3)]<-as.data.frame(lapply(sk[,c(2:3)], as.numeric)) # not necessary - b/c doing this step in the convert to granges function below

#head(sk)
#    V1       V2       V3                                V4
#1 chr1  5199000  5219000 Gorilla_beringei_beringei-Bwiruka
#2 chr1  5402000  5460000 Gorilla_beringei_beringei-Bwiruka

# individuals for GBB & GBG
#sk_ids_gbb<-unique(sk$V4)[1:12]
#sk_ids_gbg<-unique(sk$V4)[13:21]

#sk<-subset(sk, V1 != "chrX")

#-----------------------------------------------------------------------------------------------------------------------


# convert skov outliers to granges objects per individual
#skov_proc_proportion<-function(lids) {

# split skov data -> skov fragments per individual
#sk_per_id<-list()
#for (ind in (1:length(lids))) {
#sk_per_id[[ind]]<-subset(sk, V4 == lids[ind])
#}

# convert skov per id also to granges 
#skov_converttoranges_function<-function(nput) {
#sk<-sk_per_id[[nput]]  
#skranges<-GRanges(seqnames=sk[,1],ranges=IRanges(start=as.numeric(sk[,2]),end=as.numeric(sk[,3]),names=sk[,1]),strand=rep("*",length(sk[,1])))
#return(skranges)
#}

#sk_allids<-list()
#for (ind in (1:length(lids))) {
#sk_allids[[ind]]<-skov_converttoranges_function(ind)
#}

#reduce_sk_allids<-list()
#for (ind in (1:length(lids))) {
#reduce_sk_allids[[ind]]<-reduce(sk_allids[[ind]])
#}
# perhaps already reduced? in the sk_allids step?
#return(reduce_sk_allids)
#}

#-----------------------------------------------------------------------------------------------------------------------

#sk_regions_gbb<-skov_proc_proportion(sk_ids_gbb) # works
#sk_regions_gbg<-skov_proc_proportion(sk_ids_gbg)

# 2.1) filter strict skov regions by length - 40kb cutoff
#sk$V5<-sk$V3-sk$V2
#sk<-sk[sk$V5 >= 40000, ]
#sk_40kbregions_gbb<-skov_proc_proportion(sk_ids_gbb) # works
#sk_40kbregions_gbg<-skov_proc_proportion(sk_ids_gbg)
#-----------------------------------------------------------------------------------------------------------------------
# there are X chr in the skov outliers but not in the s* -> not a fair comparison

#rm_xchr_skov<-function(fG){
#aut_sk_regions_gbb<-list()
#for (ind in (1:length(fG))) {
#test_sk<-fG[[ind]]  
#aut_sk_regions_gbb[[ind]]<-test_sk[seqnames(test_sk) != "chrX"] }
#return(aut_sk_regions_gbb)
#}


#autosomes_sk_gbb<-rm_xchr_skov(sk_regions_gbb)
#autosomes_sk_gbg<-rm_xchr_skov(sk_regions_gbg)


#autosomes_40sk_gbb<-rm_xchr_skov(sk_40kbregions_gbb)
#autosomes_40sk_gbg<-rm_xchr_skov(sk_40kbregions_gbg)

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

# overlap skov outliers with S* outliers

# Wed 20 Jul 2022 17:35:51 CEST
#calc_prop_fun<-function(fG,fB){
#fG<-reduce(fG);fB<-reduce(fB) # amendment - Wed 23 Feb 2022 17:28:27 CET
#pos.overlaps<-reduce(subsetByOverlaps(fG,fB)) ## this may be more useful output? (gives ranges which overlap ie the exact positions) (yes - these are the unique overlapping regions)
## shoudl calc proportion overlapping both ways (ie of A how much overlaps with B + vice versa)
## & calc propotion of only unique overlapping regions
#prop.bp.x<-sum(reduce(subsetByOverlaps(fG,fB))@ranges@width)/sum(fG@ranges@width) # proportion of overlapping bp
#prop.frag.x<-length(reduce(subsetByOverlaps(fG,fB)))/length(fG) # proportion of overlapping fragments
#prop.bp.y<-sum(reduce(subsetByOverlaps(fG,fB))@ranges@width)/sum(fB@ranges@width) # proportion of overlapping bp
#prop.frag.y<-length(reduce(subsetByOverlaps(fG,fB)))/length(fB) # proportion of overlapping fragments
## shoudl combine the following into 1 output? b/c also need to output for vice versa B in A
#out.x<-(cbind(prop.bp.x, prop.frag.x))
#out.y<-(cbind(prop.bp.y, prop.frag.y))
#return(pos.overlaps)
#return(rbind(out.x,out.y))
#}
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
# process output of overlaps : regions overlapping, length dist, proportion (overlapping bp & fragments) - order of query matters
#process_tooverlap_fun<-function(sstar_subs,skov_subs){
# pairwise comparison of outlier windows : diff number of ranges obtained depending on which is query & which labelled as source
# ie query s* outliers by skov bed files & vice versa
#ranges_1<-reduce(subsetByOverlaps(sstar_subs,skov_subs)) 
#ranges_2<-reduce(subsetByOverlaps(skov_subs,sstar_subs))
# length distribution of overlapping regions
#lengths_1<-summary(reduce(subsetByOverlaps(sstar_subs,skov_subs))@ranges@width)/1000
#lengths_2<-summary(reduce(subsetByOverlaps(skov_subs,sstar_subs))@ranges@width)/1000
# proportion overlapping
#col1 = proportion of overlapping bp for unique overlapping regions
#col2 = proportion of overlapping fragments for unique overlapping regions
#prop_1<-calc_prop_fun(sstar_subs,skov_subs)
#prop_2<-calc_prop_fun(skov_subs,sstar_subs)
#return(list(ranges_1,ranges_2,lengths_1,lengths_2,prop_1,prop_2))
# only output ranges here? - RUN THROUGH INTERACTIVELY TMRW *******
#return(list(ranges_2,lengths_2,prop_2))
#}

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

# SEE IF THIS WORKS & GIVES EQUIVALENT RESULTS TO PREVIOUS FN ******** yes, but only outputs the ranges (not the proportions or the lengths)
# or only output ranges_2
#intersect_function<-function(sstar_subs,skov_subs) {
# 1) compare GBB - 95% CI for S*
#sstar_subs #sk_regions_gb
#out_gbb_95<-list()
#for (ind in (1:length(sstar_subs))) {
#out_gbb_95[[ind]]<-process_tooverlap_fun(sstar_subs[[ind]],skov_subs[[ind]])[[1]]
#}
#return(out_gbb_95)
#}

#ov_gbb_99<-intersect_function(sstar_GBB_99,autosomes_sk_gbb)
#ov_gbg_99<-intersect_function(sstar_GBG_99,autosomes_sk_gbg)
#ov_gbb40_99<-intersect_function(sstar_GBB_99,autosomes_40sk_gbb)
#ov_gbg40_99<-intersect_function(sstar_GBG_99,autosomes_40sk_gbg)

#-----------------------------------------------------------------------------------------------------------------------


load(file="/scratch/devel/hpawar/admix/overlap.s*.skov/overlap_objects/ov_gbb_99",verbose=T)
load(file="/scratch/devel/hpawar/admix/overlap.s*.skov/overlap_objects/ov_gbg_99",verbose=T)
load(file="/scratch/devel/hpawar/admix/overlap.s*.skov/overlap_objects/ov_gbb40_99",verbose=T)
load(file="/scratch/devel/hpawar/admix/overlap.s*.skov/overlap_objects/ov_gbg40_99",verbose=T)

#-----------------------------------------------------------------------------------------------------------------------
library(phangorn)  # check if installed in this version of R>.- no - will need to install an older version, or move to a newer R version - that may be preferable
library(ape) 
library('pegas')
#-----------------------------------------------------------------------------------------------------------------------

# granges of introg regions for 1 individual -> bed file 
  # contains introg regions across all autosomes for the specified individual
#granges_tobed<-function(scen,ind,SPE,spe){

#df <- data.frame(seqnames=seqnames(scen[[ind]]),
#  start=start(scen[[ind]])-1,
#  end=end(scen[[ind]]))

# convert to same structure of sk
#df[,1]<-as.character(df[,1])
#df[,2]<-as.integer(df[,2])

#a=ind
#tmp<-paste("/scratch/devel/hpawar/admix/overlap.s*.skov/regionsperid/",SPE,"/",spe,".",a,".tmp.bed",sep="")

#write.table(df,tmp,sep="\t",row.names=F,col.names=F,quote=F) # should write this out only once

#}

#-----------------------------------------------------------------------------------------------------------------------

#bed file -> intersect with vcf, extract biallelic snps only -> output chr, pos, gts

intersect_introgbed_regions<-function(ind,SPE,spe,chrom){
a=ind
tmp<-paste("/scratch/devel/hpawar/admix/overlap.s*.skov/regionsperid/",SPE,"/",spe,".",a,".tmp.bed",sep="")

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

save(rootedNJ, file=paste("/scratch/devel/hpawar/admix/overlap.s*.skov/trees/",SPE,"/",spe,".id.",ind,".rootedNJploutg.Robj.tree",sep=""))

# drop the outgroup before doing the bootstrapping
rootedNJ <- drop.tip(rootedNJ, "Homo_sapiens")
#bss = bootstrap.pml(fitJC, bs=10, optNni=TRUE, control = pml.control(trace = 0)) # this works
bss = bootstrap.pml(fitJC, bs=100, optNni=TRUE, control = pml.control(trace = 0)) # generate 100 bs replicates instead of 10


#bss = bootstrap.pml(fitJC, bs=100, optNni=TRUE, control = pml.control(trace = 0)) # this step was causing recursive errors in the cluster

# Add bootstraps to PHYLO object
#rootedNJ$node.label<-boots

#[hpawar@cnb4 regionsperid]$ mkdir -p /scratch/devel/hpawar/admix/overlap.s*.skov/trees/GBB
#[hpawar@cnb4 regionsperid]$ mkdir -p /scratch/devel/hpawar/admix/overlap.s*.skov/trees/GBG

out<-list(rootedNJ,bss)

#write.tree(rootedNJ, file=paste("/scratch/devel/hpawar/admix/overlap.s*.skov/trees/",SPE,"/",spe,".id.",ind,".tree",sep="")  )
# gives segmenation faults:
#test<-read.tree('/scratch/devel/hpawar/admix/overlap.s*.skov/trees/GBB/gbb.id.1.tree')
# *** caught segfault ***
#address 0x7ffd8126b000, cause 'memory not mapped'
#Segmentation fault

save(out, file=paste("/scratch/devel/hpawar/admix/overlap.s*.skov/trees/",SPE,"/",spe,".id.",ind,".Robj.tree",sep=""))

}



generate_tree_perid(il,SPEC,spec)

#generate_tree_perid(ov_gbb_99,1,"GBB","gbb") # test


# apply generate_tree_perid to all individuals of a given scenario
#nj_allids<-function(scen,SPE,spe) {
## run for all individuals - may be better to split into separate jobs??
#for (ind in (1:length(scen))) {
#generate_tree_perid(scen,ind,SPE,spe)}
#}

# where need to specify il
#generate_tree_perid(ov_gbb_99,il,"GBB","gbb")
#generate_tree_perid(ov_gbb40_99,il,"GBB","gbb40")
#generate_tree_perid(ov_gbg_99,il,"GBG","gbg")
#generate_tree_perid(ov_gbg40_99,il,"GBG","gbg40")

#nj_allids(ov_gbb_99,"GBB","gbb") # try for the first scenario only
#nj_allids(ov_gbb40_99,"GBB","gbb.40")
#nj_allids(ov_gbg_99,"GBG","gbg")
#nj_allids(ov_gbg40_99,"GBG","gbg.40")

# see if this works..

