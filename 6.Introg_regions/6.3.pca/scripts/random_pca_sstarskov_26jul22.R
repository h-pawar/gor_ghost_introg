# Tue 26 Jul 2022 11:54:57 CEST
# random regions of equivalent length distribution to overlap of (99% s* - strict skov) & generate pcas of these random regions

# amending from /scratch/devel/hpawar/admix/sstar/scripts/simul_postabc/downstream/random_pca_skov_3jun22.R

#-----------------------------------------------------------------------------------------------------------------------
#module load gcc/6.3.0 openssl/1.0.2q  R/4.0.1  BCFTOOLS/1.14
require(data.table)
options(stringsAsFactors=F)
options(scipen=100)
library(adegenet)
#library(mgcv)
library(GenomicRanges)
#library(tidyr)
#library(ggbio)
#library(pheatmap)
library(ggplot2)
#library(ggpubr)
library(valr) # need to generate random regions
library(tidyr)
library(tidyverse)

#-----------------------------------------------------------------------------------------------------------------------

# directly read in the intersect regions - Wed 27 Jul 2022 10:15:02 CEST

# need to read in the empirical data in order to generate equivalent random regions
# (1) overlap s*-skov regions - convert to granges -> obtain reduced granges -> bed -> intersect with vcf
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
 # starperind.gbb[[ind]]<-indval[which(indval[,1]>indval[,4]),,drop=F]
 # }

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

# these random regions are all fine
# could output here also - the error was at the end of the run script -> perhaps unnecessary to output the random regions here

#ran_regions_out<-list(sk_ran_gbb,sk_ran_gbg,sk_ran_gbb_40,sk_ran_gbg_40)
#save(ran_regions_out,file="/scratch/devel/hpawar/admix/skov/GBB/onerep_randomregions_3jun22")

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
# check overlap b/n these 2 objects (ie between the empirical & the random regions for individual 1)

# check this section **

#calc_prop_fun<-function(fG,fB){
#fG<-reduce(fG);fB<-reduce(fB) # amendment - Wed 23 Feb 2022 17:28:27 CET
#pos.overlaps<-reduce(subsetByOverlaps(fG,fB)) ## this may be more useful output? (gives ranges which overlap ie the exact positions) (yes - these are the unique overlapping regions)
#prop.bp.x<-sum(reduce(subsetByOverlaps(fG,fB))@ranges@width)/sum(fG@ranges@width) # proportion of overlapping bp
#return(prop.bp.x)
#return(rbind(out.x,out.y))
#}

# run this for each scenario as a sanity check - that the regions are random
#calc_checkran<-function(fG,ran_fG){
#check_ran_regions<-list()
#for (ind in (1:length(fG))) {
#check_ran_regions[[ind]]<-calc_prop_fun(fG[[ind]],ran_fG[[ind]])}
#return(check_ran_regions)
#}


#calc_checkran(ov_gbb_99,sk_ran_gbb)
#calc_checkran(ov_gbg_99,sk_ran_gbg)
#calc_checkran(ov_gbb40_99,sk_ran_gbb_40)
#calc_checkran(ov_gbg40_99,sk_ran_gbg_40)

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

#[hpawar@login2 ~]$ mkdir -p /scratch/devel/hpawar/admix/overlap.s*.skov/regionsperid/random/GBB
#[hpawar@login2 ~]$ mkdir -p /scratch/devel/hpawar/admix/overlap.s*.skov/regionsperid/random/GBG

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
tmp<-paste("/scratch/devel/hpawar/admix/overlap.s*.skov/regionsperid/random/",SPE,"/",spe,".",a,".random.tmp.bed",sep="")

write.table(df,tmp,sep="\t",row.names=F,col.names=F,quote=F) # should write this out only once

}
# 2/6/22 - random_granges_tobed - function & output looks fine - 

#>  head(read.table(tmp))
#     V1        V2        V3
#1  chr1  67179803  67180804
#2 chr20  14127100  14128101
#3  chr2 235120469 235122470

#-----------------------------------------------------------------------------------------------------------------------

#bed file -> intersect with vcf, extract biallelic snps only -> output chr, pos, gts

random_intersect_introgbed_regions<-function(scen,ind,SPE,spe,chrom){
a=ind
tmp<-paste("/scratch/devel/hpawar/admix/overlap.s*.skov/regionsperid/random/",SPE,"/",spe,".",a,".random.tmp.bed",sep="")

chrom=chrom
# filter by biallelic snps only, query by the putative introg regions for id 1 & output chr, pos & gts
testsnps<-system(paste("bcftools view -m2 -M2 -v snps -R ",tmp," /scratch/devel/mkuhlwilm/arch/subsets/3_gorilla_",chrom,".vcf.gz | bcftools query -f '%CHROM %POS [%GT ]\n' ",sep=""),intern=T) 
testsnps<-do.call(rbind,strsplit(testsnps,split=" "))  ## ie list merged vertically into df

return(testsnps)
}

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

# Wed 18 May 2022 12:34:45 CEST - works but memory heavy - need to send as jobs or at least mnsh -c 4

# generate bed file of introg regions for 1 individual
# intersect bed with all autosomal vcfs - extract chr, pos, gts
# convert gts to matrix of 0,1,2
# calculate pca
# output PC 1,2 & % variance explained by each

random_test_gts_topca<-function(scen,ind,SPE,spe) {

random_granges_tobed(scen,ind,SPE,spe)

# extract biallelic GTs for introg regions across all autosomes for the first individual
# loop over all chr - to obtain testsnps df for all chr
autosome_sk_oneid<-list()
for (chr in (1:22)) {  
autosome_sk_oneid[[chr]]<-random_intersect_introgbed_regions(scen,ind,SPE,spe,chr)
}

aut_oneid<-do.call(rbind,autosome_sk_oneid) 

aut_oneid[which(aut_oneid=="0|0")]<-0
aut_oneid[which(aut_oneid=="0|1")]<-1
aut_oneid[which(aut_oneid=="1|1")]<-2

#return(aut_oneid)
gts<-aut_oneid

# convert gts to matrix
pcsub<-matrix(as.numeric(gts[,c(3:ncol(gts))]),nrow=nrow(gts)) # set to matrix

# calculate pca
pca_tes <- dudi.pca(t(pcsub),nf=30,scannf=F)
#scannf = FALSE,   # Hide scree plot: % of variance explained by each principal component # http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/119-pca-in-r-using-ade4-quick-scripts/
# nf = 30           # Number of components kept in the results

#pca_tes$eig # eigenvalues 
#(pca_tes$li) li #the row coordinates i.e. the principal components

# % of variance explained by each component:
(k <- 100 * pca_tes$eig/sum(pca_tes$eig))

#-----------------------------------------------------------------------------------------------------------------------
# samples & their populations
# for ids
#samples.in.vcf<-system(paste("bcftools query -l /scratch/devel/mkuhlwilm/arch/subsets/3_gorilla_22.vcf.gz",sep=""),intern=T) # -l, --list-samples: list sample names and exit

# add population identifier
#identifiers<-cbind(samples.in.vcf, c(rep("GBB",12), rep("GBG",9), rep("GGD",1),rep("GGG",27)))
#colnames(identifiers)<-c("ids","pop")
#-----------------------------------------------------------------------------------------------------------------------

# take the PCs 1-2, with the population identifiers
#pc1pc2<-as.data.frame(cbind(pca_tes$li[,c(1:2)],identifiers[,2]))
# may be best to output all PCs & all ks Thu 19 May 2022 09:31:18 CEST
pc1pc2<-as.data.frame(cbind(pca_tes$li,identifiers[,2]))
#colnames(pc1pc2)<-c("pc1","pc2","pop")

x<-as.numeric(format(round(k[[1]], 2), nsmall = 2) )
xaxis=paste("PC1 (",x,"%)",sep="")
y<-as.numeric(format(round(k[[2]], 2), nsmall = 2) )
yaxis=paste("PC2 (",y,"%)",sep="")

x3<-as.numeric(format(round(k[[3]], 2), nsmall = 2) )
x3axis=paste("PC3 (",x3,"%)",sep="")
y4<-as.numeric(format(round(k[[4]], 2), nsmall = 2) )
y4axis=paste("PC4 (",y4,"%)",sep="")

a=ind

# write this out
pca_out<-paste("/scratch/devel/hpawar/admix/overlap.s*.skov/regionsperid/random/",SPE,"/",spe,".",a,".random.pca",sep="")
#tes<-list(pc1pc2,xaxis,yaxis)
tes<-list(pc1pc2,k,xaxis,yaxis,x3axis,y4axis)

save(tes,file=pca_out)

#return(list(pc1pc2,xaxis,yaxis))

return(list(pc1pc2,k,xaxis,yaxis,x3axis,y4axis))

}
#-----------------------------------------------------------------------------------------------------------------------

# is not outputting the plots - will need to save the .random.pca objects, output of random_test_gts_topca then subsequently read in & plot

# run random_test_gts_topca function for all individuals in the population & output the PCs 1-2 to R objects (then subsequently plot these)

random_pca_allids<-function(scen,SPE,spe) {
out<-list()
for (ind in (1:length(scen))) {
#out[[ind]]<- random_test_gts_topca(sk_regions_gbb,ind,SPE,spe)} 
# here the sk_regions_gbb must be the problem *** & same issue in the empirical - will need to rerun for the other 3 scenarios ** (gbg, gbb40, gbg40)
out[[ind]]<- random_test_gts_topca(scen,ind,SPE,spe)}
#return(out)
}

# apply pca_all_ids function to all scenarios - GBB, GBG, with & without 40kb cutoff
random_pca_allids(sk_ran_gbb,"GBB","gbb")
random_pca_allids(sk_ran_gbg,"GBG","gbg")

random_pca_allids(sk_ran_gbb_40,"GBB","gbb.40")
random_pca_allids(sk_ran_gbg_40,"GBG","gbg.40")

