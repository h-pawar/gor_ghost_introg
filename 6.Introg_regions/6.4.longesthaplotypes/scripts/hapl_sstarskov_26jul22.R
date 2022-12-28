# Tue 26 Jul 2022 12:23:10 CEST
# generate haplotype networks for intersect regions of length >=100kb
  # first identify how many such regions exist

# amending from - /Volumes/"Ultra USB 3.0"/IBE/further.analysis.feb.2020/gorillas/abc/simul_postabc/sstar.longhapl.easterns.29apr22.R

# RAN INTERACTIVELY - until line 472
# errors in the plotting - need to send as a separate script 

#-----------------------------------------------------------------------------------------------------------------------

# 0) read in empirical data & perform intersect
# A) identify & write out longest windows in easterns
# B) generate haplotype networks for these windows

#-----------------------------------------------------------------------------------------------------------------------
## module load gcc/6.3.0 openssl/1.0.2q  R/4.0.1  BCFTOOLS/1.12
require(data.table)
options(stringsAsFactors=F)
library(tidyverse)
options(scipen=100)
library(adegenet)
library(mgcv)
library(GenomicRanges)
#BiocManager::install(c("GenomicRanges", "plyranges", "HelloRangesData"))
#library(plyranges) # may not be necessary here
#library(GenomicFeatures) # if need to read in the gtf
library(tidyr)
#library(ggbio)
library(ggplot2)

#-----------------------------------------------------------------------------------------------------------------------
# overlap s* 99% outliers with strict skov outliers


# glm generated for the S* CIs calc from simulated data (null simulations - gor demog, no ghost)
load(file=paste("/scratch/devel/hpawar/admix/sstar/simul_postabc/stepwisesegs/final_8apr22/check_22apr22/glm",sep=""),verbose=T)
  # issue this needs to be in module load gcc/6.3.0 openssl/1.0.2q  R/4.0.1  : installed adegenet in newer r version

#-----------------------------------------------------------------------------------------------------------------------
# function to read in sstar outlier data per scenario & per CI & split this into s* windows per target individual

proc_proportion<-function(nput, ci) {
  # for nput=1 # GBG

# read in the empirical data
chroms=1:23 #there is no chrom 23.star file -> this is effectively 1:22 (reading in the autosomal data) 
# will need to specify either GB or GG
#cn1<-list(c("GB","GG"))

  # update to per subspecies comparisons
# 1) WL outgroup, EL ingroup (GGG outg, GBG ing) = "GBG"
# 2) WL outgroup, EM ingroup (GGG outg, GBB ing) = "GBB"
# 3) EL outgroup, WL ingroup (GBG outg, GGG ing) = "GGG"
# 4) EL outgroup, WC ingroup (GBG outg, GGD ing) = "GGD"


cn1<-list(c("GBG","GBB","GGG","GGD"))


# 1)  GBG, CI 3, 99.5%

# for nput=1 # GBG

#nput=1

starout<-list()
for (chrom in (chroms)) {
starout[[chrom]]<-read.table(paste("/scratch/devel/hpawar/admix/sstar/out.subsp.comp/",cn1[[1]][[nput]],"_",chrom,".star",sep=""),sep="\t",header=T)[,c(9,4,10,52,12,13,7,1,2)]
}
# read in data for s* applied to chr 9 for target individuals plus newly processed tumani (whose chr 9 had sparse data issues)
starout[[9]]<-read.table(paste("/scratch/devel/hpawar/admix/sstar/out.subsp.comp/",cn1[[1]][[nput]],"_newT_9.star",sep=""),sep="\t",header=T)[,c(9,4,10,52,12,13,7,1,2)] 

staroutgbb<-do.call(rbind,starout)
# ensure all tumani fragments are named consistently
staroutgbb$ind_id[which(staroutgbb$ind_id=="Gbg-Tumani")]<-"Gorilla_beringei_graueri-Tumani"  

#get values per individual - for each comparison
staroutgbb[,8]<-paste(staroutgbb[,8],staroutgbb[,9])
indiv.gbb<-unique(staroutgbb[,7]) # col 7 = ind_id


starperind.gbb<-list()
for (ind in (1:length(indiv.gbb))) {
  starperind.gbb[[ind]]<-list()
  # only use windows where there is enough data (33,333 bp out of 40,000, with at least 30 segregating sites)
  allstars<- staroutgbb[which(staroutgbb[,7]==indiv.gbb[ind] & staroutgbb[,4]>=33333 &staroutgbb[,2]>=30),c(1,8,2)]
  allstars[,3]<-as.numeric(jitter(allstars[,3]))
  # create a data frame of the same segregating sites:
  newdatA=data.frame(sS=allstars[,3])
  # predict S* vals (given the segregating sites) at given ci
  out.pred<-predict(mods[[nput]][[ci]],newdata=newdatA,type="response")
  indval<-cbind(allstars,out.pred)
  # which windows lie outside the expectation for the ci (3) 
  starperind.gbb[[ind]]<-indval[which(indval[,1]>indval[,4]),,drop=F]
  }

## you can see that each individual has a different number of significant regions:
iva<-c()
for (ind in 1:length(indiv.gbb)) { iva[ind]<-nrow(starperind.gbb[[ind]]) }

#iva
#[1] 900 842 898 876 883 800 864 818 876
#> length(iva)
#[1] 9

converttoranges_function<-function(nput) {
test<-starperind.gbb[[nput]][,c(1,2)]
test1<-separate(test, chrom, into = c("chr", "winstart"), sep = " (?=[^ ]+$)")
test1$winstart<-as.numeric(test1$winstart)
test1$winend<-test1$winstart + 40000
fG<-GRanges(seqnames=test1[,2],ranges=IRanges(start=as.numeric(test1[,3]),end=as.numeric(test1[,4]),names=test1[,2]),sstar=test1[,1],strand=rep("*",length(test1[,1])))
return(fG)
}


allids<-list()
for (ind in (1:length(indiv.gbb))) {
allids[[ind]]<-converttoranges_function(ind)
}

reduceallids<-list()
for (ind in (1:length(indiv.gbb))) {
reduceallids[[ind]]<-reduce(allids[[ind]])
}

return(reduceallids)
}

#-----------------------------------------------------------------------------------------------------------------------
# 1) read in sstar outlier data, per scenario, per CI (splits per individual & converts to reduced granges objects)

# read in sstar outlier data for easterns, 99% CI
# GBG, CI 2 (99%)
sstar_GBG_99<-proc_proportion(1,2)

# GBB, CI 2 (99%)
sstar_GBB_99<-proc_proportion(2,2)

#-----------------------------------------------------------------------------------------------------------------------

# (1) read in skov strict outliers & using function below split per inidivdual & generate reduced granges objects

# skov HMM bed file with the putative introgressed fragments ( ~2% with strict cutoff)
sk<-read.table("/scratch/devel/mkuhlwilm/arch/skov/gor_fragments_strict.bed")
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
sk_ids_gbb<-unique(sk$V4)[1:12]
sk_ids_gbg<-unique(sk$V4)[13:21]

sk<-subset(sk, V1 != "chrX")

#-----------------------------------------------------------------------------------------------------------------------


# convert skov outliers to granges objects per individual
skov_proc_proportion<-function(lids) {

# split skov data -> skov fragments per individual
sk_per_id<-list()
for (ind in (1:length(lids))) {
sk_per_id[[ind]]<-subset(sk, V4 == lids[ind])
}

# convert skov per id also to granges 
skov_converttoranges_function<-function(nput) {
sk<-sk_per_id[[nput]]  
skranges<-GRanges(seqnames=sk[,1],ranges=IRanges(start=as.numeric(sk[,2]),end=as.numeric(sk[,3]),names=sk[,1]),strand=rep("*",length(sk[,1])))
return(skranges)
}

sk_allids<-list()
for (ind in (1:length(lids))) {
sk_allids[[ind]]<-skov_converttoranges_function(ind)
}

reduce_sk_allids<-list()
for (ind in (1:length(lids))) {
reduce_sk_allids[[ind]]<-reduce(sk_allids[[ind]])
}
# perhaps already reduced? in the sk_allids step?
return(reduce_sk_allids)
}

#-----------------------------------------------------------------------------------------------------------------------

sk_regions_gbb<-skov_proc_proportion(sk_ids_gbb) # works
sk_regions_gbg<-skov_proc_proportion(sk_ids_gbg)

# 2.1) filter strict skov regions by length - 40kb cutoff
sk$V5<-sk$V3-sk$V2
sk<-sk[sk$V5 >= 40000, ]
sk_40kbregions_gbb<-skov_proc_proportion(sk_ids_gbb) # works
sk_40kbregions_gbg<-skov_proc_proportion(sk_ids_gbg)
#-----------------------------------------------------------------------------------------------------------------------
# there are X chr in the skov outliers but not in the s* -> not a fair comparison

rm_xchr_skov<-function(fG){
aut_sk_regions_gbb<-list()
for (ind in (1:length(fG))) {
test_sk<-fG[[ind]]  
aut_sk_regions_gbb[[ind]]<-test_sk[seqnames(test_sk) != "chrX"] }
return(aut_sk_regions_gbb)
}


autosomes_sk_gbb<-rm_xchr_skov(sk_regions_gbb)
autosomes_sk_gbg<-rm_xchr_skov(sk_regions_gbg)


autosomes_40sk_gbb<-rm_xchr_skov(sk_40kbregions_gbb)
autosomes_40sk_gbg<-rm_xchr_skov(sk_40kbregions_gbg)

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

# overlap skov outliers with S* outliers

# Wed 20 Jul 2022 17:35:51 CEST
calc_prop_fun<-function(fG,fB){
fG<-reduce(fG);fB<-reduce(fB) # amendment - Wed 23 Feb 2022 17:28:27 CET
pos.overlaps<-reduce(subsetByOverlaps(fG,fB)) ## this may be more useful output? (gives ranges which overlap ie the exact positions) (yes - these are the unique overlapping regions)
## shoudl calc proportion overlapping both ways (ie of A how much overlaps with B + vice versa)
## & calc propotion of only unique overlapping regions
prop.bp.x<-sum(reduce(subsetByOverlaps(fG,fB))@ranges@width)/sum(fG@ranges@width) # proportion of overlapping bp
prop.frag.x<-length(reduce(subsetByOverlaps(fG,fB)))/length(fG) # proportion of overlapping fragments
prop.bp.y<-sum(reduce(subsetByOverlaps(fG,fB))@ranges@width)/sum(fB@ranges@width) # proportion of overlapping bp
prop.frag.y<-length(reduce(subsetByOverlaps(fG,fB)))/length(fB) # proportion of overlapping fragments
## shoudl combine the following into 1 output? b/c also need to output for vice versa B in A
out.x<-(cbind(prop.bp.x, prop.frag.x))
out.y<-(cbind(prop.bp.y, prop.frag.y))
#return(pos.overlaps)
return(rbind(out.x,out.y))
}
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
# process output of overlaps : regions overlapping, length dist, proportion (overlapping bp & fragments) - order of query matters
process_tooverlap_fun<-function(sstar_subs,skov_subs){
# pairwise comparison of outlier windows : diff number of ranges obtained depending on which is query & which labelled as source
# ie query s* outliers by skov bed files & vice versa
#ranges_1<-reduce(subsetByOverlaps(sstar_subs,skov_subs)) 
ranges_2<-reduce(subsetByOverlaps(skov_subs,sstar_subs))
# length distribution of overlapping regions
#lengths_1<-summary(reduce(subsetByOverlaps(sstar_subs,skov_subs))@ranges@width)/1000
lengths_2<-summary(reduce(subsetByOverlaps(skov_subs,sstar_subs))@ranges@width)/1000
# proportion overlapping
#col1 = proportion of overlapping bp for unique overlapping regions
#col2 = proportion of overlapping fragments for unique overlapping regions
#prop_1<-calc_prop_fun(sstar_subs,skov_subs)
prop_2<-calc_prop_fun(skov_subs,sstar_subs)
#return(list(ranges_1,ranges_2,lengths_1,lengths_2,prop_1,prop_2))
# only output ranges here? - RUN THROUGH INTERACTIVELY TMRW *******
return(list(ranges_2,lengths_2,prop_2))
}

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

# SEE IF THIS WORKS & GIVES EQUIVALENT RESULTS TO PREVIOUS FN ******** yes, but only outputs the ranges (not the proportions or the lengths)
# or only output ranges_2
intersect_function<-function(sstar_subs,skov_subs) {
# 1) compare GBB - 95% CI for S*
#sstar_subs #sk_regions_gb
out_gbb_95<-list()
for (ind in (1:length(sstar_subs))) {
out_gbb_95[[ind]]<-process_tooverlap_fun(sstar_subs[[ind]],skov_subs[[ind]])[[1]]
}
return(out_gbb_95)
}

ov_gbb_99<-intersect_function(sstar_GBB_99,autosomes_sk_gbb)
ov_gbg_99<-intersect_function(sstar_GBG_99,autosomes_sk_gbg)
ov_gbb40_99<-intersect_function(sstar_GBB_99,autosomes_40sk_gbb)
ov_gbg40_99<-intersect_function(sstar_GBG_99,autosomes_40sk_gbg)

#-----------------------------------------------------------------------------------------------------------------------
# write these out

save(ov_gbb_99,file="/scratch/devel/hpawar/admix/overlap.s*.skov/overlap_objects/ov_gbb_99")
save(ov_gbg_99,file="/scratch/devel/hpawar/admix/overlap.s*.skov/overlap_objects/ov_gbg_99")
save(ov_gbb40_99,file="/scratch/devel/hpawar/admix/overlap.s*.skov/overlap_objects/ov_gbb40_99")
save(ov_gbg40_99,file="/scratch/devel/hpawar/admix/overlap.s*.skov/overlap_objects/ov_gbg40_99")

#-----------------------------------------------------------------------------------------------------------------------

# check length of each intersect object
#longestwindowsperid<-list()
#for (ind in (1:length(indiv.gbb))) {
#longestwindowsperid[[ind]]<-reduceallids[[ind]][which(width(reduceallids[[ind]])>100000),]
#}


#extract_longwin_fun<-function(scen){
#longestwindowsperid<-list()
#for (ind in (1:length(scen))) {
#longestwindowsperid[[ind]]<-scen[[ind]][which(width(scen[[ind]])>100000),]
#}
#return(longestwindowsperid)
#}

#extract_longwin_fun(ov_gbb_99)
  # i do obtain quite some regions, perhaps i should take a higher cutoff?

#  [[12]]
#GRanges object with 133 ranges and 0 metadata columns:
#        seqnames              ranges strand
#           <Rle>           <IRanges>  <Rle>
#    [1]     chr1   87184000-87293000      *
#    [2]     chr1 114713000-114827000      *
#    [3]     chr1 119334000-119480000      *
#    [4]     chr1 173707000-173906000      *
 #   [5]     chr1 215410000-215569000      *
 #   ...      ...                 ...    ...
 # [129]    chr19   31178000-31292000      *
 # [130]    chr20   31681000-31911000      *
 # [131]    chr20   36920000-37083000      *
 # [132]    chr20   59987000-60137000      *
 # [133]    chr22   23204000-23314000      *
 # -------
 # seqinfo: 22 sequences from an unspecified genome; no seqlengths



#extract_longwin_fun<-function(scen){
#longestwindowsperid<-list()
#for (ind in (1:length(scen))) {
#longestwindowsperid[[ind]]<-scen[[ind]][which(width(scen[[ind]])>250000),]
#}
#return(longestwindowsperid)
#}

#extract_longwin_fun(ov_gbb_99)  # much fewer ranges obtained, perhaps could perform for both length cutoffs *
  # yes do this 
#[[12]]
#GRanges object with 16 ranges and 0 metadata columns:
 #      seqnames              ranges strand
 #         <Rle>           <IRanges>  <Rle>
 #  [1]     chr2 169290000-169582000      *
 #  [2]     chr3 119406000-119931000      *
 #  [3]     chr3 167973000-168246000      *
#   [4]     chr5     8202000-8555000      *
#   [5]     chr5   74653000-74972000      *
#   ...      ...                 ...    ...
#  [12]    chr10 121532000-121973000      *
#  [13]    chr12   27235000-27514000      *
#  [14]    chr15   85833000-86146000      *
#  [15]    chr16   12070000-12514000      *
#  [16]    chr17   25738000-26011000      *
#  -------
#  seqinfo: 22 sequences from an unspecified genome; no seqlengths

  # perhaps only run this analysis for s* - skov (without 40kb length cutoff - b/c will just be a subset of the variation)
# write for both length cutoffs (100,000 bp & 250,000 bp)

# ie run this function for 
#-----------------------------------------------------------------------------------------------------------------------


extract_longwin_fun<-function(scen,x){
longestwindowsperid<-list()
for (ind in (1:length(scen))) {
longestwindowsperid[[ind]]<-scen[[ind]][which(width(scen[[ind]])>x),]
}
return(longestwindowsperid)
}


# mkdir -p /scratch/devel/hpawar/admix/overlap.s*.skov/hapl/tmp

l_gbb_100k<-extract_longwin_fun(ov_gbb_99,100000)
l_gbb_250k<-extract_longwin_fun(ov_gbb_99,250000)
l_gbg_100k<-extract_longwin_fun(ov_gbg_99,100000)
l_gbg_250k<-extract_longwin_fun(ov_gbg_99,250000)


# output this here
longwindows<-list(l_gbb_100k,l_gbb_250k,l_gbg_100k,l_gbg_250k)

save(longwindows,file=paste("/scratch/devel/hpawar/admix/overlap.s*.skov/hapl/longwindows",sep="")) 

# this has been output - Wed 27 Jul 2022 09:29:26 CEST

q()
# RUN THIS SCRIPT UNTIL HERE, THEN RUN /Volumes/"Ultra USB 3.0"/IBE/further.analysis.feb.2020/gorillas/abc/simul_postabc/hapl_sstarskov_26oct22.R   

#-----------------------------------------------------------------------------------------------------------------------


# AMEND BELOW


# B) haplotype networks

# 1) for GBG
# process regions for the haplotypes
proc_hapl<-function(input,spe) {

cn1<-list(c("GBG","GBB","GGG","GGD"))

# longestwindowsperid for this scenario

GBG<-input


# test for overlap b/n the individuals for these longest regions
common<- Reduce(intersect, GBG)

# convert 'GBG' object of regions to a GRangesList object 
myGRangesList<-GRangesList(GBG)

non_overlapping_region_counts<-function(x){
reduced <- reduce(unlist(myGRangesList))
consensusIDs <- paste0("consensus_", seq(1, length(reduced)))
mcols(reduced) <- do.call(cbind, lapply(myGRangesList, function(x) (reduced %over% x) + 0))
reducedConsensus <- reduced
mcols(reducedConsensus) <- cbind(as.data.frame(mcols(reducedConsensus)), consensusIDs)
consensusIDs <- paste0("consensus_", seq(1, length(reducedConsensus)))
return(reducedConsensus)
}

# has this worked?? # then take the rowsums? 
overlaps_GBG<-non_overlapping_region_counts()

rowSums(as.matrix(mcols(overlaps_GBG)[,1:length(GBG)]))

length(rowSums(as.matrix(mcols(overlaps_GBG)[,1:length(GBG)])))



extract_regions_fun<-function(){

regions<-as.data.frame(overlaps_GBG)

rout <- strsplit(as.character(regions$seqnames),'chr') 
#do.call(rbind, rout)
regions$chr<-as.numeric(do.call(rbind, rout)[,2])

# generate vcfs for these regions 

for (i in (1:nrow(regions))) {
out1<-paste("/scratch/devel/hpawar/admix/overlap.s*.skov/hapl/tmp/",spe,"/",spe,"_",regions$seqnames[i],"_",regions$start[i],"_",regions$end[i],".overlaps.vcf",sep="")
system(paste("bcftools view -r ",regions$seqnames[i],":",regions$start[i],"-",regions$end[i]," /scratch/devel/mkuhlwilm/arch/subsets/3_gorilla_",regions$chr[i],".vcf.gz > ",out1," ",sep=""),intern=T)
}

}

extract_regions_fun()

return(list(overlaps_GBG,rowSums(as.matrix(mcols(overlaps_GBG)[,1:length(GBG)])),length(rowSums(as.matrix(mcols(overlaps_GBG)[,1:length(GBG)])))
))
}




proc_hapl(l_gbb_100k,"gbb100k")
proc_hapl(l_gbb_250k,"gbb250k")
proc_hapl(l_gbg_100k,"gbg100k")
proc_hapl(l_gbg_250k,"gbg250k")

# set up sep dir - done
# this has been output

q()

#------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------

# having errors in this step - go through interactively - Wed 27 Jul 2022 10:48:09 CEST

# now read in & plot this - changing R versions **
#module unload R/4.0.1
#module load gcc/6.3.0 R/3.4.2 BCFTOOLS/1.12

require(data.table)
options(stringsAsFactors=F)
library(tidyr)
library('pegas')
library(vcfR)
options("scipen"=100)


# go over all regions

# generate the haplotypes
generate_hapl<-function(nput) {

cn1<-list(c("GBG","GBB","GGG","GGD"))

dir=paste("/scratch/devel/hpawar/admix/overlap.s*.skov/hapl/tmp/",nput,"/",sep="")

vcffiles <- paste(dir, list.files(path =dir, pattern=".vcf"), sep = "")

# write as function & plot for these 8

hhpnet<-list()
hind.hap<-list()
for(i in (1:length(vcffiles))) {
locs<-read.vcfR(file=vcffiles[[i]])
blocs<-vcfR2DNAbin(locs,consensus=T,extract.haps=F)
rownames(blocs)<-c(rep("GBB",12),rep("GBG",9),rep("GGD",1),rep("GGG",27))
hploc<-haplotype(blocs)
hhpnet[[i]]<-try(haploNet(hploc))
hind.hap[[i]]<-with(
  stack(setNames(attr(hploc, "index"), rownames(hploc))),
  table(hap=ind, pop=rownames(blocs)[values])
)
}


# now obtain the regions
format_region_fun<-function(i) {
t<-strsplit(as.character(vcffiles[[i]]), 'GGG_')
t1<- strsplit(as.character(t[[1]][2]), '_')
t2<-strsplit(as.character(t1[[1]]), '.overlaps')
#region<-unlist(t2)
region<-unlist(t2)[1:3]
return(region)
}


hregions<-list()
for(i in (1:length(hhpnet))) {
hregions[[i]]<-format_region_fun(i)
}



pdf(paste("/scratch/devel/hpawar/admix/overlap.s*.skov/hapl/tmp/",nput,"/plot/",nput,"_overlaps.hapl.pdf",sep=""))
for(i in (1:length(hhpnet))) {
region=hregions[[i]]
plot(hhpnet[[i]], size=attr(hhpnet[[i]], "freq"), scale.ratio = 2, cex = 0.8, pie=hind.hap[[i]],show.mutation=3,main=paste(region,sep="-"))
legend("bottomright", c("GBB","GBG","GGD","GGG"), col=rainbow(ncol(hind.hap[[i]])), pch=20) ## needs to be changed
}


dev.off()
}


#------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------

# may need to run this interactively 

generate_hapl("gbb100k")
generate_hapl("gbb250k")
generate_hapl("gbg100k")
generate_hapl("gbg250k")

#------------------------------------------------------------------------------------------------------------------------

q()

generate_hapl("gbb250k")
Error in plot.window(...) : need finite 'ylim' values
In addition: Warning messages:
1: In haplotype.DNAbin(blocs) :
  some sequences of different lengths were assigned to the same haplotype
2: In haplotype.DNAbin(blocs) :
  some sequences were not assigned to the same haplotype because of ambiguities
3: In xy.coords(x, y, xlabel, ylabel, log) : NAs introduced by coercion
4: In min(x) : no non-missing arguments to min; returning Inf
5: In max(x) : no non-missing arguments to max; returning -Inf
6: In plot.window(...) : "scale.ratio" is not a graphical parameter
7: In plot.window(...) : "pie" is not a graphical parameter
8: In plot.window(...) : "show.mutation" is not a graphical parameter

/scratch/devel/hpawar/admix/overlap.s\*.skov/hapl/tmp/gbb250k/plot/gbb250k_overlaps.hapl.pdf
#------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------


# outputs to screen of running proc_hapl function

proc_hapl(l_gbb_100k,"gbb100k")
[[1]]
GRanges object with 628 ranges and 13 metadata columns:
        seqnames            ranges strand |        V1        V2        V3
           <Rle>         <IRanges>  <Rle> | <numeric> <numeric> <numeric>
    [1]     chr1   9785000-9915000      * |         0         1         0
    [2]     chr1 14012000-14127000      * |         0         1         0
    [3]     chr1 24676000-24778000      * |         0         0         0
    [4]     chr1 31173000-31278000      * |         0         0         0
    [5]     chr1 34780000-34960000      * |         1         0         0
    ...      ...               ...    ... .       ...       ...       ...
  [624]    chr22 36832000-36978000      * |         1         0         1
  [625]    chr22 42543000-42762000      * |         0         0         0
  [626]    chr22 43251000-43404000      * |         0         0         0
  [627]    chr22 43898000-44307000      * |         0         1         0
  [628]    chr22 47337000-47536000      * |         0         0         0
               V4        V5        V6        V7        V8        V9       V10
        <numeric> <numeric> <numeric> <numeric> <numeric> <numeric> <numeric>
    [1]         0         0         0         0         0         0         0
    [2]         0         0         0         0         0         0         0
    [3]         0         0         0         0         0         0         1
    [4]         0         0         0         0         0         1         1
    [5]         0         0         0         0         0         0         0
    ...       ...       ...       ...       ...       ...       ...       ...
  [624]         0         0         0         1         0         1         0
  [625]         1         0         0         0         0         0         0
  [626]         0         0         1         0         0         0         0
  [627]         0         0         1         0         0         1         0
  [628]         0         0         1         0         0         0         0
              V11       V12  consensusIDs
        <numeric> <numeric>   <character>
    [1]         0         0   consensus_1
    [2]         0         0   consensus_2
    [3]         0         0   consensus_3
    [4]         1         0   consensus_4
    [5]         0         0   consensus_5
    ...       ...       ...           ...
  [624]         0         0 consensus_624
  [625]         0         0 consensus_625
  [626]         0         0 consensus_626
  [627]         1         0 consensus_627
  [628]         0         0 consensus_628
  -------
  seqinfo: 22 sequences from an unspecified genome; no seqlengths

[[2]]
  [1]  1  1  1  3  1  1  4  1  2  1  3  4  3  6  1  1  1  2  1  5  1  2  1  5  2
 [26]  1  2  2  2  1  1  4  1  7  1  1  2  1  1  2  1  3  1  1  3  2  1  1  1  3
 [51]  2  2  3  3  1  3  1  2  3  1  2  2  1  3  2  1  1  2  1  1  2  5  4  5  5
 [76]  2  1  1  1  3  1  1  1  4  2  1  2  4  2  5  4  2  5  2  2  2  1  1  3  3
[101]  1  3  1  4  2  1  1  4  2  2  1  2  1  4  1  3  1  1  2  2  4  4  3  2  2
[126]  1  3  2 11  7  2  1  2  1  1  5  1  6  3  3  1  5  1  1  3  4  3  6  1  5
[151]  2  2  2  3  3  3  3  5  1  1  3  2  4  1  3  2  1  1  6  1  1  1  3  1  1
[176]  4  4  2  3  1  5  1  1  1  1  4  4  1  1  3  4  3  7  1  5  1  2  3  1  1
[201]  1  1  2  1  1  3  1  1  1  1  1  5  1  8  1  1  4  4  1  2  5  2  7  1  4
[226]  2  1  1  1  5  1  1  1  1  1  2  2  1  4  1  2  2  4  1  2  2  2  1  1  3
[251]  2  4  1  1  1  2  1  2  1  1  8  3  3  1  3  2  1  7  1  1  3  3  1  1  2
[276]  2  2  3 10  4  2  4  1  1  3  2  4  4  3  2  3  1  3  5  1  1  9  1  4  1
[301]  3  2  4  1  2  1  1  1  1  3  8  1 11  1  3  1  1  4  2  1  1  3  1  2 10
[326]  2  2  2  1  4  6  1  2  6  2  1  2  1  1  1  1  1  1  1  2  1  2  3  1  5
[351]  1  1  4  1  8  1  1  1  1  6  5  1  6  2  4  3  1  1  2  3  2  1  1  1  1
[376]  4  3  3  2  4  1  5  2  5  4  1  2  4  1  2  1  1  1  2  1  6  1  8  5  2
[401]  3  1  7  2  1  1  7  1  6  1  2  2 10  1  5  1  1  1  1  5  1  2  2  2  3
[426]  1  1  1  1  4  3  1  1  1  4  1  1  5  4  1  1  2  1  6  2  1  1  4  1  1
[451]  2  1  2  4  2  1  2  1  1 10  2  1  1  2  2  3  1  7  1  1  1  1  1  8  2
[476]  2  3  4  1  1  2  1  1  5  1  1  3  1  3  5  8  3  1  7  2  4  3  2  6  3
[501]  3  5  1 10  3  1  1  1  1  1  6  1  2  1  1  2  2  2  5  1  1  1  1  1  3
[526]  4  1  3  2  2  6  2  2  1  1  1  1  4  2  2  5  1  2  2  2  2  1  5  1  1
[551]  3  1  5  1  1  4  6  6  1  9  1  1  2  7  1  2  1  6  1  2  5  2  6  1  3
[576]  4  1  3  1  6  1  4  7  9  1  1  2  2  1  6  1  1  4  8  1  6  5  1  1  2
[601]  5  1  3  5  3  5  1  1  3  1  1  3  1  1  1  1  1  2  5  4  2  2  3  4  1
[626]  1  4  1

[[3]]
[1] 628


> proc_hapl(l_gbb_250k,"gbb250k")
[[1]]
GRanges object with 103 ranges and 13 metadata columns:
        seqnames              ranges strand |        V1        V2        V3
           <Rle>           <IRanges>  <Rle> | <numeric> <numeric> <numeric>
    [1]     chr1   36173000-36436000      * |         0         0         1
    [2]     chr1   86909000-87293000      * |         0         1         0
    [3]     chr1 239693000-239993000      * |         0         0         0
    [4]     chr2   17501000-18017000      * |         0         0         1
    [5]     chr2   32220000-33064000      * |         0         0         0
    ...      ...                 ...    ... .       ...       ...       ...
   [99]    chr20   39703000-39973000      * |         0         0         0
  [100]    chr21   27185000-27489000      * |         0         0         0
  [101]    chr22   17869000-18411000      * |         0         0         1
  [102]    chr22   35786000-36046000      * |         0         0         1
  [103]    chr22   43898000-44307000      * |         0         1         0
               V4        V5        V6        V7        V8        V9       V10
        <numeric> <numeric> <numeric> <numeric> <numeric> <numeric> <numeric>
    [1]         0         0         0         0         0         0         0
    [2]         0         0         0         0         1         0         0
    [3]         0         1         0         0         0         0         0
    [4]         0         0         0         0         0         0         0
    [5]         1         1         0         0         0         0         0
    ...       ...       ...       ...       ...       ...       ...       ...
   [99]         0         1         0         1         1         0         0
  [100]         0         1         0         0         0         0         0
  [101]         0         0         0         0         0         1         0
  [102]         0         0         0         0         0         1         0
  [103]         0         0         1         0         0         1         0
              V11       V12  consensusIDs
        <numeric> <numeric>   <character>
    [1]         0         0   consensus_1
    [2]         0         0   consensus_2
    [3]         0         0   consensus_3
    [4]         0         0   consensus_4
    [5]         0         0   consensus_5
    ...       ...       ...           ...
   [99]         0         0  consensus_99
  [100]         0         0 consensus_100
  [101]         0         0 consensus_101
  [102]         1         0 consensus_102
  [103]         1         0 consensus_103
  -------
  seqinfo: 22 sequences from an unspecified genome; no seqlengths

[[2]]
  [1] 1 2 1 1 2 2 1 1 1 4 4 2 1 4 1 1 2 5 3 1 2 1 3 1 1 3 5 1 4 2 2 1 2 3 1 2 1
 [38] 2 3 3 2 1 3 2 1 1 1 2 1 2 1 4 1 2 6 2 1 1 1 2 3 1 1 1 4 1 1 3 1 1 1 2 2 1
 [75] 4 3 1 1 7 4 4 3 1 2 5 1 4 3 6 2 1 1 1 2 3 1 1 5 3 1 2 3 4

[[3]]
[1] 103

> proc_hapl(l_gbg_100k,"gbg100k")
[[1]]
GRanges object with 490 ranges and 10 metadata columns:
        seqnames            ranges strand |        V1        V2        V3
           <Rle>         <IRanges>  <Rle> | <numeric> <numeric> <numeric>
    [1]     chr1   9783000-9917000      * |         0         0         0
    [2]     chr1 25027000-25132000      * |         0         0         0
    [3]     chr1 33927000-34048000      * |         0         0         0
    [4]     chr1 36137000-36261000      * |         1         1         0
    [5]     chr1 37294000-37460000      * |         1         0         0
    ...      ...               ...    ... .       ...       ...       ...
  [486]    chr21 44167000-44358000      * |         0         0         0
  [487]    chr22 26004000-26182000      * |         0         0         0
  [488]    chr22 27581000-27693000      * |         0         0         1
  [489]    chr22 36843000-36978000      * |         0         0         0
  [490]    chr22 43906000-44304000      * |         1         1         0
               V4        V5        V6        V7        V8        V9
        <numeric> <numeric> <numeric> <numeric> <numeric> <numeric>
    [1]         0         1         0         0         0         0
    [2]         0         0         0         0         0         1
    [3]         0         0         0         0         0         1
    [4]         0         0         0         1         0         1
    [5]         0         0         0         0         0         0
    ...       ...       ...       ...       ...       ...       ...
  [486]         0         0         1         0         0         0
  [487]         0         0         1         0         0         0
  [488]         0         0         0         0         0         0
  [489]         0         0         1         0         1         0
  [490]         0         0         1         1         0         0
         consensusIDs
          <character>
    [1]   consensus_1
    [2]   consensus_2
    [3]   consensus_3
    [4]   consensus_4
    [5]   consensus_5
    ...           ...
  [486] consensus_486
  [487] consensus_487
  [488] consensus_488
  [489] consensus_489
  [490] consensus_490
  -------
  seqinfo: 22 sequences from an unspecified genome; no seqlengths

[[2]]
  [1] 1 1 1 4 1 2 1 1 1 1 2 8 2 2 2 2 1 2 1 1 2 1 1 3 3 1 1 4 1 1 1 4 1 1 1 1 1
 [38] 1 1 1 1 1 1 1 3 2 1 1 1 3 2 2 1 3 1 3 1 1 4 8 1 1 2 1 1 1 2 1 1 1 5 1 1 1
 [75] 1 1 1 1 2 1 2 1 1 3 1 1 2 2 2 1 2 2 2 2 1 1 4 2 9 1 1 6 1 1 1 2 2 1 1 1 1
[112] 1 1 4 1 1 1 1 1 1 1 1 2 1 7 2 2 1 1 1 1 4 1 2 1 2 4 2 6 3 3 1 2 3 2 4 1 6
[149] 1 1 1 2 2 2 2 1 1 1 5 3 1 1 1 1 1 1 3 1 2 1 1 2 1 1 1 3 1 1 1 1 1 1 1 1 1
[186] 1 1 1 1 3 1 1 1 1 1 2 3 5 3 6 1 2 1 1 3 1 2 1 1 3 1 1 1 2 1 4 9 1 5 1 1 7
[223] 1 1 1 1 1 3 2 6 1 1 4 2 3 2 3 1 1 2 4 2 3 1 1 2 2 1 5 1 1 8 1 1 9 1 1 1 1
[260] 2 1 5 1 2 1 3 5 2 1 1 2 4 1 6 1 5 3 1 3 1 1 1 2 1 7 2 6 1 2 1 3 3 2 1 3 1
[297] 1 2 3 2 1 3 2 3 1 4 2 2 4 2 4 2 1 1 1 4 4 1 1 1 2 4 1 5 2 1 1 9 1 7 1 5 1
[334] 4 1 4 2 1 3 8 1 1 5 1 1 4 1 1 4 2 3 2 2 2 3 1 2 1 1 1 2 6 2 2 5 1 3 2 1 1
[371] 1 1 1 1 1 2 1 1 2 3 2 1 3 1 1 4 8 1 2 1 1 3 1 1 2 5 1 1 1 1 1 8 5 1 2 3 3
[408] 1 2 2 2 2 2 1 1 1 1 2 1 2 5 3 1 3 2 2 3 8 1 1 1 3 2 2 2 3 1 1 5 2 1 2 1 4
[445] 1 1 1 1 2 1 1 2 1 1 2 1 2 2 1 1 1 1 5 1 2 3 3 7 2 2 5 2 2 2 1 2 2 2 3 5 4
[482] 1 4 1 2 1 1 1 2 4

[[3]]
[1] 490

> proc_hapl(l_gbg_250k,"gbg250k")
[[1]]
GRanges object with 91 ranges and 10 metadata columns:
       seqnames              ranges strand |        V1        V2        V3
          <Rle>           <IRanges>  <Rle> | <numeric> <numeric> <numeric>
   [1]     chr1   48152000-48407000      * |         0         0         0
   [2]     chr1   77372000-77714000      * |         1         0         0
   [3]     chr1   86933000-87380000      * |         0         0         1
   [4]     chr1   88060000-88364000      * |         0         1         0
   [5]     chr1 184558000-184939000      * |         0         0         0
   ...      ...                 ...    ... .       ...       ...       ...
  [87]    chr18   54697000-55007000      * |         0         0         0
  [88]    chr19   31880000-32156000      * |         0         0         0
  [89]    chr20   25191000-25565000      * |         0         1         0
  [90]    chr20   37384000-37674000      * |         0         0         0
  [91]    chr22   43906000-44304000      * |         1         1         0
              V4        V5        V6        V7        V8        V9 consensusIDs
       <numeric> <numeric> <numeric> <numeric> <numeric> <numeric>  <character>
   [1]         0         0         0         1         0         1  consensus_1
   [2]         0         0         0         0         0         0  consensus_2
   [3]         0         0         0         0         0         1  consensus_3
   [4]         0         0         0         0         0         0  consensus_4
   [5]         0         0         0         0         0         1  consensus_5
   ...       ...       ...       ...       ...       ...       ...          ...
  [87]         1         0         0         0         0         0 consensus_87
  [88]         1         0         0         1         0         0 consensus_88
  [89]         1         0         1         0         0         0 consensus_89
  [90]         0         0         0         0         0         1 consensus_90
  [91]         0         0         1         1         0         0 consensus_91
  -------
  seqinfo: 22 sequences from an unspecified genome; no seqlengths

[[2]]
 [1] 2 1 2 1 1 1 2 3 6 1 1 1 2 7 1 1 1 1 2 1 1 1 1 2 4 1 2 1 1 1 1 1 1 1 1 1 1 1
[39] 2 1 6 1 3 1 1 2 3 3 4 1 1 2 1 1 1 2 2 1 2 1 4 1 2 4 2 1 4 1 3 1 1 1 2 3 2 1
[77] 1 1 1 2 1 1 4 1 1 1 1 2 3 1 4

[[3]]
[1] 91
