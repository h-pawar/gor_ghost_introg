#!/usr/bin/r

# Mon 18 Jul 2022 14:50:02 CEST
# calculate gene density for the intersect regions (intersect b/n S* 99% CI & >=40kb strict Skov)

#-----------------------------------------------------------------------------------------------------------------------
# definition of gene density from kuhlwilm 2019 suppl -
# 'genic content, i.e. the fraction of sites that are protein- coding within archaic windows'

#-----------------------------------------------------------------------------------------------------------------------

# amending below from /scratch/devel/hpawar/admix/sstar/scripts/simul_postabc/downstream/gene.density.skov.12may22.R
# called by /scratch/devel/hpawar/admix/sstar/scripts/simul_postabc/downstream/gene.density.skov.12may22.arr
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

## module load gcc/6.3.0 openssl/1.0.2q R/4.0.1
require(data.table)
options(stringsAsFactors=F)
options(scipen=100)
library(mgcv)
library(GenomicRanges)
library(tidyr)
library(ggbio)
#library(pheatmap)
#library(ggplot2)
#library(ggpubr)
library(valr)
library(GenomicFeatures)


#-----------------------------------------------------------------------------------------------------------------------
# glm generated for the S* CIs calc from simulated data (null simulations - gor demog, no ghost)
load(file=paste("/scratch/devel/hpawar/admix/sstar/simul_postabc/stepwisesegs/final_8apr22/check_22apr22/glm",sep=""),verbose=T)

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
# GBG, CI 1 (99%)
sstar_GBG_99<-proc_proportion(1,2)

# GBB, CI 1 (99%)
sstar_GBB_99<-proc_proportion(2,2)
#-----------------------------------------------------------------------------------------------------------------------

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


#intersect_function(sstar_GBB_99,autosomes_sk_gbb)
#[12]]
#GRanges object with 259 ranges and 0 metadata columns:
#        seqnames              ranges strand
#           <Rle>           <IRanges>  <Rle>
#    [1]     chr1   31236000-31276000      *
#    [2]     chr1   87184000-87293000      *
#    [3]     chr1 114713000-114827000      *
#    [4]     chr1 119334000-119480000      *
#    [5]     chr1 146595000-146646000      *
#    ...      ...                 ...    ...
#  [255]    chr20   31681000-31911000      *
#  [256]    chr20   36920000-37083000      *
#  [257]    chr20   51864000-51928000      *
#  [258]    chr20   59987000-60137000      *
#  [259]    chr22   23204000-23314000      *
#  -------
#  seqinfo: 22 sequences from an unspecified genome; no seqlengths


#-----------------------------------------------------------------------------------------------------------------------

#intersect(autosome.hum.genes, sk_regions_gbb[[1]], ignore.strand=TRUE)
#GRanges object with 46 ranges and 0 metadata columns:
#       seqnames              ranges strand
#          <Rle>           <IRanges>  <Rle>
#   [1]     chr1   29815000-29823405      *
#   [2]     chr1   94490000-94536000      *
#   [3]     chr1 159808000-159825137      *
#   [4]     chr1 180966000-180992047      *
#   [5]     chr2   11696000-11744000      *
#   ...      ...                 ...    ...
#  [42]    chr18   63107000-63114748      *
#  [43]    chr20   56532182-56534716      *
#  [44]    chr22   18124172-18124285      *
#  [45]    chr22   36956917-36958962      *
#  [46]    chr22   46877226-46877657      *
#  -------
#  seqinfo: 265 sequences from an unspecified genome; no seqlengths

#  GenomicRanges::intersect(autosome.hum.genes, sk_regions_gbb[[1]], ignore.strand=TRUE) # is equivalent to the above
 

#sum(reduce(GenomicRanges::intersect(autosome.hum.genes, sk_regions_gbb[[1]], ignore.strand=TRUE))@ranges@width)
#[1] 1088576

#sum((GenomicRanges::intersect(autosome.hum.genes, sk_regions_gbb[[1]], ignore.strand=TRUE))@ranges@width) # is equivalent with or without the reduce()
#[1] 1088576

# sum(sk_regions_gbb[[1]]@ranges@width)
#[1] 54356675

#>  1088576/54356675 # realistic?
#[1] 0.02002654

#pintersect(autosome.hum.genes, sk_regions_gbb[[1]], ignore.strand=TRUE)
#Error in .pintersect_GRanges_GRanges(x, y, drop.nohit.ranges = drop.nohit.ranges,  : 
#  'y' must have the length of 'x' or length 1

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

# (2) read in gtf of gene coordinates from hg19 reference (protein coding bp)
	# initially just use the human annotation (rather than filtering by orthologs)

# following - https://www.biostars.org/p/169171/ to convert gtf file -> granges object

# need 2 objects : 1 with coordinates of genes, 2nd object of outlier regions

# generate object- with coordinates of genes

# use makeTxDbFromGFF to create the TxDB object 
gtf="/scratch/project/devel/mkuhlwilm/chimpdiv/index/Homo_sapiens.GRCh37.75_c.gtf"

test<-makeTxDbFromGFF(gtf, format='gtf')
#Import genomic features from the file as a GRanges object ... 
#OK
#Prepare the 'metadata' data frame ... OK
#Make the TxDb object ... 
#OK
#Warning message:
#In .get_cds_IDX(mcols0$type, mcols0$phase) :
#  The "phase" metadata column contains non-NA values for features of type
#  stop_codon. This information was ignored.

#str(test)
#Reference class 'TxDb' [package "GenomicFeatures"] with 5 fields

human.genes = genes(test)

# only retain autosomal chr
autosome.hum.genes<-human.genes[seqnames(human.genes) == c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22")]
#Warning message:","chr17","chr18","chr19","chr20","chr21","chr22")]
#In e1 == Rle(e2) :
#  longer object length is not a multiple of shorter object length

# data structure - contains gene coor
#str(human.genes)
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

# (3) Intersect protein coding bp (autosome.hum.genes) with the putative introg regions - empirical data (skov, s*)
# want to perform the intersect & output these ranges as well as the counts

intersect_withgenes_function<-function(empirical_id) {
intersect_id<-list()
counts_id<-list()
for (ind in (1:length(empirical_id))) {
intersect_id[[ind]]<-intersect(autosome.hum.genes, empirical_id[[ind]], ignore.strand=TRUE)
counts_id[[ind]]<-cbind(sum(reduce(intersect_id[[ind]])@ranges@width), sum(empirical_id[[ind]]@ranges@width))
}
return(list(intersect_id,counts_id))
}

#test<-intersect_withgenes_function(sk_regions_gbb)
# outputs : list 1 = the intersect ranges (protein coding bp)
#[[1]][[12]]
#GRanges object with 62 ranges and 0 metadata columns:
#       seqnames              ranges strand
#          <Rle>           <IRanges>  <Rle>
#   [1]     chr1   29815000-29823405      *
# outputs : list 2: the counts of protein coding bp, the total bp of the target regions
#[2]]
#[[2]][[1]]
#        [,1]     [,2]
#[1,] 1088576 54356675

#[[2]][[2]]
#        [,1]     [,2]
#[1,] 1812576 65301760

#scen_sk_gbb<-list(sk_regions_gbb,sk_40kbregions_gbb)
#scen_sk_gbg<-list(sk_regions_gbg,sk_40kbregions_gbg)

# changing scen from skov regions only (to the overlapping regions)

scen_gbg<-list(ov_gbg_99,ov_gbg40_99)
scen_gbb<-list(ov_gbb_99,ov_gbb40_99)


overall_calc_overlap<-function(scen){
overall_overlap<-list()
overall_df<-list()
for (ind in (1:length(scen))) {
overall_overlap[[ind]]<-intersect_withgenes_function(scen[[ind]])
overall_df[[ind]]<-do.call(rbind, overall_overlap[[ind]][[2]])}
return(overall_df)
}

prot_gbb<-overall_calc_overlap(scen_gbb)
prot_gbg<-overall_calc_overlap(scen_gbg)

empirical_pointestimate<-function(scen){
out_pt<-list()
for (ind in (1:length(scen))) {
out_pt[[ind]]<-cbind(mean(scen[[ind]][,1]/scen[[ind]][,2]), sd(scen[[ind]][,1]/scen[[ind]][,2]))}
return(out_pt)
}


emp_pt_gbb<-empirical_pointestimate(prot_gbb)
emp_pt_gbg<-empirical_pointestimate(prot_gbg)



#-----------------------------------------------------------------------------------------------------------------------

q()

# BELOW HAS BEEN SUPERCEDED BY  /scratch/devel/hpawar/admix/sstar/scripts/simul_postabc/downstream/random.gene.density.s*skov.overlap.21oct22.R 
# WAS REWRITTEN TO GENERATE RANDOM REGIONS OF SUFFICIENT CALLABILITY FOR GENE DENSITY ANALYSIS

# (4) generate random regions of equal length distribution as the empirical data for each individual

#-----------------------------------------------------------------------------------------------------------------------
# lengths of hg19
hg19<-fread('/home/devel/marcmont/scratch/Harvi/geneflow/test/test.reconst.ids/ora_1_1/proportion/hg19.autosomes.chrom.sizes.txt')

hg19.lengths<-apply(hg19[,2],2,function(x) as.numeric(as.character(x)))

# hg19 chr in the correct order
hg19<-hg19[c(1:18,20,19,22,21),]
colnames(hg19)<-c("chrom","size") 
# or whether need to read in the genome itself?
  # path to hg19 ref


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

# write as a function as well - done
generate_random<-function(fG) {
random_reg_x<-list()
for (ind in (1:length(fG))) {
random_reg_x[[ind]]<-test_proc_randomregions(fG,ind) }
return(random_reg_x)
}

#-----------------------------------------------------------------------------------------------------------------------

#sk_tes_gbb<-generate_random(sk_regions_gbb) # generates for 1 rep (random regions for each id)
# check if intersect_withgenes_function works with random regions ***
#intersect_withgenes_function(sk_tes_gbb)

#do.call(rbind, intersect_withgenes_function(sk_tes_gbb)[[2]])
#do.call(rbind, intersect_withgenes_function(sk_tes_gbb)[[2]])
#         [,1]     [,2]
# [1,] 1318844 54356675
# [2,] 1589943 65301760
# [3,] 1947838 65794751

 # add another column of the divisions

# x<-do.call(rbind, intersect_withgenes_function(sk_tes_gbb)[[2]])
# y<-x[,1]/x[,2]

# testing - fine
#fB=sk_regions_gbb
#randomreg_100iter<-list()
#x<-list()
#y<-list()
#for (i in (1:10)) {
#randomreg_100iter[[i]]<-generate_random(fB)
#x[[i]]<-do.call(rbind, intersect_withgenes_function(randomreg_100iter[[i]])[[2]])
#y[[i]]<-x[[i]][,1]/x[[i]][,2]
#}

#> y
#[[1]]
# [1] 0.03216963 0.02920228 0.02506177 0.03087500 0.01879587 0.02578618
# [7] 0.03570058 0.02993782 0.03537771 0.03625615 0.03602203 0.02433592
# works fine

#-----------------------------------------------------------------------------------------------------------------------

# only returning the proportions (protein coding bp / total bp for all random regions per individual)
random_proteincoding<-function(fB){
randomreg_100iter<-list()
x<-list()
y<-list()
for (i in (1:100)) {
randomreg_100iter[[i]]<-generate_random(fB)
x[[i]]<-do.call(rbind, intersect_withgenes_function(randomreg_100iter[[i]])[[2]])
y[[i]]<-x[[i]][,1]/x[[i]][,2]
}
return(y)
}


#scen_sk_gbb<-list(sk_regions_gbb,sk_40kbregions_gbb)
#scen_sk_gbg<-list(sk_regions_gbg,sk_40kbregions_gbg)

intersect_allrandom<-function(scen){
int_random<-list()
for (ind in (1:length(scen))) {
int_random[[ind]]<-random_proteincoding(scen[[ind]])}
return(int_random)
}

ran_gbb<-intersect_allrandom(scen_gbb)
ran_gbg<-intersect_allrandom(scen_gbg)

# output these in order to plot 
save(ran_gbb,file=paste("/scratch/devel/hpawar/admix/overlap.s*.skov/random_s*_sk_overlap_gbb_gene.density"))
save(ran_gbg,file=paste("/scratch/devel/hpawar/admix/overlap.s*.skov/random_s*_sk_overlap_gbg_gene.density"))

q()

#-----------------------------------------------------------------------------------------------------------------------

# (5) calculate p-val

#p_val=sum(abs(meansdf_random_gbg_95[,1]) >= abs(original))/(n)

# is the n here the number of reps or number of reps * number of individuals?

# testing p-val
#-----------------------------------------------------------------------------------------------------------------------
