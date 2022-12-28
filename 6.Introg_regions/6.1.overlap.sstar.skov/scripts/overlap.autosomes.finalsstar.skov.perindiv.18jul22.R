# Mon 18 Jul 2022 09:54:24 CEST
# overlap S* windows of 99% CI with strict skov hmm autosomal outliers - calculate per eastern individual

# following -   /Volumes/"Ultra USB 3.0"/IBE/further.analysis.feb.2020/gorillas/abc/simul_postabc/overlap.autosomes.sstar.skov.perindiv.3may22.R

#-----------------------------------------------------------------------------------------------------------------------
## module load gcc/6.3.0 openssl/1.0.2q R/4.0.1
require(data.table)
options(stringsAsFactors=F)
options(scipen=100)
library(mgcv)
library(GenomicRanges)
library(tidyr)
library(ggbio)
library(pheatmap)
library(ggplot2)
library(ggpubr)

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
# 2) read in skov strict outliers & using function below split per inidivdual & generate reduced granges objects
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
#-----------------------------------------------------------------------------------------------------------------------

# calculate proportion of bp overlapping in intersect

# Wed 23 Feb 2022 17:16:30 CET
# MK
#could it be that you are underestimating the overlap? 
#In the end, many regions are observed in more than one individual, but you count the total number/width across, right? 
#But then, you take the reduce() function for the overlap. 
#Using individual-based estimates will help with that, obviously, but for the overall estimate, you may just add
#fG<-reduce(fG);fB<-reduce(fB)
#at the beginning of the calc_prof_fun function.

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
return(list(ranges_2,lengths_2,prop_2))
}

#-----------------------------------------------------------------------------------------------------------------------

# perform intersect b/n s* & skov outliers for all ids
intersect_function<-function(sstar_subs,skov_subs) {
# 1) compare GBB - 95% CI for S*
#sstar_subs #sk_regions_gb
out_gbb_95<-list()
for (ind in (1:length(sstar_subs))) {
out_gbb_95[[ind]]<-process_tooverlap_fun(sstar_subs[[ind]],skov_subs[[ind]])
}

# for the proportion of overlapping bp 
# extract first rows of 3rd df
 #out_gbb_95[[1]][[3]][1,]
 # prop.bp.x prop.frag.x 
 # 0.6345073   0.4444444 

prop_overlap_gbb_95<-list()
for (ind in (1:length(sstar_subs))) {
prop_overlap_gbb_95[[ind]]<-out_gbb_95[[ind]][[3]][1,]
}

prop_df<-do.call(rbind, prop_overlap_gbb_95)

# length dists in kb
 #(out_gbb_95[[1]][[1]]@ranges@width)/1000
length_overlap_gbb_95<-list()
for (ind in (1:length(sstar_subs))) {
length_overlap_gbb_95[[ind]]<-out_gbb_95[[ind]][[2]]
}

length_df<- do.call(rbind,length_overlap_gbb_95)

return(list(out_gbb_95,prop_df,length_df))
}
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------


intersect_function(sstar_GBB_99,autosomes_sk_gbb)[[2]]
intersect_function(sstar_GBG_99,autosomes_sk_gbg)[[2]]

intersect_function(sstar_GBB_99,autosomes_40sk_gbb)[[2]]
intersect_function(sstar_GBG_99,autosomes_40sk_gbg)[[2]]

#-----------------------------------------------------------------------------------------------------------------------


#> intersect_function(sstar_GBB_99,autosomes_sk_gbb)[[2]]
#      prop.bp.x prop.frag.x
# [1,] 0.5303163   0.3437037
# [2,] 0.5305414   0.3776316
# [3,] 0.5216871   0.3595206
# [4,] 0.5293359   0.3304965
# [5,] 0.4862865   0.3093525
# [6,] 0.5119994   0.3455979
# [7,] 0.5213613   0.3430079
# [8,] 0.5053765   0.3269476
# [9,] 0.5330628   0.3602151
#[10,] 0.4462756   0.2991304
#[11,] 0.4958976   0.3514774
#[12,] 0.5308228   0.3592233

as.data.frame(intersect_function(sstar_GBB_99,autosomes_sk_gbb)[[2]][,1])
as.data.frame(intersect_function(sstar_GBG_99,autosomes_sk_gbg)[[2]][,1])
as.data.frame(intersect_function(sstar_GBB_99,autosomes_40sk_gbb)[[2]][,1])
as.data.frame(intersect_function(sstar_GBG_99,autosomes_40sk_gbg)[[2]][,1])



> as.data.frame(intersect_function(sstar_GBB_99,autosomes_sk_gbb)[[2]][,1])
   intersect_function(sstar_GBB_99, autosomes_sk_gbb)[[2]][, 1]
1                                                     0.5303163
2                                                     0.5305414
3                                                     0.5216871
4                                                     0.5293359
5                                                     0.4862865
6                                                     0.5119994
7                                                     0.5213613
8                                                     0.5053765
9                                                     0.5330628
10                                                    0.4462756
11                                                    0.4958976
12                                                    0.5308228
> as.data.frame(intersect_function(sstar_GBG_99,autosomes_sk_gbg)[[2]][,1])
  intersect_function(sstar_GBG_99, autosomes_sk_gbg)[[2]][, 1]
1                                                    0.4112736
2                                                    0.4290258
3                                                    0.4247600
4                                                    0.4286217
5                                                    0.4679428
6                                                    0.4091474
7                                                    0.4255768
8                                                    0.3983371
9                                                    0.3793612
> as.data.frame(intersect_function(sstar_GBB_99,autosomes_40sk_gbb)[[2]][,1])
   intersect_function(sstar_GBB_99, autosomes_40sk_gbb)[[2]][, 1]
1                                                       0.5714986
2                                                       0.5625861
3                                                       0.5504886
4                                                       0.5694366
5                                                       0.5182644
6                                                       0.5429275
7                                                       0.5629811
8                                                       0.5396514
9                                                       0.5648412
10                                                      0.4775200
11                                                      0.5278887
12                                                      0.5660033
> as.data.frame(intersect_function(sstar_GBG_99,autosomes_40sk_gbg)[[2]][,1])
  intersect_function(sstar_GBG_99, autosomes_40sk_gbg)[[2]][, 1]
1                                                      0.4293668
2                                                      0.4545590
3                                                      0.4458455
4                                                      0.4532607
5                                                      0.4957247
6                                                      0.4384527
7                                                      0.4560856
8                                                      0.4283009
9                                                      0.3995148


# perhaps for comparison - could also calculate for the other ci intervals of S*

#-----------------------------------------------------------------------------------------------------------------------

# GBG, CI 1 (95%)
sstar_GBG_95<-proc_proportion(1,1)
# GBB, CI 1 (95%)
sstar_GBB_95<-proc_proportion(2,1)

# GBG, CI 2 (99.5%)
sstar_GBG_99.5<-proc_proportion(1,3)
# GBB, CI 2 (99%)
sstar_GBB_99.5<-proc_proportion(2,3)


scen_GBB<-list(sstar_GBB_95,sstar_GBB_99,sstar_GBB_99.5)
scen_GBG<-list(sstar_GBG_95,sstar_GBG_99,sstar_GBG_99.5)


re_calc_overlap<-function(scen,skov_version){
overlap40kb_gbb<-list()
for (ind in (1:length(scen))) {
overlap40kb_gbb[[ind]]<-intersect_function(scen[[ind]],skov_version)[[2]]}
return(overlap40kb_gbb)
}


oprop_gbb<-re_calc_overlap(scen_GBB,autosomes_sk_gbb)
oprop_gbg<-re_calc_overlap(scen_GBG,autosomes_sk_gbg)

oprop_40gbb<-re_calc_overlap(scen_GBB,autosomes_40sk_gbb)
oprop_40gbg<-re_calc_overlap(scen_GBG,autosomes_40sk_gbg)




#hold_outputs<-list(oprop_gbb,oprop_gbg,oprop_40gbb,oprop_40gbg)

#proc_out_overlap<-function(){
#out_f<-list()
#for (ind in (1:length(hold_outputs))) {
#out_f[[ind]]<-cbind(as.data.frame(hold_outputs[[ind]][[1]][,1]),as.data.frame(hold_outputs[[ind]][[2]][,1]),as.data.frame(hold_outputs[[ind]][[3]][,1]))}
#return(out_f)
#}

#proc_out_overlap()

#-----------------------------------------------------------------------------------------------------------------------

# results added to ABC-Gor excel

# proportion of overlapping bp b/n S* - skov fragments for GBB at each S* CI (for each individual)
> cbind(as.data.frame(oprop_gbb[[1]][,1]),as.data.frame(oprop_gbb[[2]][,1]),as.data.frame(oprop_gbb[[3]][,1]))
   oprop_gbb[[1]][, 1] oprop_gbb[[2]][, 1] oprop_gbb[[3]][, 1]
1            0.6558773           0.5303163           0.4416788
2            0.6609680           0.5305414           0.4248159
3            0.6412267           0.5216871           0.4242316
4            0.6633056           0.5293359           0.4365472
5            0.6184662           0.4862865           0.4174412
6            0.6579004           0.5119994           0.4382001
7            0.6587137           0.5213613           0.4167493
8            0.6264497           0.5053765           0.4116871
9            0.6563864           0.5330628           0.4446639
10           0.5902518           0.4462756           0.3567607
11           0.6309916           0.4958976           0.4179757
12           0.6413670           0.5308228           0.4255281


# proportion of overlapping bp b/n S* - skov fragments for GBG at each S* CI
> cbind(as.data.frame(oprop_gbg[[1]][,1]),as.data.frame(oprop_gbg[[2]][,1]),as.data.frame(oprop_gbg[[3]][,1]))
  oprop_gbg[[1]][, 1] oprop_gbg[[2]][, 1] oprop_gbg[[3]][, 1]
1           0.5523482           0.4112736           0.3299910
2           0.5515809           0.4290258           0.3383781
3           0.5689909           0.4247600           0.3498864
4           0.5639518           0.4286217           0.3633259
5           0.6026423           0.4679428           0.3558730
6           0.5664158           0.4091474           0.3219495
7           0.5761456           0.4255768           0.3397743
8           0.5248341           0.3983371           0.3369652
9           0.5078236           0.3793612           0.3052611



# proportion of overlapping bp b/n S* - skov fragments >= 40kb for GBB at each S* CI
>  cbind(as.data.frame(oprop_40gbb[[1]][,1]),as.data.frame(oprop_40gbb[[2]][,1]),as.data.frame(oprop_40gbb[[3]][,1]))
   oprop_40gbb[[1]][, 1] oprop_40gbb[[2]][, 1] oprop_40gbb[[3]][, 1]
1              0.7018215             0.5714986             0.4765317
2              0.6984304             0.5625861             0.4505572
3              0.6746059             0.5504886             0.4507133
4              0.7063196             0.5694366             0.4699466
5              0.6522715             0.5182644             0.4468485
6              0.6925561             0.5429275             0.4677435
7              0.7065064             0.5629811             0.4497369
8              0.6641537             0.5396514             0.4402349
9              0.6917164             0.5648412             0.4726999
10             0.6291400             0.4775200             0.3844159
11             0.6669144             0.5278887             0.4454471
12             0.6790633             0.5660033             0.4556758


# proportion of overlapping bp b/n S* - skov fragments >= 40kb for GBG at each S* CI
>  cbind(as.data.frame(oprop_40gbg[[1]][,1]),as.data.frame(oprop_40gbg[[2]][,1]),as.data.frame(oprop_40gbg[[3]][,1]))
  oprop_40gbg[[1]][, 1] oprop_40gbg[[2]][, 1] oprop_40gbg[[3]][, 1]
1             0.5744700             0.4293668             0.3450991
2             0.5799716             0.4545590             0.3589075
3             0.5961032             0.4458455             0.3682384
4             0.5959788             0.4532607             0.3852635
5             0.6348049             0.4957247             0.3789234
6             0.6056224             0.4384527             0.3449453
7             0.6152972             0.4560856             0.3643673
8             0.5621184             0.4283009             0.3643832
9             0.5339404             0.3995148             0.3220402
