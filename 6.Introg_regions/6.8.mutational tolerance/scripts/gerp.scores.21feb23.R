
# Tue 21 Feb 2023 12:13:04 CET
## module load gcc/6.3.0 openssl/1.0.2q R/4.0.1  BCFTOOLS/1.12
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


# have already subset gerp score files into high & low impact files (step 1)
# step 2 - filter gerp mutations of each category by those polymorphic in eastern gorillas
# step 3 - overlap gerp mutations of each category with introgressed regions

#-----------------------------------------------------------------------------------------------------------------------

# 1) subset gerp score files by high & low impact files  - this section has already been run & gerp scores subset into high/low scores per chrom

#proc_fun<-function(a){

#gerp22<-read.table(paste("/scratch/devel/shan/gerp/chr",a,".scores",sep=""))

# filter by high impact sites
#gerp_hi<-gerp22[which(gerp22$V4>4),] # high impact gerp sites
#colnames(gerp_hi)<-NULL
#write.table(gerp_hi,file=paste("/scratch/devel/hpawar/admix/overlap.s*.skov/revisions_8feb23/gerp/chr",a,".high.scores",sep=""),sep="\t",row.names=F,col.names=F,quote=F)


#-2<gerp<2
#gerp_low<-gerp22[which(-2<gerp22$V4 & gerp22$V4<2),] # low impact gerp sites
#colnames(gerp_low)<-NULL

#write.table(gerp_low,file=paste("/scratch/devel/hpawar/admix/overlap.s*.skov/revisions_8feb23/gerp/chr",a,".low.scores",sep=""),sep="\t",row.names=F,col.names=F,quote=F)

#}


#for (ind in (1:22)) {
#print(ind)
#proc_fun(ind)
#}

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

# this section only needs to be run once *

# 2) filter gerp files by mutations segregating in easterns
    # ie remove sites fixed alternate across all gorillas (mutations arising in gorilla lineage after divergence from human ref)
    # & remove sites fixed reference across all easterns 



process_muts<-function(chrom,typ){


tmp=paste("/scratch/devel/hpawar/admix/overlap.s*.skov/revisions_8feb23/gerp/chr",chrom,".",typ,".scores",sep="")

# filter by tmp = score file rather than by chrom
# ie -> extract gerp muts of interest which are polymorphic in E & intersect this with the introg regions
testsnps<-system(paste("bcftools view -m2 -M2 -v snps -R ",tmp," /scratch/devel/mkuhlwilm/arch/subsets/3_gorilla_",chrom,".vcf.gz | bcftools query -f '%CHROM %POS [%GT ]\n' ",sep=""),intern=T) 
testsnps<-do.call(rbind,strsplit(testsnps,split=" "))  ## ie list merged vertically into df


#2.1) filter out those where all gorillas are 1/1
# fixed differences since the split from humans are not informative for our purpose, we are interested in sites segregating within gorillas

testsnps[which(testsnps=="0|0")]<-0
testsnps[which(testsnps=="1|1")]<-2
testsnps[which(testsnps=="0|1")]<-1

# check number of samples
y<-as.data.frame(testsnps[,c(3:51)],ncol=49)
y[,(1:49)] <- sapply(y[,(1:49)],as.numeric)

counts<-rowSums(y)

# filter these out: which(counts==(2*49))
# which(counts==(2*49)) # rows where all easterns are ref alternate - mutations fi>xed in gorilla lineage since divergence from human ref


#2.2) filter where all easterns are 0/0

# ie extract equivalent cols -> keep only segregating sites

e_y<-y[,(1:21)]
ecounts<-rowSums(e_y)

# which(ecounts==0)  # filter these out

# combine the 2 conditions: filter out all gor = 1/1 & all E = 0/0 

keep_snps<-testsnps[(which(counts!=(2*49) & ecounts!=0)),]

#-----------------------------------------------------------------------------------------------------------------------

# convert these polymorphic in E snps with x gerp scores to granges objects
keep_ranges<-GRanges(seqnames=keep_snps[,1],ranges=IRanges(start=as.numeric(keep_snps[,2])-1,end=as.numeric(keep_snps[,2]),names=keep_snps[,1]),strand=rep("*",length(keep_snps[,1])))


# then annotate these snps with their corresponding gerp scores

gerp_hi<-read.table(tmp)

gerp_ranges<-GRanges(seqnames=gerp_hi[,1],ranges=IRanges(start=as.numeric(gerp_hi[,2]),end=as.numeric(gerp_hi[,3]),names=gerp_hi[,1]),strand=rep("*",length(gerp_hi[,1])), gscore=(gerp_hi[,4]))

# perform overlap b/n the polymorphic snps in E & the gerp scores

y<- findOverlaps(gerp_ranges,keep_ranges)
 
y1<-unique(as.data.frame(y)[,1])


#gerp_ranges[y1] # return this object - snps annotated with gerp scores


return(gerp_ranges[y1])

}



#-----------------------------------------------------------------------------------------------------------------------



high_polym_gerp_scores<-list()
for (ind in (1:22)) {
print(ind)
high_polym_gerp_scores[[ind]]<-process_muts(ind,'high')
}


save(high_polym_gerp_scores,file=paste("/scratch/devel/hpawar/admix/overlap.s*.skov/revisions_8feb23/gerp/high_polym_gerp_scores"))

# then output this object

low_polym_gerp_scores<-list()
for (ind in (1:22)) {
print(ind)
low_polym_gerp_scores[[ind]]<-process_muts(ind,'low')
}

save(low_polym_gerp_scores,file=paste("/scratch/devel/hpawar/admix/overlap.s*.skov/revisions_8feb23/gerp/low_polym_gerp_scores"))

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------


# 3) intersect putative introgressed regions with the high/low gerp objects


# read in high gerp ranges
load(file=paste("/scratch/devel/hpawar/admix/overlap.s*.skov/revisions_8feb23/gerp/high_polym_gerp_scores"),verbose=T)

# read in low gerp ranges
load(file=paste("/scratch/devel/hpawar/admix/overlap.s*.skov/revisions_8feb23/gerp/low_polym_gerp_scores"),verbose=T)


# read in introgressed regions
load(file="/scratch/devel/hpawar/admix/overlap.s*.skov/overlap_objects/ov_gbb_99",verbose=T)
load(file="/scratch/devel/hpawar/admix/overlap.s*.skov/overlap_objects/ov_gbg_99",verbose=T)


high_polym_gerp_scoresList<-GRangesList(high_polym_gerp_scores) 
# -> do not reduce the object, to keep the individual snps & their corresponding gerp scores
all_high_polym_gerp_scores<-(unlist(high_polym_gerp_scoresList))
#GRanges object with 2637 ranges and 1 metadata column:

low_polym_gerp_scoresList<-GRangesList(low_polym_gerp_scores) 
# -> do not reduce the object, to keep the individual snps & their corresponding gerp scores
all_low_polym_gerp_scores<-(unlist(low_polym_gerp_scoresList))
#GRanges object with 139182 ranges and 1 metadata column:


# 1) gerp bp in introgressed regions

# how many mutations in a certain GERP category are observed in the introg regions
    # consider high (>4) & low (-2<x<2) impact sites
        # these gerp mutations are polymorphic in eastern gorillas 

# calculate the number of gerp bp in introgressed regions: sum of gerp mutations by introgressed region length # same protocol as for genic bp & regulatory bp

#-----------------------------------------------------------------------------------------------------------------------
gerp_bp_perind<-function(scen,gerp_rang){

# gives hits per ind - ie how many snps in the introg regions per individual
    # how many high impact gerp mutations which are polymorphic in E gorillas
hold<-list()
for (ind in (1:length(scen))) {
hold[[ind]]<-findOverlaps(gerp_rang,scen[[ind]])
}

# then take hits per ind / width of introg regions per ind : ie high gerp bp in introgressed regions 
count_per_region<-list()
for (ind in (1:length(scen))) {
# snp counts by introgressed region length
count_per_region[[ind]]<-length(hold[[ind]])/sum(width(scen[[ind]])-1)
}

# pop-level estimate 
out<- cbind(mean( unlist(count_per_region)),sd( unlist(count_per_region)))

return(out)
}


mg_highgerp<-gerp_bp_perind(ov_gbb_99,all_high_polym_gerp_scores)
el_highgerp<-gerp_bp_perind(ov_gbg_99,all_high_polym_gerp_scores)



mg_lowgerp<-gerp_bp_perind(ov_gbb_99,all_low_polym_gerp_scores)
el_lowgerp<-gerp_bp_perind(ov_gbg_99,all_low_polym_gerp_scores)



> mg_highgerp
               [,1]            [,2]
[1,] 0.000001279447 0.0000002159444
> el_highgerp
               [,1]            [,2]
[1,] 0.000001629846 0.0000002650428

> mg_lowgerp
              [,1]           [,2]
[1,] 0.00007010927 0.000001576169
> el_lowgerp
              [,1]           [,2]
[1,] 0.00006721701 0.000002494419

# these are the proportions of gerp mutations in the introgressed regions (mean across individuals for both populations & sd)

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

# this section has been superceded by the next **
    # ie only output - of number of gerp sites identified - what % are low/high gerp scores

# gerp bp in introgressed regions


# 2) low:high ratio of sites 


# of raw counts?
# or of number of gerp sites identified - what % are low/high gerp scores

# check with MK what is the desired output here - for low:hgih ratio
# then write function accordingly


# once low impact regions generated could perhaps generate both




#ratiogerp_perind<-function(scen,gerp_rang1,gerp_rang2){

# how many high impact gerp mutations 
#hold1<-list()
#for (ind in (1:length(scen))) {
#hold1[[ind]]<-length(findOverlaps(gerp_rang1,scen[[ind]]))
#}

# how many low impact gerp mutations
#hold2<-list()
#for (ind in (1:length(scen))) {
#hold2[[ind]]<-length(findOverlaps(gerp_rang2,scen[[ind]]))
#}

# ratio of low:high counts: # of raw counts
#ratio_counts<-list()
#for (ind in (1:length(scen))) {
#ratio_counts[[ind]]<-hold2[[ind]]/hold1[[ind]]}

# of number of gerp sites identified - what % are low/high gerp scores
#p_counts<-list()
#for (ind in (1:length(scen))) {
#p_counts[[ind]]<-hold2[[ind]]/(hold1[[ind]]+hold2[[ind]])}

#pout<-list(unlist(ratio_counts),unlist(p_counts))
#return(pout)}



#mg_highlowrat<-ratiogerp_perind(ov_gbb_99,all_high_polym_gerp_scores,all_low_polym_gerp_scores)
#el_highlowrat<-ratiogerp_perind(ov_gbg_99,all_high_polym_gerp_scores,all_low_polym_gerp_scores)


#-----------------------------------------------------------------------------------------------------------------------
#> mg_highlowrat
#[[1]]
# [1] 71.57143 42.47368 44.25926 60.11429 64.20000 65.33333 56.90000 54.41463
# [9] 59.70000 56.72000 45.61905 52.46512

#[[2]]
# [1] 0.9862205 0.9769976 0.9779051 0.9836372 0.9846626 0.9849246 0.9827288
# [8] 0.9819542 0.9835255 0.9826750 0.9785495 0.9812962

#> el_highlowrat
#[[1]]
#[1] 38.02128 43.60526 36.46154 38.25532 37.72340 66.42308 37.78571 43.79487
#[9] 40.96610

#[[2]]
#[1] 0.9743730 0.9775811 0.9733060 0.9745257 0.9741758 0.9851683 0.9742173
#[8] 0.9776760 0.9761712


# pop-level : ratio of low:high counts: # of raw counts


#> mean(mg_highlowrat[[1]])
#[1] 56.14757
#> mean(el_highlowrat[[1]])
#[1] 42.55962
# eg ~40-55X more low impact than gerp impact sites

#> mean(mg_highlowrat[[2]])
#[1] 0.9820897
#> mean(el_highlowrat[[2]])
#[1] 0.9763549
# of identified gerp sites 98% are low impact


#-----------------------------------------------------------------------------------------------------------------------

# Thu 23 Feb 2023 10:00:16 CET
#for the ratio of low:high impact sites, what would be more informative here,

#- the ratio of raw counts? eg 10X more low gerp mutations than high gerp mutations in the introgressed regions
#- or of the number of gerp sites identified - what % are low/high gerp scores? eg 90% of gerp mutations in introgressed regions of one individual are low gerp mutations


# MK
#I think both values would be similar, maybe % of high gerp would be most straightforward. 




ratiogerp_perind<-function(scen,gerp_rang1,gerp_rang2){

# how many high impact gerp mutations 
hold1<-list()
for (ind in (1:length(scen))) {
hold1[[ind]]<-length(findOverlaps(gerp_rang1,scen[[ind]]))
}

# how many low impact gerp mutations
hold2<-list()
for (ind in (1:length(scen))) {
hold2[[ind]]<-length(findOverlaps(gerp_rang2,scen[[ind]]))
}

# ratio of low:high counts: # of raw counts
#ratio_counts<-list()
#for (ind in (1:length(scen))) {
#ratio_counts[[ind]]<-hold2[[ind]]/hold1[[ind]]}

# of number of gerp sites identified - what % are low/high gerp scores
p_counts<-list()
for (ind in (1:length(scen))) {
p_counts[[ind]]<-hold1[[ind]]/(hold1[[ind]]+hold2[[ind]])}

#pout<-list(unlist(ratio_counts),unlist(p_counts))

#return(pout)

return(unlist(p_counts))
}

mg_highlowrat<-ratiogerp_perind(ov_gbb_99,all_high_polym_gerp_scores,all_low_polym_gerp_scores)
el_highlowrat<-ratiogerp_perind(ov_gbg_99,all_high_polym_gerp_scores,all_low_polym_gerp_scores)


mean(mg_highlowrat)

mean(el_highlowrat)

#> mean(mg_highlowrat)
#[1] 0.01791026
#> 
#> mean(el_highlowrat)
#[1] 0.02364506

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

# what i estimate:

For the proportion of gerp bp in introgressed regions (across MG, then across EL)


> mg_highgerp[,1]

[1] 0.000001279447

> el_highgerp[,1]

[1] 0.000001629846


Of the number of gerp sites identified - what proportion are high gerp scores? 


> mean(mg_highlowrat)

[1] 0.01791026

> mean(el_highlowrat)

[1] 0.02364506


# Thu 23 Feb 2023 10:12:31 CET
# MK
#ok, so the ELG have slightly higher proportion of deleterious alleles there than MG. That is what we want to look at. How does it compare to random regions?


q()

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

# write function to be more efficient ***

tesratiogerp_perind<-function(scen,gerp_rang1,gerp_rang2){


# then take hits per ind / width of introg regions per ind : ie high gerp bp in introgressed regions 
        # how many high impact gerp mutations which are polymorphic in E gorillas
hold<-list()
for (ind in (1:length(scen))) {
hold[[ind]]<-length(findOverlaps(gerp_rang1,scen[[ind]]))/sum(width(scen[[ind]])-1)
}


# how many high impact gerp mutations 
hold1<-list()
for (ind in (1:length(scen))) {
hold1[[ind]]<-length(findOverlaps(gerp_rang1,scen[[ind]]))
}

# how many low impact gerp mutations
hold2<-list()
for (ind in (1:length(scen))) {
hold2[[ind]]<-length(findOverlaps(gerp_rang2,scen[[ind]]))
}


# of number of gerp sites identified - what % are low/high gerp scores
p_counts<-list()
for (ind in (1:length(scen))) {
p_counts[[ind]]<-hold1[[ind]]/(hold1[[ind]]+hold2[[ind]])}

pout<-list(unlist(hold),unlist(p_counts))

return(pout)

}


test<-tesratiogerp_perind(ov_gbb_99,all_high_polym_gerp_scores,all_low_polym_gerp_scores)


 test
[[1]]
 [1] 0.0000009713453 0.0000016463520 0.0000015782552 0.0000011671335
 [5] 0.0000010662930 0.0000010951903 0.0000012826268 0.0000012374743
 [9] 0.0000011539682 0.0000012637113 0.0000015663459 0.0000013246665

[[2]]
 [1] 0.01377953 0.02300242 0.02209493 0.01636279 0.01533742 0.01507538
 [7] 0.01727116 0.01804577 0.01647446 0.01732502 0.02145046 0.01870378


> mean(test[[1]])
[1] 0.000001279447
> mean(test[[2]])
[1] 0.01791026

# yes this combines the 2 functions in a more efficient way -> process this way for the random reps

# output this, then take means later 
 cbind(mean(test[[1]]), mean(test[[2]]))
               [,1]       [,2]
[1,] 0.000001279447 0.01791026

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

# MK - Wed  1 Mar 2023 16:19:29 CET
#so does that mean there are way more deleterious mutations in the introgressed regions than expected (in EL)?
#That is quite surprising! The extent of it at least... I guess these may be low numbers anyway,
#could you also plot with the raw numbers instead of ratios?



tesratio_highgerp_perind<-function(scen,gerp_rang1){

# how many high impact gerp mutations which are polymorphic in E gorillas

hold<-list()
h_ratio<-list()
for (ind in (1:length(scen))) {
# output raw counts, as well as ratios
hold[[ind]]<-length(findOverlaps(gerp_rang1,scen[[ind]]))
# then take hits per ind / width of introg regions per ind : ie high gerp bp in introgressed regions 
h_ratio[[ind]]<-length(findOverlaps(gerp_rang1,scen[[ind]]))/sum(width(scen[[ind]])-1)
}

pout<-list(unlist(hold),unlist(h_ratio))
return(pout)
}
#-----------------------------------------------------------------------------------------------------------------------

# output for introgressed regions (as well as random regions)

mg_high<-tesratio_highgerp_perind(ov_gbb_99,all_high_polym_gerp_scores)
el_high<-tesratio_highgerp_perind(ov_gbg_99,all_high_polym_gerp_scores)

> mg_high
[[1]]
 [1] 28 57 54 35 35 36 40 41 40 25 42 43

[[2]]
 [1] 0.0000009713453 0.0000016463520 0.0000015782552 0.0000011671335
 [5] 0.0000010662930 0.0000010951903 0.0000012826268 0.0000012374743
 [9] 0.0000011539682 0.0000012637113 0.0000015663459 0.0000013246665

> el_high
[[1]]
[1] 47 38 52 47 47 26 42 39 59

[[2]]
[1] 0.0000017644630 0.0000014882118 0.0000017461970 0.0000018167762
[5] 0.0000018125024 0.0000009814284 0.0000017387704 0.0000015976404
[9] 0.0000017226277


> mean(mg_high[[1]])
[1] 39.66667
> mean(el_high[[1]])
[1] 44.11111


> mean(mg_high[[2]])
[1] 0.000001279447
> mean(el_high[[2]])
[1] 0.000001629846


#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

# Mon  6 Mar 2023 09:43:52 CET
# RH proportion plot more similar to genic content - makes more sense
#-this panel B should be the one to include

# output proportion of gerp sites (already plotted)
# output also low:high ratio (& for random regions)

# i think go back to this function - yes
ratiogerp_perind<-function(scen,gerp_rang1,gerp_rang2){

# how many high impact gerp mutations 
hold1<-list()
for (ind in (1:length(scen))) {
hold1[[ind]]<-length(findOverlaps(gerp_rang1,scen[[ind]]))
}

# how many low impact gerp mutations
hold2<-list()
for (ind in (1:length(scen))) {
hold2[[ind]]<-length(findOverlaps(gerp_rang2,scen[[ind]]))
}

# ratio of low:high counts: # of raw counts
ratio_counts<-list()
for (ind in (1:length(scen))) {
ratio_counts[[ind]]<-hold2[[ind]]/hold1[[ind]]}

# proportion
# of number of gerp sites identified - what % are low/high gerp scores
p_counts<-list()
for (ind in (1:length(scen))) {
p_counts[[ind]]<-hold1[[ind]]/(hold1[[ind]]+hold2[[ind]])}

pout<-list(unlist(ratio_counts),unlist(p_counts))


return(pout)

}

mg_highlowrat<-ratiogerp_perind(ov_gbb_99,all_high_polym_gerp_scores,all_low_polym_gerp_scores)
el_highlowrat<-ratiogerp_perind(ov_gbg_99,all_high_polym_gerp_scores,all_low_polym_gerp_scores)

> mg_highlowrat
[[1]]
 [1] 71.57143 42.47368 44.25926 60.11429 64.20000 65.33333 56.90000 54.41463
 [9] 59.70000 56.72000 45.61905 52.46512

[[2]]
 [1] 0.01377953 0.02300242 0.02209493 0.01636279 0.01533742 0.01507538
 [7] 0.01727116 0.01804577 0.01647446 0.01732502 0.02145046 0.01870378


> el_highlowrat
[[1]]
[1] 38.02128 43.60526 36.46154 38.25532 37.72340 66.42308 37.78571 43.79487
[9] 40.96610

[[2]]
[1] 0.02562704 0.02241888 0.02669405 0.02547425 0.02582418 0.01483172 0.02578269
[8] 0.02232398 0.02382876



mean(mg_highlowrat[[1]])
mean(el_highlowrat[[1]])

mean(mg_highlowrat[[2]])
 mean(el_highlowrat[[2]])

> mean(mg_highlowrat[[1]])
[1] 56.14757
> mean(el_highlowrat[[1]])
[1] 42.55962
> 
> mean(mg_highlowrat[[2]])
[1] 0.01791026
>  mean(el_highlowrat[[2]])
[1] 0.02364506