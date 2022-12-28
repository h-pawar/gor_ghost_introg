# Fri 30 Sep 2022 16:08:34 CEST
# checking random.private.pca.sstarskov.30sep22.R
	# ie are the regions truly random - seems yes

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

library(phangorn)
library(ape) 
library('pegas')

#-----------------------------------------------------------------------------------------------------------------------
load(file="/scratch/devel/hpawar/admix/overlap.s*.skov/overlap_objects/ov_gbb_99",verbose=T)
load(file="/scratch/devel/hpawar/admix/overlap.s*.skov/overlap_objects/ov_gbg_99",verbose=T)

#-----------------------------------------------------------------------------------------------------------------------

#scen=ov_gbb_99
#lids=12

#/scratch/devel/hpawar/admix/overlap.s*.skov/regionsperid/
#[hpawar@login1 ~]$ mkdir -p /scratch/devel/hpawar/admix/overlap.s*.skov/privateregionsperid/GBB
#[hpawar@login1 ~]$ mkdir -p /scratch/devel/hpawar/admix/overlap.s*.skov/privateregionsperid/GBG

#[hpawar@login1 ~]$ mkdir -p /scratch/devel/hpawar/admix/overlap.s*.skov/privateregionsperid/random/GBB
#[hpawar@login1 ~]$ mkdir -p /scratch/devel/hpawar/admix/overlap.s*.skov/privateregionsperid/random/GBG
#-----------------------------------------------------------------------------------------------------------------------


# extract private regions
	# intersect regions (sstar-skov) which are private to each individual

extract_private_regions<-function(scen,lids,SPE,spe){

myGRangesList<-GRangesList(scen) # convert to GRangesList object of length equal to number of individuals

# 1) calculate frequency of introgressed regions across individuals (how many regions in 1..N individuals)
	# by calling non_overlapping_region_counts function

non_overlapping_region_counts<-function(x){
reduced <- reduce(unlist(myGRangesList))
consensusIDs <- paste0("consensus_", seq(1, length(reduced)))
mcols(reduced) <- do.call(cbind, lapply(myGRangesList, function(x) (reduced %over% x) + 0))
reducedConsensus <- reduced
mcols(reducedConsensus) <- cbind(as.data.frame(mcols(reducedConsensus)), consensusIDs)
consensusIDs <- paste0("consensus_", seq(1, length(reducedConsensus)))
return(reducedConsensus)
}

testcount<-non_overlapping_region_counts()
counts<-(as.data.frame((mcols(testcount))[,1:lids]))
# freq of introg regions across the pop
freq<-rowSums(counts)
freq1<-as.data.frame(freq)


# 2) extract which freq=1 # correspond to private fragments found within one individual only
private_rows<- which(freq1$freq==1)
sbs<-testcount[private_rows,]


#return(sbs)

# 3) for each individual, extract their private regions
private_regions_perid<-function(id){
# subset the sbs object to the chr, start, end pos, width, strand & counts of 0 or 1 for a given ind (where 0 = ind does not carry fragment, 1 = ind has this fragment as a private fragment)  
id_1<-sbs[,id]
# convert granges to df 
id_1<-as.data.frame(id_1)
id_1_s<-id_1[,c(1:3,6)]
colnames(id_1_s)<-c("seqnames","start","end","id")
# fragments found only within this individual
test<-id_1_s[which(id_1_s$id==1),]
return(test)
}

# loop over private_regions_perid function, for all individuals of this population -> to obtain list object of private regions for each id of pop
private_allids<-list()
for (ind in (1:lids)) {
private_allids[[ind]]<-private_regions_perid(ind)
}

return(private_allids)
# function to convert private regions (unique to 1 individual) from df to bed format to bed file
#private_ind_df_tobed<-function(ind,SPE,spe){
#df<-(private_allids[[ind]][,c(1:3)])
#df$start<-df$start-1
# convert to same structure of sk
#df[,1]<-as.character(df[,1])
#df[,2]<-as.integer(df[,2])


#a=ind
#tmp<-paste0("/scratch/devel/hpawar/admix/overlap.s*.skov/privateregionsperid/",SPE,"/",spe,".",a,".private.tmp.bed",sep="")
#write.table(df,tmp,sep="\t",row.names=F,col.names=F,quote=F) # should write this out only once
#}

# 5) generate bed files of the private regions per individual (then intersect this with the vcfs)
#for (ind in (1:lids)) {
#private_ind_df_tobed(ind,SPE,spe) }

}
#-----------------------------------------------------------------------------------------------------------------------


priv_gbb<-extract_private_regions(ov_gbb_99,12,"GBB","gbb")
priv_gbg<-extract_private_regions(ov_gbg_99,9,"GBG","gbg")
#-----------------------------------------------------------------------------------------------------------------------

# convert back to granges objects

#str(priv_gbb)
#List of 12
# $ :'data.frame':	61 obs. of  4 variables:
#  ..$ seqnames: Factor w/ 22 levels "chr1","chr2",..: 1 1 1 1 1 1 2 2 2 3 ...
#  ..$ start   : int [1:61] 5402000 34780000 38054000 94490000 102303000 112379000 71030000 1524930

df_granges_fun<-function(scen,id){
out_rang<-GRanges(seqnames=scen[[id]][,1],ranges=IRanges(start=as.numeric(scen[[id]][,2]),end=as.numeric(scen[[id]][,3]),names=scen[[id]][,1]),strand=rep("*",length(scen[[id]][,1])))
return(out_rang)
}

lids=12
private_gbb_gr<-list()
for (ind in (1:lids)) {
private_gbb_gr[[ind]]<-df_granges_fun(priv_gbb,ind)
}


#-----------------------------------------------------------------------------------------------------------------------

# gives output

calc_prop_fun<-function(fG,fB){
fG<-reduce(fG);fB<-reduce(fB) # amendment - Wed 23 Feb 2022 17:28:27 CET
pos.overlaps<-reduce(subsetByOverlaps(fG,fB)) ## this may be more useful output? (gives ranges which overlap ie the exact positions) (yes - these are the unique overlapping regions)
prop.bp.x<-sum(reduce(subsetByOverlaps(fG,fB))@ranges@width)/sum(fG@ranges@width) # proportion of overlapping bp
return(prop.bp.x)
return(rbind(out.x,out.y))
}

# run this for each scenario as a sanity check - that the regions are random
calc_checkran<-function(fG,ran_fG){
check_ran_regions<-list()
for (ind in (1:length(fG))) {
check_ran_regions[[ind]]<-calc_prop_fun(fG[[ind]],ran_fG[[ind]])}
return(check_ran_regions)
}



#calc_checkran(priv_gbb,ran_regions_out[[1]])
#calc_checkran(priv_gbg,ran_regions_out[[2]])

#> calc_checkran(priv_gbb,ran_regions_out[[1]])
#Error in (function (classes, fdef, mtable)  : 
#  unable to find an inherited method for function ‘reduce’ for signature ‘"data.frame"’
#> calc_checkran(priv_gbg,ran_regions_out[[2]])
#Error in (function (classes, fdef, mtable)  : 
#  unable to find an inherited method for function ‘reduce’ for signature ‘"data.frame"’

  # ie need to convert back to granges **

#-----------------------------------------------------------------------------------------------------------------------
# load in the random regions
load(file="/scratch/devel/hpawar/admix/overlap.s*.skov/privateregionsperid/random/onerep_randomregions_30sep22",verbose=T)
#-----------------------------------------------------------------------------------------------------------------------

# check overlap b/n gbb private regions (empirical vs random)
  calc_checkran(private_gbb_gr,ran_regions_out[[1]])
[[1]]
[1] 0.01878734

[[2]]
[1] 0

[[3]]
[1] 0

[[4]]
[1] 0

[[5]]
[1] 0

[[6]]
[1] 0

[[7]]
[1] 0.02697376

[[8]]
[1] 0

[[9]]
[1] 0

[[10]]
[1] 0

[[11]]
[1] 0

[[12]]
[1] 0

# yes minimal overlap obtained


lids=9
private_gbg_gr<-list()
for (ind in (1:lids)) {
private_gbg_gr[[ind]]<-df_granges_fun(priv_gbg,ind)
}

# check overlap b/n gbg private regions (empirical vs random)

calc_checkran(private_gbg_gr,ran_regions_out[[2]])

[[1]]
[1] 0

[[2]]
[1] 0

[[3]]
[1] 0

[[4]]
[1] 0

[[5]]
[1] 0

[[6]]
[1] 0

[[7]]
[1] 0

[[8]]
[1] 0

[[9]]
[1] 0.01286974
