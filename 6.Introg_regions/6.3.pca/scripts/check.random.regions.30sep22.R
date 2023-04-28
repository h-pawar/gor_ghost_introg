# Fri 30 Sep 2022 16:08:34 CEST
# checking random.private.pca.sstarskov.30sep22.R
	# ie are the regions truly random - yes

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
load(file="/scratch/devel/hpawar/admix/overlap.s*.skov/overlap_objects/ov_gbb_99",verbose=T)
load(file="/scratch/devel/hpawar/admix/overlap.s*.skov/overlap_objects/ov_gbg_99",verbose=T)

#-----------------------------------------------------------------------------------------------------------------------
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


# 3) for each individual, extract their private regions
private_regions_perid<-function(id){
# subset the sbs object to the chr, start, end pos, width, strand & counts of 0 or 1 for a given ind (where 0 = ind does not carry fragment, 1 = ind has this fragment as a private fragment)  
id_1<-sbs[,id]

id_1<-as.data.frame(id_1)
id_1_s<-id_1[,c(1:3,6)]
colnames(id_1_s)<-c("seqnames","start","end","id")
# fragments found only within this individual
test<-id_1_s[which(id_1_s$id==1),]
return(test)
}


private_allids<-list()
for (ind in (1:lids)) {
private_allids[[ind]]<-private_regions_perid(ind)
}

return(private_allids)

}
#-----------------------------------------------------------------------------------------------------------------------


priv_gbb<-extract_private_regions(ov_gbb_99,12,"GBB","gbb")
priv_gbg<-extract_private_regions(ov_gbg_99,9,"GBG","gbg")
#-----------------------------------------------------------------------------------------------------------------------

# convert back to granges objects

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


calc_prop_fun<-function(fG,fB){
fG<-reduce(fG);fB<-reduce(fB) 
pos.overlaps<-reduce(subsetByOverlaps(fG,fB))
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



#-----------------------------------------------------------------------------------------------------------------------
# load in the random regions
load(file="/scratch/devel/hpawar/admix/overlap.s*.skov/privateregionsperid/random/onerep_randomregions_30sep22",verbose=T)
#-----------------------------------------------------------------------------------------------------------------------


lids=9
private_gbg_gr<-list()
for (ind in (1:lids)) {
private_gbg_gr[[ind]]<-df_granges_fun(priv_gbg,ind)
}

# check overlap b/n gbg private regions (empirical vs random) # yes minimal overlap obtained
