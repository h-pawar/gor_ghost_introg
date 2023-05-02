# Fri 30 Sep 2022 16:08:34 CEST
# extract private regions
	# intersect regions (sstar-skov) which are private to each individual

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


extract_private_regions<-function(scen,lids,SPE,spe){

myGRangesList<-GRangesList(scen)

# 1) calculate frequency of introgressed regions across individuals (how many regions in 1..N individuals)
	
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
freq<-rowSums(counts)
freq1<-as.data.frame(freq)


# 2) extract which freq=1 # correspond to private fragments found within one individual only
private_rows<- which(freq1$freq==1)
sbs<-testcount[private_rows,]


# 3) for each individual, extract their private regions
private_regions_perid<-function(id){
id_1<-sbs[,id]
id_1<-as.data.frame(id_1)
id_1_s<-id_1[,c(1:3,6)]
colnames(id_1_s)<-c("seqnames","start","end","id")
test<-id_1_s[which(id_1_s$id==1),]
return(test)
}

private_allids<-list()
for (ind in (1:lids)) {
private_allids[[ind]]<-private_regions_perid(ind)
}


private_ind_df_tobed<-function(ind,SPE,spe){
df<-(private_allids[[ind]][,c(1:3)])
df$start<-df$start-1
df[,1]<-as.character(df[,1])
df[,2]<-as.integer(df[,2])


a=ind
tmp<-paste0("/scratch/devel/hpawar/admix/overlap.s*.skov/privateregionsperid/",SPE,"/",spe,".",a,".private.tmp.bed",sep="")
write.table(df,tmp,sep="\t",row.names=F,col.names=F,quote=F) 
}


for (ind in (1:lids)) {
private_ind_df_tobed(ind,SPE,spe) }

}
#-----------------------------------------------------------------------------------------------------------------------


extract_private_regions(ov_gbb_99,12,"GBB","gbb")
extract_private_regions(ov_gbg_99,9,"GBG","gbg")

#-----------------------------------------------------------------------------------------------------------------------
