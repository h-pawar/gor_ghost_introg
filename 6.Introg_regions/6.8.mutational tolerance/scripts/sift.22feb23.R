#Wed 22 Feb 2023 10:37:07 CET
## module load gcc/6.3.0 openssl/1.0.2q R/4.0.1  BCFTOOLS/1.12
require(data.table)
options(stringsAsFactors=F)
library(tidyverse)
library(GenomicRanges)
library(GenomicFeatures) # if need to read in the gtf

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

# 1) subset vep files by missense mts with sift scores -> into tolerated & deleterious impact files 

#-----------------------------------------------------------------------------------------------------------------------
proc_fun<-function(chrom){

vep<-read.table(paste("/scratch/devel/shan/gorillas/vep/chr",chrom,".txt",sep=""))

nonsyn <- subset(vep, grepl("missense_variant", V7))

nonsynsift<-nonsyn[grepl("SIFT",nonsyn$V14),] 

out <- strsplit(nonsynsift[,14], "SIFT", fixed = TRUE)
o<-do.call(rbind,out)


out1 <- strsplit(o[,2], ")", fixed = TRUE)
o1<-do.call(rbind,out1)
hold<-as.data.frame(o1[,c(1:2)])

hold[,1]<-(gsub("=", "", hold[,1]))
hold1 <- strsplit(hold[,1], "(", fixed = TRUE)

hold2<-as.data.frame(do.call(rbind,hold1))
hold2[,2]<-as.numeric(hold2[,2])


pos <- strsplit(nonsynsift[,2], ":", fixed = TRUE)
pos1<-do.call(rbind,pos)

chr22sift<-(cbind(pos1, nonsynsift[,c(7,11)],hold2))

chr22sift$chr<-paste('chr',chrom,sep='')
chr22sift$end<-(as.numeric(chr22sift[,2])+1)


chr22siftdf<-chr22sift[,c(7,2,8,6,5)]

chr22tol<-chr22siftdf[which(chr22siftdf[,5]=="tolerated" | (chr22siftdf[,5]=="tolerated_low_confidence")),]
chr22del<-chr22siftdf[which(chr22siftdf[,5]=="deleterious" | (chr22siftdf[,5]=="deleterious_low_confidence")),]


colnames(chr22tol)<-NULL
colnames(chr22del)<-NULL

write.table(chr22tol,file=paste("/scratch/devel/hpawar/admix/overlap.s*.skov/revisions_8feb23/sift/chr",chrom,".tol.txt",sep=""),sep="\t",row.names=F,col.names=F,quote=F)

write.table(chr22del,file=paste("/scratch/devel/hpawar/admix/overlap.s*.skov/revisions_8feb23/sift/chr",chrom,".del.txt",sep=""),sep="\t",row.names=F,col.names=F,quote=F)


}


for (ind in (1:22)) {
print(ind)
proc_fun(ind)
}


#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------


# extract gerp muts of interest which are polymorphic in E & intersect this with the introg regions

process_muts<-function(chrom,typ){

tmp=paste("/scratch/devel/hpawar/admix/overlap.s*.skov/revisions_8feb23/sift/chr",chrom,".",typ,".txt",sep="")

testsnps<-system(paste("bcftools view -m2 -M2 -v snps -R ",tmp," /scratch/devel/mkuhlwilm/arch/subsets/3_gorilla_",chrom,".vcf.gz | bcftools query -f '%CHROM %POS [%GT ]\n' ",sep=""),intern=T) 
testsnps<-do.call(rbind,strsplit(testsnps,split=" ")) 

testsnps[which(testsnps=="0|0")]<-0
testsnps[which(testsnps=="1|1")]<-2
testsnps[which(testsnps=="0|1")]<-1

y<-as.data.frame(testsnps[,c(3:51)],ncol=49)
y[,(1:49)] <- sapply(y[,(1:49)],as.numeric)

counts<-rowSums(y)

e_y<-y[,(1:21)]
ecounts<-rowSums(e_y)

# filter 1) sites where all gorillas are 1/1 & 2) where all easterns are 0/0
keep_snps<-testsnps[(which(counts!=(2*49) & ecounts!=0)),]

#-----------------------------------------------------------------------------------------------------------------------


keep_ranges<-GRanges(seqnames=keep_snps[,1],ranges=IRanges(start=as.numeric(keep_snps[,2])-1,end=as.numeric(keep_snps[,2]),names=keep_snps[,1]),strand=rep("*",length(keep_snps[,1])))

gerp_hi<-read.table(tmp)

gerp_ranges<-GRanges(seqnames=gerp_hi[,1],ranges=IRanges(start=as.numeric(gerp_hi[,2]),end=as.numeric(gerp_hi[,3]),names=gerp_hi[,1]),strand=rep("*",length(gerp_hi[,1])), siftscore=(gerp_hi[,4]),siftid=(gerp_hi[,5]))


y<- findOverlaps(gerp_ranges,keep_ranges)
 
y1<-unique(as.data.frame(y)[,1])

return(gerp_ranges[y1])

}



#-----------------------------------------------------------------------------------------------------------------------


del_sift<-list()
for (ind in (1:22)) {
print(ind)
del_sift[[ind]]<-process_muts(ind,'del')
}


save(del_sift,file=paste("/scratch/devel/hpawar/admix/overlap.s*.skov/revisions_8feb23/sift/del_sift"))


tol_sift<-list()
for (ind in (1:22)) {
print(ind)
tol_sift[[ind]]<-process_muts(ind,'tol')
}

save(tol_sift,file=paste("/scratch/devel/hpawar/admix/overlap.s*.skov/revisions_8feb23/sift/tol_sift"))

#-----------------------------------------------------------------------------------------------------------------------


# 3) intersect putative introgressed regions with the deleterious/tolerant sift scores


# read in sift scores
load(file=paste("/scratch/devel/hpawar/admix/overlap.s*.skov/revisions_8feb23/sift/del_sift"),verbose=T)
load(file=paste("/scratch/devel/hpawar/admix/overlap.s*.skov/revisions_8feb23/sift/tol_sift"),verbose=T)

# read in introgressed regions
load(file="/scratch/devel/hpawar/admix/overlap.s*.skov/overlap_objects/ov_gbb_99",verbose=T)
load(file="/scratch/devel/hpawar/admix/overlap.s*.skov/overlap_objects/ov_gbg_99",verbose=T)


del_sift_scoresList<-GRangesList(del_sift) 
all_del_sift<-(unlist(del_sift_scoresList))

tol_sift_scoresList<-GRangesList(tol_sift) 
all_tol_sift<-(unlist(tol_sift_scoresList))

#-----------------------------------------------------------------------------------------------------------------------


tesratiogerp_perind<-function(scen,gerp_rang1,gerp_rang2){
 
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


mg_sift<-tesratiogerp_perind(ov_gbb_99,all_del_sift,all_tol_sift)

el_sift<-tesratiogerp_perind(ov_gbg_99,all_del_sift,all_tol_sift)



cbind(mean(mg_sift[[1]]), mean(mg_sift[[2]]))

cbind(mean(el_sift[[1]]), mean(el_sift[[2]]))


#-----------------------------------------------------------------------------------------------------------------------

q()


