
# Tue 21 Feb 2023 12:13:04 CET
## module load gcc/6.3.0 openssl/1.0.2q R/4.0.1  BCFTOOLS/1.12
require(data.table)
options(stringsAsFactors=F)
options(scipen=100)
library(mgcv)
library(GenomicRanges)
library(tidyr)
library(ggbio)
library(valr)
library(GenomicFeatures)


# step 1 - subset gerp score files into high & low impact files
# step 2 - filter gerp mutations of each category by those polymorphic in eastern gorillas
# step 3 - overlap gerp mutations of each category with introgressed regions

#-----------------------------------------------------------------------------------------------------------------------

# 1) subset gerp score files by high & low impact files  
proc_fun<-function(a){
gerp22<-read.table(paste("/scratch/devel/shan/gerp/chr",a,".scores",sep=""))

# filter by high impact sites
gerp_hi<-gerp22[which(gerp22$V4>4),] # high impact gerp sites
colnames(gerp_hi)<-NULL
write.table(gerp_hi,file=paste("/scratch/devel/hpawar/admix/overlap.s*.skov/revisions_8feb23/gerp/chr",a,".high.scores",sep=""),sep="\t",row.names=F,col.names=F,quote=F)


#-2<gerp<2
gerp_low<-gerp22[which(-2<gerp22$V4 & gerp22$V4<2),] # low impact gerp sites
colnames(gerp_low)<-NULL

write.table(gerp_low,file=paste("/scratch/devel/hpawar/admix/overlap.s*.skov/revisions_8feb23/gerp/chr",a,".low.scores",sep=""),sep="\t",row.names=F,col.names=F,quote=F)

}


for (ind in (1:22)) {
print(ind)
proc_fun(ind)
}

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

# this section only needs to be run once *

# 2) filter gerp files by mutations segregating in easterns
    # ie remove sites fixed alternate across all gorillas (mutations arising in gorilla lineage after divergence from human ref)
    # & remove sites fixed reference across all easterns 



process_muts<-function(chrom,typ){


tmp=paste("/scratch/devel/hpawar/admix/overlap.s*.skov/revisions_8feb23/gerp/chr",chrom,".",typ,".scores",sep="")

testsnps<-system(paste("bcftools view -m2 -M2 -v snps -R ",tmp," /scratch/devel/mkuhlwilm/arch/subsets/3_gorilla_",chrom,".vcf.gz | bcftools query -f '%CHROM %POS [%GT ]\n' ",sep=""),intern=T) 
testsnps<-do.call(rbind,strsplit(testsnps,split=" "))  ## ie list merged vertically into df

testsnps[which(testsnps=="0|0")]<-0
testsnps[which(testsnps=="1|1")]<-2
testsnps[which(testsnps=="0|1")]<-1

y<-as.data.frame(testsnps[,c(3:51)],ncol=49)
y[,(1:49)] <- sapply(y[,(1:49)],as.numeric)

counts<-rowSums(y)

e_y<-y[,(1:21)]
ecounts<-rowSums(e_y)

# filter out all gor = 1/1 & all E = 0/0 

keep_snps<-testsnps[(which(counts!=(2*49) & ecounts!=0)),]

#-----------------------------------------------------------------------------------------------------------------------

keep_ranges<-GRanges(seqnames=keep_snps[,1],ranges=IRanges(start=as.numeric(keep_snps[,2])-1,end=as.numeric(keep_snps[,2]),names=keep_snps[,1]),strand=rep("*",length(keep_snps[,1])))

gerp_hi<-read.table(tmp)

gerp_ranges<-GRanges(seqnames=gerp_hi[,1],ranges=IRanges(start=as.numeric(gerp_hi[,2]),end=as.numeric(gerp_hi[,3]),names=gerp_hi[,1]),strand=rep("*",length(gerp_hi[,1])), gscore=(gerp_hi[,4]))

y<- findOverlaps(gerp_ranges,keep_ranges)
 
y1<-unique(as.data.frame(y)[,1])

return(gerp_ranges[y1])

}



#-----------------------------------------------------------------------------------------------------------------------



high_polym_gerp_scores<-list()
for (ind in (1:22)) {
print(ind)
high_polym_gerp_scores[[ind]]<-process_muts(ind,'high')
}


save(high_polym_gerp_scores,file=paste("/scratch/devel/hpawar/admix/overlap.s*.skov/revisions_8feb23/gerp/high_polym_gerp_scores"))

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
all_high_polym_gerp_scores<-(unlist(high_polym_gerp_scoresList))

low_polym_gerp_scoresList<-GRangesList(low_polym_gerp_scores) 
all_low_polym_gerp_scores<-(unlist(low_polym_gerp_scoresList))

#-----------------------------------------------------------------------------------------------------------------------


# Mon  6 Mar 2023 09:43:52 CET

# output proportion of gerp sites 
# output also low:high ratio (& for random regions)


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


mean(mg_highlowrat[[1]])
mean(el_highlowrat[[1]])

mean(mg_highlowrat[[2]])
 mean(el_highlowrat[[2]])
