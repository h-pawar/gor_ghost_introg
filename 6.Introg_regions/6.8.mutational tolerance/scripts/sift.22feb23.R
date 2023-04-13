#sift/polyphen -> extract from the vep files

#Wed 22 Feb 2023 10:37:07 CET
## module load gcc/6.3.0 openssl/1.0.2q R/4.0.1  BCFTOOLS/1.12
require(data.table)
options(stringsAsFactors=F)
library(tidyverse)
library(GenomicRanges)
library(GenomicFeatures) # if need to read in the gtf

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

# Fri 24 Feb 2023 14:24:11 CET

# 1) subset vep files by missense mts with sift scores -> into tolerated & deleterious impact files 

#unique(nonsyn$V7)
#[1] "missense_variant"                       
#[2] "missense_variant,splice_region_variant" 
#[3] "missense_variant,NMD_transcript_variant"

# vep files = 1-based

#-----------------------------------------------------------------------------------------------------------------------
proc_fun<-function(chrom){

# read in vep info
vep<-read.table(paste("/scratch/devel/shan/gorillas/vep/chr",chrom,".txt",sep=""))

#To get SIFT scores:
# filter by missense variants
nonsyn <- subset(vep, grepl("missense_variant", V7))

nonsynsift<-nonsyn[grepl("SIFT",nonsyn$V14),] 

# split by SIFT
out <- strsplit(nonsynsift[,14], "SIFT", fixed = TRUE)
o<-do.call(rbind,out)

# split second col by ;
out1 <- strsplit(o[,2], ")", fixed = TRUE)
o1<-do.call(rbind,out1)
hold<-as.data.frame(o1[,c(1:2)])

hold[,1]<-(gsub("=", "", hold[,1]))
hold1 <- strsplit(hold[,1], "(", fixed = TRUE)

hold2<-as.data.frame(do.call(rbind,hold1))
hold2[,2]<-as.numeric(hold2[,2])


pos <- strsplit(nonsynsift[,2], ":", fixed = TRUE)
pos1<-do.call(rbind,pos)

# convert this to granges & query for polymorphic snps then query for introgressed regions
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

# ran interactively for chr 21 & 22 -> fine 


# running function interactively - check if works

for (ind in (1:20)) {
print(ind)
proc_fun(ind)
}


#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------


# next step -> query the vcfs by these sift txt files
# filter out positions where all gor are 1/1 & where all easterns are 0/0



# chrom=22 # testing - ran fine
# typ="del"

process_muts<-function(chrom,typ){


tmp=paste("/scratch/devel/hpawar/admix/overlap.s*.skov/revisions_8feb23/sift/chr",chrom,".",typ,".txt",sep="")


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

gerp_ranges<-GRanges(seqnames=gerp_hi[,1],ranges=IRanges(start=as.numeric(gerp_hi[,2]),end=as.numeric(gerp_hi[,3]),names=gerp_hi[,1]),strand=rep("*",length(gerp_hi[,1])), siftscore=(gerp_hi[,4]),siftid=(gerp_hi[,5]))

# perform overlap b/n the polymorphic snps in E & the gerp scores

y<- findOverlaps(gerp_ranges,keep_ranges)
 
y1<-unique(as.data.frame(y)[,1])


#gerp_ranges[y1] # return this object - snps annotated with gerp scores


return(gerp_ranges[y1])

}



#-----------------------------------------------------------------------------------------------------------------------


del_sift<-list()
for (ind in (1:22)) {
print(ind)
del_sift[[ind]]<-process_muts(ind,'del')
}


save(del_sift,file=paste("/scratch/devel/hpawar/admix/overlap.s*.skov/revisions_8feb23/sift/del_sift"))

# then output this object

tol_sift<-list()
for (ind in (1:22)) {
print(ind)
tol_sift[[ind]]<-process_muts(ind,'tol')
}

save(tol_sift,file=paste("/scratch/devel/hpawar/admix/overlap.s*.skov/revisions_8feb23/sift/tol_sift"))

#-----------------------------------------------------------------------------------------------------------------------

# write next step - intersect with introgressed regions

# 3) intersect putative introgressed regions with the deleterious/tolerant sift scores


# read in sift scores
load(file=paste("/scratch/devel/hpawar/admix/overlap.s*.skov/revisions_8feb23/sift/del_sift"),verbose=T)
load(file=paste("/scratch/devel/hpawar/admix/overlap.s*.skov/revisions_8feb23/sift/tol_sift"),verbose=T)



# read in introgressed regions
load(file="/scratch/devel/hpawar/admix/overlap.s*.skov/overlap_objects/ov_gbb_99",verbose=T)
load(file="/scratch/devel/hpawar/admix/overlap.s*.skov/overlap_objects/ov_gbg_99",verbose=T)


del_sift_scoresList<-GRangesList(del_sift) 
all_del_sift<-(unlist(del_sift_scoresList))
#GRanges object with 3312 ranges and 2 metadata columns:

tol_sift_scoresList<-GRangesList(tol_sift) 
all_tol_sift<-(unlist(tol_sift_scoresList))
#GRanges object with 7534 ranges and 2 metadata columns:

#-----------------------------------------------------------------------------------------------------------------------


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


mg_sift<-tesratiogerp_perind(ov_gbb_99,all_del_sift,all_tol_sift)

el_sift<-tesratiogerp_perind(ov_gbg_99,all_del_sift,all_tol_sift)



> cbind(mean(mg_sift[[1]]), mean(mg_sift[[2]]))
             [,1]      [,2]
[1,] 1.393881e-06 0.3115521
> 
> cbind(mean(el_sift[[1]]), mean(el_sift[[2]]))
             [,1]      [,2]
[1,] 1.553546e-06 0.3107029


# second val is surprisingly high.. (compared to the gerp vals)
# perhaps fewer tolerated mutations identified than expected???
        # or fewer kept after filtering out those 1/1 in all gorillas..

#-----------------------------------------------------------------------------------------------------------------------


> mg_sift
[[1]]
 [1] 1.387636e-06 1.270868e-06 1.285986e-06 1.233827e-06 1.035827e-06
 [6] 1.277722e-06 1.571218e-06 1.478933e-06 1.038571e-06 1.364808e-06
[11] 1.901991e-06 1.879178e-06

[[2]]
 [1] 0.3053435 0.3034483 0.2820513 0.3109244 0.2556391 0.3181818 0.3266667
 [8] 0.3101266 0.2950820 0.2842105 0.3750000 0.3719512

> el_sift
[[1]]
[1] 1.914630e-06 1.644866e-06 1.410390e-06 1.352918e-06 1.079789e-06
[6] 1.623132e-06 1.945767e-06 1.433780e-06 1.576642e-06

[[2]]
[1] 0.3493151 0.3529412 0.2709677 0.3097345 0.2666667 0.2925170 0.3154362
[8] 0.3211009 0.3176471

#-----------------------------------------------------------------------------------------------------------------------

q()


