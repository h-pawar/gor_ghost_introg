#Fri 24 Feb 2023 14:51:09 CET
## module load gcc/6.3.0 openssl/1.0.2q R/4.0.1  BCFTOOLS/1.12
require(data.table)
options(stringsAsFactors=F)
library(tidyverse)
library(GenomicRanges)
library(GenomicFeatures) # if need to read in the gtf

#-----------------------------------------------------------------------------------------------------------------------


proc_polyphen_fun<-function(chrom){

# read in vep info
vep<-read.table(paste("/scratch/devel/shan/gorillas/vep/chr",chrom,".txt",sep=""))

#To get poly scores:
# filter by missense variants
nonsyn <- subset(vep, grepl("missense_variant", V7))

nonsynpoly<-nonsyn[grepl("PolyPhen",nonsyn$V14),] 

# split by poly
out <- strsplit(nonsynpoly[,14], "PolyPhen", fixed = TRUE)
o<-do.call(rbind,out)

# split second col by ;
out1 <- strsplit(o[,2], ")", fixed = TRUE)
o1<-do.call(rbind,out1)
hold<-as.data.frame(o1[,c(1:2)])

hold[,1]<-(gsub("=", "", hold[,1]))
hold1 <- strsplit(hold[,1], "(", fixed = TRUE)

hold2<-as.data.frame(do.call(rbind,hold1))
hold2[,2]<-as.numeric(hold2[,2])


pos <- strsplit(nonsynpoly[,2], ":", fixed = TRUE)
pos1<-do.call(rbind,pos)

# convert this to granges & query for polymorphic snps then query for introgressed regions
chr22poly<-(cbind(pos1, nonsynpoly[,c(7,11)],hold2))

chr22poly$chr<-paste('chr',chrom,sep='')
chr22poly$end<-(as.numeric(chr22poly[,2])+1)


chr22polydf<-chr22poly[,c(7,2,8,6,5)]


#> unique(chr22polydf$V1)
#[1] "benign"            "possibly_damaging" "probably_damaging"
#[4] "unknown" 

# ok filter this into benign / possibly&probably damaging


chr22ben<-chr22polydf[which(chr22polydf[,5]=="benign"),]
chr22dam<-chr22polydf[which(chr22polydf[,5]=="possibly_damaging" | (chr22polydf[,5]=="probably_damaging")),]



colnames(chr22ben)<-NULL
colnames(chr22dam)<-NULL

write.table(chr22ben,file=paste("/scratch/devel/hpawar/admix/overlap.s*.skov/revisions_8feb23/polyphen/chr",chrom,".benign.txt",sep=""),sep="\t",row.names=F,col.names=F,quote=F)

write.table(chr22dam,file=paste("/scratch/devel/hpawar/admix/overlap.s*.skov/revisions_8feb23/polyphen/chr",chrom,".damage.txt",sep=""),sep="\t",row.names=F,col.names=F,quote=F)

}

# ran interactively for chr 21 & 22 -> fine 


# running function interactively - check if works

for (ind in (1:21)) {
print(ind)
proc_polyphen_fun(ind)
}


#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------


# next step -> query the vcfs by these sift txt files
# filter out positions where all gor are 1/1 & where all easterns are 0/0



# chrom=22 # testing - ran fine
# typ="damage"

process_muts<-function(chrom,typ){


tmp=paste("/scratch/devel/hpawar/admix/overlap.s*.skov/revisions_8feb23/polyphen/chr",chrom,".",typ,".txt",sep="")


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

gerp_ranges<-GRanges(seqnames=gerp_hi[,1],ranges=IRanges(start=as.numeric(gerp_hi[,2]),end=as.numeric(gerp_hi[,3]),names=gerp_hi[,1]),strand=rep("*",length(gerp_hi[,1])), polyphenscore=(gerp_hi[,4]),polyphenid=(gerp_hi[,5]))

# perform overlap b/n the polymorphic snps in E & the gerp scores

y<- findOverlaps(gerp_ranges,keep_ranges)
 
y1<-unique(as.data.frame(y)[,1])


#gerp_ranges[y1] # return this object - snps annotated with gerp scores


return(gerp_ranges[y1])

}


#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------



damage_polyphen<-list()
for (ind in (1:22)) {
print(ind)
damage_polyphen[[ind]]<-process_muts(ind,'damage')
}


save(damage_polyphen,file=paste("/scratch/devel/hpawar/admix/overlap.s*.skov/revisions_8feb23/polyphen/damage_polyphen"))

# then output this object

benign_polyphen<-list()
for (ind in (1:22)) {
print(ind)
benign_polyphen[[ind]]<-process_muts(ind,'benign')
}

save(benign_polyphen,file=paste("/scratch/devel/hpawar/admix/overlap.s*.skov/revisions_8feb23/polyphen/benign_polyphen"))

#-----------------------------------------------------------------------------------------------------------------------

# 3) intersect putative introgressed regions with the polyphen scores

# read in introgressed regions
load(file="/scratch/devel/hpawar/admix/overlap.s*.skov/overlap_objects/ov_gbb_99",verbose=T)
load(file="/scratch/devel/hpawar/admix/overlap.s*.skov/overlap_objects/ov_gbg_99",verbose=T)



# read in polyphen scores
load(file=paste("/scratch/devel/hpawar/admix/overlap.s*.skov/revisions_8feb23/polyphen/damage_polyphen"),verbose=T)
load(file=paste("/scratch/devel/hpawar/admix/overlap.s*.skov/revisions_8feb23/polyphen/benign_polyphen"),verbose=T)


damage_polyphen_scoresList<-GRangesList(damage_polyphen) 
all_damage_polyphen<-(unlist(damage_polyphen_scoresList))

benign_polyphen_scoresList<-GRangesList(benign_polyphen) 
all_benign_polyphen<-(unlist(benign_polyphen_scoresList))


#all_damage_polyphen
#GRanges object with 2994 ranges and 2 metadata columns:

#> all_benign_polyphen
#GRanges object with 8472 ranges and 2 metadata columns:

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


mg_polyphen<-tesratiogerp_perind(ov_gbb_99,all_damage_polyphen,all_benign_polyphen)

el_polyphen<-tesratiogerp_perind(ov_gbg_99,all_damage_polyphen,all_benign_polyphen)


cbind(mean(mg_polyphen[[1]]), mean(mg_polyphen[[2]]))
cbind(mean(el_polyphen[[1]]), mean(el_polyphen[[2]]))


#-----------------------------------------------------------------------------------------------------------------------

> cbind(mean(mg_polyphen[[1]]), mean(mg_polyphen[[2]]))
             [,1]      [,2]
[1,] 1.193226e-06 0.2506758
> cbind(mean(el_polyphen[[1]]), mean(el_polyphen[[2]]))
             [,1]      [,2]
[1,] 1.385589e-06 0.2743534