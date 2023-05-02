#-----------------------------------------------------------------------------------------------------------------------
#module load gcc/6.3.0 openssl/1.0.2q  R/4.0.1  BCFTOOLS/1.14
require(data.table)
options(stringsAsFactors=F)
library(tidyverse)
options(scipen=100)
library(adegenet)
library(GenomicRanges)
library(tidyr)
#-----------------------------------------------------------------------------------------------------------------------

args = commandArgs(trailingOnly=TRUE)

il=args[1] 
SPEC=args[2]
spec=args[3] 

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------


load(file="/scratch/devel/hpawar/admix/overlap.s*.skov/overlap_objects/ov_gbb_99",verbose=T)
load(file="/scratch/devel/hpawar/admix/overlap.s*.skov/overlap_objects/ov_gbg_99",verbose=T)
load(file="/scratch/devel/hpawar/admix/overlap.s*.skov/overlap_objects/ov_gbb40_99",verbose=T)
load(file="/scratch/devel/hpawar/admix/overlap.s*.skov/overlap_objects/ov_gbg40_99",verbose=T)

#-----------------------------------------------------------------------------------------------------------------------
library(phangorn) 
library(ape) 
library('pegas')
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
intersect_introgbed_regions<-function(ind,SPE,spe,chrom){
a=ind
tmp<-paste("/scratch/devel/hpawar/admix/overlap.s*.skov/regionsperid/",SPE,"/",spe,".",a,".tmp.bed",sep="")

chrom=chrom
# filter by biallelic snps only, query by the putative introg regions for id 1 & output chr, pos, ref allele, alt allele & gts
#testsnps<-system(paste("bcftools view -m2 -M2 -v snps -R ",tmp," /scratch/devel/mkuhlwilm/arch/subsets/3_gorilla_",chrom,".vcf.gz | bcftools query -f '%CHROM %POS [%GT ]\n' ",sep=""),intern=T) 
testsnps<-system(paste("bcftools view -m2 -M2 -v snps -R ",tmp," /scratch/devel/mkuhlwilm/arch/subsets/3_gorilla_",chrom,".vcf.gz | bcftools query -f '%CHROM %POS %REF %ALT [%GT ]\n' ",sep=""),intern=T) 
testsnps<-do.call(rbind,strsplit(testsnps,split=" "))  ## ie list merged vertically into df

return(testsnps)
}

#-----------------------------------------------------------------------------------------------------------------------

samples.in.vcf<-system(paste("bcftools query -l /scratch/devel/mkuhlwilm/arch/subsets/3_gorilla_22.vcf.gz",sep=""),intern=T) # -l, --list-samples: list sample names and exit

identifiers<-cbind(samples.in.vcf, c(rep("GBB",12), rep("GBG",9), rep("GGD",1),rep("GGG",27)))
colnames(identifiers)<-c("ids","pop")
#-----------------------------------------------------------------------------------------------------------------------

test_gts_topca<-function(ind,SPE,spe) {

autosome_sk_oneid<-list()
for (chr in (1:22)) {  
autosome_sk_oneid[[chr]]<-intersect_introgbed_regions(ind,SPE,spe,chr)
}

aut_oneid<-do.call(rbind,autosome_sk_oneid) 

# randomly sample heterozygotes for the GTs per site - equal probability 0 or 1 - 0.5
for(id in (5:ncol(aut_oneid))) {
aut_oneid[,id][which(aut_oneid[,id]=="0|1")]<-sample(c(0,1),size=length(which(aut_oneid[,id]=="0|1")),replace=T,prob=c(0.5,0.5))
}

aut_oneid[which(aut_oneid=="0|0")]<-0
aut_oneid[which(aut_oneid=="1|1")]<-1

return(aut_oneid)
}



generate_tree_perid<-function(ind,SPE,spe) {

sbset<-test_gts_topca(ind,SPE,spe)

# remove monomorphic sites 
ssb<-matrix(as.numeric(sbset[,-c(1:4)]),ncol=ncol(sbset)-4)
sr<-rowSums(ssb)
r2<-which(sr==0)
sbset2<-sbset[-r2,]

r1<-which(sr[-r2]==49) # sample sites where all gor ids are fixed derived (relative to the human ref)
sbset2<-sbset2[-sample(r1,size=round(length(r1)*.9)),]

#-----------------------------------------------------------------------------------------------------------------------

hold_rows<-list()
for (i in (1:nrow(sbset2))) { 
tes<-sbset2[i,]
ref<-tes[3]
alt<-tes[4]
tes[tes == 0] <- paste0(ref)
tes[tes == 1] <- paste0(alt)
hold_rows[[i]]<-tes
}
#-----------------------------------------------------------------------------------------------------------------------

gts<-do.call(rbind,hold_rows)

gts_2<-t(gts[,c(3,5:ncol(gts))])

sample<-c('Homo_sapiens',samples.in.vcf)
rownames(gts_2)<-sample

db_2<-as.DNAbin(gts_2)

dm_2 = dist.dna(db_2,pairwise.deletion=F,model="K80")
outg<-c('Homo_sapiens')

treeNJ = NJ(dm_2);  fit = pml(treeNJ, data=phyDat(gts_2)); fitJC = optim.pml(fit, TRUE)

#-----------------------------------------------------------------------------------------------------------------------
rootedNJ <- ladderize(root(treeNJ, outgroup=outg, resolve.root=TRUE))
#-----------------------------------------------------------------------------------------------------------------------

save(rootedNJ, file=paste("/scratch/devel/hpawar/admix/overlap.s*.skov/trees/",SPE,"/",spe,".id.",ind,".rootedNJploutg.Robj.tree",sep=""))

rootedNJ <- drop.tip(rootedNJ, "Homo_sapiens")

bss = bootstrap.pml(fitJC, bs=100, optNni=TRUE, control = pml.control(trace = 0)) 

out<-list(rootedNJ,bss)

save(out, file=paste("/scratch/devel/hpawar/admix/overlap.s*.skov/trees/",SPE,"/",spe,".id.",ind,".Robj.tree",sep=""))

}



generate_tree_perid(il,SPEC,spec)


