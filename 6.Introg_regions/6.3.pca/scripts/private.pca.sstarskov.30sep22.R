# Fri 30 Sep 2022 17:03:35 CEST
# generate pcas of private regions (intersect sstar-skov regions which are found only in one id of the pop)

#-----------------------------------------------------------------------------------------------------------------------
#module load gcc/6.3.0 openssl/1.0.2q  R/4.0.1  BCFTOOLS/1.14
require(data.table)
options(stringsAsFactors=F)
library(tidyverse)
options(scipen=100)
library(adegenet)
library(mgcv)
library(GenomicRanges)
library(tidyr)
library(ggplot2)

#-----------------------------------------------------------------------------------------------------------------------

intersect_introgbed_regions<-function(ind,SPE,spe,chrom){
a=ind
tmp<-paste0("/scratch/devel/hpawar/admix/overlap.s*.skov/privateregionsperid/",SPE,"/",spe,".",a,".private.tmp.bed",sep="")

chrom=chrom
testsnps<-system(paste("bcftools view -m2 -M2 -v snps -R ",tmp," /scratch/devel/mkuhlwilm/arch/subsets/3_gorilla_",chrom,".vcf.gz | bcftools query -f '%CHROM %POS [%GT ]\n' ",sep=""),intern=T)
testsnps<-do.call(rbind,strsplit(testsnps,split=" "))

return(testsnps)
}
#-----------------------------------------------------------------------------------------------------------------------

samples.in.vcf<-system(paste("bcftools query -l /scratch/devel/mkuhlwilm/arch/subsets/3_gorilla_22.vcf.gz",sep=""),intern=T) # -l, --list-samples: list sample names and exit

identifiers<-cbind(samples.in.vcf, c(rep("GBB",12), rep("GBG",9), rep("GGD",1),rep("GGG",27)))
colnames(identifiers)<-c("ids","pop")
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

test_gts_topca<-function(ind,SPE,spe) {

autosome_sk_oneid<-list()
for (chr in (1:22)) {
autosome_sk_oneid[[chr]]<-intersect_introgbed_regions(ind,SPE,spe,chr)
}

aut_oneid<-do.call(rbind,autosome_sk_oneid)

aut_oneid[which(aut_oneid=="0|0")]<-0
aut_oneid[which(aut_oneid=="0|1")]<-1
aut_oneid[which(aut_oneid=="1|1")]<-2

gts<-aut_oneid

pcsub<-matrix(as.numeric(gts[,c(3:ncol(gts))]),nrow=nrow(gts))

pca_tes <- dudi.pca(t(pcsub),nf=30,scannf=F)

# % of variance explained by each component:
(k <- 100 * pca_tes$eig/sum(pca_tes$eig))

pc1pc2<-as.data.frame(cbind(pca_tes$li,identifiers[,2]))

x<-as.numeric(format(round(k[[1]], 2), nsmall = 2) )
xaxis=paste("PC1 (",x,"%)",sep="")
y<-as.numeric(format(round(k[[2]], 2), nsmall = 2) )
yaxis=paste("PC2 (",y,"%)",sep="")

x3<-as.numeric(format(round(k[[3]], 2), nsmall = 2) )
x3axis=paste("PC3 (",x3,"%)",sep="")
y4<-as.numeric(format(round(k[[4]], 2), nsmall = 2) )
y4axis=paste("PC4 (",y4,"%)",sep="")

a=ind

pca_out<-paste("/scratch/devel/hpawar/admix/overlap.s*.skov/privateregionsperid/",SPE,"/pca/",spe,".",a,".pca",sep="")
tes<-list(pc1pc2,k,xaxis,yaxis,x3axis,y4axis)

save(tes,file=pca_out)

return(list(pc1pc2,k,xaxis,yaxis,x3axis,y4axis))

}
#-----------------------------------------------------------------------------------------------------------------------

pca_allids<-function(lids,SPE,spe) {
out<-list()
for (ind in (1:lids)) {
out[[ind]]<- test_gts_topca(ind,SPE,spe)}
}


pca_allids(12,"GBB","gbb")
pca_allids(9,"GBG","gbg")



#-----------------------------------------------------------------------------------------------------------------------


