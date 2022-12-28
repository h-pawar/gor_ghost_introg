# Fri 30 Sep 2022 17:03:35 CEST
# generate pcas of private regions (intersect sstar-skov regions which are found only in one id of the pop)
# follows from /Volumes/"Ultra USB 3.0"/IBE/further.analysis.feb.2020/gorillas/abc/simul_postabc/private.intersectregions.30sep22.R


# amending the below from 
#/scratch/devel/hpawar/admix/sstar/scripts/simul_postabc/downstream/pca_sstarskov_26jul22.R
#/scratch/devel/hpawar/admix/sstar/scripts/simul_postabc/downstream/pca_sstarskov_26jul22.arr


#-----------------------------------------------------------------------------------------------------------------------
#module load gcc/6.3.0 openssl/1.0.2q  R/4.0.1  BCFTOOLS/1.14
require(data.table)
options(stringsAsFactors=F)
library(tidyverse)
options(scipen=100)
library(adegenet)
library(mgcv)
library(GenomicRanges)
#BiocManager::install(c("GenomicRanges", "plyranges", "HelloRangesData"))
#library(plyranges) # may not be necessary here
#library(GenomicFeatures) # if need to read in the gtf
library(tidyr)
#library(ggbio)
library(ggplot2)

#-----------------------------------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------------------------------

#bed file -> intersect with vcf, extract biallelic snps only -> output chr, pos, gts

intersect_introgbed_regions<-function(ind,SPE,spe,chrom){
a=ind
#tmp<-paste("/scratch/devel/hpawar/admix/overlap.s*.skov/regionsperid/",SPE,"/",spe,".",a,".tmp.bed",sep="")

tmp<-paste0("/scratch/devel/hpawar/admix/overlap.s*.skov/privateregionsperid/",SPE,"/",spe,".",a,".private.tmp.bed",sep="")

chrom=chrom
# filter by biallelic snps only, query by the putative introg regions for id 1 & output chr, pos & gts
testsnps<-system(paste("bcftools view -m2 -M2 -v snps -R ",tmp," /scratch/devel/mkuhlwilm/arch/subsets/3_gorilla_",chrom,".vcf.gz | bcftools query -f '%CHROM %POS [%GT ]\n' ",sep=""),intern=T)
testsnps<-do.call(rbind,strsplit(testsnps,split=" "))  ## ie list merged vertically into df

return(testsnps)
}
#-----------------------------------------------------------------------------------------------------------------------
# samples & their populations
# for ids
samples.in.vcf<-system(paste("bcftools query -l /scratch/devel/mkuhlwilm/arch/subsets/3_gorilla_22.vcf.gz",sep=""),intern=T) # -l, --list-samples: list sample names and exit

# add population identifier
identifiers<-cbind(samples.in.vcf, c(rep("GBB",12), rep("GBG",9), rep("GGD",1),rep("GGG",27)))
colnames(identifiers)<-c("ids","pop")
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

# Wed 18 May 2022 12:34:45 CEST - works but memory heavy - need to send as jobs or at least mnsh -c 4

# generate bed file of introg regions for 1 individual
# intersect bed with all autosomal vcfs - extract chr, pos, gts
# convert gts to matrix of 0,1,2
# calculate pca
# output PC 1,2 & % variance explained by each

test_gts_topca<-function(ind,SPE,spe) {

# extract biallelic GTs for introg regions across all autosomes for the first individual
# loop over all chr - to obtain testsnps df for all chr
autosome_sk_oneid<-list()
for (chr in (1:22)) {
autosome_sk_oneid[[chr]]<-intersect_introgbed_regions(ind,SPE,spe,chr)
}

aut_oneid<-do.call(rbind,autosome_sk_oneid)

aut_oneid[which(aut_oneid=="0|0")]<-0
aut_oneid[which(aut_oneid=="0|1")]<-1
aut_oneid[which(aut_oneid=="1|1")]<-2

#return(aut_oneid)
gts<-aut_oneid

# convert gts to matrix
pcsub<-matrix(as.numeric(gts[,c(3:ncol(gts))]),nrow=nrow(gts)) # set to matrix

# calculate pca
pca_tes <- dudi.pca(t(pcsub),nf=30,scannf=F)
#scannf = FALSE,   # Hide scree plot: % of variance explained by each principal component # http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/119-pca-in-r-using-ade4-quick-scripts/
# nf = 30           # Number of components kept in the results

#pca_tes$eig # eigenvalues 
#(pca_tes$li) li #the row coordinates i.e. the principal components
# % of variance explained by each component:
(k <- 100 * pca_tes$eig/sum(pca_tes$eig))

# take the PCs 1-2, with the population identifiers
#pc1pc2<-as.data.frame(cbind(pca_tes$li[,c(1:2)],identifiers[,2]))
# may be best to output all PCs & all ks Thu 19 May 2022 09:31:18 CEST
pc1pc2<-as.data.frame(cbind(pca_tes$li,identifiers[,2]))
#colnames(pc1pc2)<-c("pc1","pc2","pop")

x<-as.numeric(format(round(k[[1]], 2), nsmall = 2) )
xaxis=paste("PC1 (",x,"%)",sep="")
y<-as.numeric(format(round(k[[2]], 2), nsmall = 2) )
yaxis=paste("PC2 (",y,"%)",sep="")

x3<-as.numeric(format(round(k[[3]], 2), nsmall = 2) )
x3axis=paste("PC3 (",x3,"%)",sep="")
y4<-as.numeric(format(round(k[[4]], 2), nsmall = 2) )
y4axis=paste("PC4 (",y4,"%)",sep="")

a=ind

# write this out
pca_out<-paste("/scratch/devel/hpawar/admix/overlap.s*.skov/privateregionsperid/",SPE,"/pca/",spe,".",a,".pca",sep="")
tes<-list(pc1pc2,k,xaxis,yaxis,x3axis,y4axis)

save(tes,file=pca_out)

#return(list(pc1pc2,xaxis,yaxis))

return(list(pc1pc2,k,xaxis,yaxis,x3axis,y4axis))

}
#-----------------------------------------------------------------------------------------------------------------------

#[hpawar@login1 ~]$ mkdir -p /scratch/devel/hpawar/admix/overlap.s*.skov/privateregionsperid/GBG/pca
#[hpawar@login1 ~]$ mkdir -p /scratch/devel/hpawar/admix/overlap.s*.skov/privateregionsperid/GBB/pca
#[hpawar@login1 ~]$ mkdir -p /scratch/devel/hpawar/admix/overlap.s*.skov/privateregionsperid/GBG/nj
#[hpawar@login1 ~]$ mkdir -p /scratch/devel/hpawar/admix/overlap.s*.skov/privateregionsperid/GBB/nj



# is not outputting the plots - will need to save the .pca objects, output of test_gts_topca then subsequently read in & plot

# run test_gts_topca function for all individuals in the population & output the PCs 1-2 to R objects (then subsequently plot these)

pca_allids<-function(lids,SPE,spe) {
out<-list()
for (ind in (1:lids)) {
out[[ind]]<- test_gts_topca(ind,SPE,spe)}
#pca_out<-paste("/scratch/devel/hpawar/admix/overlap.s*.skov/regionsperid/",SPE,"/",spe,".allids.pca",sep="")
#save(out,file=pca_out)
#return(out)
}


pca_allids(12,"GBB","gbb")
pca_allids(9,"GBG","gbg")



#-----------------------------------------------------------------------------------------------------------------------


