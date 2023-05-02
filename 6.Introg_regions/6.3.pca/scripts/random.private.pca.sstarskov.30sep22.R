#Fri 30 Sep 2022 17:23:41 CEST
# generate random regions of equivalent length distribution as the private regions (unique to one gor id per pop)


#-----------------------------------------------------------------------------------------------------------------------
#module load gcc/6.3.0 openssl/1.0.2q  R/4.0.1  BCFTOOLS/1.14
require(data.table)
options(stringsAsFactors=F)
options(scipen=100)
library(adegenet)
library(GenomicRanges)
library(ggplot2)
library(valr)
library(tidyr)
library(tidyverse)

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

load(file="/scratch/devel/hpawar/admix/overlap.s*.skov/overlap_objects/ov_gbb_99",verbose=T)
load(file="/scratch/devel/hpawar/admix/overlap.s*.skov/overlap_objects/ov_gbg_99",verbose=T)
#-----------------------------------------------------------------------------------------------------------------------

extract_private_regions<-function(scen,lids){

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

return(private_allids)
}



sbs_gbb<-extract_private_regions(ov_gbb_99,12)
sbs_gbg<-extract_private_regions(ov_gbg_99,9)
 

df_granges_fun<-function(scen,id){
out_rang<-GRanges(seqnames=scen[[id]][,1],ranges=IRanges(start=as.numeric(scen[[id]][,2]),end=as.numeric(scen[[id]][,3]),names=scen[[id]][,1]),strand=rep("*",length(scen[[id]][,1])))
return(out_rang)
}

lids=12
private_gbb_gr<-list()
for (ind in (1:lids)) {
private_gbb_gr[[ind]]<-df_granges_fun(sbs_gbb,ind)
}


private_gbg_gr<-list()
for (ind in (1:9)) {
private_gbg_gr[[ind]]<-df_granges_fun(sbs_gbg,ind)
}


#-----------------------------------------------------------------------------------------------------------------------

samples.in.vcf<-system(paste("bcftools query -l /scratch/devel/mkuhlwilm/arch/subsets/3_gorilla_22.vcf.gz",sep=""),intern=T) # -l, --list-samples: list sample names and exit

identifiers<-cbind(samples.in.vcf, c(rep("GBB",12), rep("GBG",9), rep("GGD",1),rep("GGG",27)))
colnames(identifiers)<-c("ids","pop")
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

hg19<-fread('/home/devel/marcmont/scratch/Harvi/geneflow/test/test.reconst.ids/ora_1_1/proportion/hg19.autosomes.chrom.sizes.txt')

hg19.lengths<-apply(hg19[,2],2,function(x) as.numeric(as.character(x)))
hg19<-hg19[c(1:18,20,19,22,21),]
colnames(hg19)<-c("chrom","size") 

#-----------------------------------------------------------------------------------------------------------------------

calc_winlengths<-function(fG){
win_lengths<-list()
for (ind in (1:length(fG))) {
win_lengths[[ind]]<-table(fG[[ind]]@ranges@width-1)}
return(win_lengths)
}
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

# generate random regions of a given length distribution (matching that of the fG input for each individual li)

test_proc_randomregions<-function(fG, li) {
y<-as.data.frame(calc_winlengths(fG)[[li]])
y$Var1<-as.numeric(as.character(y$Var1))

test_randomreg<-list()
for (ind in (1:nrow(y))) {
test_randomreg[[ind]]<-bed_random(hg19, length = y[ind,1], n =  y[ind,2]) 
}

convert_tib_df<-function(lin) {
x<- test_randomreg[[lin]]
x<-as.data.frame(x)
return(x)
}

all_windowcategories<-list()
for (ind in (1:nrow(y))) {
all_windowcategories[[ind]]<-convert_tib_df(ind) }

x1<-do.call("rbind", all_windowcategories)
test_s<-GRanges(seqnames=x1[,1],ranges=IRanges(start=as.numeric(x1[,2]),end=as.numeric(x1[,3]),names=x1[,1]),strand=rep("*",length(x1[,1])))

return(test_s)
}

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

generate_random<-function(fG) {
random_reg_x<-list()
for (ind in (1:length(fG))) {
random_reg_x[[ind]]<-test_proc_randomregions(fG,ind) }
return(random_reg_x)
}

#-----------------------------------------------------------------------------------------------------------------------

sk_ran_gbb<-generate_random(private_gbb_gr) # generates for 1 rep (random regions for each id)
sk_ran_gbg<-generate_random(private_gbg_gr)


ran_regions_out<-list(sk_ran_gbb,sk_ran_gbg)
save(ran_regions_out,file="/scratch/devel/hpawar/admix/overlap.s*.skov/privateregionsperid/random/onerep_randomregions_30sep22")

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------


random_granges_tobed<-function(scen,ind,SPE,spe){

df <- data.frame(seqnames=seqnames(scen[[ind]]),
  start=start(scen[[ind]])-1,
  end=end(scen[[ind]]))

df[,1]<-as.character(df[,1])
df[,2]<-as.integer(df[,2])

a=ind
tmp<-paste("/scratch/devel/hpawar/admix/overlap.s*.skov/privateregionsperid/random/",SPE,"/",spe,".",a,".random.tmp.bed",sep="")

write.table(df,tmp,sep="\t",row.names=F,col.names=F,quote=F) # should write this out only once

}


#-----------------------------------------------------------------------------------------------------------------------

random_intersect_introgbed_regions<-function(scen,ind,SPE,spe,chrom){
a=ind
tmp<-paste(" /scratch/devel/hpawar/admix/overlap.s*.skov/privateregionsperid/random/",SPE,"/",spe,".",a,".random.tmp.bed",sep="")

chrom=chrom
testsnps<-system(paste("bcftools view -m2 -M2 -v snps -R ",tmp," /scratch/devel/mkuhlwilm/arch/subsets/3_gorilla_",chrom,".vcf.gz | bcftools query -f '%CHROM %POS [%GT ]\n' ",sep=""),intern=T) 
testsnps<-do.call(rbind,strsplit(testsnps,split=" "))

return(testsnps)
}

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

random_test_gts_topca<-function(scen,ind,SPE,spe) {

random_granges_tobed(scen,ind,SPE,spe)

autosome_sk_oneid<-list()
for (chr in (1:22)) {  
autosome_sk_oneid[[chr]]<-random_intersect_introgbed_regions(scen,ind,SPE,spe,chr)
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


pca_out<-paste("/scratch/devel/hpawar/admix/overlap.s*.skov/privateregionsperid/random/",SPE,"/pca/",spe,".",a,".random.pca",sep="")
tes<-list(pc1pc2,k,xaxis,yaxis,x3axis,y4axis)

save(tes,file=pca_out)


return(list(pc1pc2,k,xaxis,yaxis,x3axis,y4axis))

}
#-----------------------------------------------------------------------------------------------------------------------


random_pca_allids<-function(scen,SPE,spe) {
out<-list()
for (ind in (1:length(scen))) {
out[[ind]]<- random_test_gts_topca(scen,ind,SPE,spe)}
}

random_pca_allids(sk_ran_gbb,"GBB","gbb")
random_pca_allids(sk_ran_gbg,"GBG","gbg")


