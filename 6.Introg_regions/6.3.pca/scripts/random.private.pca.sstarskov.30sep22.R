#Fri 30 Sep 2022 17:23:41 CEST
# generate random regions of equivalent length distribution as the private regions (unique to one gor id per pop)


# amendign below from 
#/Volumes/"Ultra USB 3.0"/IBE/further.analysis.feb.2020/gorillas/abc/simul_postabc/random_pca_sstarskov_26jul22.R # ran interactively


#-----------------------------------------------------------------------------------------------------------------------
#module load gcc/6.3.0 openssl/1.0.2q  R/4.0.1  BCFTOOLS/1.14
require(data.table)
options(stringsAsFactors=F)
options(scipen=100)
library(adegenet)
#library(mgcv)
library(GenomicRanges)
#library(tidyr)
#library(ggbio)
#library(pheatmap)
library(ggplot2)
#library(ggpubr)
library(valr) # need to generate random regions
library(tidyr)
library(tidyverse)

#-----------------------------------------------------------------------------------------------------------------------

# directly read in the intersect regions - Wed 27 Jul 2022 10:15:02 CEST

# need to read in the empirical data in order to generate equivalent random regions
# (1) overlap s*-skov regions - convert to granges -> obtain reduced granges -> bed -> intersect with vcf
# glm generated for the S* CIs calc from simulated data (null simulations - gor demog, no ghost)
#load(file=paste("/scratch/devel/hpawar/admix/sstar/simul_postabc/stepwisesegs/final_8apr22/check_22apr22/glm",sep=""),verbose=T)
  # issue this needs to be in module load gcc/6.3.0 openssl/1.0.2q  R/4.0.1  : installed adegenet in newer r version

#-----------------------------------------------------------------------------------------------------------------------

# read in the bed files & convert to granges - or reperform steps in private.intersectregions.30sep22.R

# ie needs output of sbs for next steps 
#-----------------------------------------------------------------------------------------------------------------------

load(file="/scratch/devel/hpawar/admix/overlap.s*.skov/overlap_objects/ov_gbb_99",verbose=T)
load(file="/scratch/devel/hpawar/admix/overlap.s*.skov/overlap_objects/ov_gbg_99",verbose=T)
#-----------------------------------------------------------------------------------------------------------------------

extract_private_regions<-function(scen,lids){

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
# convert granges to df 
id_1<-as.data.frame(id_1)
id_1_s<-id_1[,c(1:3,6)]
colnames(id_1_s)<-c("seqnames","start","end","id")
# fragments found only within this individual
test<-id_1_s[which(id_1_s$id==1),]
return(test)
}

# loop over private_regions_perid function, for all individuals of this population -> to obtain list object of private regions for each id of pop
private_allids<-list()
for (ind in (1:lids)) {
private_allids[[ind]]<-private_regions_perid(ind)
}

return(private_allids)
}



sbs_gbb<-extract_private_regions(ov_gbb_99,12)
sbs_gbg<-extract_private_regions(ov_gbg_99,9)
 
# gives list of dataframes

# convert back from df to granges

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
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
# samples & their populations
# for ids
samples.in.vcf<-system(paste("bcftools query -l /scratch/devel/mkuhlwilm/arch/subsets/3_gorilla_22.vcf.gz",sep=""),intern=T) # -l, --list-samples: list sample names and exit

# add population identifier
identifiers<-cbind(samples.in.vcf, c(rep("GBB",12), rep("GBG",9), rep("GGD",1),rep("GGG",27)))
colnames(identifiers)<-c("ids","pop")
#-----------------------------------------------------------------------------------------------------------------------


# A) generate random regions of equal length distribution as the empirical data for each individual
# (B) intersect w vcfs -> matrix of gts -> pca

#-----------------------------------------------------------------------------------------------------------------------
# lengths of hg19
hg19<-fread('/home/devel/marcmont/scratch/Harvi/geneflow/test/test.reconst.ids/ora_1_1/proportion/hg19.autosomes.chrom.sizes.txt')

hg19.lengths<-apply(hg19[,2],2,function(x) as.numeric(as.character(x)))

# hg19 chr in the correct order
hg19<-hg19[c(1:18,20,19,22,21),]
colnames(hg19)<-c("chrom","size") 

#-----------------------------------------------------------------------------------------------------------------------

# gives lengths in bp
calc_winlengths<-function(fG){
win_lengths<-list()
for (ind in (1:length(fG))) {
win_lengths[[ind]]<-table(fG[[ind]]@ranges@width-1)}
return(win_lengths)
}
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

# generate random regions of a given length distribution (matching that of the fG input for each individual li)
# run this function per individual
test_proc_randomregions<-function(fG, li) {
# for one individual (li)
y<-as.data.frame(calc_winlengths(fG)[[li]])
# convert from factor -> character -> numeric
y$Var1<-as.numeric(as.character(y$Var1))

# write as a function - go through each of lengths & n for this individual
test_randomreg<-list()
for (ind in (1:nrow(y))) {
#test_randomreg[[ind]]<-bed_random(hg19, length = y[ind,1], n =  y[ind,2], seed = 10104) 
# run without specifying the seed
test_randomreg[[ind]]<-bed_random(hg19, length = y[ind,1], n =  y[ind,2]) 
}


convert_tib_df<-function(lin) {
x<- test_randomreg[[lin]]
x<-as.data.frame(x)
# convert to granges - then this will be what will be being intersected - later 
#test_s<-GRanges(seqnames=x[,1],ranges=IRanges(start=as.numeric(x[,2]),end=as.numeric(x[,3]),names=x[,1]),strand=rep("*",length(x[,1])))
#return(test_s)
return(x)
}

all_windowcategories<-list()
for (ind in (1:nrow(y))) {
all_windowcategories[[ind]]<-convert_tib_df(ind) }

#do.call("rbind", list(dataframes to merge))
x1<-do.call("rbind", all_windowcategories)
# convert to granges - then this will be what will be being intersected
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
sk_ran_gbg<-generate_random(private_gbg_gr) # generates for 1 rep (random regions for each id)


# could output here also - the error was at the end of the run script -> perhaps unnecessary to output the random regions here

ran_regions_out<-list(sk_ran_gbb,sk_ran_gbg)
save(ran_regions_out,file="/scratch/devel/hpawar/admix/overlap.s*.skov/privateregionsperid/random/onerep_randomregions_30sep22")


#[hpawar@login1 ~]$ mkdir -p /scratch/devel/hpawar/admix/overlap.s*.skov/privateregionsperid/random/GBG/pca
#[hpawar@login1 ~]$ mkdir -p /scratch/devel/hpawar/admix/overlap.s*.skov/privateregionsperid/random/GBG/nj
#[hpawar@login1 ~]$ mkdir -p /scratch/devel/hpawar/admix/overlap.s*.skov/privateregionsperid/random/GBB/pca
#[hpawar@login1 ~]$ mkdir -p /scratch/devel/hpawar/admix/overlap.s*.skov/privateregionsperid/random/GBB/nj

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
# check overlap b/n these 2 objects (ie between the empirical & the random regions for individual 1)

# check this section **

#calc_prop_fun<-function(fG,fB){
#fG<-reduce(fG);fB<-reduce(fB) # amendment - Wed 23 Feb 2022 17:28:27 CET
#pos.overlaps<-reduce(subsetByOverlaps(fG,fB)) ## this may be more useful output? (gives ranges which overlap ie the exact positions) (yes - these are the unique overlapping regions)
#prop.bp.x<-sum(reduce(subsetByOverlaps(fG,fB))@ranges@width)/sum(fG@ranges@width) # proportion of overlapping bp
#return(prop.bp.x)
#return(rbind(out.x,out.y))
#}

# run this for each scenario as a sanity check - that the regions are random
#calc_checkran<-function(fG,ran_fG){
#check_ran_regions<-list()
#for (ind in (1:length(fG))) {
#check_ran_regions[[ind]]<-calc_prop_fun(fG[[ind]],ran_fG[[ind]])}
#return(check_ran_regions)
#}


#calc_checkran(ov_gbb_99,sk_ran_gbb)
#calc_checkran(ov_gbg_99,sk_ran_gbg)
#calc_checkran(ov_gbb40_99,sk_ran_gbb_40)
#calc_checkran(ov_gbg40_99,sk_ran_gbg_40)

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------


# functions with paths outputting to the random data directories
# granges of introg regions for 1 individual -> bed file 
  # contains introg regions across all autosomes for the specified individual
random_granges_tobed<-function(scen,ind,SPE,spe){

df <- data.frame(seqnames=seqnames(scen[[ind]]),
  start=start(scen[[ind]])-1,
  end=end(scen[[ind]]))

# convert to same structure of sk
df[,1]<-as.character(df[,1])
df[,2]<-as.integer(df[,2])

a=ind
tmp<-paste("/scratch/devel/hpawar/admix/overlap.s*.skov/privateregionsperid/random/",SPE,"/",spe,".",a,".random.tmp.bed",sep="")

write.table(df,tmp,sep="\t",row.names=F,col.names=F,quote=F) # should write this out only once

}

#random_granges_tobed(sk_ran_gbb,1,"GBB","gbb")


#-----------------------------------------------------------------------------------------------------------------------

#bed file -> intersect with vcf, extract biallelic snps only -> output chr, pos, gts

random_intersect_introgbed_regions<-function(scen,ind,SPE,spe,chrom){
a=ind
tmp<-paste(" /scratch/devel/hpawar/admix/overlap.s*.skov/privateregionsperid/random/",SPE,"/",spe,".",a,".random.tmp.bed",sep="")

chrom=chrom
# filter by biallelic snps only, query by the putative introg regions for id 1 & output chr, pos & gts
testsnps<-system(paste("bcftools view -m2 -M2 -v snps -R ",tmp," /scratch/devel/mkuhlwilm/arch/subsets/3_gorilla_",chrom,".vcf.gz | bcftools query -f '%CHROM %POS [%GT ]\n' ",sep=""),intern=T) 
testsnps<-do.call(rbind,strsplit(testsnps,split=" "))  ## ie list merged vertically into df

return(testsnps)
}

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

# Wed 18 May 2022 12:34:45 CEST - works but memory heavy - need to send as jobs or at least mnsh -c 4

# generate bed file of introg regions for 1 individual
# intersect bed with all autosomal vcfs - extract chr, pos, gts
# convert gts to matrix of 0,1,2
# calculate pca
# output PC 1,2 & % variance explained by each

random_test_gts_topca<-function(scen,ind,SPE,spe) {

random_granges_tobed(scen,ind,SPE,spe)

# extract biallelic GTs for introg regions across all autosomes for the first individual
# loop over all chr - to obtain testsnps df for all chr
autosome_sk_oneid<-list()
for (chr in (1:22)) {  
autosome_sk_oneid[[chr]]<-random_intersect_introgbed_regions(scen,ind,SPE,spe,chr)
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

#-----------------------------------------------------------------------------------------------------------------------
# samples & their populations
# for ids
#samples.in.vcf<-system(paste("bcftools query -l /scratch/devel/mkuhlwilm/arch/subsets/3_gorilla_22.vcf.gz",sep=""),intern=T) # -l, --list-samples: list sample names and exit

# add population identifier
#identifiers<-cbind(samples.in.vcf, c(rep("GBB",12), rep("GBG",9), rep("GGD",1),rep("GGG",27)))
#colnames(identifiers)<-c("ids","pop")
#-----------------------------------------------------------------------------------------------------------------------

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
pca_out<-paste("/scratch/devel/hpawar/admix/overlap.s*.skov/privateregionsperid/random/",SPE,"/pca/",spe,".",a,".random.pca",sep="")
#tes<-list(pc1pc2,xaxis,yaxis)
tes<-list(pc1pc2,k,xaxis,yaxis,x3axis,y4axis)

save(tes,file=pca_out)

#return(list(pc1pc2,xaxis,yaxis))

return(list(pc1pc2,k,xaxis,yaxis,x3axis,y4axis))

}
#-----------------------------------------------------------------------------------------------------------------------

# is not outputting the plots - will need to save the .random.pca objects, output of random_test_gts_topca then subsequently read in & plot

# run random_test_gts_topca function for all individuals in the population & output the PCs 1-2 to R objects (then subsequently plot these)

random_pca_allids<-function(scen,SPE,spe) {
out<-list()
for (ind in (1:length(scen))) {
out[[ind]]<- random_test_gts_topca(scen,ind,SPE,spe)}
#return(out)
}

# apply pca_all_ids function to all scenarios - GBB, GBG, with & without 40kb cutoff
random_pca_allids(sk_ran_gbb,"GBB","gbb")
random_pca_allids(sk_ran_gbg,"GBG","gbg")


