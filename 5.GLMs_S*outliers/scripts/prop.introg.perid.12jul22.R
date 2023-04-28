# Tue 12 Jul 2022 15:26:07 CEST
# calculate proportion of introgressed windows per individual for final S* dataset
#  outlier windows/ total windows with data per individual

#------------------------------------------------------------------------------------------------------------------------
## module load gcc/6.3.0 openssl/1.0.2q PYTHON/2.7.17 R/4.0.1 tabix
require(data.table)
options(stringsAsFactors=F)
options(scipen=100)
library(mgcv)
library(GenomicRanges)
library(tidyr)
library(ggbio)
library(pheatmap)
library(ggplot2)
library(ggpubr)
library("viridis")  
#-----------------------------------------------------------------------------------------------------------------------
# glm generated for the S* CIs calc from simulated data (null simulations - gor demog, no ghost)
load(file=paste("/scratch/devel/hpawar/admix/sstar/simul_postabc/stepwisesegs/final_8apr22/check_22apr22/glm",sep=""),verbose=T)
#Loading objects:
#  mods

# mods = output of glm per CI per scenario
  # 4 scenarios -> for each 3 CIs

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

# 1) calculate proportion of introgressed windows per individual

#------------------------------------------------------------------------------------------------------------------------


proc_proportion<-function(nput,ci) {

# read in the empirical data
chroms=1:23

  # update to per subspecies comparisons
# 1) WL outgroup, EL ingroup (GGG outg, GBG ing) = "GBG"
# 2) WL outgroup, EM ingroup (GGG outg, GBB ing) = "GBB"
# 3) EL outgroup, WL ingroup (GBG outg, GGG ing) = "GGG"
# 4) EL outgroup, WC ingroup (GBG outg, GGD ing) = "GGD"


cn1<-list(c("GBG","GBB","GGG","GGD"))


# 1)  GBG, CI 3, 99.5%

# for nput=1 # GBG

starout<-list()
for (chrom in (chroms)) {
starout[[chrom]]<-read.table(paste("/scratch/devel/hpawar/admix/sstar/out.subsp.comp/",cn1[[1]][[nput]],"_",chrom,".star",sep=""),sep="\t",header=T)[,c(9,4,10,52,12,13,7,1,2)]
}
# read in data for s* applied to empirical data
starout[[9]]<-read.table(paste("/scratch/devel/hpawar/admix/sstar/out.subsp.comp/",cn1[[1]][[nput]],"_newT_9.star",sep=""),sep="\t",header=T)[,c(9,4,10,52,12,13,7,1,2)] 

staroutgbb<-do.call(rbind,starout)
staroutgbb$ind_id[which(staroutgbb$ind_id=="Gbg-Tumani")]<-"Gorilla_beringei_graueri-Tumani"  

#get values per individual - for each comparison
staroutgbb[,8]<-paste(staroutgbb[,8],staroutgbb[,9])
indiv.gbb<-unique(staroutgbb[,7]) # col 7 = ind_id

starperind.gbb<-list()
total_windows_perid<-list()
for (ind in (1:length(indiv.gbb))) {
  starperind.gbb[[ind]]<-list()
  # only use windows where there is enough data (33,333 bp out of 40,000, with at least 30 segregating sites)
  allstars<- staroutgbb[which(staroutgbb[,7]==indiv.gbb[ind] & staroutgbb[,4]>=33333 &staroutgbb[,2]>=30),c(1,8,2)]
  allstars[,3]<-as.numeric(jitter(allstars[,3]))
  # create a data frame of the same segregating sites:
  newdatA=data.frame(sS=allstars[,3])
  # predict S* vals (given the segregating sites) at given ci
  out.pred<-predict(mods[[nput]][[ci]],newdata=newdatA,type="response")
  indval<-cbind(allstars,out.pred)
  # which windows lie outside the expectation for the ci (3) 
  starperind.gbb[[ind]]<-indval[which(indval[,1]>indval[,4]),,drop=F]
  # output total windows per individual
  total_windows_perid[[ind]]<-indval 
  }

## you can see that each individual has a different number of significant regions:
iva<-c()
for (ind in 1:length(indiv.gbb)) { iva[ind]<-nrow(starperind.gbb[[ind]]) }


# proportion inferred by sstar to be introgressed per individual
# outlier windows / total windows


calc_prop<-list()
for (ind in (1:length(indiv.gbb))) {
  calc_prop[[ind]]<-iva[[ind]]/nrow(total_windows_perid[[ind]])}

out_prop<-unlist(calc_prop)

return((list(iva,out_prop)))

}


# apply for the 3 CIs for the 2 eastern scenarios
calc_gbg<-list()
for (i in 1:3) {
  calc_gbg[[i]]<-proc_proportion(1,i) }


calc_gbb<-list()
for (i in 1:3) {
  calc_gbb[[i]]<-proc_proportion(2,i)}
#-----------------------------------------------------------------------------------------------------------------------
as.data.frame((calc_gbg))[,c(2,4,6)]
as.data.frame((calc_gbb))[,c(2,4,6)]
#-----------------------------------------------------------------------------------------------------------------------
