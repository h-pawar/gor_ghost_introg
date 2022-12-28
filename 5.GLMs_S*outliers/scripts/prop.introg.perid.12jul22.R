# Tue 12 Jul 2022 15:26:07 CEST
# calculate proportion of introgressed windows per individual for final S* dataset

# following procedure of  /Volumes/"Ultra USB 3.0"/IBE/further.analysis.feb.2020/gorillas/abc/simul_postabc/sstarperindividual.28apr22.R
# & /Volumes/"Ultra USB 3.0"/IBE/further.analysis.feb.2020/gorillas/abc/simul_postabc/prop.introg.perid.13jun22.R # for calculating introgression proportion per individual
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
# Fri 10 Jun 2022 16:14:31 CEST
# MK
#you should calculate the outlier windows/ total windows with data per individual. 
#The outlier/genome size is skewed, as you have large regions without data. 
#------------------------------------------------------------------------------------------------------------------------


proc_proportion<-function(nput,ci) {
  # for nput=1 # GBG

# read in the empirical data
chroms=1:23
# will need to specify either GB or GG
#cn1<-list(c("GB","GG"))

  # update to per subspecies comparisons
# 1) WL outgroup, EL ingroup (GGG outg, GBG ing) = "GBG"
# 2) WL outgroup, EM ingroup (GGG outg, GBB ing) = "GBB"
# 3) EL outgroup, WL ingroup (GBG outg, GGG ing) = "GGG"
# 4) EL outgroup, WC ingroup (GBG outg, GGD ing) = "GGD"


cn1<-list(c("GBG","GBB","GGG","GGD"))


# 1)  GBG, CI 3, 99.5%

# for nput=1 # GBG

#nput=1

starout<-list()
for (chrom in (chroms)) {
starout[[chrom]]<-read.table(paste("/scratch/devel/hpawar/admix/sstar/out.subsp.comp/",cn1[[1]][[nput]],"_",chrom,".star",sep=""),sep="\t",header=T)[,c(9,4,10,52,12,13,7,1,2)]
}
# read in data for s* applied to chr 9 for target individuals plus newly processed tumani (whose chr 9 had sparse data issues)
starout[[9]]<-read.table(paste("/scratch/devel/hpawar/admix/sstar/out.subsp.comp/",cn1[[1]][[nput]],"_newT_9.star",sep=""),sep="\t",header=T)[,c(9,4,10,52,12,13,7,1,2)] 

staroutgbb<-do.call(rbind,starout)

# ensure all tumani fragments are named consistently
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

#iva
#[1] 900 842 898 876 883 800 864 818 876
#> length(iva)
#[1] 9


# proportion inferred by sstar to be introgressed per individual
# outlier windows / total windows
# equivalent to length(starperind.gbb)/nrow(allstars): from the glms script

# this should work
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
> #-----------------------------------------------------------------------------------------------------------------------
> as.data.frame((calc_gbg))[,c(2,4,6)]
  c.0.0306722050554683..0.0289499012208095..0.0305329010688415..
1                                                     0.03067221
2                                                     0.02894990
3                                                     0.03053290
4                                                     0.03003901
5                                                     0.03033028
6                                                     0.02792412
7                                                     0.03026696
8                                                     0.02867129
9                                                     0.02910187
  c.0.0170457423636087..0.0160832784560053..0.0169444303733347..
1                                                     0.01704574
2                                                     0.01608328
3                                                     0.01694443
4                                                     0.01679246
5                                                     0.01705841
6                                                     0.01526012
7                                                     0.01636189
8                                                     0.01581733
9                                                     0.01656451
  c.0.0115875588875943..0.0109543589483815..0.0117268628742212..
1                                                     0.01158756
2                                                     0.01095436
3                                                     0.01172686
4                                                     0.01134694
5                                                     0.01138493
6                                                     0.01027050
7                                                     0.01113165
8                                                     0.01063776
9                                                     0.01125829
> as.data.frame((calc_gbb))[,c(2,4,6)]
   c.0.046364713893022..0.0389503946737952..0.0376264090989332..
1                                                     0.04636471
2                                                     0.03895039
3                                                     0.03762641
4                                                     0.04176229
5                                                     0.04075354
6                                                     0.03954304
7                                                     0.04340151
8                                                     0.04163619
9                                                     0.04024916
10                                                    0.03439841
11                                                    0.03823166
12                                                    0.03930346
   c.0.0269966963407561..0.0234660681411243..0.0221673013391168..
1                                                      0.02699670
2                                                      0.02346607
3                                                      0.02216730
4                                                      0.02472701
5                                                      0.02383174
6                                                      0.02297430
7                                                      0.02539530
8                                                      0.02462613
9                                                      0.02336519
10                                                     0.01997327
11                                                     0.02248254
12                                                     0.02345346
   c.0.0186240637530578..0.0160265301490429..0.0153834514412529..
1                                                      0.01862406
2                                                      0.01602653
3                                                      0.01538345
4                                                      0.01670744
5                                                      0.01669483
6                                                      0.01573651
7                                                      0.01692180
8                                                      0.01663178
9                                                      0.01586261
10                                                     0.01346682
11                                                     0.01519431
12                                                     0.01579956

# slightly lower proportions inferred (as introgressed per individual) than previously??
# makes sense, if now have a slightly higher denominator (more regions inferred on chr 9)?