#!/usr/bin/env
#module load R/4.0.1

#Fri 21 Jan 2022 17:33:15 GMT
# 3) correlate each of the empirical summary stat measures with each other

# rerunning - generateabc.6oct.R
# to obtain matrix of summary stats for empirical data - then perform correlations
# which of the summary stats are correlated with each other?

# MK: For the summary stats, you could just calculate pairwise correlations and plot that as a heatmap, that would be the simplest solution.


#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
# Wed  6 Oct 2021 09:42:10 CEST
# rerunning abc analysis after recalculating "fixedsites" and "fixedsitesperid" in the empirical data
# needed to  exclude differences to human (1/1 across all individuals) in the empirical data
#-----------------------------------------------------------------------------------------------------------------------

# Wed 29 Sep 2021 09:36:54 CEST
# PATHS - Mon 27 Sep 2021 11:15:44 CEST
# new batch of simulated data (700 simn reps) & summary stats were generated with 
    # /scratch/devel/hpawar/admix/abc/simul/scripts/test.abc.model.arr - which calls
    # /scratch/devel/hpawar/admix/abc/simul/scripts/test.abc.model.v5.R # new version supercedes test.abc.model.v4.R (used to generate prev 1-2000 simns)

# empirical data calculate 
#- number of population-wise fixed sites and the number of population-wise segregating sites; 
# - the fixed sites per individual
    #/scratch/devel/hpawar/admix/abc/simul/scripts/segs.emp.27sep.R  # called by 
    #/scratch/devel/hpawar/admix/abc/simul/scripts/het.emp.26aug.arr

# 3) filtered by the informative windows - & only retain summary statistics in sufficiently informative windows (>0.5)
# ran interactively - /Volumes/"Ultra USB 3.0"/IBE/further.analysis.feb.2020/gorillas/abc/filter.segs.by.clbl.27sep.R
    # note this is equivalent to /Volumes/"Ultra USB 3.0"/IBE/further.analysis.feb.2020/gorillas/abc/filter.het.by.clbl.26aug.R
    # but with paths amended for new statistics

# process new stats & normalise by data coverage
# /Volumes/"Ultra USB 3.0"/IBE/further.analysis.feb.2020/gorillas/abc/process.seg.emp.27sep.R

# 4) format data for abc analysis - run the following interactively
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
library(abc)
options(scipen=100)
require(data.table)
options(stringsAsFactors=F)

# read in for empirical data
files <- paste("/scratch/devel/hpawar/admix/abc/emp.data/window.info/subset/", list.files(path = "/scratch/devel/hpawar/admix/abc/emp.data/window.info/subset/", pattern="_informativewindows"), sep = "")

genomewide=list()
for (i in 1:length(files)){
load(file=files[[i]],verbose=T)
genomewide[[i]]<-stats_inform1
}

#replace NaN values with 0
find.nas_fun1<-function(x) {
for (j in 1:length(genomewide[[x]])){
    # convert window info from character -> numeric - so do not lose this info when replacing non-finite numbers with 0
genomewide[[x]][[j]][[1]]<-as.numeric(genomewide[[x]][[j]][[1]])
}
out_1<-list()
for (j in 1:length(genomewide[[x]])){
out_1[[j]]<-lapply(genomewide[[x]][[j]], function(y) replace(y, !is.finite(y), 0))
}
return(out_1)
}

#MK: If they are NaN (not a number), they are due to division by zero (I guess actually often 0/0) and might be corrected to 0 instead. 

# run for all 22 chr
finalwindows=list()
for (i in 1:length(genomewide)){
finalwindows[[i]]<-find.nas_fun1(i)
}

#> str(genomewide_out)
#List of 22
# $ :List of 1735
#  ..$ :List of 6
#-----------------------------------------------------------------------------------------------------------------------

format_stats_fun<-function(x,y) {
otest=list()
for (i in 1:length(finalwindows[[x]])){
otest[[i]]<-t(as.data.frame(finalwindows[[x]][[i]][[y]]))
}
outest<-data.frame(matrix(unlist(otest), nrow=length(otest), byrow=TRUE),stringsAsFactors=FALSE)
return(outest)
}


# could write this second section into a function also - yes => applicable for all summary stats
stats_asdf_fun<-function(x){
stats_o=list()
for (i in 1:22){
stats_o[[i]]<-format_stats_fun(i,x)
}
ptest<-do.call(rbind, stats_o)
return(ptest)
}

# summary stats for all retained windows (genome-wide) - for fst, pi, tajimasd
fst_df<-stats_asdf_fun(4)
pi_df<-stats_asdf_fun(5)
tajima_df<-stats_asdf_fun(6)

inf_fsts=list()
for (i in 1:ncol(fst_df)){
inf_fsts[[i]]<-which(fst_df[,i] > 1)
}

#sort(unlist(inf_fsts))
# [1]   551   741   984  1820  2339  2781  3202  3216  4047  4047  4461  5887
#[13]  6161  7976  9234  9234  9888 10136 10235 10304 10743 12926 12973 13803
#[25] 14435 15478 16690 16899 17010 17010 17369 17369 17402 18040 18275 19526
#[37] 19781 19860 20050 20353 20745 20878 20977 22584 22734 22779 22842 23572
#[49] 24226 24383 24901 24901 25196 25276 27006 27006 27104 28020 28474 28488
#[61] 28661 29120 29614 29624 29645 30666 30822 31228

#fst_df[551,]
#           X1        X2        X3              X4        X5        X6
#551 0.8502885 0.1604438 0.1859062 269122666808322 0.4949833 0.1756261
 
# fst_df[741,]
#            X1        X2        X3        X4              X5       X6
#741 0.09075206 0.2516275 0.2162011 0.5186349 299325369861235 0.133672

#length(sort(unlist(inf_fsts)))
#[1] 68

#-----------------------------------------------------------------------------------------------------------------------
# me
#1. for the empirical data after replacing NaNs and Inf with 0 there remains 68 windows with 
#fsts > 1 - can I remove these windows? As it looks like something has gone wrong there (eg below).

#MK
#good that you checked that! Indeed, this FST seems to be an error.
#-----------------------------------------------------------------------------------------------------------------------

fst_df1<-fst_df[-c(sort(unlist(inf_fsts))),]

#sapply(fst_df1, mean) # more normal vals
#   X1        X2        X3        X4        X5        X6 
#0.1344107 0.1852550 0.1766149 0.3601883 0.3375883 0.2005140 

# remove the same rows from the other dfs **

pi_df1<-pi_df[-c(sort(unlist(inf_fsts))),]
tajima_df1<-tajima_df[-c(sort(unlist(inf_fsts))),]

# nrow(tajima_df1)
#[1] 31552

#-----------------------------------------------------------------------------------------------------------------------

# 2) read in heterozygosity files 
    # generated using
        # /scratch/devel/hpawar/admix/abc/simul/scripts/het.emp.31aug.R, 
          # which was rerun Wed  1 Sep 2021 15:08:55 CEST - recalculating seg sites - removing positions fixed between all gor subspecies
        #called by - /scratch/devel/hpawar/admix/abc/simul/scripts/het.emp.26aug.arr
        #  #/Volumes/"Ultra USB 3.0"/IBE/further.analysis.feb.2020/gorillas/abc/filter.het.by.26aug.R - ran interactively
# only the informative windows

hetfiles <- paste("/scratch/devel/hpawar/admix/abc/emp.data/test/", list.files(path = "/scratch/devel/hpawar/admix/abc/emp.data/test/", pattern="_informativewindows"), sep = "")


hetgenomewide=list()
for (i in 1:length(hetfiles)){
load(file=hetfiles[[i]],verbose=T)
hetgenomewide[[i]]<-stats_inform1
}


# 2.1) calc heterozygosity per individual

format_hetperid_fun<-function(x,y) {
htest=list()
for (i in 1:length(hetgenomewide[[x]])){
htest[[i]]<-t(as.data.frame(hetgenomewide[[x]][[i]][[2]][[y]]))
}
hutest<-data.frame(matrix(unlist(htest), nrow=length(htest), byrow=TRUE),stringsAsFactors=FALSE)
return(hutest)
}


hetperid_asdf_fun<-function(y){
stats_h=list()
for (i in 1:22){
stats_h[[i]]<-format_hetperid_fun(i,y)
}
ptest<-do.call(rbind, stats_h)
return(ptest)
}

hetperid_wl<-hetperid_asdf_fun(1)
hetperid_wc<-hetperid_asdf_fun(2)
hetperid_el<-hetperid_asdf_fun(3)
hetperid_em<-hetperid_asdf_fun(4)

# removing windows remvoed for the other summary stats b/c infinite fsts
hetperid_wl1<-hetperid_wl[-c(sort(unlist(inf_fsts))),]
hetperid_wc1<-hetperid_wc[-c(sort(unlist(inf_fsts))),]
hetperid_el1<-hetperid_el[-c(sort(unlist(inf_fsts))),]
hetperid_em1<-hetperid_em[-c(sort(unlist(inf_fsts))),]


# 3) sum of all heterozygous sites divided by the sum of all sites with data
# *1000 for het per kb

# estimation of number of sites with data - go to sites.wdata.31aug.R 
    # /Volumes/"Ultra USB 3.0"/IBE/further.analysis.feb.2020/gorillas/abc/sites.wdata.31aug.R


# get in the required format
#"het_mu_WL", "het_mu_WC", "het_mu_EL", "het_mu_EM",
#"het_sd_WL", "het_sd_WC", "het_sd_EL", "het_sd_EM",

het_mu1<-c(
mean((colSums(hetperid_wl1)/767730079))*1000,
sum(hetperid_wc1/767730079)*1000,
mean((colSums(hetperid_el1)/767730079))*1000,
mean((colSums(hetperid_em1)/767730079))*1000
    )

het_sd1<-c(
sd((colSums(hetperid_wl1)/767730079))*1000,
sd((colSums(hetperid_el1)/767730079))*1000,
sd((colSums(hetperid_em1)/767730079))*1000
    )

#-----------------------------------------------------------------------------------------------------------------------

# from process.seg.emp.27sep.R for new stats 
#- number of population-wise fixed sites and the number of population-wise segregating sites;
# - the fixed sites per individual: normalising this by number of sites with data 

# Tue 28 Sep 2021 11:46:21 CEST
# 2) read in files w newly calculated stats
#- number of population-wise fixed sites and the number of population-wise segregating sites; 
# - the fixed sites per individual
# generated using
    #/scratch/devel/hpawar/admix/abc/simul/scripts/segs.emp.27sep.R  # called by 
    #/scratch/devel/hpawar/admix/abc/simul/scripts/het.emp.26aug.arr
#then ran interactively - /Volumes/"Ultra USB 3.0"/IBE/further.analysis.feb.2020/gorillas/abc/filter.segs.by.clbl.27sep.R
    # note this is equivalent to /Volumes/"Ultra USB 3.0"/IBE/further.analysis.feb.2020/gorillas/abc/filter.het.by.clbl.26aug.R
    # but with paths amended for new statistics

segsfiles <- paste("/scratch/devel/hpawar/admix/abc/emp.data/test/segs/", list.files(path = "/scratch/devel/hpawar/admix/abc/emp.data/test/segs/", pattern="_informativewindows"), sep = "")

segsgenomewide=list()
for (i in 1:length(segsfiles)){
load(file=segsfiles[[i]],verbose=T)
segsgenomewide[[i]]<-stats_inform1
}

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
# A) number of population-wise fixed sites

format_fixedsitesperid_fun<-function(x,y,z) {
htest=list()
for (i in 1:length(segsgenomewide[[x]])){
htest[[i]]<-t(as.data.frame(segsgenomewide[[x]][[i]][[2]][z,][[y]]))
}
hutest<-data.frame(matrix(unlist(htest), nrow=length(htest), byrow=TRUE),stringsAsFactors=FALSE)
return(hutest)
}


fixedsitesperid_asdf_fun<-function(y,z){
stats_h=list()
for (i in 1:22){
stats_h[[i]]<-format_fixedsitesperid_fun(i,y,z)
}
ptest<-do.call(rbind, stats_h)
return(ptest)
}



fixedsitesperid_wl<-fixedsitesperid_asdf_fun(1,1)
# remove stats from windows with infinite fsts
fixedsitesperid_wl1<-fixedsitesperid_wl[-c(sort(unlist(inf_fsts))),]
fixedsitesperid_wc<-fixedsitesperid_asdf_fun(2,1)
fixedsitesperid_wc1<-fixedsitesperid_wc[-c(sort(unlist(inf_fsts))),]
fixedsitesperid_el<-fixedsitesperid_asdf_fun(3,1)
fixedsitesperid_el1<-fixedsitesperid_el[-c(sort(unlist(inf_fsts))),]
fixedsitesperid_em<-fixedsitesperid_asdf_fun(4,1)
fixedsitesperid_em1<-fixedsitesperid_em[-c(sort(unlist(inf_fsts))),]


# perhaps should normalise by data coverage here also?

pop_fixedsites<-rbind(sum(fixedsitesperid_wl1/767730079)*1000,
sum(fixedsitesperid_wc1/767730079)*1000,
sum(fixedsitesperid_el1/767730079)*1000,
sum(fixedsitesperid_em1/767730079)*1000)


# A) number of population-wise fixed sites
# WL, WC, EL, EM
#           [,1]
#[1,] 0.07708829
#[2,] 0.72709409
#[3,] 0.53358597
#[4,] 0.48170967
#-----------------------------------------------------------------------------------------------------------------------
# B) number of population-wise segregating sites
# z=2 for B

segsitesperid_wl<-fixedsitesperid_asdf_fun(1,2)
segsitesperid_wl1<-segsitesperid_wl[-c(sort(unlist(inf_fsts))),]
segsitesperid_wc<-fixedsitesperid_asdf_fun(2,2)
segsitesperid_wc1<-segsitesperid_wc[-c(sort(unlist(inf_fsts))),]
segsitesperid_el<-fixedsitesperid_asdf_fun(3,2)
segsitesperid_el1<-segsitesperid_el[-c(sort(unlist(inf_fsts))),]
segsitesperid_em<-fixedsitesperid_asdf_fun(4,2)
segsitesperid_em1<-segsitesperid_em[-c(sort(unlist(inf_fsts))),]

pop_segsites<-rbind(sum(segsitesperid_wl1/767730079)*1000,
sum(segsitesperid_wc1/767730079)*1000,
sum(segsitesperid_el1/767730079)*1000,
sum(segsitesperid_em1/767730079)*1000)

# B) number of population-wise segregating sites
# WL, WC, EL, EM
#    [,1]
#[1,] 3.9409085
#[2,] 0.6936383
#[3,] 1.2544721
#[4,] 1.4698317

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

# 2.2) fixed sites per individual
#now fixed sites are output in the same way as the het per id
#-----------------------------------------------------------------------------------------------------------------------
# 2.2) calc fixed sites per individual

# this calculation is equivalent to processing het sites per id in het.emp.31aug.R
format_fixedperid_fun<-function(x,y) {
htest=list()
for (i in 1:length(segsgenomewide[[x]])){
htest[[i]]<-t(as.data.frame(segsgenomewide[[x]][[i]][[3]][[y]]))
}
hutest<-data.frame(matrix(unlist(htest), nrow=length(htest), byrow=TRUE),stringsAsFactors=FALSE)
return(hutest)
}


fixedperid_asdf_fun<-function(y){
stats_h=list()
for (i in 1:22){
stats_h[[i]]<-format_fixedperid_fun(i,y)
}
ptest<-do.call(rbind, stats_h)
return(ptest)
}

fixedperid_wl<-fixedperid_asdf_fun(1)
fixedperid_wc<-fixedperid_asdf_fun(2)
fixedperid_el<-fixedperid_asdf_fun(3)
fixedperid_em<-fixedperid_asdf_fun(4)

fixedperid_wl1<-fixedperid_wl[-c(sort(unlist(inf_fsts))),]
fixedperid_wc1<-fixedperid_wc[-c(sort(unlist(inf_fsts))),]
fixedperid_el1<-fixedperid_el[-c(sort(unlist(inf_fsts))),]
fixedperid_em1<-fixedperid_em[-c(sort(unlist(inf_fsts))),]

# 3) sum of all fixed sites divided by the sum of all sites with data
# *1000 for fixed sites per kb

all_fixedperid<-rbind(
(cbind(mean((colSums(fixedperid_wl1)/767730079))*1000,sd((colSums(fixedperid_wl1)/767730079))*1000)),
(cbind(sum(fixedperid_wc1/767730079)*1000,0)),
(cbind(mean((colSums(fixedperid_el1)/767730079))*1000,sd((colSums(fixedperid_el1)/767730079))*1000)),
(cbind(mean((colSums(fixedperid_em1)/767730079))*1000,sd((colSums(fixedperid_em1)/767730079))*1000)))

# all_fixedperid
#          [,1]       [,2]
#[1,] 0.6275814 0.03072747
#[2,] 0.7270941 0.00000000
#[3,] 0.8805440 0.01569050
#[4,] 0.8756855 0.01364289
#-----------------------------------------------------------------------------------------------------------------------

# now combine all into 1 vector with identifers for each

# order of parameters

#het_mu1 # mu het - per kbp & normalised by data coverage of empirical data
#het_sd1 # sd het
#pop_fixedsites # population fixed sites  
#pop_segsites # population segregating sites
#all_fixedperid # fixed sites per individual
#fst.out, # fst
#fout.pi, # pi
#fout.taj # tajima's d


# add in het_sd_WC as 0 here - so that same number of rows as simulated data & can filter after   
target=c(het_mu1,
     het_sd1[1],0,het_sd1[2:3],
pop_fixedsites,
pop_segsites,
all_fixedperid,
sapply(fst_df1, mean),
sapply(pi_df1, mean),
sapply(pi_df1, sd),
sapply(tajima_df1, mean),
sapply(tajima_df1, sd))


names(target)<-c(
"het_mu_WL", "het_mu_WC", "het_mu_EL", "het_mu_EM",
"het_sd_WL", "het_sd_WC","het_sd_EL", "het_sd_EM",
"fixedsites_WL", "fixedsites_WC", "fixedsites_EL", "fixedsites_EM",
"segsites_WL", "segsites_WC", "segsites_EL", "segsites_EM",
"fixedperid_mu_WL", "fixedperid_mu_WC", "fixedperid_mu_EL", "fixedperid_mu_EM",
"fixedperid_sd_WL", "fixedperid_sd_WC", "fixedperid_sd_EL", "fixedperid_sd_EM",
"fst_WL.WC","fst_WL.EL","fst_WL.EM","fst_WC.EL","fst_WC.EM","fst_EL.EM",
"pi_mu_WL", "pi_mu_WC", "pi_mu_EL", "pi_mu_EM",
"pi_sd_WL", "pi_sd_WC", "pi_sd_EL", "pi_sd_EM",
 "tajima_mu_WL", "tajima_mu_EL", "tajima_mu_EM",
  "tajima_sd_WL", "tajima_sd_EL", "tajima_sd_EM")

#length(target)
#[1] 44

target
#       het_mu_WL        het_mu_WC        het_mu_EL        het_mu_EM 
#      0.91084593       0.70084007       0.41157419       0.42950219 
#       het_sd_WL        het_sd_WC        het_sd_EL        het_sd_EM 
#      0.06199505       0.00000000       0.03066836       0.03339686 
#   fixedsites_WL    fixedsites_WC    fixedsites_EL    fixedsites_EM 
#      0.07708829       0.72709409       0.53358597       0.48170967 
#     segsites_WL      segsites_WC      segsites_EL      segsites_EM 
#      3.94090851       0.69363832       1.25447215       1.46983169 
#fixedperid_mu_WL fixedperid_mu_WC fixedperid_mu_EL fixedperid_mu_EM 
#      0.62758139       0.72709409       0.88054400       0.87568548 
#fixedperid_sd_WL fixedperid_sd_WC fixedperid_sd_EL fixedperid_sd_EM 
#      0.03072747       0.00000000       0.01569050       0.01364289 
#       fst_WL.WC        fst_WL.EL        fst_WL.EM        fst_WC.EL 
#      0.13441066       0.18525497       0.17661490       0.36018834 
#       fst_WC.EM        fst_EL.EM         pi_mu_WL         pi_mu_WC 
#      0.33758827       0.20051395       0.06813919       0.00000000 
#        pi_mu_EL         pi_mu_EM         pi_sd_WL         pi_sd_WC 
#      0.03044540       0.03674017       0.02216627       0.00000000 
#        pi_sd_EL         pi_sd_EM     tajima_mu_WL     tajima_mu_EL 
#      0.02362419       0.02474919       0.09748612       0.08333107 
#    tajima_mu_EM     tajima_sd_WL     tajima_sd_EL     tajima_sd_EM 
#      0.33797093       0.46499798       1.00589942       0.95881797 
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

# Fri 21 Jan 2022 17:39:46 GMT
# generate matrix before taking the means/sds * so have full dists, rather than single point values

#-----------------------------------------------------------------------------------------------------------------------

#fst_df1, pi_df1, tajima_df1
# (colSums(fixedperid_wl1)/767730079))*1000 , chekc segsitesperid_em before or after standardisation??

# perhaps try for fst-pi first?
#cor.test(x=test$n_snps, y=test$s_star, method = 'spearman')

# 6 cols
# head(fst_df1)
#          X1         X2         X3        X4        X5
#1 0.07276079 0.09839044 0.08428689 0.1689794 0.1480510
#2 0.15155204 0.12226314 0.11647895 0.3058824 0.2921739

#head(pi_df1)
#          X1 X2          X3         X4
#1 0.07429066  0 0.069850552 0.05608719
#2 0.02442411  0 0.018565047 0.02801888

#cor.test(x=fst_df1$X1, y=pi_df1$X1, method = 'spearman')

#    Spearmans rank correlation rho

#data:  fst_df1$X1 and pi_df1$X1
#S = 3602971088073, p-value < 0.00000000000000022
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#      rho 
#0.3117736 

#Warning message:
#In cor.test.default(x = fst_df1$X1, y = pi_df1$X1, method = "spearman") :
#  Cannot compute exact p-value with ties

# ie loop throguh this

#http://www.sthda.com/english/wiki/ggplot2-quick-correlation-matrix-heatmap-r-software-and-data-visualization
#Compute the correlation matrix
#Correlation matrix can be created using the R function cor() :

#cormat <- round(cor(mydata),2)
#head(cormat)  
#-----------------------------------------------------------------------------------------------------------------------

colnames(fst_df1)<-c("fst_WL.WC","fst_WL.EL","fst_WL.EM","fst_WC.EL","fst_WC.EM","fst_EL.EM")
colnames(pi_df1)<-c("pi_WL", "pi_WC", "pi_EL", "pi_EM")
colnames(tajima_df1)<-c("tajima_WL", "tajima_EL", "tajima_EM")

# target[c(6,22,32,36)]
#       het_sd_WC fixedperid_sd_WC         pi_mu_WC         pi_sd_WC 
#               0                0                0                0 
              # parameters to remove

pi_df2<-pi_df1[,-c(2)] # remove pi_WC


test<-cbind(fst_df1, pi_df2, tajima_df1)
cormat <- round(cor(test),2)
#head(cormat)  
 (cormat)  
          fst_WL.WC fst_WL.EL fst_WL.EM fst_WC.EL fst_WC.EM fst_EL.EM pi_WL
fst_WL.WC      1.00      0.15      0.16      0.63      0.64      0.06  0.30
fst_WL.EL      0.15      1.00      0.63      0.39      0.22      0.19  0.20
fst_WL.EM      0.16      0.63      1.00      0.22      0.40      0.15  0.18
fst_WC.EL      0.63      0.39      0.22      1.00      0.78      0.12  0.08
fst_WC.EM      0.64      0.22      0.40      0.78      1.00      0.11  0.10
fst_EL.EM      0.06      0.19      0.15      0.12      0.11      1.00  0.08
pi_WL          0.30      0.20      0.18      0.08      0.10      0.08  1.00
pi_EL          0.08     -0.18      0.03     -0.14     -0.01     -0.08  0.25
pi_EM          0.10      0.08     -0.16      0.00     -0.11      0.02  0.31
tajima_WL      0.33      0.27      0.26      0.06      0.07      0.07  0.66
tajima_EL      0.03     -0.12      0.00     -0.03      0.01     -0.08  0.04
tajima_EM      0.02      0.02     -0.11     -0.01     -0.01      0.01  0.10
          pi_EL pi_EM tajima_WL tajima_EL tajima_EM
fst_WL.WC  0.08  0.10      0.33      0.03      0.02
fst_WL.EL -0.18  0.08      0.27     -0.12      0.02
fst_WL.EM  0.03 -0.16      0.26      0.00     -0.11
fst_WC.EL -0.14  0.00      0.06     -0.03     -0.01
fst_WC.EM -0.01 -0.11      0.07      0.01     -0.01
fst_EL.EM -0.08  0.02      0.07     -0.08      0.01
pi_WL      0.25  0.31      0.66      0.04      0.10
pi_EL      1.00  0.41      0.17      0.41      0.10
pi_EM      0.41  1.00      0.21      0.08      0.50
tajima_WL  0.17  0.21      1.00      0.03      0.07
tajima_EL  0.41  0.08      0.03      1.00      0.04
tajima_EM  0.10  0.50      0.07      0.04      1.00


# ask re the normalisations?
# first try for fst-pi-tajima

# order of parameters
# take vals of het rahter than mus & sds?

#het_mu1 # mu het - per kbp & normalised by data coverage of empirical data
#het_sd1 # sd het
#pop_fixedsites # population fixed sites  
#pop_segsites # population segregating sites
#all_fixedperid # fixed sites per individual
#fst.out, # fst
#fout.pi, # pi
#fout.taj # tajima's d

#-----------------------------------------------------------------------------------------------------------------------
#head(hetperid_wl1)
#  X1  X2 X3 X4 X5 X6 X7 X8 X9 X10 X11 X12 X13 X14 X15 X16 X17 X18 X19 X20 X21
#1 28  88 38 28 26 29 90 29 31  33  32  33  15  29  16  34   6  24  31   6   4
#2 10  10 14  2 10 12 14  6  8   2  10  11   2  10   2   9   0   3   7   2   0

#head((colSums(hetperid_wl1)/767730079))*1000)

#((colSums(hetperid_wl1)/767730079))*1000
#       X1        X2        X3        X4        X5        X6        X7        X8 
#0.9710926 0.9466569 0.9617469 0.9427675 0.9925598 0.9669193 0.9378921 0.9934403 
#       X9       X10       X11       X12       X13       X14       X15       X16 
#0.9331144 0.8775793 0.9295689 0.9172013 0.8990595 0.9164471 0.8688080 0.9266590 
#      X17       X18       X19       X20       X21       X22       X23       X24 
#0.9414416 0.7258137 0.9175438 0.8851014 0.7531657 0.8303582 0.9271826 0.9321232 
#      X25       X26       X27 
#0.8737954 0.9286975 0.8961040 

#head(((hetperid_wl1)/767730079))*1000
#              X1            X2            X3             X4            X5
#1 0.000036471151 0.00011462362 0.00004949656 0.000036471151 0.00003386607
#2 0.000013025411 0.00001302541 0.00001823558 0.000002605082 0.00001302541


# str((rowSums(hetperid_wl1)/767730079))*1000
# Named num [1:31552] 0.000001196 0.000000233 0.000000513 0.000000831 0.000001344 ...
# - attr(*, "names")= chr [1:31552] "1" "2" "3" "4" ...
#numeric(0)
#-----------------------------------------------------------------------------------------------------------------------

# normalise the heterozygosities
# fpr wl heterozygosity per id? is this right? - take population level
hetperid_wl2<-((rowSums(hetperid_wl1)/767730079))*1000
#hetperid_wc2<-((rowSums(hetperid_wc1)/767730079))*1000
#Error in rowSums(hetperid_wc1) : 
#  'x' must be an array of at least two dimensions
hetperid_wc2<-(((hetperid_wc1)/767730079))*1000
hetperid_el2<-((rowSums(hetperid_el1)/767730079))*1000
hetperid_em2<-((rowSums(hetperid_em1)/767730079))*1000



fixedsitesperid_wl2<-(((fixedsitesperid_wl1)/767730079))*1000
fixedsitesperid_wc2<-(((fixedsitesperid_wc1)/767730079))*1000
fixedsitesperid_el2<-(((fixedsitesperid_el1)/767730079))*1000
fixedsitesperid_em2<-(((fixedsitesperid_em1)/767730079))*1000


segsitesperid_wl2<-(((segsitesperid_wl1)/767730079))*1000
segsitesperid_wc2<-(((segsitesperid_wc1)/767730079))*1000
segsitesperid_el2<-(((segsitesperid_el1)/767730079))*1000
segsitesperid_em2<-(((segsitesperid_em1)/767730079))*1000


fixedperid_wl2<-((rowSums(fixedperid_wl1)/767730079))*1000
fixedperid_wc2<-(((fixedperid_wc1)/767730079))*1000
fixedperid_el2<-((rowSums(fixedperid_el1)/767730079))*1000
fixedperid_em2<-((rowSums(fixedperid_em1)/767730079))*1000

#het_mu1 # mu het - per kbp & normalised by data coverage of empirical data
#het_sd1 # sd het
#pop_fixedsites # population fixed sites  
#pop_segsites # population segregating sites
#all_fixedperid # fixed sites per individual
#fst.out, # fst
#fout.pi, # pi
#fout.taj # tajima's d


test1<-cbind(hetperid_wl2,hetperid_wc2,hetperid_el2,hetperid_em2,
fixedsitesperid_wl2,fixedsitesperid_wc2,fixedsitesperid_el2,fixedsitesperid_em2,
segsitesperid_wl2,segsitesperid_wc2,segsitesperid_el2,segsitesperid_em2,
fixedperid_wl2,fixedperid_wc2,fixedperid_el2,fixedperid_em2
    )

colnames(test1)<-c("het_wl","het_wc","het_el","het_em","fixedsites_wl","fixedsites_wc","fixedsites_el","fixedsites_em", "segsites_wl","segsites_wc","segsites_el","segsites_em","fixedperid_wl","fixedperid_wc","fixedperid_el","fixedperid_em")

test2<-cbind(test1,test)

cormat2 <- round(cor(test2),2)
head(cormat2)  

(cormat2)
              het_wl het_wc het_el het_em fixedsites_wl fixedsites_wc
het_wl          1.00   0.70   0.77   0.79          0.13          0.48
het_wc          0.70   1.00   0.58   0.60          0.09         -0.09
het_el          0.77   0.58   1.00   0.82          0.17          0.25
het_em          0.79   0.60   0.82   1.00          0.18          0.26
fixedsites_wl   0.13   0.09   0.17   0.18          1.00          0.34
fixedsites_wc   0.48  -0.09   0.25   0.26          0.34          1.00
fixedsites_el   0.44   0.27  -0.03   0.14          0.16          0.55
fixedsites_em   0.43   0.26   0.11  -0.01          0.11          0.54
segsites_wl     0.86   0.57   0.56   0.57          0.20          0.68
segsites_wc     0.70   1.00   0.58   0.60          0.09         -0.09
segsites_el     0.73   0.51   0.85   0.70          0.28          0.43
segsites_em     0.75   0.53   0.70   0.82          0.33          0.48
fixedperid_wl   0.55   0.34   0.31   0.32          0.47          0.79
fixedperid_wc   0.48  -0.09   0.25   0.26          0.34          1.00
fixedperid_el   0.60   0.36   0.17   0.29          0.32          0.75
fixedperid_em   0.60   0.36   0.28   0.20          0.32          0.75
fst_WL.WC       0.10  -0.40   0.01   0.01          0.02          0.49
fst_WL.EL       0.05   0.02  -0.15  -0.03          0.13          0.13
fst_WL.EM       0.04   0.01  -0.04  -0.13          0.14          0.13
fst_WC.EL       0.00  -0.45  -0.12  -0.05          0.02          0.43
fst_WC.EM       0.01  -0.45  -0.05  -0.09          0.02          0.44
fst_EL.EM       0.03   0.01  -0.08  -0.02          0.05          0.09
pi_WL           0.43   0.24   0.15   0.17         -0.05          0.38
pi_EL           0.16   0.09   0.41   0.21          0.17          0.17
pi_EM           0.18   0.11   0.21   0.36          0.22          0.21
tajima_WL       0.30   0.17   0.11   0.12          0.01          0.20
tajima_EL       0.03   0.02   0.19   0.05          0.03          0.04
tajima_EM       0.08   0.05   0.07   0.17          0.05          0.08
              fixedsites_el fixedsites_em segsites_wl segsites_wc segsites_el
het_wl                 0.44          0.43        0.86        0.70        0.73
het_wc                 0.27          0.26        0.57        1.00        0.51
het_el                -0.03          0.11        0.56        0.58        0.85
het_em                 0.14         -0.01        0.57        0.60        0.70
fixedsites_wl          0.16          0.11        0.20        0.09        0.28
fixedsites_wc          0.55          0.54        0.68       -0.09        0.43
fixedsites_el          1.00          0.68        0.61        0.27       -0.04
fixedsites_em          0.68          1.00        0.60        0.26        0.20
segsites_wl            0.61          0.60        1.00        0.57        0.68
segsites_wc            0.27          0.26        0.57        1.00        0.51
segsites_el           -0.04          0.20        0.68        0.51        1.00
segsites_em            0.27          0.02        0.73        0.53        0.78
fixedperid_wl          0.67          0.65        0.80        0.34        0.52
fixedperid_wc          0.55          0.54        0.68       -0.09        0.43
fixedperid_el          0.83          0.71        0.82        0.36        0.39
fixedperid_em          0.73          0.79        0.83        0.36        0.48
fst_WL.WC              0.10          0.10        0.10       -0.40        0.06
fst_WL.EL              0.31          0.14        0.06        0.02       -0.12
fst_WL.EM              0.16          0.29        0.05        0.01        0.01
fst_WC.EL              0.22          0.09        0.05       -0.45       -0.13
fst_WC.EM              0.10          0.20        0.05       -0.45       -0.02
fst_EL.EM              0.14          0.07        0.08        0.01       -0.03
pi_WL                  0.38          0.36        0.47        0.24        0.26
pi_EL                 -0.32         -0.04        0.15        0.09        0.58
pi_EM                  0.01         -0.28        0.20        0.11        0.33
tajima_WL              0.19          0.18        0.23        0.17        0.18
tajima_EL             -0.10         -0.01        0.04        0.02        0.15
tajima_EM              0.02         -0.08        0.10        0.05        0.11
              segsites_em fixedperid_wl fixedperid_wc fixedperid_el
het_wl               0.75          0.55          0.48          0.60
het_wc               0.53          0.34         -0.09          0.36
het_el               0.70          0.31          0.25          0.17
het_em               0.82          0.32          0.26          0.29
fixedsites_wl        0.33          0.47          0.34          0.32
fixedsites_wc        0.48          0.79          1.00          0.75
fixedsites_el        0.27          0.67          0.55          0.83
fixedsites_em        0.02          0.65          0.54          0.71
segsites_wl          0.73          0.80          0.68          0.82
segsites_wc          0.53          0.34         -0.09          0.36
segsites_el          0.78          0.52          0.43          0.39
segsites_em          1.00          0.59          0.48          0.54
fixedperid_wl        0.59          1.00          0.79          0.90
fixedperid_wc        0.48          0.79          1.00          0.75
fixedperid_el        0.54          0.90          0.75          1.00
fixedperid_em        0.49          0.91          0.75          0.93
fst_WL.WC            0.07          0.12          0.49          0.15
fst_WL.EL            0.04          0.17          0.13          0.27
fst_WL.EM           -0.10          0.16          0.13          0.18
fst_WC.EL           -0.01          0.08          0.43          0.15
fst_WC.EM           -0.10          0.09          0.44          0.11
fst_EL.EM            0.03          0.11          0.09          0.15
pi_WL                0.29          0.42          0.38          0.49
pi_EL                0.31          0.20          0.17          0.04
pi_EM                0.58          0.27          0.21          0.21
tajima_WL            0.20          0.20          0.20          0.27
tajima_EL            0.07          0.05          0.04         -0.05
tajima_EM            0.19          0.09          0.08          0.09
              fixedperid_em fst_WL.WC fst_WL.EL fst_WL.EM fst_WC.EL fst_WC.EM
het_wl                 0.60      0.10      0.05      0.04      0.00      0.01
het_wc                 0.36     -0.40      0.02      0.01     -0.45     -0.45
het_el                 0.28      0.01     -0.15     -0.04     -0.12     -0.05
het_em                 0.20      0.01     -0.03     -0.13     -0.05     -0.09
fixedsites_wl          0.32      0.02      0.13      0.14      0.02      0.02
fixedsites_wc          0.75      0.49      0.13      0.13      0.43      0.44
fixedsites_el          0.73      0.10      0.31      0.16      0.22      0.10
fixedsites_em          0.79      0.10      0.14      0.29      0.09      0.20
segsites_wl            0.83      0.10      0.06      0.05      0.05      0.05
segsites_wc            0.36     -0.40      0.02      0.01     -0.45     -0.45
segsites_el            0.48      0.06     -0.12      0.01     -0.13     -0.02
segsites_em            0.49      0.07      0.04     -0.10     -0.01     -0.10
fixedperid_wl          0.91      0.12      0.17      0.16      0.08      0.09
fixedperid_wc          0.75      0.49      0.13      0.13      0.43      0.44
fixedperid_el          0.93      0.15      0.27      0.18      0.15      0.11
fixedperid_em          1.00      0.15      0.19      0.25      0.10      0.14
fst_WL.WC              0.15      1.00      0.15      0.16      0.63      0.64
fst_WL.EL              0.19      0.15      1.00      0.63      0.39      0.22
fst_WL.EM              0.25      0.16      0.63      1.00      0.22      0.40
fst_WC.EL              0.10      0.63      0.39      0.22      1.00      0.78
fst_WC.EM              0.14      0.64      0.22      0.40      0.78      1.00
fst_EL.EM              0.12      0.06      0.19      0.15      0.12      0.11
pi_WL                  0.48      0.30      0.20      0.18      0.08      0.10
pi_EL                  0.15      0.08     -0.18      0.03     -0.14     -0.01
pi_EM                  0.13      0.10      0.08     -0.16      0.00     -0.11
tajima_WL              0.27      0.33      0.27      0.26      0.06      0.07
tajima_EL              0.03      0.03     -0.12      0.00     -0.03      0.01
tajima_EM              0.03      0.02      0.02     -0.11     -0.01     -0.01
              fst_EL.EM pi_WL pi_EL pi_EM tajima_WL tajima_EL tajima_EM
het_wl             0.03  0.43  0.16  0.18      0.30      0.03      0.08
het_wc             0.01  0.24  0.09  0.11      0.17      0.02      0.05
het_el            -0.08  0.15  0.41  0.21      0.11      0.19      0.07
het_em            -0.02  0.17  0.21  0.36      0.12      0.05      0.17
fixedsites_wl      0.05 -0.05  0.17  0.22      0.01      0.03      0.05
fixedsites_wc      0.09  0.38  0.17  0.21      0.20      0.04      0.08
fixedsites_el      0.14  0.38 -0.32  0.01      0.19     -0.10      0.02
fixedsites_em      0.07  0.36 -0.04 -0.28      0.18     -0.01     -0.08
segsites_wl        0.08  0.47  0.15  0.20      0.23      0.04      0.10
segsites_wc        0.01  0.24  0.09  0.11      0.17      0.02      0.05
segsites_el       -0.03  0.26  0.58  0.33      0.18      0.15      0.11
segsites_em        0.03  0.29  0.31  0.58      0.20      0.07      0.19
fixedperid_wl      0.11  0.42  0.20  0.27      0.20      0.05      0.09
fixedperid_wc      0.09  0.38  0.17  0.21      0.20      0.04      0.08
fixedperid_el      0.15  0.49  0.04  0.21      0.27     -0.05      0.09
fixedperid_em      0.12  0.48  0.15  0.13      0.27      0.03      0.03
fst_WL.WC          0.06  0.30  0.08  0.10      0.33      0.03      0.02
fst_WL.EL          0.19  0.20 -0.18  0.08      0.27     -0.12      0.02
fst_WL.EM          0.15  0.18  0.03 -0.16      0.26      0.00     -0.11
fst_WC.EL          0.12  0.08 -0.14  0.00      0.06     -0.03     -0.01
fst_WC.EM          0.11  0.10 -0.01 -0.11      0.07      0.01     -0.01
fst_EL.EM          1.00  0.08 -0.08  0.02      0.07     -0.08      0.01
pi_WL              0.08  1.00  0.25  0.31      0.66      0.04      0.10
pi_EL             -0.08  0.25  1.00  0.41      0.17      0.41      0.10
pi_EM              0.02  0.31  0.41  1.00      0.21      0.08      0.50
tajima_WL          0.07  0.66  0.17  0.21      1.00      0.03      0.07
tajima_EL         -0.08  0.04  0.41  0.08      0.03      1.00      0.04
tajima_EM          0.01  0.10  0.10  0.50      0.07      0.04      1.00
#-----------------------------------------------------------------------------------------------------------------------
library(reshape2)
melted_cormat2 <- melt(cormat2)
head(melted_cormat2)
head(melted_cormat2)
           Var1   Var2 value
1        het_wl het_wl  1.00
2        het_wc het_wl  0.70
3        het_el het_wl  0.77
4        het_em het_wl  0.79
5 fixedsites_wl het_wl  0.13
6 fixedsites_wc het_wl  0.48
library(ggplot2)

# Get lower triangle of the correlation matrix
  get_lower_tri<-function(cormat){
    cormat[upper.tri(cormat)] <- NA
    return(cormat)
  }
  # Get upper triangle of the correlation matrix
  get_upper_tri <- function(cormat){
    cormat[lower.tri(cormat)]<- NA
    return(cormat)
  }

upper_tri <- get_upper_tri(cormat2)
upper_tri

#http://www.sthda.com/english/wiki/ggplot2-quick-correlation-matrix-heatmap-r-software-and-data-visualization
#Melt the correlation data and drop the rows with NA values :

# Melt the correlation matrix
#library(reshape2)
melted_cormat2 <- melt(upper_tri, na.rm = TRUE)
# Heatmap
#library(ggplot2)
ggplot(data = melted_cormat2, aes(Var2, Var1, fill = value))+
 geom_tile(color = "white")+
 scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
   midpoint = 0, limit = c(-1,1), space = "Lab", 
   name="Pearson\nCorrelation") +
  theme_minimal()+ 
 theme(axis.text.x = element_text(angle = 45, vjust = 1, 
    size = 12, hjust = 1))+
 coord_fixed()

reorder_cormat <- function(cormat){
# Use correlation between variables as distance
dd <- as.dist((1-cormat)/2)
hc <- hclust(dd)
cormat <-cormat[hc$order, hc$order]
}

# Reorder the correlation matrix
cormat2 <- reorder_cormat(cormat2)
upper_tri <- get_upper_tri(cormat2)
# Melt the correlation matrix
melted_cormat2 <- melt(upper_tri, na.rm = TRUE)
# Create a ggheatmap
ggheatmap <- ggplot(melted_cormat2, aes(Var2, Var1, fill = value))+
 geom_tile(color = "white")+
 scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
   midpoint = 0, limit = c(-1,1), space = "Lab", 
    name="Pearson\nCorrelation") +
  theme_minimal()+ # minimal theme
 theme(axis.text.x = element_text(angle = 45, vjust = 1, 
    size = 12, hjust = 1))+
 coord_fixed()
# Print the heatmap
print(ggheatmap)


pdf("/scratch/devel/hpawar/admix/abc/simul/test/modelcomp/out/corrl.emp.summarystats.24jan22.pdf") 
print(ggheatmap)
dev.off()

#scp -r hpawar@172.16.10.20:/scratch/devel/hpawar/admix/abc/simul/test/modelcomp/out/corrl.emp.summarystats.24jan22.pdf  /Users/harvi/Downloads/gorilla_abc/modelchoice
#g6n.tm4D2L73


#-----------------------------------------------------------------------------------------------------------------------
# correlations between the summary statistics & how these were calculated. I normalised the heterozygosity, pop-wise fixed sites, pop-wise seg sites & fixed sites per individual by data coverage. 

#Mon 24 Jan 2022 11:39:30 GMT
# MK
#that looks interesting. Basically, the fixed sites per pop and ind are
#quite correlated, also with segsites, while hets, fst, tajima and pi are
#less correlated within and across? I would suggest to change the order,
#i.e. the same type of statistics next to each other.
#-----------------------------------------------------------------------------------------------------------------------

# => replot

library(reshape2)
melted_cormat2 <- melt(cormat2)
head(melted_cormat2)
library(ggplot2)

# Get lower triangle of the correlation matrix
  get_lower_tri<-function(cormat){
    cormat[upper.tri(cormat)] <- NA
    return(cormat)
  }
  # Get upper triangle of the correlation matrix
  get_upper_tri <- function(cormat){
    cormat[lower.tri(cormat)]<- NA
    return(cormat)
  }

upper_tri <- get_upper_tri(cormat2)
upper_tri

# Melt the correlation matrix

melted_cormat2 <- melt(upper_tri, na.rm = TRUE)
# Heatmap

ggheatmap2 <- ggplot(data = melted_cormat2, aes(Var2, Var1, fill = value))+
 geom_tile(color = "white")+
 scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
   midpoint = 0, limit = c(-1,1), space = "Lab", 
   name="Pearson\nCorrelation") +
  theme_minimal()+ 
 theme(axis.text.x = element_text(angle = 45, vjust = 1, 
    size = 12, hjust = 1))+
 coord_fixed()

pdf("/scratch/devel/hpawar/admix/abc/simul/test/modelcomp/out/corrl.emp.summarystats.unordered.24jan22.pdf") 
print(ggheatmap2)
dev.off()

#scp -r hpawar@172.16.10.20:/scratch/devel/hpawar/admix/abc/simul/test/modelcomp/out/corrl.emp.summarystats.unordered.24jan22.pdf  /Users/harvi/Downloads/gorilla_abc/modelchoice
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

# Tue 27 Sep 2022 12:23:06 CEST
# replot for the paper - reduce size of x axis test to be more legible
# remove axes labels
# => replot

# change these to be consistent - ie all lower case
rownames(cormat2)[17:28]<-tolower(rownames(cormat2)[17:28])
colnames(cormat2)[17:28]<-tolower(colnames(cormat2)[17:28])


library(reshape2)
melted_cormat2 <- melt(cormat2)
head(melted_cormat2)
library(ggplot2)

# Get lower triangle of the correlation matrix
  get_lower_tri<-function(cormat){
    cormat[upper.tri(cormat)] <- NA
    return(cormat)
  }
  # Get upper triangle of the correlation matrix
  get_upper_tri <- function(cormat){
    cormat[lower.tri(cormat)]<- NA
    return(cormat)
  }

upper_tri <- get_upper_tri(cormat2)
upper_tri

# Melt the correlation matrix

melted_cormat2 <- melt(upper_tri, na.rm = TRUE)
# Heatmap

ggheatmap2 <- ggplot(data = melted_cormat2, aes(Var2, Var1, fill = value))+
 geom_tile(color = "white")+
 scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
   midpoint = 0, limit = c(-1,1), space = "Lab", 
   name="Pearson\nCorrelation") +
  theme_classic() + 
 theme(axis.text.x = element_text(angle = 45, vjust = 1, 
    size = 10, hjust = 1))+
 coord_fixed()

pdf("/scratch/devel/hpawar/admix/abc/simul/test/modelcomp/out/corrl.emp.summarystats.unordered.24jan22.replot.pdf") 
print(ggheatmap2)
dev.off()

#scp -r hpawar@172.16.10.20:/scratch/devel/hpawar/admix/abc/simul/test/modelcomp/out/corrl.emp.summarystats.unordered.24jan22.replot.pdf  /Users/harvi/Downloads/gorilla_abc/modelchoice



