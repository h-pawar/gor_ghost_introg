#!/usr/bin/env
#module load R/4.0.1
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
#MK: If they are NaN (not a number), they are due to division by zero (I guess actually often 0/0) and might be corrected to 0 instead. 
find.nas_fun1<-function(x) {
for (j in 1:length(genomewide[[x]])){
genomewide[[x]][[j]][[1]]<-as.numeric(genomewide[[x]][[j]][[1]])
}
out_1<-list()
for (j in 1:length(genomewide[[x]])){
out_1[[j]]<-lapply(genomewide[[x]][[j]], function(y) replace(y, !is.finite(y), 0))
}
return(out_1)
}

finalwindows=list()
for (i in 1:length(genomewide)){
finalwindows[[i]]<-find.nas_fun1(i)
}

#-----------------------------------------------------------------------------------------------------------------------

format_stats_fun<-function(x,y) {
otest=list()
for (i in 1:length(finalwindows[[x]])){
otest[[i]]<-t(as.data.frame(finalwindows[[x]][[i]][[y]]))
}
outest<-data.frame(matrix(unlist(otest), nrow=length(otest), byrow=TRUE),stringsAsFactors=FALSE)
return(outest)
}


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



#-----------------------------------------------------------------------------------------------------------------------


fst_df1<-fst_df[-c(sort(unlist(inf_fsts))),]

# remove the same rows from the other dfs (those with infinite fst vals) **

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


# normalise by data coverage

pop_fixedsites<-rbind(sum(fixedsitesperid_wl1/767730079)*1000,
sum(fixedsitesperid_wc1/767730079)*1000,
sum(fixedsitesperid_el1/767730079)*1000,
sum(fixedsitesperid_em1/767730079)*1000)


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

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

# 2.2) fixed sites per individual

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
#-----------------------------------------------------------------------------------------------------------------------


# simulated data 

simufiles <- paste("/scratch/devel/hpawar/admix/abc/simul/test/11sep21/", list.files(path = "/scratch/devel/hpawar/admix/abc/simul/test/11sep21/", pattern="abc_sim"), sep = "")

simns=list()
simns_param=list()
for (i in 1:length(simufiles)){
load(file=simufiles[[i]],verbose=T)
simns[[i]]<-simuresults
simns_param[[i]]<-inputsets
}


#> head(str(simns))
#List of 700 # 700  simns -> 51 iteration each of 5000 windows -> with list of 6 summary stats
# $ :List of 51
#  ..$ V13 :List of 6

# simulations which failed 
probl<-grep("Error", simns)

problsimns=list()
for (i in 1:length(probl)){
problsimns[[i]]<-grep("Error", simns[[probl[i]]])
}


#Thu 30 Sep 2021 11:59:08 CEST
#STREAMLINED VERSION OF FUNCTIONS 
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
# individual functions


# fst functions
simu_fst_stats_fun<-function(x) {
test=list()
for (i in 1:length(simns[[x]])){
test[[i]]<-(simns[[x]][[i]][[4]])
}
simu_fst<-data.frame(matrix(unlist(test), nrow=length(test), byrow=TRUE),stringsAsFactors=FALSE)
}

out_fst_fun1<-function(x0,x) {
stats_s1=list()
for (i in x0:x){
stats_s1[[i]]<-simu_fst_stats_fun(i)
}
return(stats_s1)
}

format_problsimns_fun<-function(x) {
y<-probl[[x]]
hold_problout=list()
for(i in 1:length(simns[[y]]) ){
  if( (i %in% problsimns[[x]]) ){ next }
 hold_problout[[i]]<-(simns[[y]][[i]][[4]]) 
}
t<-hold_problout[lengths(hold_problout) != 0]
test_11<-data.frame(matrix(unlist(t), nrow=length(t), byrow=TRUE),stringsAsFactors=FALSE)
return(test_11)
}

#-----------------------------------------------------------------------------------------------------------------------
# pi/tajima functions
simu_pt_stats_fun<-function(x,y,z) {
test=list()
for (i in 1:length(simns[[x]])){
test[[i]]<-(simns[[x]][[i]][[y]][,z])
}
simu_pi<-data.frame(matrix(unlist(test), nrow=length(test), byrow=TRUE),stringsAsFactors=FALSE)
}

out_pi_fun1<-function(x0,x,a,b) {
stats_s1=list()
for (i in x0:x){
stats_s1[[i]]<-simu_pt_stats_fun(i,a,b)
}
return(stats_s1)
}

pi_format_problsimns_fun<-function(x,a,z) {
y<-probl[[x]]
hold_problout=list()
for(i in 1:length(simns[[y]]) ){
  if( (i %in% problsimns[[x]]) ){ next }
 hold_problout[[i]]<-(simns[[y]][[i]][[a]][,z]) 
}
t<-hold_problout[lengths(hold_problout) != 0]
test_11<-data.frame(matrix(unlist(t), nrow=length(t), byrow=TRUE),stringsAsFactors=FALSE)
return(test_11)
}

# a=5 for pi, a=6 for tajima

# pi mu, pi sd, tajima mu, tajima sd
# pi mu 5,1
# pi sd 5,2
# taj mu 6,1
# taj sd 6,2

#-----------------------------------------------------------------------------------------------------------------------

# heterozygosity/fixed sites per id/popn fixed sites/popn seg sites functions

simu_het_stats_fun<-function(x,y,z) {
test=list()
for (i in 1:length(simns[[x]])){
test[[i]]<-(simns[[x]][[i]][[y]][z,])
}
simu_het<-data.frame(matrix(unlist(test), nrow=length(test), byrow=TRUE),stringsAsFactors=FALSE)
}

out_het_fun1<-function(x0,x,a,b) {
stats_s1=list()
for (i in x0:x){
stats_s1[[i]]<-simu_het_stats_fun(i,a,b)
}
return(stats_s1)
}

het_format_problsimns_fun<-function(x,a,z) {
y<-probl[[x]]
hold_problout=list()
for(i in 1:length(simns[[y]]) ){
  if( (i %in% problsimns[[x]]) ){ next }
 hold_problout[[i]]<-(simns[[y]][[i]][[a]][z,]) 
}
t<-hold_problout[lengths(hold_problout) != 0]
test_11<-data.frame(matrix(unlist(t), nrow=length(t), byrow=TRUE),stringsAsFactors=FALSE)
return(test_11)
}

# het mu 1,1
# het sd 1,2

# ifix mu 3,1
# ifix sd 3,2

# popfix #2,1
# popseg #2,2


#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

fst_process_fun<-function(afun,afun1) {

out_1_10<-afun(1,10)
#11
out_12_45<-afun(12,45)
#46
out_47_102<-afun(47,102)
#103
out_104_241<-afun(104,241)
#242
out_243_246<-afun(243,246)
#247
out_248_293<-afun(248,293)
#294,295
out_296_426<-afun(296,426)
#427
out_428_465<-afun(428,465)
#466
out_467_475<-afun(467,475)
#476
out_477_482<-afun(477,482)
#483
out_484_509<-afun(484,509)
#510
out_511_563<-afun(511,563)
#564
out_565_672<-afun(565,672)
#673
out_674_700<-afun(674,700)



# apply for all 14 simn reps which failed
out_rest=list()
for (i in 1:length(probl)){
out_rest[[i]]<-afun1(i)
}


out_12_45<-out_12_45[lengths(out_12_45) != 0]
out_47_102<-out_47_102[lengths(out_47_102) != 0]
out_104_241<-out_104_241[lengths(out_104_241) != 0]
out_243_246<-out_243_246[lengths(out_243_246) != 0]
out_248_293<-out_248_293[lengths(out_248_293) != 0]
out_296_426<-out_296_426[lengths(out_296_426) != 0]
out_428_465<-out_428_465[lengths(out_428_465) != 0]
out_467_475<-out_467_475[lengths(out_467_475) != 0]
out_477_482<-out_477_482[lengths(out_477_482) != 0]
out_484_509<-out_484_509[lengths(out_484_509) != 0]
out_511_563<-out_511_563[lengths(out_511_563) != 0]
out_565_672<-out_565_672[lengths(out_565_672) != 0]
out_674_700<-out_674_700[lengths(out_674_700) != 0]


test2<-do.call(c, list(out_1_10,list(out_rest[[1]]),
out_12_45,
list(out_rest[[2]]),
out_47_102,
list(out_rest[[3]]),
out_104_241,
list(out_rest[[4]]),
out_243_246,
list(out_rest[[5]]),
out_248_293,
list(out_rest[[6]]),
list(out_rest[[7]]),
out_296_426,
list(out_rest[[8]]),
out_428_465,
list(out_rest[[9]]),
out_467_475,
list(out_rest[[10]]),
out_477_482,
list(out_rest[[11]]),
out_484_509,
list(out_rest[[12]]),
out_511_563,
list(out_rest[[13]]),
out_565_672,
list(out_rest[[14]]),
out_674_700))

out_1_700<-as.data.frame(do.call(rbind,test2)) # out df for the 700 simns

return(out_1_700)
}
#-----------------------------------------------------------------------------------------------------------------------
comb_process_fun<-function(afun,afun1,a,b) {

out_1_10<-afun(1,10,a,b)
#11
out_12_45<-afun(12,45,a,b)
#46
out_47_102<-afun(47,102,a,b)
#103
out_104_241<-afun(104,241,a,b)
#242
out_243_246<-afun(243,246,a,b)
#247
out_248_293<-afun(248,293,a,b)
#294,295
out_296_426<-afun(296,426,a,b)
#427
out_428_465<-afun(428,465,a,b)
#466
out_467_475<-afun(467,475,a,b)
#476
out_477_482<-afun(477,482,a,b)
#483
out_484_509<-afun(484,509,a,b)
#510
out_511_563<-afun(511,563,a,b)
#564
out_565_672<-afun(565,672,a,b)
#673
out_674_700<-afun(674,700,a,b)



# apply for all 14 simn reps which failed
out_rest=list()
for (i in 1:length(probl)){
out_rest[[i]]<-afun1(i,a,b)
}


out_12_45<-out_12_45[lengths(out_12_45) != 0]
out_47_102<-out_47_102[lengths(out_47_102) != 0]
out_104_241<-out_104_241[lengths(out_104_241) != 0]
out_243_246<-out_243_246[lengths(out_243_246) != 0]
out_248_293<-out_248_293[lengths(out_248_293) != 0]
out_296_426<-out_296_426[lengths(out_296_426) != 0]
out_428_465<-out_428_465[lengths(out_428_465) != 0]
out_467_475<-out_467_475[lengths(out_467_475) != 0]
out_477_482<-out_477_482[lengths(out_477_482) != 0]
out_484_509<-out_484_509[lengths(out_484_509) != 0]
out_511_563<-out_511_563[lengths(out_511_563) != 0]
out_565_672<-out_565_672[lengths(out_565_672) != 0]
out_674_700<-out_674_700[lengths(out_674_700) != 0]


test2<-do.call(c, list(out_1_10,list(out_rest[[1]]),
out_12_45,
list(out_rest[[2]]),
out_47_102,
list(out_rest[[3]]),
out_104_241,
list(out_rest[[4]]),
out_243_246,
list(out_rest[[5]]),
out_248_293,
list(out_rest[[6]]),
list(out_rest[[7]]),
out_296_426,
list(out_rest[[8]]),
out_428_465,
list(out_rest[[9]]),
out_467_475,
list(out_rest[[10]]),
out_477_482,
list(out_rest[[11]]),
out_484_509,
list(out_rest[[12]]),
out_511_563,
list(out_rest[[13]]),
out_565_672,
list(out_rest[[14]]),
out_674_700))

out_1_700<-as.data.frame(do.call(rbind,test2)) # out df for the 700 simns

return(out_1_700)
}

#-----------------------------------------------------------------------------------------------------------------------

# apply as follows:

# fst
s_fst<-fst_process_fun(out_fst_fun1,format_problsimns_fun)

# pi_mu
s_pi_mu<-comb_process_fun(out_pi_fun1,pi_format_problsimns_fun,5,1)

# pi_sd
s_pi_sd<-comb_process_fun(out_pi_fun1,pi_format_problsimns_fun,5,2)

# tajima_mu
s_tajima_mu<-comb_process_fun(out_pi_fun1,pi_format_problsimns_fun,6,1)
s_tajima_mu[is.na(s_tajima_mu)] = 0 

# pi_sd
s_tajima_sd<-comb_process_fun(out_pi_fun1,pi_format_problsimns_fun,6,2)
s_tajima_sd[is.na(s_tajima_sd)] = 0 



# het mu
s_het_mu<-comb_process_fun(out_het_fun1,het_format_problsimns_fun,1,1)

# het sd
s_het_sd<-comb_process_fun(out_het_fun1,het_format_problsimns_fun,1,2)
s_het_sd[is.na(s_het_sd)] = 0

# fixed sites per id mu
s_ifix_mu<-comb_process_fun(out_het_fun1,het_format_problsimns_fun,3,1)

# fixed sites per id sd
s_ifix_sd<-comb_process_fun(out_het_fun1,het_format_problsimns_fun,3,2)
s_ifix_sd[is.na(s_ifix_sd)] = 0

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

#  popn-wise fixed sites
s_popfix<-comb_process_fun(out_het_fun1,het_format_problsimns_fun,2,1)
s_popfix_perkb<-((s_popfix/(2500*40000))*1000)

#  popn-wise seg sites
s_popseg<-comb_process_fun(out_het_fun1,het_format_problsimns_fun,2,2)
s_popseg_perkb<-((s_popseg/(2500*40000))*1000)

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

sumstat<-cbind(s_het_mu,s_het_sd,
s_popfix_perkb,s_popseg_perkb,
s_ifix_mu,s_ifix_sd,
s_fst,
s_pi_mu,
s_pi_sd,
s_tajima_mu,
s_tajima_sd)


sumstat<-data.matrix(sumstat)
# matching simulation vectors, each of them in the same order as summary stats from the real data (the "sumstat")

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

# format input parameters used to generate the simulations generating these summary stats

parameters_fun<-function(x) {
pn=list()
for (i in 1:length(simns_param[[x]])){
pn[[i]]<-simns_param[[x]][[i]]
}
return(pn)
}

pn_out=list()
for (i in 1:length(simns_param)){
pn_out[[i]]<-parameters_fun(i)
}


#  remove those iterations of simulations which failed & did not output summary stats

for (i in 1:length(probl)){
pn_out[[probl[[i]]]]<-pn_out[[probl[[i]]]][-c(problsimns[[i]])]
}

pn_format=list()
for (i in 1:length(pn_out)){
pn_format[[i]]<-data.frame(matrix(unlist(pn_out[[i]]), nrow=length(pn_out[[i]]), byrow=TRUE),stringsAsFactors=FALSE)
}

param_df<-as.data.frame(do.call(rbind, pn_format))

nrow(param_df) #[1] 35543
ncol(param_df) #[1] 19

# nrow(param_df)
#[1] 101991

#ncol(param_df)
#[1] 23

param<-data.matrix(param_df)
#parameter vector for each simulation vector, i.e. the input values you used to create each simulation (the "param"

paraname=c("w_lowl_t0","w_cros_t0","e_lowl_t0","e_moun_t0","e_lowl_t1","e_lowl_t2","e_moun_t3","e_moun_t3.1","t4","e_anc_t4","t5","w_lowl_t5","admix_w_e_t6","admix_e_w_t6","t7","w_anc_t7","t8","gor_anc","id")
  
colnames(param)<-paraname
head(param)

# remove the id parameter - ie now parameter 19
param<-param[,c(1:18)]

#-----------------------------------------------------------------------------------------------------------------------

# rm uninformative parameters
# het_sd_WC (6)
#fixedperid_sd_WC (22)
# pi_mu_WC (32)
# pi_sd_WC (36)

target1<-target[-c(6,22,32,36)]

# target1
# het_mu_WL        het_mu_WC        het_mu_EL        het_mu_EM 
#      0.91084593       0.70084007       0.41157419       0.42950219 
#       het_sd_WL        het_sd_EL        het_sd_EM    fixedsites_WL 
#      0.06199505       0.03066836       0.03339686       0.07708829 
#   fixedsites_WC    fixedsites_EL    fixedsites_EM      segsites_WL 
#      0.72709409       0.53358597       0.48170967       3.94090851 
#     segsites_WC      segsites_EL      segsites_EM fixedperid_mu_WL 
#      0.69363832       1.25447215       1.46983169       0.62758139 
#fixedperid_mu_WC fixedperid_mu_EL fixedperid_mu_EM fixedperid_sd_WL 
#      0.72709409       0.88054400       0.87568548       0.03072747 
#fixedperid_sd_EL fixedperid_sd_EM        fst_WL.WC        fst_WL.EL 
#      0.01569050       0.01364289       0.13441066       0.18525497 
#       fst_WL.EM        fst_WC.EL        fst_WC.EM        fst_EL.EM 
#      0.17661490       0.36018834       0.33758827       0.20051395 
#        pi_mu_WL         pi_mu_EL         pi_mu_EM         pi_sd_WL 
#      0.06813919       0.03044540       0.03674017       0.02216627 
#        pi_sd_EL         pi_sd_EM     tajima_mu_WL     tajima_mu_EL 
#      0.02362419       0.02474919       0.09748612       0.08333107 
#    tajima_mu_EM     tajima_sd_WL     tajima_sd_EL     tajima_sd_EM 
#      0.33797093       0.46499798       1.00589942       0.95881797 

# length(target1) #[1] 40

# remove same cols from the simulated vals
sumstat1<-sumstat[,-c(6,22,32,36)]

# ncol(sumstat1) #[1] 40

save(param,file=paste("/scratch/devel/hpawar/admix/abc/results/test/11sep21simns/segrecalc/param_6oct"))
save(sumstat1,file=paste("/scratch/devel/hpawar/admix/abc/results/test/11sep21simns/segrecalc/sumstat1_6oct"))
save(target1,file=paste("/scratch/devel/hpawar/admix/abc/results/test/11sep21simns/segrecalc/target1_6oct"))
