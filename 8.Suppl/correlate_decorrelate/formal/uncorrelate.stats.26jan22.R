#module load R/4.0.1

#install.packages('pls')
# install.packages('MASS')
library("pls"); library("MASS");


# ~35,700 ms simulations for null demog

# processing from generateabc.6oct.R
options(scipen=100)
require(data.table)
options(stringsAsFactors=F)

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

# read in the simulations & process to a matrix

# simulated data 

# new batch of simulations. - output has path
#/scratch/devel/hpawar/admix/abc/simul/test/11sep21/abc_sim*

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

#head(simns[[1]][[1]])

# simulations which failed 
probl<-grep("Error", simns)

problsimns=list()
for (i in 1:length(probl)){
problsimns[[i]]<-grep("Error", simns[[probl[i]]])
}


#Thu 30 Sep 2021 11:59:08 CEST
#STREAMLINED VERSION OF FUNCTIONS -TEST WHEN CLUSTER IS BACK UP - works
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
# MK - b) Another question concerns the normalization of the "full_segsites" object part for the simulated values. 
#As far as I see, you save the raw numbers of fixed sites, and then process interactively.
#But why is the normalization (s_popseg/(700*40000))*1000). Shouldnt it be (s_popseg/(2500*40000))*1000). 
#If you have 2500 simulations, that should be the factor, or am I missing something here? 
#I just checked how this would influence the segsites:

## ~Median of normalized value   
#> c(949874,144163, 356160, 381739)/(2500*40000)*1000
#[1] 9.49874 1.44163 3.56160 3.81739
## empirical values
#segsites_WL      segsites_WC      segsites_EL      segsites_EM
#3.94090851       0.69363832       1.25447215       1.46983169
#Thats not so far off anymore.

# =>  I was mixing up number of simulation reps (the 700) with number of windows generated (the 2500). It should be 2500.
# => change factor in division from 700 -> 2500
#-----------------------------------------------------------------------------------------------------------------------

#  popn-wise fixed sites
s_popfix<-comb_process_fun(out_het_fun1,het_format_problsimns_fun,2,1)
s_popfix_perkb<-((s_popfix/(2500*40000))*1000)

#  popn-wise seg sites
s_popseg<-comb_process_fun(out_het_fun1,het_format_problsimns_fun,2,2)
s_popseg_perkb<-((s_popseg/(2500*40000))*1000)

# when cluster is working again - test if this works ** yes it does

#-----------------------------------------------------------------------------------------------------------------------
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

#ncol(sumstat) # [1] 44
#nrow(sumstat) # [1] 35543

sumstat<-data.matrix(sumstat)
# matching simulation vectors, each of them in the same order as summary stats from the real data (the "sumstat")

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

# format input parameters used to generate the simulations generating these summary stats

#head(simns_param)
#length(simns_param[[1]])
#[1] 51

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

#torm<- c(3,  9, 15, 21, 27, 33, 39, 45, 51)
#pn_out[[862]]<-pn_out[[862]][-c(torm)]

# head(pn_out[[probl[[1]]]])
#pn_out[[11]]<-pn_out[[11]][-c(problsimns[[1]])]
#pn_out[[probl[[1]]]]<-pn_out[[probl[[1]]]][-c(problsimns[[1]])]


for (i in 1:length(probl)){
pn_out[[probl[[i]]]]<-pn_out[[probl[[i]]]][-c(problsimns[[i]])]
}

# length(pn_out[[11]])
#[1] 42
# length(pn_out[[483]])
#[1] 35
#length(pn_out)
#[1] 700


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

# check paranames - which have been fixed
paraname=c("w_lowl_t0","w_cros_t0","e_lowl_t0","e_moun_t0","e_lowl_t1","e_lowl_t2","e_moun_t3","e_moun_t3.1","t4","e_anc_t4","t5","w_lowl_t5","admix_w_e_t6","admix_e_w_t6","t7","w_anc_t7","t8","gor_anc","id")
  
colnames(param)<-paraname
head(param)

# remove the id parameter - ie now parameter 19
param<-param[,c(1:18)]

#-----------------------------------------------------------------------------------------------------------------------
# remove same cols from the simulated vals
sumstat1<-sumstat[,-c(6,22,32,36)]

#-----------------------------------------------------------------------------------------------------------------------
# str(sumstat1)
# num [1:35543, 1:40] 3.11 2.79 2.08 3.32 1.14 ...
# - attr(*, "dimnames")=List of 2
#  ..$ : NULL
#  ..$ : chr [1:40] "X1" "X2" "X3" "X4" ...


# str(param)
# num [1:35543, 1:18] 97.6 32.3 77.7 50.3 52.7 ...
# - attr(*, "dimnames")=List of 2
#  ..$ : NULL
#  ..$ : chr [1:18] "w_lowl_t0" "w_cros_t0" "e_lowl_t0" "e_moun_t0" ...

# box cox transf - so that it follows normal dist

# try this -
#transform statistics via boxcox
#for(i in 1:length(sumstat1)){
#d<-cbind(sumstat1[,i], param);
#mylm<-lm(as.formula(d), data=d);
#myboxcox<-boxcox(mylm, lambda=seq(-20,100,1/10), interp=T, eps=1/50);
#lambda<-c(lambda, myboxcox$x[myboxcox$y==max(myboxcox$y)]);
#myGM<-c(myGM, mean(exp(log(stat[,i]))));
#}

#Error in formula.default(object, env = baseenv()) : invalid formula

# perhaps need to go through the abc toolbox & go through the - 
#example files are provided in #the subfolder exampleFiles.

#head(sumstat1[,1])
#[1] 3.111270 2.789345 2.081357 3.323228 1.135773 2.993201

#> head(d)
#              w_lowl_t0 w_cros_t0  e_lowl_t0 e_moun_t0  e_lowl_t1 e_lowl_t2
#[1,] 3.111270  97.62212 13.575184  0.4760922  8.095278  4.0822297 14.408516
#[2,] 2.789345  32.32839 17.651724 18.0069939  6.395847  0.6318954 24.965467
#[3,] 2.081357  77.70069 17.269579  8.9075054  2.079896 14.3745150 17.306423
#[4,] 3.323228  50.29358 14.867581  8.0067277  3.308055  6.4509626 22.000146
#[5,] 1.135773  52.66651  5.971104 10.3816306 17.430553  3.2258722  9.993514
#[6,] 2.993201  84.14413  4.094294 15.7704928  1.259206 14.8080218 19.159895
#     e_moun_t3 e_moun_t3.1        t4  e_anc_t4        t5 w_lowl_t5 admix_w_e_t6
#[1,] 1.4252746    7.797648 0.1660231 27.869704 0.4238028  17.15203     35.19760
#[2,] 1.8327422   13.121675 0.1598330  2.676607 0.5642337  13.62836     15.98596
#[3,] 1.4892850    2.687306 0.2424358 24.153495 0.5937883  20.38505     17.99188
#[4,] 4.8182251   11.703166 0.1524235  4.309535 0.5210987  23.46663     68.90734
#[5,] 4.6947332    7.895361 0.2211435  5.277714 0.5474891  12.82225     54.85539
#[6,] 0.1810329   19.912959 0.2053703 15.554648 0.5574161  33.51546     69.88029
#     admix_e_w_t6        t7 w_anc_t7        t8  gor_anc
#[1,]    43.424713 2.0520676 88.39160  9.923953 67.38882
#[2,]     1.467758 4.2312926 77.11491  9.419132 89.50979
#[3,]    54.456309 0.7896675 61.95525  9.405389 35.38746
#[4,]    47.979495 2.2901779 92.95131 14.393179 67.81708
#[5,]    15.084633 5.8375626 50.99454 13.769801 27.42961
#[6,]    62.567449 5.0597520 98.40054 11.807262 61.99037
#> mylm<-lm(as.formula(d), data=d);
#Error in formula.default(object, env = baseenv()) : invalid formula

# this is the problem b/c multiple y cols (for each of the parameters..)

# abcsampler: perform the Box-Cox transformation prior to the linear transformation on the statistics.  http://cmpg.unibe.ch/software/ABCtoolbox/ABCtoolbox_manual.pdf


# whether the param - shoudl instead be the summary stats names? no doesnt seem so
# paramnames<-c(
#"het_mu_WL", "het_mu_WC", "het_mu_EL", "het_mu_EM",
#"het_sd_WL", "het_sd_WC","het_sd_EL", "het_sd_EM",
#"fixedsites_WL", "fixedsites_WC", "fixedsites_EL", "fixedsites_EM",
#"segsites_WL", "segsites_WC", "segsites_EL", "segsites_EM",
#"fixedperid_mu_WL", "fixedperid_mu_WC", "fixedperid_mu_EL", "fixedperid_mu_EM",
#"fixedperid_sd_WL", "fixedperid_sd_WC", "fixedperid_sd_EL", "fixedperid_sd_EM",
#"fst_WL.WC","fst_WL.EL","fst_WL.EM","fst_WC.EL","fst_WC.EM","fst_EL.EM",
#"pi_mu_WL", "pi_mu_WC", "pi_mu_EL", "pi_mu_EM",
#"pi_sd_WL", "pi_sd_WC", "pi_sd_EL", "pi_sd_EM",
# "tajima_mu_WL", "tajima_mu_EL", "tajima_mu_EM",
#  "tajima_sd_WL", "tajima_sd_EL", "tajima_sd_EM")

#parmanames1<-paramnames[-c(6,22,32,36)]


#for(i in 1:length(sumstat1)){
#d<-cbind(sumstat1[,i], parmanames1);
#mylm<-lm(as.formula(d), data=d);
#myboxcox<-boxcox(mylm, lambda=seq(-20,100,1/10), interp=T, eps=1/50);
#lambda<-c(lambda, myboxcox$x[myboxcox$y==max(myboxcox$y)]);
#myGM<-c(myGM, mean(exp(log(stat[,i]))));
#}

#Error in formula.character(object, env = baseenv()) : 
#  invalid formula structure(c("3.11127", "2.78934518518519", "2.08135666666667", "3.32322814814815", "1.13577296296296", "2.99320148148148", "1.43993259259259", "2.90965962962963", "2.95969481481481", "3.19606592592593", "1.97996592592593", "2.71406407407407", "2.91530407407407", "2.45313888888889", "1.46124962962963", "1.91106777777778", "3.8835562962963", "1.92335481481481", "1.14133518518519", "2.56094740740741", "3.47167962962963", "1.81943925925926", "1.76071444444444", "1.36255333333333", "2.52441851851852",  "2.51319814814815", "4.08436518518519", "2.13900148148148", "4.29908518518519", "1.89109407407407", "0.921878518518518", "0.835276666666667", "2.11150481481481", "3.19447296296296", "2.0802162962963", "3.44183444444444", "2.01377851851852", "3.13584074074074", "1.12241444444444", "1.51562962962963", "1.32499111111111", "1.58606851851852", "3.47022481481481", "1.57026703703704", "1.56563074074074", "2.73508962962963", "2.57429925925926", "3.08389", "3.05360074074074",
#In addition: Warning messages:
#1: In cbind(sumstat1[, i], parmanames1) :
#  number of rows of result is not a multiple of vector length (arg 2)
#2: Using formula(x) is deprecated when x is a character vector of length > 1.
#  Consider formula(paste(x, collapse = " ")) instead. 

#-----------------------------------------------------------------------------------------------------------------------

#Mon 14 Feb 2022 15:17:43 CET
# MK
#- the error message seems to be due to the input not being data.frames
#- they show a previous step to force values to be between 1 and 2, which is indeed a weird scaling factor
#- applying the below should result in a table of boxcox transformed parameters, 
#then the rest of the stuff also seems to work - at least technically, you have to see how the output looks like


#
##force stat in [1,2]
myMax<-c(); myMin<-c(); lambda<-c(); myGM<-c();
for(i in 1:ncol(sumstat1)){
  myMax<-c(myMax, max(sumstat1[,i])); myMin<-c(myMin, min(sumstat1[,i]));
  sumstat1[,i]<-1+(sumstat1[,i] -myMin[i])/(myMax[i]-myMin[ i]);
  print(c(range(sumstat1[,i]),length(which(sumstat1[,i]==NA)),mean(sumstat1[,i])))
}

#transform statistics via boxcox
for(i in 1:ncol(sumstat1)){
print(i)
d<-cbind(data.frame(sumstat1[,i]), data.frame(param));
mylm<-lm(as.formula(d), data=d);
myboxcox<-boxcox(mylm, lambda=seq(-20,100,1/10), interp=T, eps=1/50);
lambda<-c(lambda, myboxcox$x[myboxcox$y==max(myboxcox$y)]);
myGM<-c(myGM, mean(exp(log(sumstat[,i]))));
}

#Warning messages:
#1: In log(sumstat[, i]) : NaNs produced
#2: In log(sumstat[, i]) : NaNs produced

sstat<-sumstat1
#standardize the BC-stat
myBCMeans<-c(); myBCSDs<-c();
for(i in 1:ncol(sstat)){
    sstat[,i]<-(sstat[,i] ^lambda[i] - 1)/(lambda[i]*myGM[i] ^(lambda[i]-1));
    myBCSDs<-c(myBCSDs, sd(sstat[,i]));
    myBCMeans<-c(myBCMeans, mean(sstat[,i]));
    sstat[,i]<-(sstat[,i] -myBCMeans[i])/myBCSDs[i];
  }
#perform pls
 
 # myPlsr<-plsr(as.matrix(param)~as.matrix(sstat), scale=F, validation='LOO');

#Error in matrix(0, ncol = ncomp, nrow = npred) : 
#  invalid 'ncol' value (< 0)


#> ncol(param)
#[1] 18
#> nrow(param)
#[1] 35543
#> ncol(sstat)
#[1] 40
#> nrow(sstat)
#[1] 35543

# maybe remove the nan cols?
#head(sstat) # multiple cols with only NaNs
#head(sstat[,c(1:5,7:21,23:38)])

# cols 6, 22,39,40 # only contain NaN
# check which stats these correspond to **

test<-sstat[,c(1:5,7:21,23:38)]

# perform pls
myPlsr<-plsr(as.matrix(param)~as.matrix(test), scale=F, validation='LOO');
# is taking a long time to run..
  # may need to send as a job.. yes

#write pls to a file
myPlsrDataFrame<-data.frame(comp1=myPlsr$loadings[,1]);
for(i in 2:numComp){myPlsrDataFrame<-cbind(myPlsrDataFrame, myPlsr$loadings[,i]);}
write.table(cbind(names(stats), myMax, myMin, lambda, myGM, myBCMeans, myBCSDs, myPlsrDataFrame), file="/scratch/devel/hpawar/admix/abc/simul/test/modelcomp/out/PLSfile_sstat.nullabc.15feb22.txt", col.names=F,row.names=F, sep=’’, quote=F); 

#make RMSEP plot
pdf("/scratch/devel/hpawar/admix/abc/simul/test/modelcomp/out/RMSEP_sstat.nullabc.15feb22.pdf");
plot(RMSEP(myPlsr));
dev.off();

#"/scratch/devel/hpawar/admix/abc/simul/test/modelcomp/out/PLSfile_sstat.nullabc.15feb22.txt"
#"/scratch/devel/hpawar/admix/abc/simul/test/modelcomp/out/RMSEP_sstat.nullabc.15feb22.pdf"

q()
#-----------------------------------------------------------------------------------------------------------------------


# https://www.statology.org/box-cox-transformation-in-r/
#box-cox transformation is a commonly used method for transforming a non-normally distributed dataset into a more normally distributed one.
# boxcox() function from the MASS() library.

#find optimal lambda for Box-Cox transformation 
#bc <- boxcox(y ~ x)
# but here what is the x & y??

#library(MASS)

#create data
#y=c(1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 6, 7, 8)
#x=c(7, 7, 8, 3, 2, 4, 4, 6, 6, 7, 5, 3, 3, 5, 8)

#fit linear regression model
#model <- lm(y~x)

#find optimal lambda for Box-Cox transformation 
#bc <- boxcox(y ~ x)
#(lambda <- bc$x[which.max(bc$y)])

#[1] -0.4242424

#fit new linear regression model using the Box-Cox transformation
n#ew_model <- lm(((y^lambda-1)/lambda) ~ x)


#inearize the relationship
#between model parameters and statistics before
#transforming them linearly in some cases