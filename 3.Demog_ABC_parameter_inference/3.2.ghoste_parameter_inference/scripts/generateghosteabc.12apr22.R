#!/usr/bin/env
#module load R/4.0.1
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
# Mon 11 Apr 2022 16:12:50 CEST
# final ABC parameter inference for model with possible ghost gene flow into e_anc

#-----------------------------------------------------------------------------------------------------------------------

# model B)
# simulated data (700 reps) & summary stats for demographic model + ghost gene flow into e_anc
    # /scratch/devel/hpawar/admix/abc/simul/scripts/ghost.model.v1.arr - which calls
    # /scratch/devel/hpawar/admix/abc/simul/scripts/ghost.model.1mar22.R 
# targeting 700*51 = 35,700 reps

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------


# & incorporate processing of new way of calculating tajimas d from generateabc.amendtajima.R  
# amending the below from generateabc.6oct.R
# processing summary statistics for empirical data & simulated data should be identical
    # only diff should be the path to the simulated data (ghost -> e_anc simulations)

# run the below interactively when the cluster is back up

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

# new way of calculating tajimas d -
#/scratch/devel/hpawar/admix/abc/simul/scripts/recalc.tajima.21mar22.R
#/scratch/devel/hpawar/admix/abc/simul/scripts/recalc.tajima.21mar22.arr

# will also need to rerun - for reps with failed iter due to time limits - pending on cluster (Fri  1 Apr 2022 11:45:36 CEST)
# /scratch/devel/hpawar/admix/abc/simul/scripts/recalc.tajima.failedreps.R 

# /Volumes/"Ultra USB 3.0"/IBE/further.analysis.feb.2020/gorillas/abc/generateabc.amendtajima.R - final null abc parameter inference

# 4) format data for abc analysis - run the following interactively
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
library(abc)
options(scipen=100)
require(data.table)
options(stringsAsFactors=F)


# load in the empirical data **

# directly read in target1
load("/scratch/devel/hpawar/admix/abc/results/test/11sep21simns/segrecalc/target1_6oct", verbose=T)
#Loading objects:
#  target1 
#summary stats for empirical data after removing uninformative parameters # het_sd_WC (6), fixedperid_sd_WC (22), pi_mu_WC (32), pi_sd_WC (36)
#target1<-target[-c(6,22,32,36)]


# & amend paths
# mkdir -p /scratch/devel/hpawar/admix/abc/simul/test/ghost/final_1mar22/segrecalc/ # output files dir

#-----------------------------------------------------------------------------------------------------------------------

# simulated data 
simufiles <- paste("/scratch/devel/hpawar/admix/abc/simul/test/ghost/final_1mar22/", list.files(path = "/scratch/devel/hpawar/admix/abc/simul/test/ghost/final_1mar22/", pattern="ghost.abc_sim"), sep = "")



simns=list()
simns_param=list()
for (i in 1:length(simufiles)){
load(file=simufiles[[i]],verbose=T)
simns[[i]]<-simuresults
simns_param[[i]]<-inputsets
}


# check if any simulations failed - none did
# simulations which failed 
probl<-grep("Error", simns)

problsimns=list()
for (i in 1:length(probl)){
problsimns[[i]]<-grep("Error", simns[[probl[i]]])
}

#Error in simns[[probl[i]]] : 
#  attempt to select less than one element in get1index


#-----------------------------------------------------------------------------------------------------------------------
# data structure of simns 700  simns -> 51 iteration each of 5000 windows -> with list of 6 summary stats (where tajima d 6th stat is itself a list)
#> str(simns)
#List of 700
# $ :List of 51
#  ..$ V1 :List of 6
#  .. ..$ : num [1:2, 1:4] 1.522 0.016 0.972 NA 0.182 ...
#  .. ..$ : num [1:2, 1:4] 13023 900209 86539 97198 112373 ...
#  .. ..$ : num [1:2, 1:4] 0.59396 0.00858 0.86539 NA 1.26017 ...
#  .. ..$ : num [1:6, 1] 0.123 0.187 0.175 0.545 0.451 ...
#  .. ..$ : num [1:4, 1:2] 0.1342 0.0848 0.0155 0.0297 0.0217 ...
#  .. ..$ :List of 3
#  .. .. ..$ : num [1, 1:2] -0.83 0.378
#  .. .. ..$ : num [1, 1:2] 0.745 1.412
#  .. .. ..$ : num [1, 1:2] 0.333 0.994
#  ..$ V2 :List of 6
#-----------------------------------------------------------------------------------------------------------------------



#Thu 30 Sep 2021 11:59:08 CEST
#STREAMLINED VERSION OF FUNCTIONS -TEST WHEN CLUSTER IS BACK UP - works
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
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

# all simulation reps for ghost -> e_anc have generated => amend what is being output

# try 
#fst_process_fun<-function(afun,afun1) {
#out_1_700<-afun(1,700)
#return(out_1_700)
#}

# check the output format should be a df - otherwise need to specify in function

#str(s_fst)
#List of 700
# $ :'data.frame':  51 obs. of  6 variables:
    # output is fine - lists of dfs - check if this was the case before - or if had already converted to 1 df?

#fst_process_fun<-function(afun,afun1) {
#out_1_700<-afun(1,700)
#rout_1_700<-as.data.frame(do.call(rbind,out_1_700))
#return(rout_1_700)
#}

#tmp<-fst_process_fun(out_fst_fun1,format_problsimns_fun)
#str(tmp)
#'data.frame':   35700 obs. of  6 variables:
# $ X1: num  0.154 0.164 0.149 0.175 0.162 ...
    # i think this is rather the way to go

#comb_process_fun<-function(afun,afun1,a,b) {
#out_1_700<-afun(1,700,a,b)
#return(out_1_700)
#}

#-----------------------------------------------------------------------------------------------------------------------


# all simulation reps for ghost -> e_anc have generated => amend what is being output

fst_process_fun<-function(afun,afun1) {
out_1_700<-afun(1,700)
rout_1_700<-as.data.frame(do.call(rbind,out_1_700))
return(rout_1_700)
}


comb_process_fun<-function(afun,afun1,a,b) {
out_1_700<-afun(1,700,a,b)
rout_1_700<-as.data.frame(do.call(rbind,out_1_700))
return(rout_1_700)
}

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------


# apply as follows:

# fst
s_fst<-fst_process_fun(out_fst_fun1,format_problsimns_fun)

#str(s_fst)
#'data.frame':   35673 obs. of  6 variables:
# $ X1: num  0.132 0.121 0.184 0.34 0.3 ...
# $ X2: num  0.174 0.223 0.224 0.327 0.352 ...


# pi_mu
s_pi_mu<-comb_process_fun(out_pi_fun1,pi_format_problsimns_fun,5,1)

# pi_sd
s_pi_sd<-comb_process_fun(out_pi_fun1,pi_format_problsimns_fun,5,2)

# tajima_mu - updated way of calculating, process below
#s_tajima_mu<-comb_process_fun(out_pi_fun1,pi_format_problsimns_fun,6,1)
#s_tajima_mu[is.na(s_tajima_mu)] = 0 

# pi_sd
#s_tajima_sd<-comb_process_fun(out_pi_fun1,pi_format_problsimns_fun,6,2)
#s_tajima_sd[is.na(s_tajima_sd)] = 0 



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
#-----------------------------------------------------------------------------------------------------------------------
# process tajimas d stat
# bind in order means WL,EL,EM then sds

# also incorporate if any failed reps **

process_tajima_fun<-function(x) {
test=list()
for (i in 1:length(simns[[x]])){
test[[i]]<-cbind(
simns[[x]][[i]][[6]][[1]][,1],
simns[[x]][[i]][[6]][[2]][,1],
simns[[x]][[i]][[6]][[3]][,1],
simns[[x]][[i]][[6]][[1]][,2],
simns[[x]][[i]][[6]][[2]][,2],
simns[[x]][[i]][[6]][[3]][,2])
}
test_df<-data.frame(matrix(unlist(test), nrow=length(test), byrow=TRUE),stringsAsFactors=FALSE)
}

out_10=list()
for (i in 1:length(simns)){
# first way of calculating - rm nas
out_10[[i]]<-process_tajima_fun(i)
}

tajima_df<-do.call(rbind, out_10)
#-----------------------------------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------------------------------

sumstat<-cbind(s_het_mu,s_het_sd,
s_popfix_perkb,s_popseg_perkb,
s_ifix_mu,s_ifix_sd,
s_fst,
s_pi_mu,
s_pi_sd,
tajima_df)

#> ncol(sumstat)
#[1] 44
#> nrow(sumstat)
#[1] 35700

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


#for (i in 1:length(probl)){
#pn_out[[probl[[i]]]]<-pn_out[[probl[[i]]]][-c(problsimns[[i]])]
#}
#Error in probl[[i]] : subscript out of bounds


pn_format=list()
for (i in 1:length(pn_out)){
pn_format[[i]]<-data.frame(matrix(unlist(pn_out[[i]]), nrow=length(pn_out[[i]]), byrow=TRUE),stringsAsFactors=FALSE)
}

param_df<-as.data.frame(do.call(rbind, pn_format))

nrow(param_df)
#[1] 35700
ncol(param_df)
#[1] 9


param<-data.matrix(param_df)
#parameter vector for each simulation vector, i.e. the input values you used to create each simulation (the "param"

# amend parameters - these are those being sampled from priors

paraname=c("ghost_t0", "e_moun_t3.1", "t4", "t_archintrog", "archaicintrog", "admix_w_e_t6", "t9","gor_ghost_anc","id")

colnames(param)<-paraname
head(param)

# remove the id parameter - ie now parameter 9
param<-param[,c(1:8)]

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
# Q - whether to log transform any of the parameters? no b/c now all are on similar orders of magnitude
    # b/c have also normalised seg sites
#-----------------------------------------------------------------------------------------------------------------------

# rm uninformative parameters
# het_sd_WC (6)
#fixedperid_sd_WC (22)
# pi_mu_WC (32)
# pi_sd_WC (36)

#target1<-target[-c(6,22,32,36)]

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

 #str(target1)
 #Named num [1:40] 0.911 0.701 0.412 0.43 0.062 ...
 #- attr(*, "names")= chr [1:40] "het_mu_WL" "het_mu_WC" "het_mu_EL" "het_mu_EM" ...


# remove same cols from the simulated vals
sumstat1<-sumstat[,-c(6,22,32,36)]

# ncol(sumstat1) #[1] 40

#-----------------------------------------------------------------------------------------------------------------------


# remove ghost_t0 from the parameter search
param1<-param[,-c(1)]
head(param1)

#     e_moun_t3.1        t4 t_archintrog archaicintrog admix_w_e_t6       t9
#[1,]   12.203244 0.1698861    0.7710174      49.16162     47.91701 42.88590
#[2,]   19.080218 0.1926829    5.1176096      54.82516     82.82042 24.98297

#     gor_ghost_anc
#[1,]      30.22583
#[2,]      54.61848


# should also force these posteriors to be within the priors, same logit transf
# amend re priors used in /scratch/devel/hpawar/admix/abc/simul/scripts/ghost.model.1mar22.R 

#lowerp<-c(0.1,2.9191,0.1446,0.2445,0,6.7955,15,10)
#upperp<-c(100,19.1638,0.2420,14.3616,100,99.2159,50,100)
# without ghost_t0 priors as follows
lowerp<-c(2.9191,0.1446,0.2445,0,6.7955,15,10)
upperp<-c(19.1638,0.2420,14.3616,100,99.2159,50,100)

priordf<-cbind(lowerp,upperp)
prior_ranges<-data.matrix(priordf)

# Which you feed into the ABC like this:
ghostabc_12apr22<-abc(target=target1,param=param1,tol=0.005,sumstat=sumstat1,method="neuralnet",numnet=100,transf="logit",logit.bounds=prior_ranges)
# The "tol" parameter is unknown and might be tested a bit.

> ghostabc_12apr22<-abc(target=target1,param=param1,tol=0.005,sumstat=sumstat1,method="neuralnet",numnet=100,transf="logit",logit.bounds=prior_ranges)
123456789101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899100
123456789101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899100
Warning message:
All parameters are "logit" transformed. 

#-----------------------------------------------------------------------------------------------------------------------

ghostabc_12apr22[[17]][[1]]<-colnames(param1)
ghostabc_12apr22[[17]][[2]]<-names(target1)

ghostabc_12apr22
summary(ghostabc_12apr22)

save(ghostabc_12apr22,file=paste("/scratch/devel/hpawar/admix/abc/simul/test/ghost/final_1mar22/segrecalc/ghosteabc_12apr22"))
save(param1,file=paste("/scratch/devel/hpawar/admix/abc/simul/test/ghost/final_1mar22/segrecalc/ghosteabc_param_12apr22"))
save(sumstat1,file=paste("/scratch/devel/hpawar/admix/abc/simul/test/ghost/final_1mar22/segrecalc/ghosteabc_sumstat1_12apr22"))


#scp -r hpawar@172.16.10.20:/scratch/devel/hpawar/admix/abc/simul/test/ghost/final_1mar22/segrecalc/ghosteabc_12apr22 /Users/harvi/Downloads/gorilla_abc/ghost/final_12apr22
#scp -r hpawar@172.16.10.20:/scratch/devel/hpawar/admix/abc/simul/test/ghost/final_1mar22/segrecalc/ghosteabc_param_12apr22  /Users/harvi/Downloads/gorilla_abc/ghost/final_12apr22

#-----------------------------------------------------------------------------------------------------------------------

> ghostabc_12apr22
Call:
abc(target = target1, param = param1, sumstat = sumstat1, tol = 0.005, 
    method = "neuralnet", transf = "logit", logit.bounds = prior_ranges, 
    numnet = 100)
Method:
Non-linear regression via neural networks
with correction for heteroscedasticity

Parameters:
e_moun_t3.1, t4, t_archintrog, archaicintrog, admix_w_e_t6, t9, gor_ghost_anc

Statistics:
het_mu_WL, het_mu_WC, het_mu_EL, het_mu_EM, het_sd_WL, het_sd_EL, het_sd_EM, fixedsites_WL, fixedsites_WC, fixedsites_EL, fixedsites_EM, segsites_WL, segsites_WC, segsites_EL, segsites_EM, fixedperid_mu_WL, fixedperid_mu_WC, fixedperid_mu_EL, fixedperid_mu_EM, fixedperid_sd_WL, fixedperid_sd_EL, fixedperid_sd_EM, fst_WL.WC, fst_WL.EL, fst_WL.EM, fst_WC.EL, fst_WC.EM, fst_EL.EM, pi_mu_WL, pi_mu_EL, pi_mu_EM, pi_sd_WL, pi_sd_EL, pi_sd_EM, tajima_mu_WL, tajima_mu_EL, tajima_mu_EM, tajima_sd_WL, tajima_sd_EL, tajima_sd_EM

Total number of simulations 35700 

Number of accepted simulations:  179 

> summary(ghostabc_12apr22)
Call: 
abc(target = target1, param = param1, sumstat = sumstat1, tol = 0.005, 
    method = "neuralnet", transf = "logit", logit.bounds = prior_ranges, 
    numnet = 100)
Data:
 abc.out$adj.values (179 posterior samples)
Weights:
 abc.out$weights

                       e_moun_t3.1      t4 t_archintrog archaicintrog
Min.:                       2.9264  0.1784       0.2469       88.0123
Weighted 2.5 % Perc.:       2.9523  0.1872       0.2874       95.5370
Weighted Median:            3.0090  0.1980       0.5037       99.0813
Weighted Mean:              3.0271  0.1977       0.5806       98.7352
Weighted Mode:              2.9878  0.1941       0.4223       99.5246
Weighted 97.5 % Perc.:      3.1943  0.2081       1.4211       99.9014
Max.:                       4.3775  0.2201       2.3397       99.9303
                       admix_w_e_t6      t9 gor_ghost_anc
Min.:                        8.0753 29.8582       10.6156
Weighted 2.5 % Perc.:        8.8122 39.2035       18.8026
Weighted Median:            11.0722 45.0179       35.3814
Weighted Mean:              11.8285 44.9005       39.5084
Weighted Mode:              10.6170 45.3048       28.2786
Weighted 97.5 % Perc.:      17.3850 49.6024       80.6825
Max.:                       47.8210 49.9911       93.8139
#-----------------------------------------------------------------------------------------------------------------------


# on local
setwd("/Users/harvi/Downloads/gorilla_abc/ghost/final_12apr22")
require(data.table)

# A) tol 0.005
load("/Users/harvi/Downloads/gorilla_abc/ghost/final_12apr22/ghosteabc_12apr22",verbose=T)
load("/Users/harvi/Downloads/gorilla_abc/ghost/final_12apr22/ghosteabc_param_12apr22",verbose=T)
library(abc)
# histogram of the weighted posterior sample  : https://cran.r-project.org/web/packages/abc/vignettes/abcvignette.pdf
pdf("/Users/harvi/Downloads/gorilla_abc/ghost/final_12apr22/ghosteabc_posterior_12apr22.pdf") 
hist(ghostabc_12apr22)
dev.off() 
# generate the diagnostic plots of the estimation of the posterior distribution 
# a density plot of the prior distribution, a density plot of the posterior distribution estimated with and without regression- based correction, a scatter plot of the Euclidean distances as a function of the parameter values, and a normal Q-Q plot of the residuals from the regression
pdf("/Users/harvi/Downloads/gorilla_abc/ghost/final_12apr22/ghosteabc_diagnostic_12apr22.pdf") 
plot(ghostabc_12apr22, param1)
dev.off() 
