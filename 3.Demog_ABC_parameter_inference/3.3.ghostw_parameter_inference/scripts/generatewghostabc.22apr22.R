#!/usr/bin/env
#module load R/4.0.1
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

# Fri 22 Apr 2022 12:26:33 CEST
# final ABC parameter inference for model with possible ghost gene flow into e_anc

# amend below from #/Volumes/"Ultra USB 3.0"/IBE/further.analysis.feb.2020/gorillas/abc/simul_postabc/abc+ghost/generatewghostabc.22apr22.R

#-----------------------------------------------------------------------------------------------------------------------

# model C)
# /scratch/devel/hpawar/admix/abc/simul/scripts/wanc.ghost.model.1apr22.R  # called by -
#/scratch/devel/hpawar/admix/abc/simul/scripts/wanc.ghost.model.v1.arr

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
# mkdir -p /scratch/devel/hpawar/admix/abc/simul/test/ghost/final_1mar22/wanc/segrecalc/ # output files dir

#-----------------------------------------------------------------------------------------------------------------------

# simulated data  for ghostw
simufiles <- paste("/scratch/devel/hpawar/admix/abc/simul/test/ghost/final_1mar22/wanc/", list.files(path = "/scratch/devel/hpawar/admix/abc/simul/test/ghost/final_1mar22/wanc/", pattern="wanc.ghost.abc_sim"), sep = "")



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


# should also force these posteriors to be within the priors, same logit transf
# amend re priors used in /scratch/devel/hpawar/admix/abc/simul/scripts/ghost.model.1mar22.R 

#lowerp<-c(0.1,2.9191,0.1446,0.2445,0,6.7955,15,10)
#upperp<-c(100,19.1638,0.2420,14.3616,100,99.2159,50,100)
# without ghost_t0 priors as follows

# AMEND WITH PRIORS FROM /scratch/devel/hpawar/admix/abc/simul/scripts/wanc.ghost.model.1apr22.R **
# only priors which differ b/n ghoste & ghostw is t_archintrog=runif(1000,min=5.9803,max=14.3616) # lower prior here
lowerp<-c(2.9191,0.1446,5.9803,0,6.7955,15,10)
upperp<-c(19.1638,0.2420,14.3616,100,99.2159,50,100)

priordf<-cbind(lowerp,upperp)
prior_ranges<-data.matrix(priordf)

# Which you feed into the ABC like this:
wghostabc_22apr22<-abc(target=target1,param=param1,tol=0.005,sumstat=sumstat1,method="neuralnet",numnet=100,transf="logit",logit.bounds=prior_ranges)
# The "tol" parameter is unknown and might be tested a bit.

> wghostabc_22apr22<-abc(target=target1,param=param1,tol=0.005,sumstat=sumstat1,method="neuralnet",numnet=100,transf="logit",logit.bounds=prior_ranges)
123456789101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899100
123456789101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899100
Warning message:
All parameters are "logit" transformed. 

#-----------------------------------------------------------------------------------------------------------------------

wghostabc_22apr22[[17]][[1]]<-colnames(param1)
wghostabc_22apr22[[17]][[2]]<-names(target1)

wghostabc_22apr22
summary(wghostabc_22apr22)

# CHANGE PATHS **
#mkdir -p /scratch/devel/hpawar/admix/abc/simul/test/ghost/final_1mar22/wanc/segrecalc/

save(wghostabc_22apr22,file=paste("/scratch/devel/hpawar/admix/abc/simul/test/ghost/final_1mar22/wanc/segrecalc/wghostabc_22apr22"))
save(param1,file=paste("/scratch/devel/hpawar/admix/abc/simul/test/ghost/final_1mar22/wanc/segrecalc/wghostabc_param_22apr22"))
save(sumstat1,file=paste("/scratch/devel/hpawar/admix/abc/simul/test/ghost/final_1mar22/wanc/segrecalc/wghostabc_sumstat1_22apr22"))


#scp -r hpawar@172.16.10.20:/scratch/devel/hpawar/admix/abc/simul/test/ghost/final_1mar22/wanc/segrecalc/wghostabc_22apr22 /Users/harvi/Downloads/gorilla_abc/ghost/final_12apr22
#scp -r hpawar@172.16.10.20:/scratch/devel/hpawar/admix/abc/simul/test/ghost/final_1mar22/wanc/segrecalc/wghostabc_param_22apr22  /Users/harvi/Downloads/gorilla_abc/ghost/final_12apr22

#-----------------------------------------------------------------------------------------------------------------------
> wghostabc_22apr22
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

> summary(wghostabc_22apr22)
Call: 
abc(target = target1, param = param1, sumstat = sumstat1, tol = 0.005, 
    method = "neuralnet", transf = "logit", logit.bounds = prior_ranges, 
    numnet = 100)
Data:
 abc.out$adj.values (179 posterior samples)
Weights:
 abc.out$weights

                       e_moun_t3.1      t4 t_archintrog archaicintrog
Min.:                       3.6242  0.2259      13.9104        1.5083
Weighted 2.5 % Perc.:       5.7861  0.2289      13.9999        3.7328
Weighted Median:            9.8111  0.2345      14.0514        7.0962
Weighted Mean:             10.3369  0.2341      14.0488        7.9336
Weighted Mode:              9.0623  0.2347      14.0410        6.2082
Weighted 97.5 % Perc.:     17.6685  0.2373      14.0865       19.3478
Max.:                      18.9875  0.2399      14.1974       46.9077
                       admix_w_e_t6      t9 gor_ghost_anc
Min.:                       21.4342 17.5696       52.6003
Weighted 2.5 % Perc.:       31.4315 19.3692       63.5606
Weighted Median:            80.3770 24.6563       79.0217
Weighted Mean:              76.4927 25.6181       79.4316
Weighted Mode:              89.9786 22.7944       77.9864
Weighted 97.5 % Perc.:      98.3049 39.5711       93.3204
Max.:                       99.1477 45.7853       96.0408

#-----------------------------------------------------------------------------------------------------------------------


# on local
setwd("/Users/harvi/Downloads/gorilla_abc/ghost/final_12apr22")
require(data.table)

# A) tol 0.005
load("/Users/harvi/Downloads/gorilla_abc/ghost/final_12apr22/wghostabc_22apr22",verbose=T)
load("/Users/harvi/Downloads/gorilla_abc/ghost/final_12apr22/wghostabc_param_22apr22",verbose=T)
library(abc)
# histogram of the weighted posterior sample  : https://cran.r-project.org/web/packages/abc/vignettes/abcvignette.pdf
pdf("/Users/harvi/Downloads/gorilla_abc/ghost/final_12apr22/wghostabc_posterior_22apr22.pdf") 
hist(wghostabc_22apr22)
dev.off() 
# generate the diagnostic plots of the estimation of the posterior distribution 
# a density plot of the prior distribution, a density plot of the posterior distribution estimated with and without regression- based correction, a scatter plot of the Euclidean distances as a function of the parameter values, and a normal Q-Q plot of the residuals from the regression
pdf("/Users/harvi/Downloads/gorilla_abc/ghost/final_12apr22/wghostabc_diagnostic_22apr22.pdf") 
plot(wghostabc_22apr22, param1)
dev.off() 