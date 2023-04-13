#!/usr/bin/env

# Mon 27 Feb 2023 10:35:17 CET

# revised ghost models: sampling all parameters from priors 

# 10,000 reps, model comparison simulations in paths - 
# /scratch/devel/hpawar/admix/abc/simul/test/ghost/revisions_8feb23/modelcomp/ghoste/rev.mc.ghoste_sim*
# /scratch/devel/hpawar/admix/abc/simul/test/ghost/revisions_8feb23/modelcomp/ghostw/rev.mc.ghostw_sim* 



#-----------------------------------------------------------------------------------------------------------------------


#module load R/4.0.1

#Mon 25 Apr 2022 14:16:16 CEST
# have generated 10,000 reps for nullplghost & ghoste & ghostw - using final posteriors

# dir with original mc simulated data -  
# /scratch/devel/hpawar/admix/abc/simul/test/modelcomp/final_12apr22/nullplusghost/
# /scratch/devel/hpawar/admix/abc/simul/test/modelcomp/final_12apr22/ghoste/
# /scratch/devel/hpawar/admix/abc/simul/test/modelcomp/final_12apr22/ghostw

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------


# Thu 16 Dec 2021 14:38:17 CET
# have generated 50*50=2500 reps for 
    # A.1) null + ghost
    # B) ghost -> e_anc

# reps generated with -
    # removing sites fixed across all gor ids, before calculating the summary stats - number of population-wise fixed and segregating sites & the number of fixed sites per individual


# paths to scripts
# model A.1) null + ghost

# 1) null pl ghost
#/scratch/devel/hpawar/admix/abc/simul/scripts/modelcomp/mc.nullplghost.12apr22.arr
#/scratch/devel/hpawar/admix/abc/simul/scripts/modelcomp/mc.nullplghost.12apr22.R

# 2) ghost -> e_anc
#/scratch/devel/hpawar/admix/abc/simul/scripts/modelcomp/mc.ghoste.12apr22.arr
#/scratch/devel/hpawar/admix/abc/simul/scripts/modelcomp/mc.ghoste.12apr22.R


#-----------------------------------------------------------------------------------------------------------------------

# following process of /Volumes/"Ultra USB 3.0"/IBE/further.analysis.feb.2020/gorillas/abc/simul_postabc/abc+ghost/modelcomp/modelchoice.null.9dec21.R

#-----------------------------------------------------------------------------------------------------------------------
# Fri  3 Dec 2021 12:20:18 CET
# process summary statistics from 'model comparison simulations'
    # these are simulations generated taking the parameter vals as the weighted median posteriors inferred from ABC analysis
# 1) process summary stats in same way as in previous generateabc scripts
# 2) perform model choice
    # - cross validation
    # - misclassification
#-----------------------------------------------------------------------------------------------------------------------


# following procedure of
# /Volumes/"Ultra USB 3.0"/IBE/further.analysis.feb.2020/gorillas/abc/simul_postabc/abc+ghost/wanc.generateghostabc.22nov.R

# run the following interactively
#-----------------------------------------------------------------------------------------------------------------------
library(abc)
options(scipen=100)
require(data.table)
options(stringsAsFactors=F)

# 1) functions to process the summary stats
# 2) read in each of the simulated data sets (data generated under each of the 3 models)
    # A - null demography
    # B - ghost gene flow to e_anc

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


# all simulation reps for ghost -> e_anc have generated => amend what is being output

fst_process_fun<-function(afun,afun1) {
#out_1_200<-afun(1,200)
# reducing to 1,10 - Thu  9 Dec 2021 15:18:13 CET **
# now 1,50 - Thu 16 Dec 2021 14:45:06 CET **
#out_1_200<-afun(1,50)
# increase back to 200 - Wed  2 Feb 2022 10:57:03 CET
out_1_200<-afun(1,200)
rout_1_200<-as.data.frame(do.call(rbind,out_1_200))
return(rout_1_200)
}


comb_process_fun<-function(afun,afun1,a,b) {
out_1_200<-afun(1,200,a,b)
#out_1_200<-afun(1,50,a,b)
rout_1_200<-as.data.frame(do.call(rbind,out_1_200))
return(rout_1_200)
}

#-----------------------------------------------------------------------------------------------------------------------



matrix_function<-function(){
#-----------------------------------------------------------------------------------------------------------------------
# apply as follows:

# fst
s_fst<-fst_process_fun(out_fst_fun1,format_problsimns_fun)

# pi_mu
s_pi_mu<-comb_process_fun(out_pi_fun1,pi_format_problsimns_fun,5,1)

# pi_sd
s_pi_sd<-comb_process_fun(out_pi_fun1,pi_format_problsimns_fun,5,2)

# tajima_mu - updated way of calculating, process below
#s_tajima_mu<-comb_process_fun(out_pi_fun1,pi_format_problsimns_fun,6,1)
#s_tajima_mu[is.na(s_tajima_mu)] = 0 

# tajima_sd
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
# NOTE - change division here ** 
# b/c mc simulations have generated 250 windows per iter (rather than 2500)

#  popn-wise fixed sites
s_popfix<-comb_process_fun(out_het_fun1,het_format_problsimns_fun,2,1)
s_popfix_perkb<-((s_popfix/(250*40000))*1000)


#  popn-wise seg sites
s_popseg<-comb_process_fun(out_het_fun1,het_format_problsimns_fun,2,2)
s_popseg_perkb<-((s_popseg/(250*40000))*1000)

# when cluster is working again - test if this works ** yes it does

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


#ncol(sumstat) # [1] 44
# nrow(sumstat)
#[1] 35200

sumstat<-data.matrix(sumstat)
# matching simulation vectors, each of them in the same order as summary stats from the real data (the "sumstat")

return(sumstat)
}



#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

#Mon 27 Feb 2023 10:50:34 CET
# Read in revised model comparison (mc) simulations & process the summary statistics

# B) revised ghoste - sample all parameters from priors
simufiles <- paste("/scratch/devel/hpawar/admix/abc/simul/test/ghost/revisions_8feb23/modelcomp/ghoste/", list.files(path = "/scratch/devel/hpawar/admix/abc/simul/test/ghost/revisions_8feb23/modelcomp/ghoste/", pattern="rev.mc.ghoste_sim"), sep = "")


simns=list()
simns_param=list()
for (i in 1:length(simufiles)){
load(file=simufiles[[i]],verbose=T)
simns[[i]]<-simuresults
simns_param[[i]]<-inputsets
}

#str(simns)
#str(simns)
#List of 200
# $ :List of 50
#  ..$ V1 :List of 6

# simulations which failed 
#probl<-grep("Error", simns)

#problsimns=list()
#for (i in 1:length(probl)){
#problsimns[[i]]<-grep("Error", simns[[probl[i]]])
#}
#Error in simns[[probl[i]]] : 
#  attempt to select less than one element in get1index
# -> no simulations failed 

rev.ghoste_stats<-matrix_function()

# shoudl assign simns to diff vector - to retain original data for each model
rev.ghoste_simns<-simns
rev.ghoste_param<-simns_param
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
# C) revised ghostw - sample all parameters from priors
simufiles <- paste("/scratch/devel/hpawar/admix/abc/simul/test/ghost/revisions_8feb23/modelcomp/ghostw/", list.files(path = "/scratch/devel/hpawar/admix/abc/simul/test/ghost/revisions_8feb23/modelcomp/ghostw/", pattern="rev.mc.ghostw_sim"), sep = "")


simns=list()
simns_param=list()
for (i in 1:length(simufiles)){
load(file=simufiles[[i]],verbose=T)
simns[[i]]<-simuresults
simns_param[[i]]<-inputsets
}

#str(simns)

#> str(simns)
#List of 200
# $ :List of 50
#  ..$ V1 :List of 6

# simulations which failed 
#probl<-grep("Error", simns)

#problsimns=list()
#for (i in 1:length(probl)){
#problsimns[[i]]<-grep("Error", simns[[probl[i]]])
#}


#Error in simns[[probl[i]]] : 
#  attempt to select less than one element in get1index
# -> no simulations failed 

rev.ghostw_stats<-matrix_function()

# shoudl assign simns to diff vector - to retain original data for each model
rev.ghostw_simns<-simns
rev.ghostw_param<-simns_param

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

# 1) Read in model comparison (mc) simulations & process the summary statistics

# A.2) null + ghost mc simulations - after removing sites fixed in all gor ids (weighted median posterior parameter vals)
simufiles <- paste("/scratch/devel/hpawar/admix/abc/simul/test/modelcomp/final_12apr22/nullplusghost/", list.files(path = "/scratch/devel/hpawar/admix/abc/simul/test/modelcomp/final_12apr22/nullplusghost/", pattern="test.mc.nullplusghost_sim"), sep = "")

simns=list()
simns_param=list()
for (i in 1:length(simufiles)){
load(file=simufiles[[i]],verbose=T)
simns[[i]]<-simuresults
simns_param[[i]]<-inputsets
}

# check if all reps generated
#probl<-grep("Error", simns)#

#problsimns=list()
#for (i in 1:length(probl)){
#problsimns[[i]]<-grep("Error", simns[[probl[i]]])
#}
#Error in simns[[probl[i]]] : 
#  attempt to select less than one element in get1index

nullplusghost_stats<-matrix_function()

# shoudl assign simns to diff vector - to retain original data for each model
nullplusghost_simns<-simns
nullplusghost_param<-simns_param

#-----------------------------------------------------------------------------------------------------------------------
# 1) Read in model comparison (mc) simulations & process the summary statistics

# B) ghoste - setting ne of ghost pop to 25 instead of 0.1
simufiles <- paste("/scratch/devel/hpawar/admix/abc/simul/test/modelcomp/final_12apr22/ghoste/", list.files(path = "/scratch/devel/hpawar/admix/abc/simul/test/modelcomp/final_12apr22/ghoste/", pattern="test.mc.ghoste_sim"), sep = "")


simns=list()
simns_param=list()
for (i in 1:length(simufiles)){
load(file=simufiles[[i]],verbose=T)
simns[[i]]<-simuresults
simns_param[[i]]<-inputsets
}

#str(simns)
#str(simns)
#List of 199
# $ :List of 50
#  ..$ V1 :List of 6
# should be 200 - check which failed
# or just generate 1 extra rep - this may be easiest

# simulations which failed 
#probl<-grep("Error", simns)

#problsimns=list()
#for (i in 1:length(probl)){
#problsimns[[i]]<-grep("Error", simns[[probl[i]]])
#}
#Error in simns[[probl[i]]] : 
#  attempt to select less than one element in get1index
# -> no simulations failed 

orig.ghoste_stats<-matrix_function()

# shoudl assign simns to diff vector - to retain original data for each model
orig.ghoste_simns<-simns
orig.ghoste_param<-simns_param

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
# 1) Read in model comparison (mc) simulations & process the summary statistics

# C) ghostw - setting ne of ghost pop to 25 instead of 0.1
simufiles <- paste("/scratch/devel/hpawar/admix/abc/simul/test/modelcomp/final_12apr22/ghostw/", list.files(path = "/scratch/devel/hpawar/admix/abc/simul/test/modelcomp/final_12apr22/ghostw/", pattern="test.mc.ghostw_sim"), sep = "")


simns=list()
simns_param=list()
for (i in 1:length(simufiles)){
load(file=simufiles[[i]],verbose=T)
simns[[i]]<-simuresults
simns_param[[i]]<-inputsets
}

#str(simns)
#str(simns)
#List of 199
# $ :List of 50
#  ..$ V1 :List of 6
# should be 200 - check which failed
# or just generate 1 extra rep - this may be easiest

# simulations which failed 
#probl<-grep("Error", simns)

#problsimns=list()
#for (i in 1:length(probl)){
#problsimns[[i]]<-grep("Error", simns[[probl[i]]])
#}
#Error in simns[[probl[i]]] : 
#  attempt to select less than one element in get1index
# -> no simulations failed 

orig.ghostw_stats<-matrix_function()

# shoudl assign simns to diff vector - to retain original data for each model
orig.ghostw_simns<-simns
orig.ghostw_param<-simns_param

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------


# 1.2) remove uninformative stats (same as removed for previous ABC simulations) - remove from each of the models
# het_sd_WC (6) #fixedperid_sd_WC (22) # pi_mu_WC (32)# pi_sd_WC (36)

nullplusghost_stats1<-nullplusghost_stats[,-c(6,22,32,36)]
orig.ghoste_stats1<-orig.ghoste_stats[,-c(6,22,32,36)]
orig.ghostw_stats1<-orig.ghostw_stats[,-c(6,22,32,36)]
rev.ghoste_stats1<-rev.ghoste_stats[,-c(6,22,32,36)]
rev.ghostw_stats1<-rev.ghostw_stats[,-c(6,22,32,36)]



#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

# 2) perform model choice

# 1) compare : null vs revised ghostE vs revised ghostW

# 2) compare :  null vs revised ghostE vs revised ghostW vs orig ghostE vs orig ghostW

#-----------------------------------------------------------------------------------------------------------------------

# 1) compare : null vs revised ghostE vs revised ghostW

# ie combine stats into dfs

stats3mod<-rbind(as.data.frame(nullplusghost_stats1),as.data.frame(rev.ghoste_stats1),as.data.frame(rev.ghostw_stats1))

#models<-c(
#replicate(2500, "nullplusghost"),
#replicate(2500, "ghoste")
#    )

# increased to 200*50
models<-c(
replicate(10000, "nullplusghost"),
replicate(10000, "rev.ghoste"),
replicate(10000, "rev.ghostw")
    )
#-----------------------------------------------------------------------------------------------------------------------

#  2.1) cross-validation for model selection
# if ABC can distinguish between the models

#-----------------------------------------------------------------------------------------------------------------------

# this was run interactively

cv.modsel <- cv4postpr(models, stats3mod, nval=100, tol=0.05, method="neuralnet")

s <- summary(cv.modsel)

s

Confusion matrix based on 100 samples for each model.

$tol0.05
              nullplusghost rev.ghoste rev.ghostw
nullplusghost           100          0          0
rev.ghoste                0        100          0
rev.ghostw                0          0        100


Mean model posterior probabilities (neuralnet)

$tol0.05
              nullplusghost rev.ghoste rev.ghostw
nullplusghost             1          0          0
rev.ghoste                0          1          0
rev.ghostw                0          0          1

#-----------------------------------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------------------------------

# 2.2) model selection with postpr function
# calculate posterior probabilities of each demog model

# read in summary stats from empirical data
load("/scratch/devel/hpawar/admix/abc/simul/test/ghost/3nov21/segrecalc/ghosttarget1_22nov",verbose=T)
#Loading objects:
 # target1

# tolerance rate of 0.05%
testmod<-postpr(target1, models, stats3mod, tol=.05, method="neuralnet")


Warning message:
There are 3 models but only 1 for which simulations have been accepted.
No regression is performed, method is set to rejection.
Consider increasing the tolerance rate.TRUE 

# function summary prints out posterior model probabilities and ratios of model probabilities (the Bayes factors) in a user-friendly way
summary(testmod)

#-----------------------------------------------------------------------------------------------------------------------
> summary(testmod)

> summary(testmod)way
Call: 
postpr(target = target1, index = models, sumstat = stats3mod, 
    tol = 0.05, method = "neuralnet")
Data:
 postpr.out$values (1500 posterior samples)
Models a priori:
 nullplusghost, rev.ghoste, rev.ghostw
Models a posteriori:
 nullplusghost, rev.ghoste, rev.ghostw

Proportion of accepted simulations (rejection):
nullplusghost    rev.ghoste    rev.ghostw 
            0             0             1 

Bayes factors:
              nullplusghost rev.ghoste rev.ghostw
nullplusghost                                   0
rev.ghoste                                      0
rev.ghostw              Inf        Inf          1



#-----------------------------------------------------------------------------------------------------------------------

# may need to send as a a job, where nval=1000
#cv.modsel <- cv4postpr(models, stats3mod, nval=1000, tol=0.05, method="neuralnet")

#s <- summary(cv.modsel)

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------


# 2) compare :  null vs revised ghostE vs revised ghostW vs orig ghostE vs orig ghostW


stats5mod<-rbind(as.data.frame(nullplusghost_stats1),as.data.frame(orig.ghoste_stats1),as.data.frame(orig.ghostw_stats1),as.data.frame(rev.ghoste_stats1),as.data.frame(rev.ghostw_stats1))

#models<-c(
#replicate(2500, "nullplusghost"),
#replicate(2500, "ghoste")
#    )

# increased to 200*50
models_rev<-c(
replicate(10000, "nullplusghost"),
replicate(10000, "orig.ghoste"),
replicate(10000, "orig.ghostw"),
replicate(10000, "rev.ghoste"),
replicate(10000, "rev.ghostw")
    )


#  2.1) cross-validation for model selection

cv.modsel <- cv4postpr(models_rev, stats5mod, nval=100, tol=0.05, method="neuralnet")

# is hanging a bit here..

s <- summary(cv.modsel)

Confusion matrix based on 100 samples for each model.

$tol0.05
              nullplusghost orig.ghoste orig.ghostw rev.ghoste rev.ghostw
nullplusghost           100           0           0          0          0
orig.ghoste               0         100           0          0          0
orig.ghostw               0           0         100          0          0
rev.ghoste                0           0           0        100          0
rev.ghostw                0           0           0          0        100


Mean model posterior probabilities (neuralnet)

$tol0.05
              nullplusghost orig.ghoste orig.ghostw rev.ghoste rev.ghostw
nullplusghost             1           0           0          0          0
orig.ghoste               0           1           0          0          0
orig.ghostw               0           0           1          0          0
rev.ghoste                0           0           0          1          0
rev.ghostw                0           0           0          0          1




# 2.2) model selection with postpr function
# calculate posterior probabilities of each demog model

# read in summary stats from empirical data
load("/scratch/devel/hpawar/admix/abc/simul/test/ghost/3nov21/segrecalc/ghosttarget1_22nov",verbose=T)
#Loading objects:
 # target1

# tolerance rate of 0.05%
testmod_rev<-postpr(target1, models_rev, stats5mod, tol=.05, method="neuralnet")

Warning message:
There are 5 models but only 1 for which simulations have been accepted.
No regression is performed, method is set to rejection.
Consider increasing the tolerance rate.TRUE 


# function summary prints out posterior model probabilities and ratios of model probabilities (the Bayes factors) in a user-friendly way
summary(testmod_rev)


> summary(testmod_rev)
Call: 
postpr(target = target1, index = models_rev, sumstat = stats5mod, 
    tol = 0.05, method = "neuralnet")
Data:
 postpr.out$values (2500 posterior samples)
Models a priori:
 nullplusghost, orig.ghoste, orig.ghostw, rev.ghoste, rev.ghostw
Models a posteriori:
 nullplusghost, orig.ghoste, orig.ghostw, rev.ghoste, rev.ghostw

Proportion of accepted simulations (rejection):
nullplusghost   orig.ghoste   orig.ghostw    rev.ghoste    rev.ghostw 
            0             1             0             0             0 

Bayes factors:
              nullplusghost orig.ghoste orig.ghostw rev.ghoste rev.ghostw
nullplusghost                         0                                  
orig.ghoste             Inf           1         Inf        Inf        Inf
orig.ghostw                           0                                  
rev.ghoste                            0                                  
rev.ghostw                            0                                  





#-----------------------------------------------------------------------------------------------------------------------


e_stats3mod<-rbind(as.data.frame(nullplusghost_stats1),as.data.frame(orig.ghoste_stats1),as.data.frame(rev.ghoste_stats1))

#models<-c(
#replicate(2500, "nullplusghost"),
#replicate(2500, "ghoste")
#    )

# increased to 200*50
e_models_rev<-c(
replicate(10000, "nullplusghost"),
replicate(10000, "orig.ghoste"),
replicate(10000, "rev.ghoste")
    )


e_testmod_rev<-postpr(target1, e_models_rev, e_stats3mod, tol=.05, method="neuralnet")

Warning message:
There are 3 models but only 2 for which simulations have been accepted.
No regression is performed, method is set to rejection.
Consider increasing the tolerance rate.TRUE 

# function summary prints out posterior model probabilities and ratios of model probabilities (the Bayes factors) in a user-friendly way
summary(e_testmod_rev)


> summary(e_testmod_rev)
Call: 
postpr(target = target1, index = e_models_rev, sumstat = e_stats3mod, 
    tol = 0.05, method = "neuralnet")
Data:
 postpr.out$values (1500 posterior samples)
Models a priori:
 nullplusghost, orig.ghoste, rev.ghoste
Models a posteriori:
 nullplusghost, orig.ghoste, rev.ghoste

Proportion of accepted simulations (rejection):
nullplusghost   orig.ghoste    rev.ghoste 
       0.2067        0.7933        0.0000 

Bayes factors:
              nullplusghost orig.ghoste rev.ghoste
nullplusghost        1.0000      0.2605        Inf
orig.ghoste          3.8387      1.0000        Inf
rev.ghoste           0.0000      0.0000           


#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

e_stats3mod<-rbind(as.data.frame(orig.ghoste_stats1),as.data.frame(rev.ghoste_stats1))

#models<-c(
#replicate(2500, "nullplusghost"),
#replicate(2500, "ghoste")
#    )

# increased to 200*50
e_models_rev<-c(
replicate(10000, "orig.ghoste"),
replicate(10000, "rev.ghoste")
    )


e_testmod_rev<-postpr(target1, e_models_rev, e_stats3mod, tol=.05, method="neuralnet")
#final  value 0.041031 
#converged
    # convergence only reached here *

# function summary prints out posterior model probabilities and ratios of model probabilities (the Bayes factors) in a user-friendly way
summary(e_testmod_rev)

Call: 
postpr(target = target1, index = e_models_rev, sumstat = e_stats3mod, 
    tol = 0.05, method = "neuralnet")
Data:
 postpr.out$values (1000 posterior samples)
Models a priori:
 orig.ghoste, rev.ghoste
Models a posteriori:
 orig.ghoste, rev.ghoste

Proportion of accepted simulations (rejection):
orig.ghoste  rev.ghoste 
      0.999       0.001 

Bayes factors:
            orig.ghoste rev.ghoste
orig.ghoste       1.000    999.000
rev.ghoste        0.001      1.000


Posterior model probabilities (neuralnet):
orig.ghoste  rev.ghoste 
     0.9991      0.0009 

Bayes factors:
            orig.ghoste rev.ghoste
orig.ghoste      1.0000  1151.9649
rev.ghoste       0.0009     1.0000


#-----------------------------------------------------------------------------------------------------------------------

# check if these can be distinguished...
cv.modsel <- cv4postpr(e_models_rev, e_stats3mod, nval=100, tol=0.05, method="neuralnet")

s <- summary(cv.modsel)






#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

# Mon 27 Feb 2023 14:55:14 CET
# increase tolerance rate from 0.05 to 0.1


stats5mod<-rbind(as.data.frame(nullplusghost_stats1),as.data.frame(orig.ghoste_stats1),as.data.frame(orig.ghostw_stats1),as.data.frame(rev.ghoste_stats1),as.data.frame(rev.ghostw_stats1))


# increased to 200*50
models_rev<-c(
replicate(10000, "nullplusghost"),
replicate(10000, "orig.ghoste"),
replicate(10000, "orig.ghostw"),
replicate(10000, "rev.ghoste"),
replicate(10000, "rev.ghostw")
    )


load("/scratch/devel/hpawar/admix/abc/simul/test/ghost/3nov21/segrecalc/ghosttarget1_22nov",verbose=T)
#Loading objects:
 # target1

# tolerance rate of 0.1
testmod_rev<-postpr(target1, models_rev, stats5mod, tol=.1, method="neuralnet")


# function summary prints out posterior model probabilities and ratios of model probabilities (the Bayes factors) in a user-friendly way
summary(testmod_rev)

> summary(testmod_rev)ayes factors) in a user-friendly way
Call: 
postpr(target = target1, index = models_rev, sumstat = stats5mod, 
    tol = 0.1, method = "neuralnet")
Data:
 postpr.out$values (5000 posterior samples)
Models a priori:
 nullplusghost, orig.ghoste, orig.ghostw, rev.ghoste, rev.ghostw
Models a posteriori:
 nullplusghost, orig.ghoste, orig.ghostw, rev.ghoste, rev.ghostw

Proportion of accepted simulations (rejection):
nullplusghost   orig.ghoste   orig.ghostw    rev.ghoste    rev.ghostw 
       0.0012        0.9988        0.0000        0.0000        0.0000 

Bayes factors:
              nullplusghost orig.ghoste orig.ghostw rev.ghoste rev.ghostw
nullplusghost        1.0000      0.0012         Inf        Inf        Inf
orig.ghoste        832.3333      1.0000         Inf        Inf        Inf
orig.ghostw          0.0000      0.0000                                  
rev.ghoste           0.0000      0.0000                                  
rev.ghostw           0.0000      0.0000          

#-----------------------------------------------------------------------------------------------------------------------

e_stats3mod<-rbind(as.data.frame(nullplusghost_stats1),as.data.frame(orig.ghoste_stats1),as.data.frame(rev.ghoste_stats1))


# increased to 200*50
e_models_rev<-c(
replicate(10000, "nullplusghost"),
replicate(10000, "orig.ghoste"),
replicate(10000, "rev.ghoste")
    )

e_testmod_rev<-postpr(target1, e_models_rev, e_stats3mod, tol=.1, method="neuralnet")


# function summary prints out posterior model probabilities and ratios of model probabilities (the Bayes factors) in a user-friendly way
summary(e_testmod_rev)



> e_testmod_rev<-postpr(target1, e_models_rev, e_stats3mod, tol=.1, method=Warning message:
There are 3 models but only 2 for which simulations have been accepted.
No regression is performed, method is set to rejection.
Consider increasing the tolerance rate.TRUE 
> 
> 
> # function summary prints out posterior model probabilities and ratios of> summary(e_testmod_rev)es factors) in a user-friendly way
Call: 
postpr(target = target1, index = e_models_rev, sumstat = e_stats3mod, 
    tol = 0.1, method = "neuralnet")
Data:
 postpr.out$values (3000 posterior samples)
Models a priori:
 nullplusghost, orig.ghoste, rev.ghoste
Models a posteriori:
 nullplusghost, orig.ghoste, rev.ghoste

Proportion of accepted simulations (rejection):
nullplusghost   orig.ghoste    rev.ghoste 
        0.261         0.739         0.000 

Bayes factors:
              nullplusghost orig.ghoste rev.ghoste
nullplusghost        1.0000      0.3532        Inf
orig.ghoste          2.8314      1.0000        Inf
rev.ghoste           0.0000      0.0000     




#-----------------------------------------------------------------------------------------------------------------------


e_stats3mod<-rbind(as.data.frame(orig.ghoste_stats1),as.data.frame(rev.ghoste_stats1))

#models<-c(
#replicate(2500, "nullplusghost"),
#replicate(2500, "ghoste")
#    )

# increased to 200*50
e_models_rev<-c(
replicate(10000, "orig.ghoste"),
replicate(10000, "rev.ghoste")
    )


e_testmod_rev<-postpr(target1, e_models_rev, e_stats3mod, tol=.1, method="neuralnet")
#final  value 0.041031 
#converged
    # convergence only reached here *

# function summary prints out posterior model probabilities and ratios of model probabilities (the Bayes factors) in a user-friendly way
summary(e_testmod_rev)


> summary(e_testmod_rev)es factors) in a user-friendly way
Call: 
postpr(target = target1, index = e_models_rev, sumstat = e_stats3mod, 
    tol = 0.1, method = "neuralnet")
Data:
 postpr.out$values (2000 posterior samples)
Models a priori:
 orig.ghoste, rev.ghoste
Models a posteriori:
 orig.ghoste, rev.ghoste

Proportion of accepted simulations (rejection):
orig.ghoste  rev.ghoste 
     0.9965      0.0035 

Bayes factors:
            orig.ghoste rev.ghoste
orig.ghoste      1.0000   284.7143
rev.ghoste       0.0035     1.0000


Posterior model probabilities (neuralnet):
orig.ghoste  rev.ghoste 
     0.9992      0.0008 

Bayes factors:
            orig.ghoste rev.ghoste
orig.ghoste      1.0000  1244.3354
rev.ghoste       0.0008     1.0000


