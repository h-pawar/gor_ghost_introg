#!/usr/bin/env

#-----------------------------------------------------------------------------------------------------------------------
# Mon 27 Feb 2023 10:35:17 CET

# paths to model comparison simulations

# revised demog models (models D, E): 10,000 reps, model comparison simulations in paths - 
# /scratch/devel/hpawar/admix/abc/simul/test/ghost/revisions_8feb23/modelcomp/ghoste/rev.mc.ghoste_sim*
# /scratch/devel/hpawar/admix/abc/simul/test/ghost/revisions_8feb23/modelcomp/ghostw/rev.mc.ghostw_sim* 

# demog models (models A, B, C): 10,000 reps, model comparison simulations in paths - 
# /scratch/devel/hpawar/admix/abc/simul/test/modelcomp/final_12apr22/nullplusghost/
# /scratch/devel/hpawar/admix/abc/simul/test/modelcomp/final_12apr22/ghoste/
# /scratch/devel/hpawar/admix/abc/simul/test/modelcomp/final_12apr22/ghostw

#-----------------------------------------------------------------------------------------------------------------------
# 1) process summary statistics from 'model comparison simulations'
    # these are simulations generated taking the parameter vals as the weighted median posteriors inferred from ABC parameter inference
# 2) perform model choice
    # - cross validation
    # - misclassification
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
library(abc)
options(scipen=100)
require(data.table)
options(stringsAsFactors=F)

#-----------------------------------------------------------------------------------------------------------------------

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
out_1_200<-afun(1,200)
rout_1_200<-as.data.frame(do.call(rbind,out_1_200))
return(rout_1_200)
}


comb_process_fun<-function(afun,afun1,a,b) {
out_1_200<-afun(1,200,a,b)
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

# tajima - updated way of calculating, process below

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

#  popn-wise fixed sites
s_popfix<-comb_process_fun(out_het_fun1,het_format_problsimns_fun,2,1)
s_popfix_perkb<-((s_popfix/(250*40000))*1000)


#  popn-wise seg sites
s_popseg<-comb_process_fun(out_het_fun1,het_format_problsimns_fun,2,2)
s_popseg_perkb<-((s_popseg/(250*40000))*1000)

#-----------------------------------------------------------------------------------------------------------------------
# process tajimas d stat

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
out_10[[i]]<-process_tajima_fun(i)
}

tajima_df<-do.call(rbind, out_10)
#-----------------------------------------------------------------------------------------------------------------------


sumstat<-cbind(s_het_mu,s_het_sd,
s_popfix_perkb,s_popseg_perkb,
s_ifix_mu,s_ifix_sd,
s_fst,
s_pi_mu,
s_pi_sd,
tajima_df)

sumstat<-data.matrix(sumstat)
# matching simulation vectors, each of them in the same order as summary stats from the real data (the "sumstat")

return(sumstat)
}



#-----------------------------------------------------------------------------------------------------------------------

# D) revised ghoste - sample all parameters from priors
simufiles <- paste("/scratch/devel/hpawar/admix/abc/simul/test/ghost/revisions_8feb23/modelcomp/ghoste/", list.files(path = "/scratch/devel/hpawar/admix/abc/simul/test/ghost/revisions_8feb23/modelcomp/ghoste/", pattern="rev.mc.ghoste_sim"), sep = "")


simns=list()
simns_param=list()
for (i in 1:length(simufiles)){
load(file=simufiles[[i]],verbose=T)
simns[[i]]<-simuresults
simns_param[[i]]<-inputsets
}

rev.ghoste_stats<-matrix_function()

rev.ghoste_simns<-simns
rev.ghoste_param<-simns_param
#-----------------------------------------------------------------------------------------------------------------------
# E) revised ghostw - sample all parameters from priors
simufiles <- paste("/scratch/devel/hpawar/admix/abc/simul/test/ghost/revisions_8feb23/modelcomp/ghostw/", list.files(path = "/scratch/devel/hpawar/admix/abc/simul/test/ghost/revisions_8feb23/modelcomp/ghostw/", pattern="rev.mc.ghostw_sim"), sep = "")


simns=list()
simns_param=list()
for (i in 1:length(simufiles)){
load(file=simufiles[[i]],verbose=T)
simns[[i]]<-simuresults
simns_param[[i]]<-inputsets
}

rev.ghostw_stats<-matrix_function()

rev.ghostw_simns<-simns
rev.ghostw_param<-simns_param

#-----------------------------------------------------------------------------------------------------------------------

# A.1) null + non-interacting ghost mc simulations 
simufiles <- paste("/scratch/devel/hpawar/admix/abc/simul/test/modelcomp/final_12apr22/nullplusghost/", list.files(path = "/scratch/devel/hpawar/admix/abc/simul/test/modelcomp/final_12apr22/nullplusghost/", pattern="test.mc.nullplusghost_sim"), sep = "")

simns=list()
simns_param=list()
for (i in 1:length(simufiles)){
load(file=simufiles[[i]],verbose=T)
simns[[i]]<-simuresults
simns_param[[i]]<-inputsets
}

nullplusghost_stats<-matrix_function()

nullplusghost_simns<-simns
nullplusghost_param<-simns_param

#-----------------------------------------------------------------------------------------------------------------------

# B) ghoste 
simufiles <- paste("/scratch/devel/hpawar/admix/abc/simul/test/modelcomp/final_12apr22/ghoste/", list.files(path = "/scratch/devel/hpawar/admix/abc/simul/test/modelcomp/final_12apr22/ghoste/", pattern="test.mc.ghoste_sim"), sep = "")


simns=list()
simns_param=list()
for (i in 1:length(simufiles)){
load(file=simufiles[[i]],verbose=T)
simns[[i]]<-simuresults
simns_param[[i]]<-inputsets
}

orig.ghoste_stats<-matrix_function()

orig.ghoste_simns<-simns
orig.ghoste_param<-simns_param

#-----------------------------------------------------------------------------------------------------------------------

# C) ghostw
simufiles <- paste("/scratch/devel/hpawar/admix/abc/simul/test/modelcomp/final_12apr22/ghostw/", list.files(path = "/scratch/devel/hpawar/admix/abc/simul/test/modelcomp/final_12apr22/ghostw/", pattern="test.mc.ghostw_sim"), sep = "")


simns=list()
simns_param=list()
for (i in 1:length(simufiles)){
load(file=simufiles[[i]],verbose=T)
simns[[i]]<-simuresults
simns_param[[i]]<-inputsets
}

orig.ghostw_stats<-matrix_function()

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

stats3mod<-rbind(as.data.frame(nullplusghost_stats1),as.data.frame(rev.ghoste_stats1),as.data.frame(rev.ghostw_stats1))

models<-c(
replicate(10000, "nullplusghost"),
replicate(10000, "rev.ghoste"),
replicate(10000, "rev.ghostw")
    )
#-----------------------------------------------------------------------------------------------------------------------

#  2.1) cross-validation for model selection # if ABC can distinguish between the models


cv.modsel <- cv4postpr(models, stats3mod, nval=100, tol=0.05, method="neuralnet")
s <- summary(cv.modsel)
s
# cross validation with nval=1000 sent as job (in separate script)
#-----------------------------------------------------------------------------------------------------------------------

# 2.2) model selection with postpr function
# calculate posterior probabilities of each demog model

# read in summary stats from empirical data
load("/scratch/devel/hpawar/admix/abc/simul/test/ghost/3nov21/segrecalc/ghosttarget1_22nov",verbose=T)

# tolerance rate of 0.05%
testmod<-postpr(target1, models, stats3mod, tol=.05, method="neuralnet")

summary(testmod)

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------


# 2) compare :  null vs revised ghostE vs revised ghostW vs orig ghostE vs orig ghostW


stats5mod<-rbind(as.data.frame(nullplusghost_stats1),as.data.frame(orig.ghoste_stats1),as.data.frame(orig.ghostw_stats1),as.data.frame(rev.ghoste_stats1),as.data.frame(rev.ghostw_stats1))

models_rev<-c(
replicate(10000, "nullplusghost"),
replicate(10000, "orig.ghoste"),
replicate(10000, "orig.ghostw"),
replicate(10000, "rev.ghoste"),
replicate(10000, "rev.ghostw")
    )


#  2.1) cross-validation for model selection

cv.modsel <- cv4postpr(models_rev, stats5mod, nval=100, tol=0.05, method="neuralnet")


s <- summary(cv.modsel)


# 2.2) model selection with postpr function
# calculate posterior probabilities of each demog model

# read in summary stats from empirical data
load("/scratch/devel/hpawar/admix/abc/simul/test/ghost/3nov21/segrecalc/ghosttarget1_22nov",verbose=T)

# tolerance rate of 0.05%
testmod_rev<-postpr(target1, models_rev, stats5mod, tol=.05, method="neuralnet")

summary(testmod_rev)

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

# Mon 27 Feb 2023 14:55:14 CET
# increase tolerance rate from 0.05 to 0.1


stats5mod<-rbind(as.data.frame(nullplusghost_stats1),as.data.frame(orig.ghoste_stats1),as.data.frame(orig.ghostw_stats1),as.data.frame(rev.ghoste_stats1),as.data.frame(rev.ghostw_stats1))

models_rev<-c(
replicate(10000, "nullplusghost"),
replicate(10000, "orig.ghoste"),
replicate(10000, "orig.ghostw"),
replicate(10000, "rev.ghoste"),
replicate(10000, "rev.ghostw")
    )


load("/scratch/devel/hpawar/admix/abc/simul/test/ghost/3nov21/segrecalc/ghosttarget1_22nov",verbose=T)

# tolerance rate of 0.1
testmod_rev<-postpr(target1, models_rev, stats5mod, tol=.1, method="neuralnet")

summary(testmod_rev)

#-----------------------------------------------------------------------------------------------------------------------

e_stats3mod<-rbind(as.data.frame(nullplusghost_stats1),as.data.frame(orig.ghoste_stats1),as.data.frame(rev.ghoste_stats1))


e_models_rev<-c(
replicate(10000, "nullplusghost"),
replicate(10000, "orig.ghoste"),
replicate(10000, "rev.ghoste")
    )

e_testmod_rev<-postpr(target1, e_models_rev, e_stats3mod, tol=.1, method="neuralnet")

summary(e_testmod_rev)



#-----------------------------------------------------------------------------------------------------------------------

