#!/usr/bin/env
#module load R/4.0.1
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

# Fri 22 Apr 2022 12:26:33 CEST
# final ABC parameter inference for model with possible ghost gene flow into w_anc
#-----------------------------------------------------------------------------------------------------------------------

# model C)
# /scratch/devel/hpawar/admix/abc/simul/scripts/wanc.ghost.model.1apr22.R  # called by -
#/scratch/devel/hpawar/admix/abc/simul/scripts/wanc.ghost.model.v1.arr

# targeting 700*51 = 35,700 reps

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
library(abc)
options(scipen=100)
require(data.table)
options(stringsAsFactors=F)


# load in the empirical data # directly read in target1
load("/scratch/devel/hpawar/admix/abc/results/test/11sep21simns/segrecalc/target1_6oct", verbose=T)
#Loading objects:
#  target1 
#summary stats for empirical data after removing uninformative parameters # het_sd_WC (6), fixedperid_sd_WC (22), pi_mu_WC (32), pi_sd_WC (36)

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


#Thu 30 Sep 2021 11:59:08 CEST
#STREAMLINED VERSION OF FUNCTIONS
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


# all simulation reps for ghost -> w_anc have generated

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


# apply as follows:

# fst
s_fst<-fst_process_fun(out_fst_fun1,format_problsimns_fun)


# pi_mu
s_pi_mu<-comb_process_fun(out_pi_fun1,pi_format_problsimns_fun,5,1)

# pi_sd
s_pi_sd<-comb_process_fun(out_pi_fun1,pi_format_problsimns_fun,5,2)

# tajima_mu - updated way of calculating, process below

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
s_popfix_perkb<-((s_popfix/(2500*40000))*1000)


#  popn-wise seg sites
s_popseg<-comb_process_fun(out_het_fun1,het_format_problsimns_fun,2,2)
s_popseg_perkb<-((s_popseg/(2500*40000))*1000)

#-----------------------------------------------------------------------------------------------------------------------
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


#  remove those iterations of simulations which failed & did not output summary stats - none did here
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

paraname=c("ghost_t0", "e_moun_t3.1", "t4", "t_archintrog", "archaicintrog", "admix_w_e_t6", "t9","gor_ghost_anc","id")

colnames(param)<-paraname
head(param)

# remove the id parameter - ie now parameter 9
param<-param[,c(1:8)]

#-----------------------------------------------------------------------------------------------------------------------

# remove same cols from the simulated vals
sumstat1<-sumstat[,-c(6,22,32,36)]

# ncol(sumstat1) #[1] 40

#-----------------------------------------------------------------------------------------------------------------------


# remove ghost_t0 from the parameter search
param1<-param[,-c(1)]
head(param1)


# should also force these posteriors to be within the priors, same logit transf
# AMEND WITH PRIORS FROM /scratch/devel/hpawar/admix/abc/simul/scripts/wanc.ghost.model.1apr22.R **
# only priors which differ b/n ghoste & ghostw is t_archintrog=runif(1000,min=5.9803,max=14.3616) # lower prior here
lowerp<-c(2.9191,0.1446,5.9803,0,6.7955,15,10)
upperp<-c(19.1638,0.2420,14.3616,100,99.2159,50,100)

priordf<-cbind(lowerp,upperp)
prior_ranges<-data.matrix(priordf)

# Which you feed into the ABC like this:
wghostabc_22apr22<-abc(target=target1,param=param1,tol=0.005,sumstat=sumstat1,method="neuralnet",numnet=100,transf="logit",logit.bounds=prior_ranges)
# The "tol" parameter is unknown and might be tested a bit.

#> wghostabc_22apr22<-abc(target=target1,param=param1,tol=0.005,sumstat=sumstat1,method="neuralnet",numnet=100,transf="logit",logit.bounds=prior_ranges)
#123456789101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899100
#123456789101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899100
#Warning message:
#All parameters are "logit" transformed. 

#-----------------------------------------------------------------------------------------------------------------------

wghostabc_22apr22[[17]][[1]]<-colnames(param1)
wghostabc_22apr22[[17]][[2]]<-names(target1)

wghostabc_22apr22
summary(wghostabc_22apr22)

save(wghostabc_22apr22,file=paste("/scratch/devel/hpawar/admix/abc/simul/test/ghost/final_1mar22/wanc/segrecalc/wghostabc_22apr22"))
save(param1,file=paste("/scratch/devel/hpawar/admix/abc/simul/test/ghost/final_1mar22/wanc/segrecalc/wghostabc_param_22apr22"))
save(sumstat1,file=paste("/scratch/devel/hpawar/admix/abc/simul/test/ghost/final_1mar22/wanc/segrecalc/wghostabc_sumstat1_22apr22"))

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
