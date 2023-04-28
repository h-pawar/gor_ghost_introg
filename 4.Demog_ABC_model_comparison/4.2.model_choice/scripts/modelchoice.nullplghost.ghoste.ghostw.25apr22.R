#!/usr/bin/env
#module load R/4.0.1

#Mon 25 Apr 2022 14:16:16 CEST
# have generated 10,000 reps for nullplghost & ghoste & ghostw - using final posteriors

# dir with mc simulated data -  
# /scratch/devel/hpawar/admix/abc/simul/test/modelcomp/final_12apr22/nullplusghost/
# /scratch/devel/hpawar/admix/abc/simul/test/modelcomp/final_12apr22/ghoste/
# /scratch/devel/hpawar/admix/abc/simul/test/modelcomp/final_12apr22/ghostw

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

# process summary statistics from 'model comparison simulations'
    # these are simulations generated taking the parameter vals as the weighted median posteriors inferred from ABC analysis
# 1) process summary stats in same way as in previous generateabc scripts
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
# NOTE - change division here ** 
# b/c mc simulations have generated 250 windows per iter (rather than 2500)

#  popn-wise fixed sites
s_popfix<-comb_process_fun(out_het_fun1,het_format_problsimns_fun,2,1)
s_popfix_perkb<-((s_popfix/(250*40000))*1000)


#  popn-wise seg sites
s_popseg<-comb_process_fun(out_het_fun1,het_format_problsimns_fun,2,2)
s_popseg_perkb<-((s_popseg/(250*40000))*1000)


#-----------------------------------------------------------------------------------------------------------------------
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


#ncol(sumstat) # [1] 44
# nrow(sumstat)
#[1] 35200

sumstat<-data.matrix(sumstat)
# matching simulation vectors, each of them in the same order as summary stats from the real data (the "sumstat")

return(sumstat)
}



#-----------------------------------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------------------------------

# 1) Read in model comparison (mc) simulations & process the summary statistics

# A.2) null + non-interacting ghost mc simulations 
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

nullplusghost_simns<-simns
nullplusghost_param<-simns_param

#-----------------------------------------------------------------------------------------------------------------------
# 1) Read in model comparison (mc) simulations & process the summary statistics

# B) ghoste 
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

ghoste_stats<-matrix_function()

ghoste_simns<-simns
ghoste_param<-simns_param

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
# 1) Read in model comparison (mc) simulations & process the summary statistics

# C) ghostw 
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

ghostw_stats<-matrix_function()

ghostw_simns<-simns
ghostw_param<-simns_param

#-----------------------------------------------------------------------------------------------------------------------



# 1.2) remove uninformative stats (same as removed for previous ABC simulations) - remove from each of the models
# het_sd_WC (6) #fixedperid_sd_WC (22) # pi_mu_WC (32)# pi_sd_WC (36)

nullplusghost_stats1<-nullplusghost_stats[,-c(6,22,32,36)]
ghoste_stats1<-ghoste_stats[,-c(6,22,32,36)]
ghostw_stats1<-ghostw_stats[,-c(6,22,32,36)]


# 2) perform model choice

# ie combine stats into dfs

stats3mod<-rbind(as.data.frame(nullplusghost_stats1),as.data.frame(ghoste_stats1),as.data.frame(ghostw_stats1))


models<-c(
replicate(10000, "nullplusghost"),
replicate(10000, "ghoste"),
replicate(10000, "ghostw")
    )
#-----------------------------------------------------------------------------------------------------------------------

#  2.1) cross-validation for model selection
# if ABC can distinguish between the models

# nval=100
cv.modsel <- cv4postpr(models, stats3mod, nval=100, tol=0.05, method="neuralnet")
s <- summary(cv.modsel)

#-----------------------------------------------------------------------------------------------------------------------

# cross validation with nval=1000 sent as a job
#cv.modsel <- cv4postpr(models, stats3mod, nval=1000, tol=0.05, method="neuralnet")
#s <- summary(cv.modsel)
#/scratch/devel/hpawar/admix/abc/simul/scripts/modelcomp/check.crossvalidation.25apr22.R
#/scratch/devel/hpawar/admix/abc/simul/scripts/modelcomp/check.crossvalidation.25apr22.arr
#less /scratch/devel/hpawar/admix/sstar/log/cross_valid_39152694_1.out # outputs to screen from running the crossvalidation R script



#-----------------------------------------------------------------------------------------------------------------------

# 2.2) model selection with postpr function
# calculate posterior probabilities of each demog model

# read in summary stats from empirical data
load("/scratch/devel/hpawar/admix/abc/simul/test/ghost/3nov21/segrecalc/ghosttarget1_22nov",verbose=T)
#Loading objects:
 # target1

# tolerance rate of 0.05%
testmod<-postpr(target1, models, stats3mod, tol=.05, method="neuralnet")

# function summary prints out posterior model probabilities and ratios of model probabilities (the Bayes factors) in a user-friendly way
summary(testmod)

#-----------------------------------------------------------------------------------------------------------------------
 paramnames<-c(
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


parmanames1<-paramnames[-c(6,22,32,36)]

colnames(stats3mod)<-parmanames1

pdf("/scratch/devel/hpawar/admix/abc/simul/test/modelcomp/final_12apr22/modelchoice/mc.summarystats.nullplg.ghoste.ghostw.25apr22.pdf") 



boxplot(stats3mod[,"het_mu_WL"]~models, main="het_mu_WL",ylim=c(0.85,2.3))
abline(h = target1[[1]], col = "red")

boxplot(stats3mod[,"het_mu_WC"]~models, main="het_mu_WC",ylim=c(0.6,1.6))
abline(h = target1[[2]], col = "red")

boxplot(stats3mod[,"het_mu_EL"]~models, main="het_mu_EL",ylim=c(0.1,0.45))
abline(h = target1[[3]], col = "red")

boxplot(stats3mod[,"het_mu_EM"]~models, main="het_mu_EM",ylim=c(0.2,0.65))
abline(h = target1[[4]], col = "red")

boxplot(stats3mod[,"het_sd_WL"]~models, main="het_sd_WL")
abline(h = target1[[5]], col = "red")

boxplot(stats3mod[,"het_sd_EL"]~models, main="het_sd_EL")
abline(h = target1[[6]], col = "red")

boxplot(stats3mod[,"het_sd_EM"]~models, main="het_sd_EM")
abline(h = target1[[7]], col = "red")



boxplot(stats3mod[,"fixedsites_WL"]~models, main="fixedsites_WL")
abline(h = target1[[8]], col = "red")

boxplot(stats3mod[,"fixedsites_WC"]~models, main="fixedsites_WC",ylim=c(0.5,3))
abline(h = target1[[9]], col = "red")

boxplot(stats3mod[,"fixedsites_EL"]~models, main="fixedsites_EL",ylim=c(0.5,2.1))
abline(h = target1[[10]], col = "red")

boxplot(stats3mod[,"fixedsites_EM"]~models, main="fixedsites_EM",ylim=c(0.4,1.9))
abline(h = target1[[11]], col = "red")


boxplot(stats3mod[,"segsites_WL"]~models, main="segsites_WL",ylim=c(3.5,11))
abline(h = target1[[12]], col = "red")

boxplot(stats3mod[,"segsites_WC"]~models, main="segsites_WC",ylim=c(0.6,1.6))
abline(h = target1[[13]], col = "red")

boxplot(stats3mod[,"segsites_EL"]~models, main="segsites_EL",ylim=c(0.3,1.4))
abline(h = target1[[14]], col = "red")

boxplot(stats3mod[,"segsites_EM"]~models, main="segsites_EM",ylim=c(1,6.5))
abline(h = target1[[15]], col = "red")



boxplot(stats3mod[,"fixedperid_mu_WL"]~models, main="fixedperid_mu_WL")
abline(h = target1[[16]], col = "red")

boxplot(stats3mod[,"fixedperid_mu_WC"]~models, main="fixedperid_mu_WC",ylim=c(0.5,3))
abline(h = target1[[17]], col = "red")

boxplot(stats3mod[,"fixedperid_mu_EL"]~models, main="fixedperid_mu_EL",ylim=c(0.8,2.2))
abline(h = target1[[18]], col = "red")

boxplot(stats3mod[,"fixedperid_mu_EM"]~models, main="fixedperid_mu_EM",ylim=c(0.8,2.3))
abline(h = target1[[19]], col = "red")



boxplot(stats3mod[,"fixedperid_sd_WL"]~models, main="fixedperid_sd_WL")
abline(h = target1[[20]], col = "red")

boxplot(stats3mod[,"fixedperid_sd_EL"]~models, main="fixedperid_sd_EL")
abline(h = target1[[21]], col = "red")

boxplot(stats3mod[,"fixedperid_sd_EM"]~models, main="fixedperid_sd_EM")
abline(h = target1[[22]], col = "red")



boxplot(stats3mod[,"fst_WL.WC"]~models, main="fst_WL.WC")
abline(h = target1[[23]], col = "red")

boxplot(stats3mod[,"fst_WL.EL"]~models, main="fst_WL.EL")
abline(h = target1[[24]], col = "red")

boxplot(stats3mod[,"fst_WL.EM"]~models, main="fst_WL.EM",ylim=c(0.17,0.24))
abline(h = target1[[25]], col = "red")

boxplot(stats3mod[,"fst_WC.EL"]~models, main="fst_WC.EL",ylim=c(0.3,0.65))
abline(h = target1[[26]], col = "red")

boxplot(stats3mod[,"fst_WC.EM"]~models, main="fst_WC.EM",ylim=c(0.3,0.58))
abline(h = target1[[27]], col = "red")

boxplot(stats3mod[,"fst_EL.EM"]~models, main="fst_EL.EM",ylim=c(0.1,0.3))
abline(h = target1[[28]], col = "red")


boxplot(stats3mod[,"pi_mu_WL"]~models, main="pi_mu_WL",ylim=c(0.05,0.16))
abline(h = target1[[29]], col = "red")

boxplot(stats3mod[,"pi_mu_EL"]~models, main="pi_mu_EL",ylim=c(0.01,0.031))
abline(h = target1[[30]], col = "red")

boxplot(stats3mod[,"pi_mu_EM"]~models, main="pi_mu_EM",ylim=c(0.02,0.045))
abline(h = target1[[31]], col = "red")

boxplot(stats3mod[,"pi_sd_WL"]~models, main="pi_sd_WL")
abline(h = target1[[32]], col = "red")

boxplot(stats3mod[,"pi_sd_EL"]~models, main="pi_sd_EL")
abline(h = target1[[33]], col = "red")

boxplot(stats3mod[,"pi_sd_EM"]~models, main="pi_sd_EM")
abline(h = target1[[34]], col = "red")


boxplot(stats3mod[,"tajima_mu_WL"]~models, main="tajima_mu_WL",ylim=c(-0.9,0.1))
abline(h = target1[[35]], col = "red")

boxplot(stats3mod[,"tajima_mu_EL"]~models, main="tajima_mu_EL",ylim=c(0.08,1.12))
abline(h = target1[[36]], col = "red")

boxplot(stats3mod[,"tajima_mu_EM"]~models, main="tajima_mu_EM",ylim=c(0.05,0.65))
abline(h = target1[[37]], col = "red")

boxplot(stats3mod[,"tajima_sd_WL"]~models, main="tajima_sd_WL",ylim=c(0.25,0.81))
abline(h = target1[[38]], col = "red")

boxplot(stats3mod[,"tajima_sd_EL"]~models, main="tajima_sd_EL",ylim=c(1,1.61))
abline(h = target1[[39]], col = "red")

boxplot(stats3mod[,"tajima_sd_EM"]~models, main="tajima_sd_EM",ylim=c(0.7,1.3))
abline(h = target1[[40]], col = "red")
dev.off()
