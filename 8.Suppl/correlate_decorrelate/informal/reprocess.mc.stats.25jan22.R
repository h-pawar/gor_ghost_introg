#!/usr/bin/env
#module load R/4.0.1

#-----------------------------------------------------------------------------------------------------------------------
# Tue 25 Jan 2022 10:08:05 GMT
# for the model comparison simulations - amend processing of fixed & seg sites

# MK: try to simplify the correlated statistics.
# The easiest way would be to take the correlated stats and add them up, 
#   I suggest to use the rowSums of all fixedperid, segsites, and fixedsites together (so, the sum of 12 values in each row), 
#   and leave all other stats as they are, and then see how that changes the picture.

# MK
#model comparison here - since you have the stats, you dont need to make new simulations, 
#just calculate the rowSums of those that are there as one statistic next to the FSTs, Pis, etc. 
#As a very easy first step to see how this influences the comparison. 
#Then you can also look into other approaches, not sure how much effort eh PLS transformation is.

#-----------------------------------------------------------------------------------------------------------------------
#/Volumes/"Ultra USB 3.0"/IBE/further.analysis.feb.2020/gorillas/abc/simul_postabc/abc+ghost/modelcomp/modelchoice.nullplghost.22jan22.R
# ie change how processing these stats

#-----------------------------------------------------------------------------------------------------------------------
library(abc)
options(scipen=100)
require(data.table)
options(stringsAsFactors=F)

# A.2) null + ghost mc simulations - after removing sites fixed in all gor ids (weighted median posterior parameter vals)
simufiles <- paste("/scratch/devel/hpawar/admix/abc/simul/test/modelcomp/14dec21/21jan22/nullplusghost/", list.files(path = "/scratch/devel/hpawar/admix/abc/simul/test/modelcomp/14dec21/21jan22/nullplusghost/", pattern="test.mc.nullplusghost_sim"), sep = "")

simns=list()
simns_param=list()
for (i in 1:length(simufiles)){
load(file=simufiles[[i]],verbose=T)
simns[[i]]<-simuresults
simns_param[[i]]<-inputsets
}


load(file="/scratch/devel/hpawar/admix/abc/simul/test/modelcomp/14dec21/21jan22/nullplusghost/test.mc.nullplusghost_sim2",verbose=T)
simns[[2]]<-simuresults
simns_param[[2]]<-inputsets

load(file="/scratch/devel/hpawar/admix/abc/simul/test/modelcomp/14dec21/21jan22/nullplusghost/test.mc.nullplusghost_sim50",verbose=T)
simns[[50]]<-simuresults
simns_param[[50]]<-inputsets

#-----------------------------------------------------------------------------------------------------------------------

# str(simns)
#List of 50
# $ :List of 50
#  ..$ V1 :List of 6
#  .. ..$ : num [1:2, 1:4] 1.4179 0.0391 0.8223 NA 0.8014 ...
#  .. ..$ : num [1:2, 1:4] 768 74783 9014 8223 5049 ...
#  .. ..$ : num [1:2, 1:4] 0.6125 0.0199 0.9014 NA 0.9291 ...
#  .. ..$ : num [1:6, 1] 0.154 0.179 0.189 0.344 0.385 ...
#  .. ..$ : num [1:4, 1:2] 0.1288 0.0741 0.0736 0.0668 0.0239 ...
#  .. ..$ : num [1:3, 1:2] 0.251 0.36 0.487 0.601 0.592 ...
# order of stats:
#return(list(het_op,full_segsites,fix_op,fst.out,fout.pi,fout.taj))


#full_segsites,fix_op - need to take row sums of

# ifix mu 3,1
# ifix sd 3,2

# popfix #2,1
# popseg #2,2

#-----------------------------------------------------------------------------------------------------------------------
# individual functions (same as before)


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
#out_1_200<-afun(1,200)
# reducing to 1,10 - Thu  9 Dec 2021 15:18:13 CET **
# now 1,50 - Thu 16 Dec 2021 14:45:06 CET **
out_1_200<-afun(1,50)
rout_1_200<-as.data.frame(do.call(rbind,out_1_200))
return(rout_1_200)
}


comb_process_fun<-function(afun,afun1,a,b) {
#out_1_200<-afun(1,200,a,b)
out_1_200<-afun(1,50,a,b)
rout_1_200<-as.data.frame(do.call(rbind,out_1_200))
return(rout_1_200)
}

#-----------------------------------------------------------------------------------------------------------------------


# have amended matrix function - to take the correlated stats and add them up, 
# take rowSums of all fixedperid, segsites, and fixedsites together (so, the sum of 12 values in each row),


matrix_function2<-function(){
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

# target[c(6,22,32,36)]
#       het_sd_WC fixedperid_sd_WC         pi_mu_WC         pi_sd_WC 
#               0                0                0                0 
              # parameters to remove

s_ifix_sd<-s_ifix_sd[,-c(2)]

# combine :
# fixed sites per id mu, fixed sites per id sd, popn-wise fixed sites in kb,  popn-wise seg sites in kb

comb_fixseg<-cbind(s_ifix_mu,s_ifix_sd,s_popfix_perkb,s_popseg_perkb)

comb_fixseg_sum<-rowSums(comb_fixseg)
# is this what is meant?

# but now this column is an order of magnitude diff from rest - log transform?
# comb_fixseg_sum        X1        X2        X3        X4        X5         X6
#1        19.21838 0.1537853 0.1788508 0.1889740 0.3436768 0.3850433 0.05454134
#2        19.35858 0.1581377 0.1761734 0.1859381 0.3508051 0.3871871 0.05299687
#3        19.38465 0.1494828 0.1750970 0.1866308 0.3336969 0.3737495 0.05085960

# head(log(comb_fixseg_sum))
#[1] 2.955867 2.963136 2.964482 2.961887 2.959906 2.957641
# even after log transforming the rowsums of correlated are quite higher than hte rest?

# MK
#I would not log transform them, better divide by 10 to move the order of magnitude (
  #that would equal these stats per 10kbp rather than 1kbp). But otherwise, yes, that is what I meant.

comb_fixseg_sum1<-comb_fixseg_sum/10

# head(comb_fixseg_sum1)
#[1] 1.921838 1.935858 1.938465 1.933442 1.929616 1.925250

sumstat2<-cbind(s_het_mu,s_het_sd,
comb_fixseg_sum1,
s_fst,
s_pi_mu,
s_pi_sd,
s_tajima_mu,
s_tajima_sd)

#> ncol(sumstat2)
#[1] 29
#> nrow(sumstat2)
#[1] 2500

sumstat2<-data.matrix(sumstat2) 

return(sumstat2)
}

#-----------------------------------------------------------------------------------------------------------------------
# now need to calculate for the empirical data in the same way **
# or take the sums of the empirical post-hoc?
# otherwise - referring back to generateabc.6oct.R - processing of emp summary stats & first abc analysis

#-----------------------------------------------------------------------------------------------------------------------

# apply new matrix function for null
nullplusghostlongdiv_stats<-matrix_function2()

# shoudl assign simns to diff vector - to retain original data for each model
nullplusghostlongdiv_simns<-simns
nullplusghostlongdiv_param<-simns_param

#-----------------------------------------------------------------------------------------------------------------------

# B) ghoste - setting ne of ghost pop to 25 instead of 0.1
simufiles <- paste("/scratch/devel/hpawar/admix/abc/simul/test/modelcomp/14dec21/21jan22/ghoste/", list.files(path = "/scratch/devel/hpawar/admix/abc/simul/test/modelcomp/14dec21/21jan22/ghoste/", pattern="test.mc.ghoste_sim"), sep = "")


simns=list()
simns_param=list()
for (i in 1:length(simufiles)){
load(file=simufiles[[i]],verbose=T)
simns[[i]]<-simuresults
simns_param[[i]]<-inputsets
}

ghoste_stats<-matrix_function2()

# shoudl assign simns to diff vector - to retain original data for each model
ghoste_simns<-simns
ghoste_param<-simns_param

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------


# 1.2) remove uninformative stats (same as removed for previous ABC simulations) - remove from each of the models
# het_sd_WC (6) #fixedperid_sd_WC (22) (already removed) # pi_mu_WC (32)# pi_sd_WC (36)

#sumstat2<-cbind(s_het_mu,s_het_sd,
#comb_fixseg_sum1,
#s_fst,
#s_pi_mu,
#s_pi_sd,
#s_tajima_mu,
#s_tajima_sd)


 paramnames<-c(
"het_mu_WL", "het_mu_WC", "het_mu_EL", "het_mu_EM",
"het_sd_WL", "het_sd_WC","het_sd_EL", "het_sd_EM",
"combfixseg",
"fst_WL.WC","fst_WL.EL","fst_WL.EM","fst_WC.EL","fst_WC.EM","fst_EL.EM",
"pi_mu_WL", "pi_mu_WC", "pi_mu_EL", "pi_mu_EM",
"pi_sd_WL", "pi_sd_WC", "pi_sd_EL", "pi_sd_EM",
 "tajima_mu_WL", "tajima_mu_EL", "tajima_mu_EM",
  "tajima_sd_WL", "tajima_sd_EL", "tajima_sd_EM")

# remove 6, 17,21: "het_sd_WC",  "pi_mu_WC", "pi_sd_WC"
parmanames1<-paramnames[-c(6,17,21)]


nullplusghostlongdiv_stats1<-nullplusghostlongdiv_stats[,-c(6,17,21)]
ghoste_stats1<-ghoste_stats[,-c(6,17,21)]

#ncol(nullplusghostlongdiv_stats1)
#[1] 26
#-----------------------------------------------------------------------------------------------------------------------

# 2) perform model choice

# ie combine stats into dfs

#stats3mod<-rbind(as.data.frame(nullplusghost_stats1),as.data.frame(ghoste_stats1))
stats3mod<-rbind(as.data.frame(nullplusghostlongdiv_stats1),as.data.frame(ghoste_stats1))

models<-c(
replicate(2500, "nullplusghost"),
replicate(2500, "ghoste")
    )
#-----------------------------------------------------------------------------------------------------------------------

#  2.1) cross-validation for model selection
# if ABC can distinguish between the models


cv.modsel <- cv4postpr(models, stats3mod, nval=100, tol=0.05, method="neuralnet")
# There were 50 or more warnings (use warnings() to see the first 50)


s <- summary(cv.modsel)

 s <- summary(cv.modsel)
Confusion matrix based on 100 samples for each model.

$tol0.05
              ghoste nullplusghost
ghoste           100             0
nullplusghost      0           100


Mean model posterior probabilities (neuralnet)

$tol0.05
              ghoste nullplusghost
ghoste             1             0
nullplusghost      0             1
#-----------------------------------------------------------------------------------------------------------------------

# 2.2) model selection with postpr function
# calculate posterior probabilities of each demog model

# read in summary stats from empirical data
load("/scratch/devel/hpawar/admix/abc/simul/test/ghost/3nov21/segrecalc/ghosttarget1_22nov",verbose=T)
#Loading objects:
 # target1

# target1[8:22]
#   fixedsites_WL    fixedsites_WC    fixedsites_EL    fixedsites_EM 
#      0.07708829       0.72709409       0.53358597       0.48170967 
#     segsites_WL      segsites_WC      segsites_EL      segsites_EM 
#      3.94090851       0.69363832       1.25447215       1.46983169 
#fixedperid_mu_WL fixedperid_mu_WC fixedperid_mu_EL fixedperid_mu_EM 
#      0.62758139       0.72709409       0.88054400       0.87568548 
#fixedperid_sd_WL fixedperid_sd_EL fixedperid_sd_EM 
#      0.03072747       0.01569050       0.01364289 

# take sum then divide by 10 - these stats per 10kbp rather than 1kbp)
tfixsegs<-target1[8:22]

tfixsegs1<-sum(tfixsegs)/10
names(tfixsegs1)<-c("comb_fixsegs")

target2<-c(target1[1:7],tfixsegs1,target1[23:40])

# target2
#   het_mu_WL    het_mu_WC    het_mu_EL    het_mu_EM    het_sd_WL    het_sd_EL 
#  0.91084593   0.70084007   0.41157419   0.42950219   0.06199505   0.03066836 
#   het_sd_EM comb_fixsegs    fst_WL.WC    fst_WL.EL    fst_WL.EM    fst_WC.EL 
#  0.03339686   1.23492945   0.13441066   0.18525497   0.17661490   0.36018834 
#   fst_WC.EM    fst_EL.EM     pi_mu_WL     pi_mu_EL     pi_mu_EM     pi_sd_WL 
#  0.33758827   0.20051395   0.06813919   0.03044540   0.03674017   0.02216627 
#    pi_sd_EL     pi_sd_EM tajima_mu_WL tajima_mu_EL tajima_mu_EM tajima_sd_WL 
#  0.02362419   0.02474919   0.09748612   0.08333107   0.33797093   0.46499798 
#tajima_sd_EL tajima_sd_EM 
#  1.00589942   0.95881797 

# tolerance rate of 0.05%
testmod<-postpr(target2, models, stats3mod, tol=.05, method="neuralnet")
Warning message:
There are 2 models but only 1 for which simulations have been accepted.
No regression is performed, method is set to rejection.
Consider increasing the tolerance rate.TRUE 

# function summary prints out posterior model probabilities and ratios of model probabilities (the Bayes factors) in a user-friendly way
summary(testmod)
Call: 
postpr(target = target2, index = models, sumstat = stats3mod, 
    tol = 0.05, method = "neuralnet")
Data:
 postpr.out$values (250 posterior samples)
Models a priori:
 ghoste, nullplusghost
Models a posteriori:
 ghoste, nullplusghost

Proportion of accepted simulations (rejection):
       ghoste nullplusghost 
            1             0 

Bayes factors:
              ghoste nullplusghost
ghoste             1           Inf
nullplusghost      0              
#-----------------------------------------------------------------------------------------------------------------------


colnames(stats3mod)<-names(target2)


pdf("/scratch/devel/hpawar/admix/abc/simul/test/modelcomp/out/mc.summarystats.nullplghostlongdiv.ghoste.mergefixseg.25jan22.pdf") 



boxplot(stats3mod[,"het_mu_WL"]~models, main="het_mu_WL",ylim=c(0.85,2.0))
abline(h = target2[[1]], col = "red")

boxplot(stats3mod[,"het_mu_WC"]~models, main="het_mu_WC",ylim=c(0.5,1.5))
abline(h = target2[[2]], col = "red")

boxplot(stats3mod[,"het_mu_EL"]~models, main="het_mu_EL",ylim=c(0,2.5))
abline(h = target2[[3]], col = "red")

boxplot(stats3mod[,"het_mu_EM"]~models, main="het_mu_EM",ylim=c(0,2.5))
abline(h = target2[[4]], col = "red")

boxplot(stats3mod[,"het_sd_WL"]~models, main="het_sd_WL")
abline(h = target2[[5]], col = "red")

boxplot(stats3mod[,"het_sd_EL"]~models, main="het_sd_EL")
abline(h = target2[[6]], col = "red")

boxplot(stats3mod[,"het_sd_EM"]~models, main="het_sd_EM")
abline(h = target2[[7]], col = "red")


boxplot(stats3mod[,"comb_fixsegs"]~models, main="comb_fixsegs",ylim=c(1.2,4))
abline(h = target2[[8]], col = "red")



boxplot(stats3mod[,"fst_WL.WC"]~models, main="fst_WL.WC")
abline(h = target2[[9]], col = "red")

boxplot(stats3mod[,"fst_WL.EL"]~models, main="fst_WL.EL")
abline(h = target2[[10]], col = "red")

boxplot(stats3mod[,"fst_WL.EM"]~models, main="fst_WL.EM",ylim=c(0.17,0.24))
abline(h = target2[[11]], col = "red")

boxplot(stats3mod[,"fst_WC.EL"]~models, main="fst_WC.EL")
abline(h = target2[[12]], col = "red")

boxplot(stats3mod[,"fst_WC.EM"]~models, main="fst_WC.EM")
abline(h = target2[[13]], col = "red")

boxplot(stats3mod[,"fst_EL.EM"]~models, main="fst_EL.EM",ylim=c(0.04,0.21))
abline(h = target2[[14]], col = "red")


boxplot(stats3mod[,"pi_mu_WL"]~models, main="pi_mu_WL",ylim=c(0.05,0.15))
abline(h = target2[[15]], col = "red")

boxplot(stats3mod[,"pi_mu_EL"]~models, main="pi_mu_EL",ylim=c(0.03,0.14))
abline(h = target2[[16]], col = "red")

boxplot(stats3mod[,"pi_mu_EM"]~models, main="pi_mu_EM",ylim=c(0.03,0.13))
abline(h = target2[[17]], col = "red")

boxplot(stats3mod[,"pi_sd_WL"]~models, main="pi_sd_WL")
abline(h = target2[[18]], col = "red")

boxplot(stats3mod[,"pi_sd_EL"]~models, main="pi_sd_EL")
abline(h = target2[[19]], col = "red")

boxplot(stats3mod[,"pi_sd_EM"]~models, main="pi_sd_EM")
abline(h = target2[[20]], col = "red")


boxplot(stats3mod[,"tajima_mu_WL"]~models, main="tajima_mu_WL",ylim=c(0.09,0.4))
abline(h = target2[[21]], col = "red")

boxplot(stats3mod[,"tajima_mu_EL"]~models, main="tajima_mu_EL",ylim=c(0.08,0.45))
abline(h = target2[[22]], col = "red")

boxplot(stats3mod[,"tajima_mu_EM"]~models, main="tajima_mu_EM",ylim=c(0.3,0.55))
abline(h = target2[[23]], col = "red")

boxplot(stats3mod[,"tajima_sd_WL"]~models, main="tajima_sd_WL",ylim=c(0.4,0.7))
abline(h = target2[[24]], col = "red")

boxplot(stats3mod[,"tajima_sd_EL"]~models, main="tajima_sd_EL",ylim=c(0.5,1.1))
abline(h = target2[[25]], col = "red")

boxplot(stats3mod[,"tajima_sd_EM"]~models, main="tajima_sd_EM",ylim=c(0.4,1))
abline(h = target2[[26]], col = "red")


dev.off()

#scp -r hpawar@172.16.10.20:/scratch/devel/hpawar/admix/abc/simul/test/modelcomp/out/mc.summarystats.nullplghostlongdiv.ghoste.mergefixseg.25jan22.pdf  /Users/harvi/Downloads/gorilla_abc/modelchoice
#g6n.tm4D2L73


#pdf("/scratch/devel/hpawar/admix/abc/simul/test/modelcomp/out/mc.summarystats.nullplghostlongdiv.ghoste.mergefixseg.25jan22.pdf") 
#