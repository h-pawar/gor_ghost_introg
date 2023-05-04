#!/usr/bin/env
#module load R/4.0.1

#-----------------------------------------------------------------------------------------------------------------------
# Tue 25 Jan 2022 10:08:05 GMT

# simplify the correlated statistics.
# The easiest way would be to take the correlated stats and add them - rowSums of all fixedperid, segsites, and fixedsites stats
#-----------------------------------------------------------------------------------------------------------------------
library(abc)
options(scipen=100)
require(data.table)
options(stringsAsFactors=F)

# A.1) null + non-interacting ghost
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
out_1_200<-afun(1,50)
rout_1_200<-as.data.frame(do.call(rbind,out_1_200))
return(rout_1_200)
}


comb_process_fun<-function(afun,afun1,a,b) {
out_1_200<-afun(1,50,a,b)
rout_1_200<-as.data.frame(do.call(rbind,out_1_200))
return(rout_1_200)
}

#-----------------------------------------------------------------------------------------------------------------------


matrix_function2<-function(){
#-----------------------------------------------------------------------------------------------------------------------

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

#  popn-wise fixed sites
s_popfix<-comb_process_fun(out_het_fun1,het_format_problsimns_fun,2,1)
s_popfix_perkb<-((s_popfix/(250*40000))*1000)


#  popn-wise seg sites
s_popseg<-comb_process_fun(out_het_fun1,het_format_problsimns_fun,2,2)
s_popseg_perkb<-((s_popseg/(250*40000))*1000)


#-----------------------------------------------------------------------------------------------------------------------


s_ifix_sd<-s_ifix_sd[,-c(2)]


comb_fixseg<-cbind(s_ifix_mu,s_ifix_sd,s_popfix_perkb,s_popseg_perkb)

comb_fixseg_sum<-rowSums(comb_fixseg)


  # these stats per 10kbp rather than 1kbp
comb_fixseg_sum1<-comb_fixseg_sum/10


sumstat2<-cbind(s_het_mu,s_het_sd,
comb_fixseg_sum1,
s_fst,
s_pi_mu,
s_pi_sd,
s_tajima_mu,
s_tajima_sd)

sumstat2<-data.matrix(sumstat2) 

return(sumstat2)
}

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

nullplusghostlongdiv_stats<-matrix_function2()

nullplusghostlongdiv_simns<-simns
nullplusghostlongdiv_param<-simns_param

#-----------------------------------------------------------------------------------------------------------------------

# B) ghoste 
simufiles <- paste("/scratch/devel/hpawar/admix/abc/simul/test/modelcomp/14dec21/21jan22/ghoste/", list.files(path = "/scratch/devel/hpawar/admix/abc/simul/test/modelcomp/14dec21/21jan22/ghoste/", pattern="test.mc.ghoste_sim"), sep = "")


simns=list()
simns_param=list()
for (i in 1:length(simufiles)){
load(file=simufiles[[i]],verbose=T)
simns[[i]]<-simuresults
simns_param[[i]]<-inputsets
}

ghoste_stats<-matrix_function2()

ghoste_simns<-simns
ghoste_param<-simns_param

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------


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

#-----------------------------------------------------------------------------------------------------------------------

# 2) perform model choice


stats3mod<-rbind(as.data.frame(nullplusghostlongdiv_stats1),as.data.frame(ghoste_stats1))

models<-c(
replicate(2500, "nullplusghost"),
replicate(2500, "ghoste")
    )
#-----------------------------------------------------------------------------------------------------------------------

#  2.1) cross-validation for model selection
# if ABC can distinguish between the models


cv.modsel <- cv4postpr(models, stats3mod, nval=100, tol=0.05, method="neuralnet")

s <- summary(cv.modsel)

#-----------------------------------------------------------------------------------------------------------------------

# 2.2) model selection with postpr function

# read in summary stats from empirical data
load("/scratch/devel/hpawar/admix/abc/simul/test/ghost/3nov21/segrecalc/ghosttarget1_22nov",verbose=T)

# take sum then divide by 10 - these stats per 10kbp rather than 1kbp)
tfixsegs<-target1[8:22]

tfixsegs1<-sum(tfixsegs)/10
names(tfixsegs1)<-c("comb_fixsegs")

target2<-c(target1[1:7],tfixsegs1,target1[23:40])


# tolerance rate of 0.05%
testmod<-postpr(target2, models, stats3mod, tol=.05, method="neuralnet")

summary(testmod)

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
