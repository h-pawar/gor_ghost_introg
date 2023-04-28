#!/usr/bin/env
#module load R/4.0.1

# previous scripts, used to recalculate tajimas d 
#/scratch/devel/hpawar/admix/abc/simul/scripts/recalc.tajima.21mar22.R
#/scratch/devel/hpawar/admix/abc/simul/scripts/recalc.tajima.21mar22.arr
# & for reps with failed iter due to time limits
# /scratch/devel/hpawar/admix/abc/simul/scripts/recalc.tajima.failedreps.R  # also called by the recalc.tajima.21mar22.arr wrapper


# combine the recalculated tajimas d stats with those for the original stats 

#-----------------------------------------------------------------------------------------------------------------------

# 0) check 2 lists for equality - ie that the input sets are in the same order 

# order of parameters used to generate new taj d reps
#load(file=paste("/scratch/devel/hpawar/admix/abc/results/test/11sep21simns/segrecalc/null.simns_param"),verbose=T)
#taj_inputsets<-simns_param


# original order of parameters
#simufiles <- paste("/scratch/devel/hpawar/admix/abc/simul/test/11sep21/", list.files(path = "/scratch/devel/hpawar/admix/abc/simul/test/11sep21/", pattern="abc_sim"), sep = "")

#simns=list()
#simns_param=list()
#for (i in 1:length(simufiles)){
#load(file=simufiles[[i]],verbose=T)
#simns[[i]]<-simuresults
#simns_param[[i]]<-inputsets
#}


#identical(taj_inputsets, simns_param)
#[1] TRUE

# ie need to read in the reps in the simufiles in this order (rather than in i...)
#-----------------------------------------------------------------------------------------------------------------------


library(abc)
options(scipen=100)
require(data.table)
options(stringsAsFactors=F)
# read in param, sumstat1 & target1 from 6 Oct 21
load("/scratch/devel/hpawar/admix/abc/results/test/11sep21simns/segrecalc/param_6oct", verbose=T)
load("/scratch/devel/hpawar/admix/abc/results/test/11sep21simns/segrecalc/sumstat1_6oct", verbose=T)
load("/scratch/devel/hpawar/admix/abc/results/test/11sep21simns/segrecalc/target1_6oct", verbose=T)

# ncol(sumstat1)
#[1] 40

#  uninformative parameters already removed
# het_sd_WC (6)
#fixedperid_sd_WC (22)
# pi_mu_WC (32)
# pi_sd_WC (36)

# add in tajimas mu (Wl, EL, EM) & sd (WL, EL, EM) - ie replace the last 6 cols of sumstat1


#-----------------------------------------------------------------------------------------------------------------------
# now read in regenerated tajimas d stats
# remove iter of reps which lack stats for the other stats

# read in tajimas d recalculated for null demog simulations - Sun 27 Mar 2022 17:23:52 CEST
tajfiles <- paste("/scratch/devel/hpawar/admix/abc/simul/test/11sep21/tajima_21mar22/", list.files(path = "/scratch/devel/hpawar/admix/abc/simul/test/11sep21/tajima_21mar22/", pattern="taj_sim"), sep = "")
taj=list()
taj_param=list()
for (i in 1:length(tajfiles)){  
load(file=paste("/scratch/devel/hpawar/admix/abc/simul/test/11sep21/tajima_21mar22/taj_sim",i,"",sep=""),verbose=T)
taj[[i]]<-tajresults
taj_param[[i]]<-rep_param
}

 probl<-grep("Error", taj)
 problsimns=list()
 for (i in 1:length(probl)){
 problsimns[[i]]<-grep("Error", taj[[probl[i]]])
 }

#-----------------------------------------------------------------------------------------------------------------------
# check which reps were the probl reps (which were removed from previous stats)

# to do this - need to read in all prev stats - grep which had errors - then remove these from the new tajimas d reps 
simufiles <- paste("/scratch/devel/hpawar/admix/abc/simul/test/11sep21/", list.files(path = "/scratch/devel/hpawar/admix/abc/simul/test/11sep21/", pattern="abc_sim"), sep = "")

simns=list()
simns_param=list()
for (i in 1:length(simufiles)){
load(file=simufiles[[i]],verbose=T)
simns[[i]]<-simuresults
simns_param[[i]]<-inputsets
}

#-----------------------------------------------------------------------------------------------------------------------
# check again the 2 lists for equality # to check the results are in the right order (combining exactly the same reps)
# this shoudl give 'true' then can continue 

identical(taj_param, simns_param)
#[1] TRUE

#-----------------------------------------------------------------------------------------------------------------------


# simulations which failed (for the original stats)
  # ie need to remove these same reps from the tajimas d
probl<-grep("Error", simns)

problsimns=list()
for (i in 1:length(probl)){
problsimns[[i]]<-grep("Error", simns[[probl[i]]])
}

#probl
# [1]  14  74 108 191 316 320 363 364 483 518 527 533 558 606

#-----------------------------------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------------------------------
# process tajimas d 
format_probltaj_fun<-function(x) {
y<-probl[[x]]
hold_problout=list()
for(i in 1:length(taj[[y]]) ){
  if( (i %in% problsimns[[x]]) ){ next }
 hold_problout[[i]]<-(taj[[x]][[i]]) 
}
t<-hold_problout[lengths(hold_problout) != 0]
return(t)
}

out_subs=list()
for (i in 1:length(probl)){
out_subs[[i]]<-format_probltaj_fun(i)
}

taj_1<-taj
for (i in 1:length(probl)){
taj_1[[probl[[i]]]]<-out_subs[[i]]
}

# for one simn rep (1 of 700) - go over all iter, cbind the means & sds
# x = rep number, y = which way of calculating (1 or 2) 
#(either rm na, or set na to 0, where no nas among the 2500 windows, these 2 measures should be equal)

process_tajima_fun<-function(x,y) {
test=list()
for (i in 1:length(taj_1[[x]])){
test[[i]]<-cbind(
taj[[x]][[i]][[1]][[y]][,1],
taj[[x]][[i]][[2]][[y]][,1], 
taj[[x]][[i]][[3]][[y]][,1], 
taj[[x]][[i]][[1]][[y]][,2],
taj[[x]][[i]][[2]][[y]][,2],
taj[[x]][[i]][[3]][[y]][,2])
}
test_df<-data.frame(matrix(unlist(test), nrow=length(test), byrow=TRUE),stringsAsFactors=FALSE)
}

out_10=list()
for (i in 1:length(taj_1)){
# first way of calculating - rm nas
out_10[[i]]<-process_tajima_fun(i,1)
}

tajima_df<-do.call(rbind, out_10)


out_10_1=list()
for (i in 1:length(taj_1)){
# second way of calculating - set nas to 0
out_10_1[[i]]<-process_tajima_fun(i,2)
}
tajima_df_1<-do.call(rbind, out_10_1)


# WL mu, EL mu, EM mu, then the same for the sds 

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

original<-sumstat1[,c(1:34)]
sumstat2<-cbind(original, tajima_df)
sumstat2<-data.matrix(sumstat2)

lowerp<-c(3,0.1,0.1,0.1,0.01,1,0.01,0.1,0.14,0.1,0.19,5,0,0,0.1,10,1.5,10)
upperp<-c(100,20,30,20,20,25,5,20,0.25,30,0.6,50,100,100,6,100,15,100)
priordf<-cbind(lowerp,upperp)
prior_ranges<-data.matrix(priordf)


myabc_5apr22<-abc(target=target1,param=param,tol=0.005,sumstat=sumstat2,method="neuralnet",numnet=100,transf="logit",logit.bounds=prior_ranges)

#123456789101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899100
#123456789101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899100
#Warning message:
#All parameters are "logit" transformed. 

myabc_5apr22[[17]][[1]]<-colnames(param)
myabc_5apr22[[17]][[2]]<-names(target1)

myabc_5apr22
summary(myabc_5apr22)


save(myabc_5apr22,file=paste("/scratch/devel/hpawar/admix/abc/results/test/11sep21simns/segrecalc/myabc_5apr22"))
save(sumstat2,file=paste("/scratch/devel/hpawar/admix/abc/results/test/11sep21simns/segrecalc/sumstat2_5apr22"))


#-----------------------------------------------------------------------------------------------------------------------

setwd("/Users/harvi/Downloads/gorilla_abc/11sep21simns/segrecalc")
require(data.table)

# A) tol 0.005
load("/Users/harvi/Downloads/gorilla_abc/11sep21simns/segrecalc/myabc_5apr22",verbose=T)
load("/Users/harvi/Downloads/gorilla_abc/11sep21simns/segrecalc/param_6oct",verbose=T)
library(abc)
# histogram of the weighted posterior sample  : https://cran.r-project.org/web/packages/abc/vignettes/abcvignette.pdf
pdf("/Users/harvi/Downloads/gorilla_abc/11sep21simns/segrecalc/nullabc_posterior.5apr22.pdf") 
hist(myabc_5apr22)
dev.off() 
# generate the diagnostic plots of the estimation of the posterior distribution 
# a density plot of the prior distribution, a density plot of the posterior distribution estimated with and without regression- based correction, a scatter plot of the Euclidean distances as a function of the parameter values, and a normal Q-Q plot of the residuals from the regression
pdf("/Users/harvi/Downloads/gorilla_abc/11sep21simns/segrecalc/nullabc_diagnostic.5apr22.pdf") 
plot(myabc_5apr22, param)
dev.off() 
#-----------------------------------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------------------------------
# the differences b/n the results, using as input the two ways of calculating tajimas d are small 


 # null ABC parameter inference setting nas to 0
 # set nas to 0s
sumstat3<-cbind(original, tajima_df_1)
sumstat3<-data.matrix(sumstat2)

lowerp<-c(3,0.1,0.1,0.1,0.01,1,0.01,0.1,0.14,0.1,0.19,5,0,0,0.1,10,1.5,10)
upperp<-c(100,20,30,20,20,25,5,20,0.25,30,0.6,50,100,100,6,100,15,100)
priordf<-cbind(lowerp,upperp)
prior_ranges<-data.matrix(priordf)


myabc_5apr22_1<-abc(target=target1,param=param,tol=0.005,sumstat=sumstat3,method="neuralnet",numnet=100,transf="logit",logit.bounds=prior_ranges)

#123456789101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899100
#123456789101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899100
#Warning message:
#All parameters are "logit" transformed. 

myabc_5apr22_1[[17]][[1]]<-colnames(param)
myabc_5apr22_1[[17]][[2]]<-names(target1)

myabc_5apr22_1
summary(myabc_5apr22_1)


save(myabc_5apr22_1,file=paste("/scratch/devel/hpawar/admix/abc/results/test/11sep21simns/segrecalc/myabc_5apr22_nasto0"))
save(sumstat3,file=paste("/scratch/devel/hpawar/admix/abc/results/test/11sep21simns/segrecalc/sumstat3_5apr22_nasto0"))


setwd("/Users/harvi/Downloads/gorilla_abc/11sep21simns/segrecalc")
require(data.table)

# A) tol 0.005
load("/Users/harvi/Downloads/gorilla_abc/11sep21simns/segrecalc/myabc_5apr22_nasto0",verbose=T)
load("/Users/harvi/Downloads/gorilla_abc/11sep21simns/segrecalc/param_6oct",verbose=T)
library(abc)
# histogram of the weighted posterior sample  : https://cran.r-project.org/web/packages/abc/vignettes/abcvignette.pdf
pdf("/Users/harvi/Downloads/gorilla_abc/11sep21simns/segrecalc/nullabc_posterior.5apr22_nasto0.pdf") 
hist(myabc_5apr22_1)
dev.off() 
# generate the diagnostic plots of the estimation of the posterior distribution 
# a density plot of the prior distribution, a density plot of the posterior distribution estimated with and without regression- based correction, a scatter plot of the Euclidean distances as a function of the parameter values, and a normal Q-Q plot of the residuals from the regression
pdf("/Users/harvi/Downloads/gorilla_abc/11sep21simns/segrecalc/nullabc_diagnostic.5apr22_nasto0.pdf") 
plot(myabc_5apr22_1, param)
dev.off() 
#-----------------------------------------------------------------------------------------------------------------------
