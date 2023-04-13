# Fri 24 Feb 2023 10:13:26 CET

# regenerated ghostw parameter inference simulations - fixed typo in script

# ABC ghost parameter inference - sample all parameters from priors

# /scratch/devel/hpawar/admix/abc/simul/scripts/revisions_feb23/ghostw.rev.17feb23.arr
# /scratch/devel/hpawar/admix/abc/simul/scripts/revisions_feb23/ghostw.rev.17feb23.R

# paths to simulated data: /scratch/devel/hpawar/admix/abc/simul/test/ghost/revisions_8feb23/ghostw/rev_17feb23/
  # naming ghostw.rev.abc_sim700


#module load R/4.0.1
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------


# AMEND BELOW 


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




#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------


# simulated data 
simufiles <- paste("/scratch/devel/hpawar/admix/abc/simul/test/ghost/revisions_8feb23/ghostw/rev_17feb23/", list.files(path = "/scratch/devel/hpawar/admix/abc/simul/test/ghost/revisions_8feb23/ghostw/rev_17feb23/", pattern="ghostw.rev.abc_sim"), sep = "")



simns=list()
simns_param=list()
for (i in 1:length(simufiles)){
load(file=simufiles[[i]],verbose=T)
simns[[i]]<-simuresults
simns_param[[i]]<-inputsets
}


# check if any simulations failed 
# simulations which failed 
probl<-grep("Error", simns)

problsimns=list()
for (i in 1:length(probl)){
problsimns[[i]]<-grep("Error", simns[[probl[i]]])
}


#length(probl)
#[1] 129


#howmanyproblsimns=list()
#for (i in 1:length(probl)){
#howmanyproblsimns[[i]]<-length(problsimns[[i]])
#}

#-----------------------------------------------------------------------------------------------------------------------
# data structure of simns 700  simns -> 51 iteration each of 5000 windows -> with list of 6 summary stats (where tajima d 6th stat is itself a list)
#> str(simns)
#List of 700
# $ :List of 51
#  ..$ V1 :List of 6

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

#fst_process_fun<-function(afun,afun1) {
#out_1_700<-afun(1,700)
#rout_1_700<-as.data.frame(do.call(rbind,out_1_700))
#return(rout_1_700)
#}


#comb_process_fun<-function(afun,afun1,a,b) {
#out_1_700<-afun(1,700,a,b)
#rout_1_700<-as.data.frame(do.call(rbind,out_1_700))
#return(rout_1_700)
#}


#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
# rewritten tajimas - to incorporate failed reps 

#-----------------------------------------------------------------------------------------------------------------------
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


#-----------------------------------------------------------------------------------------------------------------------

out_taj_fun1<-function(x0,x) {
stats_s1=list()
for (i in x0:x){
stats_s1[[i]]<-process_tajima_fun(i)
}
return(stats_s1)
}

#-----------------------------------------------------------------------------------------------------------------------

format_problsimns_taj_fun<-function(x) {
y<-probl[[x]]
hold_problout=list()
for(i in 1:length(simns[[y]]) ){
  if( (i %in% problsimns[[x]]) ){ next }
 hold_problout[[i]]<-(simns[[y]][[i]][[6]]) 
}
t<-hold_problout[lengths(hold_problout) != 0]

probltest=list()
#for (i in 1:length(simns[[x]])){
for (i in 1:length(t)){
probltest[[i]]<-cbind(
#simns[[x]][[i]][[6]][[2]][,1],
t[[i]][[1]][,1],
t[[i]][[2]][,1],
t[[i]][[3]][,1],
t[[i]][[1]][,2],
t[[i]][[2]][,2],
t[[i]][[3]][,2])
}
probltest_df<-data.frame(matrix(unlist(probltest), nrow=length(probltest), byrow=TRUE),stringsAsFactors=FALSE)

return(probltest_df)}

#-----------------------------------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------------------------------

# probl
# [1]   2  11  21  24  33  45  51  53  55  69  74  76  79  84  97  98  99 100 107
#[20] 110 111 112 125 127 132 149 160 164 167 173 175 178 200 201 203 204 206 214
#[39] 217 218 222 244 247 252 269 274 276 279 280 284 287 297 298 315 322 326 328
#[58] 332 333 342 348 359 360 362 366 367 368 369 396 399 401 402 407 428 437 439

#-----------------------------------------------------------------------------------------------------------------------

#out_test=list()
#for (i in 1:length(simufiles)) {
  
 # skip_to_next <- FALSE
  
  # Note that print(b) fails since b doesn't exist
  
 # tryCatch(

#out_test[[i]]<-simu_fst_stats_fun(i)

#, error = function(e) { skip_to_next <<- TRUE})
  
#  if(skip_to_next) { next }     
#}

#-----------------------------------------------------------------------------------------------------------------------

all_fun<-function(afun){

out_test=list()
for (i in 1:length(simufiles)) {
  
  skip_to_next <- FALSE
  
  # Note that print(b) fails since b doesn't exist
  
  tryCatch(

out_test[[i]]<-afun(i)

, error = function(e) { skip_to_next <<- TRUE})
  
  if(skip_to_next) { next }     
}


out_test<-out_test[lengths(out_test) != 0]
o<-as.data.frame(do.call(rbind,out_test))

return(o)
}




s_fst<-all_fun(simu_fst_stats_fun)
s_taj<-all_fun(process_tajima_fun)

# works

#> nrow(s_fst)
#[1] 29121
#> nrow(s_taj)
#[1] 29121
#-----------------------------------------------------------------------------------------------------------------------

call_fun<-function(afun,a,b){

out_test=list()
for (i in 1:length(simufiles)) {
  
  skip_to_next <- FALSE
  
  # Note that print(b) fails since b doesn't exist
  
  tryCatch(

out_test[[i]]<-afun(i,a,b)

, error = function(e) { skip_to_next <<- TRUE})
  
  if(skip_to_next) { next }     
}


out_test<-out_test[lengths(out_test) != 0]
o<-as.data.frame(do.call(rbind,out_test))

return(o)
}


s_pi_mu<-call_fun(simu_pt_stats_fun,5,1)

s_pi_sd<-call_fun(simu_pt_stats_fun,5,2)

#simu_het_stats_fun(i,a,b)

s_het_mu<--call_fun(simu_het_stats_fun,1,1)

s_het_sd<--call_fun(simu_het_stats_fun,1,2)
s_het_sd[is.na(s_het_sd)] = 0


s_ifix_mu<--call_fun(simu_het_stats_fun,3,1)

s_ifix_sd<--call_fun(simu_het_stats_fun,3,2)
s_ifix_sd[is.na(s_ifix_sd)] = 0


s_popfix<--call_fun(simu_het_stats_fun,2,1)
s_popfix_perkb<-((s_popfix/(2500*40000))*1000)


s_popseg<--call_fun(simu_het_stats_fun,2,2)
s_popseg_perkb<-((s_popseg/(2500*40000))*1000)

#-----------------------------------------------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

sumstat<-cbind(s_het_mu,s_het_sd,
s_popfix_perkb,s_popseg_perkb,
s_ifix_mu,s_ifix_sd,
s_fst,
s_pi_mu,
s_pi_sd,
s_taj)

#str(sumstat)
#'data.frame': 18564 obs. of  44 variables:
# $ X1: num  -3.549 -2.096 -0.966 -2.091 -2.931 ...

# combine this with the other meaningful simulations from 15feb23 where tarchintrog < t8


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


# remove iterations with any failed rep - probl


pn_out1<-pn_out[-probl]



pn_format=list()
for (i in 1:length(pn_out1)){
pn_format[[i]]<-data.frame(matrix(unlist(pn_out1[[i]]), nrow=length(pn_out1[[i]]), byrow=TRUE),stringsAsFactors=FALSE)
}



param_df<-as.data.frame(do.call(rbind, pn_format))

# nrow(param_df)
#[1] 29121

#for (i in 1:length(probl)){
#pn_out[[probl[[i]]]]<-pn_out[[probl[[i]]]][-c(problsimns[[i]])]
#}



param<-data.matrix(param_df)
#parameter vector for each simulation vector, i.e. the input values you used to create each simulation (the "param"

# all parameters now being sampled from priors **
paraname=c("w_lowl_t0","w_cros_t0","e_lowl_t0","e_moun_t0","ghost_t0","e_lowl_t1","e_lowl_t2","e_moun_t3","e_moun_t3.1","t4","e_anc_t4","t_archintrog","archaicintrog","t5","w_lowl_t5","admix_w_e_t6","admix_e_w_t6","t7","w_anc_t7","t8","gor_anc","t9","gor_ghost_anc","id")

colnames(param)<-paraname
head(param)

# remove the id parameter - ie now parameter 24
param<-param[,c(1:23)]



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
#-----------------------------------------------------------------------------------------------------------------------


# remove ghost_t0 from the parameter search
param1<-param[,-c(5)]
#head(param1)

#-----------------------------------------------------------------------------------------------------------------------


param1<-data.matrix(param1)  
#parameter vector for each simulation vector, i.e. the input values you used to create each simulation (the "param"


#-----------------------------------------------------------------------------------------------------------------------

# should also force these posteriors to be within the priors, same logit transf

# W priors
lowerp<-c(3,0.1,0.1,0.1,0.1,0.01,1,0.01,0.1,0.14,0.1,6.0025,0,0.19,5,0,0,0.1,10,1.5,10,15.0025,10)
upperp<-c(100,20,30,20,100,20,25,5,20,0.25,30,14.9975,100,0.6,50,100,100,6,100,15,100,50,100)

# remove ghostt0 prior from this * - remove entry 5 
# without ghost_t0 priors as follows
lowerp1<-lowerp[-5]
upperp1<-upperp[-5]


priordf<-cbind(lowerp1,upperp1)
prior_ranges<-data.matrix(priordf)

# Which you feed into the ABC like this:
ghostwabcrev_24feb23<-abc(target=target1,param=param1,tol=0.005,sumstat=sumstat1,method="neuralnet",numnet=100,transf="logit",logit.bounds=prior_ranges)


# The "tol" parameter is unknown and might be tested a bit.

> ghostwabcrev_24feb23<-abc(target=target1,param=param1,tol=0.005,sumstat=sumstat1,method="neuralnet",123456789101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899100
123456789101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899100
Warning message:
All parameters are "logit" transformed. 

#-----------------------------------------------------------------------------------------------------------------------

ghostwabcrev_24feb23[[17]][[1]]<-colnames(param1)
ghostwabcrev_24feb23[[17]][[2]]<-names(target1)

ghostwabcrev_24feb23
summary(ghostwabcrev_24feb23)


# AMEND PATHS **
# mkdir -p /scratch/devel/hpawar/admix/abc/simul/test/ghost/revisions_8feb23/ghostw/segrecalc


save(ghostwabcrev_24feb23,file=paste("/scratch/devel/hpawar/admix/abc/simul/test/ghost/revisions_8feb23/ghostw/segrecalc/ghostwabcrev_24feb23"))
save(param1,file=paste("/scratch/devel/hpawar/admix/abc/simul/test/ghost/revisions_8feb23/ghostw/segrecalc/ghostwabc_param_24feb23"))
save(sumstat1,file=paste("/scratch/devel/hpawar/admix/abc/simul/test/ghost/revisions_8feb23/ghostw/segrecalc/ghostwabc_sumstat1_24feb23"))


#-----------------------------------------------------------------------------------------------------------------------

> ghostwabcrev_24feb23
Call:
abc(target = target1, param = param1, sumstat = sumstat1, tol = 0.005, 
    method = "neuralnet", transf = "logit", logit.bounds = prior_ranges, 
    numnet = 100)
Method:
Non-linear regression via neural networks
with correction for heteroscedasticity

Parameters:
w_lowl_t0, w_cros_t0, e_lowl_t0, e_moun_t0, e_lowl_t1, e_lowl_t2, e_moun_t3, e_moun_t3.1, t4, e_anc_t4, t_archintrog, archaicintrog, t5, w_lowl_t5, admix_w_e_t6, admix_e_w_t6, t7, w_anc_t7, t8, gor_anc, t9, gor_ghost_anc

Statistics:
het_mu_WL, het_mu_WC, het_mu_EL, het_mu_EM, het_sd_WL, het_sd_EL, het_sd_EM, fixedsites_WL, fixedsites_WC, fixedsites_EL, fixedsites_EM, segsites_WL, segsites_WC, segsites_EL, segsites_EM, fixedperid_mu_WL, fixedperid_mu_WC, fixedperid_mu_EL, fixedperid_mu_EM, fixedperid_sd_WL, fixedperid_sd_EL, fixedperid_sd_EM, fst_WL.WC, fst_WL.EL, fst_WL.EM, fst_WC.EL, fst_WC.EM, fst_EL.EM, pi_mu_WL, pi_mu_EL, pi_mu_EM, pi_sd_WL, pi_sd_EL, pi_sd_EM, tajima_mu_WL, tajima_mu_EL, tajima_mu_EM, tajima_sd_WL, tajima_sd_EL, tajima_sd_EM

Total number of simulations 29121 

Number of accepted simulations:  146 

> summary(ghostwabcrev_24feb23)
Call: 
abc(target = target1, param = param1, sumstat = sumstat1, tol = 0.005, 
    method = "neuralnet", transf = "logit", logit.bounds = prior_ranges, 
    numnet = 100)
Data:
 abc.out$adj.values (146 posterior samples)
Weights:
 abc.out$weights

                       w_lowl_t0 w_cros_t0 e_lowl_t0 e_moun_t0 e_lowl_t1
Min.:                    48.0090    4.3772    0.1384    0.1259    0.0439
Weighted 2.5 % Perc.:    59.4414    9.2000    2.9353    0.1990    0.1172
Weighted Median:         83.9216   15.8840   21.0246    0.4062    0.9848
Weighted Mean:           83.7779   15.4008   19.3416    1.0271    1.8585
Weighted Mode:           82.6407   17.0895   21.9791    0.3230    0.8020
Weighted 97.5 % Perc.:   99.0277   19.1504   29.3188    7.3289    6.7332
Max.:                    99.9834   19.8506   29.9473   15.0440   17.6185
                       e_lowl_t2 e_moun_t3 e_moun_t3.1      t4 e_anc_t4
Min.:                     1.1244    0.4283      2.6215  0.1404   0.2270
Weighted 2.5 % Perc.:     1.3326    1.2398      2.7852  0.1441   0.4994
Weighted Median:          2.4187    3.6617      9.8074  0.1990   1.6466
Weighted Mean:            2.7443    3.3848      9.8021  0.1970   2.7591
Weighted Mode:            1.9429    3.7723     10.9278  0.2068   1.3278
Weighted 97.5 % Perc.:    6.2619    4.7630     17.3300  0.2456  11.3819
Max.:                    13.1267    4.9563     18.0256  0.2495  27.0995
                       t_archintrog archaicintrog      t5 w_lowl_t5
Min.:                        6.0095        0.0438  0.1902   43.7373
Weighted 2.5 % Perc.:        6.0188        0.4455  0.1924   47.5977
Weighted Median:             6.0815       17.0257  0.3919   49.6499
Weighted Mean:               6.0856       24.6037  0.3905   49.5033
Weighted Mode:               6.0790        9.1757  0.2976   49.8129
Weighted 97.5 % Perc.:       6.2074       79.2994  0.5925   49.9468
Max.:                        6.3867       99.5935  0.5984   49.9769
                       admix_w_e_t6 admix_e_w_t6      t7 w_anc_t7      t8
Min.:                        4.6747       1.4095  0.1801  19.1647  3.0827
Weighted 2.5 % Perc.:        7.3217       3.7917  0.3496  72.5528  3.8799
Weighted Median:            53.3923      21.0984  0.9728  91.8882  6.1259
Weighted Mean:              54.5402      24.6831  1.1704  90.1263  6.6798
Weighted Mode:              49.8241      17.8047  0.8037  96.8961  5.5927
Weighted 97.5 % Perc.:      99.4208      65.5610  2.9163  99.9272 11.3642
Max.:                       99.8703      89.9105  4.1209  99.9852 14.3902
                       gor_anc      t9 gor_ghost_anc
Min.:                  10.0400 21.6390       10.7712
Weighted 2.5 % Perc.:  10.1183 24.5863       10.7955
Weighted Median:       11.3216 41.3720       25.4882
Weighted Mean:         13.1733 40.6067       27.4116
Weighted Mode:         10.8907 45.1715       18.4293
Weighted 97.5 % Perc.: 27.9935 49.0967       69.0760
Max.:                  57.8614 49.6542       98.0403
> 

  
#-----------------------------------------------------------------------------------------------------------------------


# histogram of the weighted posterior sample  : https://cran.r-project.org/web/packages/abc/vignettes/abcvignette.pdf
pdf("/scratch/devel/hpawar/admix/abc/simul/test/ghost/revisions_8feb23/ghostw/segrecalc/ghostwabcrev_posterior_24feb23.pdf") 

hist(ghostwabcrev_24feb23)


# generate the diagnostic plots of the estimation of the posterior distribution 
# a density plot of the prior distribution, a density plot of the posterior distribution estimated with and without regression- based correction, a scatter plot of the Euclidean distances as a function of the parameter values, and a normal Q-Q plot of the residuals from the regression
pdf("/scratch/devel/hpawar/admix/abc/simul/test/ghost/revisions_8feb23/ghostw/segrecalc/ghostwabcrev_diagnostic_24feb23.pdf") 
plot(ghostwabcrev_24feb23, param1)
dev.off() 

#~/Downloads/gor_ghost_revisions_feb23


scp -r -oHostKeyAlgorithms=+ssh-rsa  hpawar@172.16.10.21:/scratch/devel/hpawar/admix/abc/simul/test/ghost/revisions_8feb23/ghostw/segrecalc/ghostwabcrev_diagnostic_24feb23.pdf ~/Downloads/gor_ghost_revisions_feb23
scp -r -oHostKeyAlgorithms=+ssh-rsa  hpawar@172.16.10.21:/scratch/devel/hpawar/admix/abc/simul/test/ghost/revisions_8feb23/ghostw/segrecalc/ghostwabcrev_posterior_24feb23.pdf ~/Downloads/gor_ghost_revisions_feb23
