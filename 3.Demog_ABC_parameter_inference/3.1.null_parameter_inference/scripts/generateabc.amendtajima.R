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

# transfer from cluster to local & plot
#scp -r hpawar@172.16.10.20:/scratch/devel/hpawar/admix/abc/results/test/11sep21simns/segrecalc/myabc_5apr22 /Users/harvi/Downloads/gorilla_abc/11sep21simns/segrecalc
#scp -r hpawar@172.16.10.20:/scratch/devel/hpawar/admix/abc/results/test/11sep21simns/segrecalc/sumstat2_5apr22 /Users/harvi/Downloads/gorilla_abc/11sep21simns/segrecalc
#g6n.tm4D2L73



# save output

# also run with alt way of calculating tajimas d
  # & compare these 2 ways of calculating *

 > myabc_5apr22
Call:
abc(target = target1, param = param, sumstat = sumstat2, tol = 0.005, 
    method = "neuralnet", transf = "logit", logit.bounds = prior_ranges, 
    numnet = 100)
Method:
Non-linear regression via neural networks
with correction for heteroscedasticity

Parameters:
w_lowl_t0, w_cros_t0, e_lowl_t0, e_moun_t0, e_lowl_t1, e_lowl_t2, e_moun_t3, e_moun_t3.1, t4, e_anc_t4, t5, w_lowl_t5, admix_w_e_t6, admix_e_w_t6, t7, w_anc_t7, t8, gor_anc

Statistics:
het_mu_WL, het_mu_WC, het_mu_EL, het_mu_EM, het_sd_WL, het_sd_EL, het_sd_EM, fixedsites_WL, fixedsites_WC, fixedsites_EL, fixedsites_EM, segsites_WL, segsites_WC, segsites_EL, segsites_EM, fixedperid_mu_WL, fixedperid_mu_WC, fixedperid_mu_EL, fixedperid_mu_EM, fixedperid_sd_WL, fixedperid_sd_EL, fixedperid_sd_EM, fst_WL.WC, fst_WL.EL, fst_WL.EM, fst_WC.EL, fst_WC.EM, fst_EL.EM, pi_mu_WL, pi_mu_EL, pi_mu_EM, pi_sd_WL, pi_sd_EL, pi_sd_EM, tajima_mu_WL, tajima_mu_EL, tajima_mu_EM, tajima_sd_WL, tajima_sd_EL, tajima_sd_EM

Total number of simulations 35543 

Number of accepted simulations:  178 

> summary(myabc_5apr22)
Call: 
abc(target = target1, param = param, sumstat = sumstat2, tol = 0.005, 
    method = "neuralnet", transf = "logit", logit.bounds = prior_ranges, 
    numnet = 100)
Data:
 abc.out$adj.values (178 posterior samples)
Weights:
 abc.out$weights

                       w_lowl_t0 w_cros_t0 e_lowl_t0 e_moun_t0 e_lowl_t1
Min.:                    24.0772    1.4641    0.9005    0.5445    0.0875
Weighted 2.5 % Perc.:    29.7673    5.7072    1.7186    0.8017    0.1107
Weighted Median:         66.0631   13.5755   20.9691    2.5202    0.2605
Weighted Mean:           65.7924   13.4357   19.5982    3.2900    0.2975
Weighted Mode:           66.2059   11.3421   26.8231    2.0618    0.2148
Weighted 97.5 % Perc.:   94.3898   19.5507   29.8148   10.1643    0.7918
Max.:                    97.8154   19.9128   29.9728   18.0409    5.5104
                       e_lowl_t2 e_moun_t3 e_moun_t3.1      t4 e_anc_t4      t5
Min.:                     1.9154    0.0102      0.9816  0.1427   1.1005  0.2294
Weighted 2.5 % Perc.:     4.2507    0.0162      1.7916  0.1489   1.8424  0.2799
Weighted Median:         19.6131    0.1742      7.1596  0.1894   6.5868  0.5177
Weighted Mean:           18.2514    0.3471      8.4715  0.1931   8.2575  0.4944
Weighted Mode:           22.0536    0.1022      5.3544  0.1856   4.4518  0.5413
Weighted 97.5 % Perc.:   24.6268    2.0573     18.7359  0.2379  24.0880  0.5928
Max.:                    24.9982    3.2491     19.9335  0.2491  28.3913  0.5992
                       w_lowl_t5 admix_w_e_t6 admix_e_w_t6      t7 w_anc_t7
Min.:                    19.6881       2.4435       0.2311  5.5413  86.3465
Weighted 2.5 % Perc.:    36.8678       5.4291       2.6746  5.6847  88.5089
Weighted Median:         46.9481      59.6369      35.1565  5.9324  96.9035
Weighted Mean:           45.7180      56.5273      36.8365  5.9057  95.8687
Weighted Mode:           48.0186      77.4935      11.4183  5.9463  97.5276
Weighted 97.5 % Perc.:   49.4398      99.1490      87.4340  5.9862  99.5861
Max.:                    49.9735      99.6090      99.2944  5.9934  99.9108
                            t8 gor_anc
Min.:                   5.3974 10.0928
Weighted 2.5 % Perc.:   8.7031 10.4331
Weighted Median:       11.6819 14.3404
Weighted Mean:         11.6275 14.8813
Weighted Mode:         11.7384 13.7842
Weighted 97.5 % Perc.: 14.2486 22.5235
Max.:                  14.7692 55.3827 
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
# the differences b/n thse are so small, unlikely to impact? stick with first way of calculating?

# 1) rm nas in tajima d calc
> summary(tajima_df)
       X1                X2                 X3                 X4        
 Min.   :-1.6629   Min.   :-1.78775   Min.   :-1.88787   Min.   :0.2808  
 1st Qu.:-0.2081   1st Qu.:-0.05374   1st Qu.:-0.03062   1st Qu.:0.6444  
 Median : 0.1549   Median : 0.21959   Median : 0.28464   Median :0.7569  
 Mean   : 0.1028   Mean   : 0.18114   Mean   : 0.22864   Mean   :0.7979  
 3rd Qu.: 0.4414   3rd Qu.: 0.46775   3rd Qu.: 0.56231   3rd Qu.:0.9279  
 Max.   : 1.4336   Max.   : 1.14613   Max.   : 1.34788   Max.   :1.6127  
       X5               X6        
 Min.   :0.4126   Min.   :0.3641  
 1st Qu.:0.8233   1st Qu.:0.8357  
 Median :0.9145   Median :0.9357  
 Mean   :0.9316   Mean   :0.9485  
 3rd Qu.:1.0117   3rd Qu.:1.0399  
 Max.   :1.7947   Max.   :1.9308

 # 2) set nas to 0  
> summary(tajima_df_1)
       X1                X2                 X3                 X4        
 Min.   :-1.6629   Min.   :-1.78775   Min.   :-1.88636   Min.   :0.2808  
 1st Qu.:-0.2081   1st Qu.:-0.05196   1st Qu.:-0.03062   1st Qu.:0.6444  
 Median : 0.1549   Median : 0.21776   Median : 0.28431   Median :0.7569  
 Mean   : 0.1028   Mean   : 0.17906   Mean   : 0.22745   Mean   :0.7979  
 3rd Qu.: 0.4414   3rd Qu.: 0.46468   3rd Qu.: 0.55936   3rd Qu.:0.9279  
 Max.   : 1.4336   Max.   : 1.14429   Max.   : 1.34626   Max.   :1.6127  
       X5               X6        
 Min.   :0.2756   Min.   :0.3660  
 1st Qu.:0.8216   1st Qu.:0.8357  
 Median :0.9128   Median :0.9356  
 Mean   :0.9262   Mean   :0.9473  
 3rd Qu.:1.0089   3rd Qu.:1.0397  
 Max.   :1.7168   Max.   :1.9297  


 # null ABC parameter inference setting nas to 0
 # set nas to 0s
sumstat3<-cbind(original, tajima_df_1)
sumstat3<-data.matrix(sumstat2)

lowerp<-c(3,0.1,0.1,0.1,0.01,1,0.01,0.1,0.14,0.1,0.19,5,0,0,0.1,10,1.5,10)
upperp<-c(100,20,30,20,20,25,5,20,0.25,30,0.6,50,100,100,6,100,15,100)
priordf<-cbind(lowerp,upperp)
prior_ranges<-data.matrix(priordf)


myabc_5apr22_1<-abc(target=target1,param=param,tol=0.005,sumstat=sumstat3,method="neuralnet",numnet=100,transf="logit",logit.bounds=prior_ranges)

123456789101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899100
123456789101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899100
Warning message:
All parameters are "logit" transformed. 

myabc_5apr22_1[[17]][[1]]<-colnames(param)
myabc_5apr22_1[[17]][[2]]<-names(target1)

myabc_5apr22_1
summary(myabc_5apr22_1)


save(myabc_5apr22_1,file=paste("/scratch/devel/hpawar/admix/abc/results/test/11sep21simns/segrecalc/myabc_5apr22_nasto0"))
save(sumstat3,file=paste("/scratch/devel/hpawar/admix/abc/results/test/11sep21simns/segrecalc/sumstat3_5apr22_nasto0"))


> myabc_5apr22_1
Call:
abc(target = target1, param = param, sumstat = sumstat3, tol = 0.005, 
    method = "neuralnet", transf = "logit", logit.bounds = prior_ranges, 
    numnet = 100)
Method:
Non-linear regression via neural networks
with correction for heteroscedasticity

Parameters:
w_lowl_t0, w_cros_t0, e_lowl_t0, e_moun_t0, e_lowl_t1, e_lowl_t2, e_moun_t3, e_moun_t3.1, t4, e_anc_t4, t5, w_lowl_t5, admix_w_e_t6, admix_e_w_t6, t7, w_anc_t7, t8, gor_anc

Statistics:
het_mu_WL, het_mu_WC, het_mu_EL, het_mu_EM, het_sd_WL, het_sd_EL, het_sd_EM, fixedsites_WL, fixedsites_WC, fixedsites_EL, fixedsites_EM, segsites_WL, segsites_WC, segsites_EL, segsites_EM, fixedperid_mu_WL, fixedperid_mu_WC, fixedperid_mu_EL, fixedperid_mu_EM, fixedperid_sd_WL, fixedperid_sd_EL, fixedperid_sd_EM, fst_WL.WC, fst_WL.EL, fst_WL.EM, fst_WC.EL, fst_WC.EM, fst_EL.EM, pi_mu_WL, pi_mu_EL, pi_mu_EM, pi_sd_WL, pi_sd_EL, pi_sd_EM, tajima_mu_WL, tajima_mu_EL, tajima_mu_EM, tajima_sd_WL, tajima_sd_EL, tajima_sd_EM

Total number of simulations 35543 

Number of accepted simulations:  178 

> summary(myabc_5apr22_1)
Call: 
abc(target = target1, param = param, sumstat = sumstat3, tol = 0.005, 
    method = "neuralnet", transf = "logit", logit.bounds = prior_ranges, 
    numnet = 100)
Data:
 abc.out$adj.values (178 posterior samples)
Weights:
 abc.out$weights

                       w_lowl_t0 w_cros_t0 e_lowl_t0 e_moun_t0 e_lowl_t1
Min.:                    18.6785    2.4617    1.1564    0.4127    0.0689
Weighted 2.5 % Perc.:    23.7580    6.9591    2.3079    0.7191    0.0760
Weighted Median:         64.7488   14.5582   20.8943    2.1578    0.2431
Weighted Mean:           64.2682   14.3019   19.7667    2.8079    0.2939
Weighted Mode:           64.3838   12.5476   26.5448    1.6796    0.1987
Weighted 97.5 % Perc.:   95.9855   19.6304   29.7763    7.6795    0.9448
Max.:                    98.5050   19.9268   29.9695   18.4970    6.5620
                       e_lowl_t2 e_moun_t3 e_moun_t3.1      t4 e_anc_t4      t5
Min.:                     3.0980    0.0102      1.5664  0.1409   0.6859  0.2140
Weighted 2.5 % Perc.:     8.0031    0.0150      2.9191  0.1446   1.4730  0.2839
Weighted Median:         21.9824    0.1150      9.9114  0.1839   5.3253  0.5381
Weighted Mean:           20.4935    0.2137     10.5400  0.1898   7.1841  0.5104
Weighted Mode:           23.0154    0.0671      8.0002  0.1781   3.3288  0.5605
Weighted 97.5 % Perc.:   24.6869    1.4780     19.1638  0.2420  23.0154  0.5976
Max.:                    24.9932    2.0316     19.9673  0.2493  25.9505  0.5998
                       w_lowl_t5 admix_w_e_t6 admix_e_w_t6      t7 w_anc_t7
Min.:                    27.1977       3.4615       0.1966  5.8086  84.3330
Weighted 2.5 % Perc.:    42.6356       6.7955       2.4079  5.8920  91.4236
Weighted Median:         48.2880      61.0343      33.0968  5.9778  98.1352
Weighted Mean:           47.5893      57.7143      35.0435  5.9679  97.0822
Weighted Mode:           48.9181      74.8586      11.6146  5.9832  98.6423
Weighted 97.5 % Perc.:   49.6917      99.2159      85.7601  5.9965  99.8200
Max.:                    49.9827      99.3275      98.7843  5.9978  99.9794
                            t8 gor_anc
Min.:                   5.6924 10.0954
Weighted 2.5 % Perc.:   9.5956 10.4147
Weighted Median:       12.7037 14.3641
Weighted Mean:         12.5051 14.9677
Weighted Mode:         12.8369 12.9148
Weighted 97.5 % Perc.: 14.5240 23.2983
Max.:                  14.8802 57.1283



# transfer from cluster to local & plot
#scp -r hpawar@172.16.10.20:/scratch/devel/hpawar/admix/abc/results/test/11sep21simns/segrecalc/myabc_5apr22_nasto0 /Users/harvi/Downloads/gorilla_abc/11sep21simns/segrecalc
#scp -r hpawar@172.16.10.20:/scratch/devel/hpawar/admix/abc/results/test/11sep21simns/segrecalc/sumstat3_5apr22_nasto0 /Users/harvi/Downloads/gorilla_abc/11sep21simns/segrecalc
#g6n.tm4D2L73

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
