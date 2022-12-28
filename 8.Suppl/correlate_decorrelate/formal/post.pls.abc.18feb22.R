# Thu 17 Feb 2022 18:19:30 CET
# MK
#my interpretation of the discussions with Oscar was that one could try ABC inference on reduced dimensions of the PLS scores 
#(or PCA scores, I think he mentioned first), or, as you are trying, using transformed statistics. 
#For that, the simulated and real data statistics need to be treated the same way. 
#I was playing around a bit, and paste some code below that may be worth trying on the full dataset, 
#and do a side-by-side comparison to ABC inference on transformed values (if the latter works at all).

# ie Perform ABC inference using the decorrelated statistics (ie directly on the PLS components)

#-----------------------------------------------------------------------------------------------------------------------

#module load R/4.0.1

library("pls"); library("MASS");
library(abc)

# processing from generateabc.6oct.R
options(scipen=100)
require(data.table)
options(stringsAsFactors=F)

# 1) load summary stats and parameters
load("/scratch/devel/hpawar/admix/abc/results/test/11sep21simns/segrecalc/sumstat1_6oct")
sstat1<-sumstat1[,-c(6,22,32,36)]
load("/scratch/devel/hpawar/admix/abc/results/test/11sep21simns/segrecalc/param_6oct")
parm=param
load("/scratch/devel/hpawar/admix/abc/results/test/11sep21simns/segrecalc/target1_6oct")
target1<-target1[-c(6,22,32,36)]


# 2) do the transformation steps - important here is that you are including the target stats to also get these in the same format
sumstat1<-sstat1
param=parm
trget2<-target1
myMax<-c(); myMin<-c(); lambda<-c(); myGM<-c();
for(i in 1:ncol(sumstat1)){
  myMax<-c(myMax, max(sumstat1[,i])); myMin<-c(myMin, min(sumstat1[,i]));
  sumstat1[,i]<-1+(sumstat1[,i] -myMin[i])/(myMax[i]-myMin[ i]);
  trget2[i]<-1+(trget2[i] -myMin[i])/(myMax[i]-myMin[ i]);
  print(c(range(sumstat1[,i]),length(which(sumstat1[,i]==NA)),mean(sumstat1[,i])))
}

#transform statistics via boxcox
for(i in 1:ncol(sumstat1)){
print(i)
d<-cbind(data.frame(sumstat1[,i]), data.frame(param));
mylm<-lm(as.formula(d), data=d);
myboxcox<-boxcox(mylm, lambda=seq(-20,100,1/10), interp=T, eps=1/50);
lambda<-c(lambda, myboxcox$x[myboxcox$y==max(myboxcox$y)]);
myGM<-c(myGM, mean(exp(log(sumstat1[,i]))));
}

sstat<-sumstat1
#standardize the BC-stat
myBCMeans<-c(); myBCSDs<-c();
for(i in 1:ncol(sstat)){
    sstat[,i]<-(sstat[,i] ^lambda[i] - 1)/(lambda[i]*myGM[i] ^(lambda[i]-1));
    trget2[i]<-(trget2[i] ^lambda[i] - 1)/(lambda[i]*myGM[i] ^(lambda[i]-1));
    myBCSDs<-c(myBCSDs, sd(sstat[,i]));
    myBCMeans<-c(myBCMeans, mean(sstat[,i]));
    sstat[,i]<-(sstat[,i] -myBCMeans[i])/myBCSDs[i];
    trget2[i]<-(trget2[i]-myBCMeans[i])/myBCSDs[i]
  }
#perform pls -- you have already done this step with the full data!
##  myPlsr<-plsr(as.matrix(param)~as.matrix(sstat), scale=F, validation='LOO');

#  read in the pls results i already have for 10 components of 36 components??

#save(myPlsr,file=paste("/scratch/devel/hpawar/admix/abc/simul/test/modelcomp/out/pls/plsr_woutvalid_10comp.17feb22"))

load("/scratch/devel/hpawar/admix/abc/simul/test/modelcomp/out/pls/plsr_woutvalid_10comp.17feb22", verbose=T)

# 10 components - explains 95% of the variance
#> cumsum(explvar(myPlsr))
#  Comp 1   Comp 2   Comp 3   Comp 4   Comp 5   Comp 6   Comp 7   Comp 8 
#38.27876 53.68219 81.18623 85.92309 89.36204 90.93559 93.52908 94.11261 
#  Comp 9  Comp 10 
#94.94002 95.67793 

# 2) get the scores for each simulation and each dimension, and predict the scores for the real statistics
scrs<-matrix(myPlsr$scores,ncol=ncol(myPlsr$scores))  
trget3<-predict(myPlsr,newdata=matrix(trget2,nrow=1),type="scores")
# then do PCA on the first 10 dimensions
myabc2<-abc(target=trget3[1:10],param=param,tol=0.01,sumstat=scrs[,1:10],method="neuralnet",numnet=100)
Warning messages:
1: All parameters are "none" transformed. 
2: In abc(target = trget3[1:10], param = param, tol = 0.01, sumstat = scrs[,  :
  No summary statistics names are given, using S1, S2, .

> head(trget3)
       Comp 1    Comp 2   Comp 3  Comp 4    Comp 5    Comp 6    Comp 7   Comp 8
[1,] 4.282227 0.8279513 1.122744 2.18467 -2.149207 -4.167857 0.7656653 0.589975
       Comp 9   Comp 10
[1,] 2.599719 0.4639265

myabc2[[17]][[1]]<-colnames(param)
myabc2[[17]][[2]]<-names(trget3)

myabc2
summary(myabc2)
save(myabc2,file=paste("/scratch/devel/hpawar/admix/abc/simul/test/modelcomp/out/pls/abc_pls_18feb22"))
save(trget3,file=paste("/scratch/devel/hpawar/admix/abc/simul/test/modelcomp/out/pls/trget3_18feb22"))

save(param,file=paste("/scratch/devel/hpawar/admix/abc/simul/test/modelcomp/out/pls/param_pls_18feb22"))

> myabc2
Call:
abc(target = trget3[1:10], param = param, sumstat = scrs[, 1:10], 
    tol = 0.01, method = "neuralnet", numnet = 100)
Method:
Non-linear regression via neural networks
with correction for heteroscedasticity

Parameters:
w_lowl_t0, w_cros_t0, e_lowl_t0, e_moun_t0, e_lowl_t1, e_lowl_t2, e_moun_t3, e_moun_t3.1, t4, e_anc_t4, t5, w_lowl_t5, admix_w_e_t6, admix_e_w_t6, t7, w_anc_t7, t8, gor_anc

Statistics:


Total number of simulations 35543 

Number of accepted simulations:  356 

> summary(myabc2)
Call: 
abc(target = trget3[1:10], param = param, sumstat = scrs[, 1:10], 
    tol = 0.01, method = "neuralnet", numnet = 100)
Data:
 abc.out$adj.values (356 posterior samples)
Weights:
 abc.out$weights

                       w_lowl_t0 w_cros_t0 e_lowl_t0 e_moun_t0 e_lowl_t1
Min.:                    65.5522   10.7702    4.1804    0.5618    2.3660
Weighted 2.5 % Perc.:    72.4135   14.0093    5.6060    1.5016    2.7025
Weighted Median:         89.2094   20.2728   16.8704    2.9892    4.6104
Weighted Mean:           89.1134   20.2399   16.9132    4.5613    5.6189
Weighted Mode:           88.8603   21.6758   16.7661    2.5874    4.3340
Weighted 97.5 % Perc.:  105.7322   27.0298   28.3043   15.8667   14.0589
Max.:                   116.2874   32.9661   29.8577   26.0943   26.8509
                       e_lowl_t2 e_moun_t3 e_moun_t3.1       t4 e_anc_t4
Min.:                     5.8449   -1.2246      1.9362   0.1537  16.2305
Weighted 2.5 % Perc.:     6.7670   -0.8588      2.8451   0.1658  17.6658
Weighted Median:         18.3987    1.4849     10.5757   0.2179  22.9690
Weighted Mean:           17.9128    1.6603     11.2607   0.2187  23.2499
Weighted Mode:           21.9545    0.1787      8.1417   0.1867  23.0287
Weighted 97.5 % Perc.:   28.9096    4.8840     20.2628   0.2714  30.7557
Max.:                    30.7584    5.5802     23.9765   0.2906  36.7354
                             t5 w_lowl_t5 admix_w_e_t6 admix_e_w_t6       t7
Min.:                    0.2136   -4.3325     -28.1960     -28.5556   4.9132
Weighted 2.5 % Perc.:    0.2453   21.7173     -16.9018     -19.5507   5.8338
Weighted Median:         0.4769   45.6855      17.0756      39.6015   7.7345
Weighted Mean:           0.4590   44.8998      25.0986      41.0881   7.6604
Weighted Mode:           0.5377   52.9118       3.1457      38.3964   7.9567
Weighted 97.5 % Perc.:   0.6301   66.7046     113.9061     102.2018   9.2722
Max.:                    0.6526   70.0695     165.7739     146.6777   9.9762
                       w_anc_t7       t8  gor_anc
Min.:                   60.5097   7.5132   8.0778
Weighted 2.5 % Perc.:   70.3390   9.3448  20.1600
Weighted Median:        86.2332  12.8414  35.0938
Weighted Mean:          85.1152  12.7708  35.5836
Weighted Mode:          88.2899  12.9198  33.7290
Weighted 97.5 % Perc.:  96.5090  16.0228  48.5139
Max.:                  102.1845  17.6026  82.6521

# transfer to local & plot

#scp -r hpawar@172.16.10.20:/scratch/devel/hpawar/admix/abc/simul/test/modelcomp/out/pls/abc_pls_18feb22  /Users/harvi/Downloads/gorilla_abc/modelchoice/pls
#scp -r hpawar@172.16.10.20:/scratch/devel/hpawar/admix/abc/simul/test/modelcomp/out/pls/trget3_18feb22  /Users/harvi/Downloads/gorilla_abc/modelchoice/pls

# need to output the param file **
#scp -r hpawar@172.16.10.20:/scratch/devel/hpawar/admix/abc/simul/test/modelcomp/out/pls/param_pls_18feb22 /Users/harvi/Downloads/gorilla_abc/modelchoice/pls


#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
# plotting these results

# then plot on local

#scp -r hpawar@172.16.10.20:/scratch/devel/hpawar/admix/abc/simul/test/modelcomp/out/pls/abc_pls_18feb22  /Users/harvi/Downloads/gorilla_abc/modelchoice/pls
#scp -r hpawar@172.16.10.20:/scratch/devel/hpawar/admix/abc/simul/test/modelcomp/out/pls/trget3_18feb22  /Users/harvi/Downloads/gorilla_abc/modelchoice/pls

# on local
setwd("/Users/harvi/Downloads/gorilla_abc/modelchoice/pls")
require(data.table)

# A) tol 0.005
load("/Users/harvi/Downloads/gorilla_abc/modelchoice/pls/abc_pls_18feb22",verbose=T)
load("/Users/harvi/Downloads/gorilla_abc/modelchoice/pls/param_pls_18feb22",verbose=T)
library(abc)
# histogram of the weighted posterior sample  : https://cran.r-project.org/web/packages/abc/vignettes/abcvignette.pdf
pdf("/Users/harvi/Downloads/gorilla_abc/modelchoice/pls/abc_pls_posterior.18feb22.pdf") 
hist(myabc2)
dev.off() 
# generate the diagnostic plots of the estimation of the posterior distribution 
# a density plot of the prior distribution, a density plot of the posterior distribution estimated with and without regression- based correction, a scatter plot of the Euclidean distances as a function of the parameter values, and a normal Q-Q plot of the residuals from the regression
pdf("/Users/harvi/Downloads/gorilla_abc/modelchoice/pls/abc_pls_diagnostic.18feb22.pdf") 
plot(myabc2, param)
dev.off() 

