#Fri 18 Feb 2022 14:15:58 CET

# 1) Calculate how many PLS components explain 99% of the variance
# 2) Perform ABC inference using these 17 components (rather than the optimal 10 components - as performed in post.pls.abc.18feb22.R)

#MK
#Another idea: Could you try the same with the PLS scores that explain 99% of the variance?
#cumsum(explvar(myPlsr)) will give a string with the explained variance; 
#in your file you only have the first 10 components, which explain ~95% (which already seems a good cutoff). 
#But it would be interesting to see the ABC results when going to 99%.
#-----------------------------------------------------------------------------------------------------------------------
#module load R/4.0.1

library("pls"); library("MASS");

# processing from generateabc.6oct.R
options(scipen=100)
require(data.table)
options(stringsAsFactors=F)
load("/scratch/devel/hpawar/admix/abc/simul/test/modelcomp/out/pls/plsr_woutvalid_10comp.17feb22", verbose=T)
cumsum(explvar(myPlsr))

> cumsum(explvar(myPlsr))
  Comp 1   Comp 2   Comp 3   Comp 4   Comp 5   Comp 6   Comp 7   Comp 8 
38.27876 53.68219 81.18623 85.92309 89.36204 90.93559 93.52908 94.11261 
  Comp 9  Comp 10 
94.94002 95.67793 

# rerun pls increasing the number of components, until hit 99% explained? 
# but isn't the point dimension reduction?

#-----------------------------------------------------------------------------------------------------------------------

load(file=paste("/scratch/devel/hpawar/admix/abc/simul/test/modelcomp/out/pls/sstat_16feb22"), verbose=T)
# loads 'test' which represents sstat
load(file=paste("/scratch/devel/hpawar/admix/abc/simul/test/modelcomp/out/pls/param_16feb22"), verbose=T)
# loads parameters
load(file=paste("/scratch/devel/hpawar/admix/abc/simul/test/modelcomp/out/pls/BCvals_16feb22"), verbose=T)
# loads myBCvals

myBCSDs<-myBCvals[[1]]
myBCMeans<-myBCvals[[2]]

load(file=paste("/scratch/devel/hpawar/admix/abc/simul/test/modelcomp/out/pls/maxminvals_16feb22"), verbose=T)
# loads maxminvals

myMax<-maxminvals[[1]]
myMin<-maxminvals[[2]]

load(file=paste("/scratch/devel/hpawar/admix/abc/simul/test/modelcomp/out/pls/lambdagmvals_16feb22"), verbose=T)
# loads lambdagmvals
lambda<-lambdagmvals[[1]]
myGM<-lambdagmvals[[2]]


# removing same cols as removed from test (so all of length 36)
myBCSDs1<-myBCSDs[c(1:5,7:21,23:38)]
myBCMeans1<-myBCMeans[c(1:5,7:21,23:38)]

# remove also from myMax, myMin
myMax1<-myMax[c(1:5,7:21,23:38)]
myMin1<-myMin[c(1:5,7:21,23:38)]

# with 10 pls components (seems to be optimal number from /Volumes/"Ultra USB 3.0"/IBE/further.analysis.feb.2020/gorillas/abc/simul_postabc/abc+ghost/modelcomp/uncorrelate.stats.26jan22.R)
numComp=12

# perform pls # without validation
myPlsr<-plsr(as.matrix(param)~as.matrix(test), scale=F,ncomp=numComp)
 cumsum(explvar(myPlsr))
 Comp 1   Comp 2   Comp 3   Comp 4   Comp 5   Comp 6   Comp 7   Comp 8 
38.27876 53.68219 81.18623 85.92309 89.36204 90.93559 93.52908 94.11261 
  Comp 9  Comp 10  Comp 11  Comp 12 
94.94002 95.67793 96.19082 96.92482 

numComp=15

# perform pls # without validation
myPlsr<-plsr(as.matrix(param)~as.matrix(test), scale=F,ncomp=numComp)
 cumsum(explvar(myPlsr))
   cumsum(explvar(myPlsr))
  Comp 1   Comp 2   Comp 3   Comp 4   Comp 5   Comp 6   Comp 7   Comp 8 
38.27876 53.68219 81.18623 85.92309 89.36204 90.93559 93.52908 94.11261 
  Comp 9  Comp 10  Comp 11  Comp 12  Comp 13  Comp 14  Comp 15 
94.94002 95.67793 96.19082 96.92482 97.57554 98.16050 98.48517 

numComp=20
# perform pls # without validation
myPlsr<-plsr(as.matrix(param)~as.matrix(test), scale=F,ncomp=numComp)
 cumsum(explvar(myPlsr))
 Comp 1   Comp 2   Comp 3   Comp 4   Comp 5   Comp 6   Comp 7   Comp 8 
38.27876 53.68219 81.18623 85.92309 89.36204 90.93559 93.52908 94.11261 
  Comp 9  Comp 10  Comp 11  Comp 12  Comp 13  Comp 14  Comp 15  Comp 16 
94.94002 95.67793 96.19082 96.92482 97.57554 98.16050 98.48517 98.73232 
 Comp 17  Comp 18  Comp 19  Comp 20 
98.91655 99.17081 99.36651 99.55468

#-----------------------------------------------------------------------------------------------------------------------
# run abc inference with 17 pls components?

# me
#to get to 99% I have to increase the PLS components to 17.
#Is it worth running the ABC on this? Im not sure, if this is again just increasing the dimensionality? 

# MK
#yes, I would try it. It is still reducing the number of dimensions by more than half, 
#so if we are mostly concerned about removing the correlations, that should be dealt with here. 
#Meanwhile, we still retain most of the information. This is just to explore, so I still think we should focus on the test with 10 components.
#-----------------------------------------------------------------------------------------------------------------------

numComp=17
# perform pls # without validation, but with 17 components (which explain 98.9% of the variance)
myPlsr<-plsr(as.matrix(param)~as.matrix(test), scale=F,ncomp=numComp)

#-----------------------------------------------------------------------------------------------------------------------
# prep for abc - following same approach as post.pls.abc.18feb22.R

library(abc)

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

# 2) get the scores for each simulation and each dimension, and predict the scores for the real statistics
scrs<-matrix(myPlsr$scores,ncol=ncol(myPlsr$scores))  
trget3<-predict(myPlsr,newdata=matrix(trget2,nrow=1),type="scores")
# then do PCA on the first 10 dimensions
myabc2_17<-abc(target=trget3[1:10],param=param,tol=0.01,sumstat=scrs[,1:10],method="neuralnet",numnet=100)


myabc2_17
summary(myabc2_17)
save(myabc2_17,file=paste("/scratch/devel/hpawar/admix/abc/simul/test/modelcomp/out/pls/abc_17plscomp_18feb22"))
save(trget3,file=paste("/scratch/devel/hpawar/admix/abc/simul/test/modelcomp/out/pls/trget3_17plscomp_18feb22"))

save(param,file=paste("/scratch/devel/hpawar/admix/abc/simul/test/modelcomp/out/pls/param_17plscomp_18feb22"))


> myabc2_17<-abc(target=trget3[1:10],param=param,tol=0.01,sumstat=scrs[,1:10],method="neuralnet",numnet=100)
123456789101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899100
123456789101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899100
Warning messages:
1: All parameters are "none" transformed. 
2: In abc(target = trget3[1:10], param = param, tol = 0.01, sumstat = scrs[,  :
  No summary statistics names are given, using S1, S2, ...

  > myabc2_17
Call:
abc(target = trget3[1:10], param = param, sumstat = scrs[, 1:10], 
    tol = 0.01, method = "neuralnet", numnet = 100)
Method:
Non-linear regression via neural networks
with correction for heteroscedasticity

Parameters:
w_lowl_t0, w_cros_t0, e_lowl_t0, e_moun_t0, e_lowl_t1, e_lowl_t2, e_moun_t3, e_moun_t3.1, t4, e_anc_t4, t5, w_lowl_t5, admix_w_e_t6, admix_e_w_t6, t7, w_anc_t7, t8, gor_anc

Statistics:
S1, S2, S3, S4, S5, S6, S7, S8, S9, S10

Total number of simulations 35543 

Number of accepted simulations:  356 

> summary(myabc2_17)
Call: 
abc(target = trget3[1:10], param = param, sumstat = scrs[, 1:10], 
    tol = 0.01, method = "neuralnet", numnet = 100)
Data:
 abc.out$adj.values (356 posterior samples)
Weights:
 abc.out$weights

                       w_lowl_t0 w_cros_t0 e_lowl_t0 e_moun_t0 e_lowl_t1
Min.:                    60.7877   11.7725    3.6647    0.7833    3.0091
Weighted 2.5 % Perc.:    69.6086   15.8585    5.4786    1.5169    3.1980
Weighted Median:         90.1757   21.5939   16.7151    2.9977    5.0889
Weighted Mean:           90.0941   21.5206   16.6991    4.5018    5.9511
Weighted Mode:           88.6675   23.0211   16.5648    2.5157    4.7036
Weighted 97.5 % Perc.:  110.3457   27.4912   28.6469   15.8308   14.5810
Max.:                   120.4449   32.3652   29.9810   25.6524   26.8977
                       e_lowl_t2 e_moun_t3 e_moun_t3.1       t4 e_anc_t4
Min.:                     5.9334   -1.1956      1.2885   0.1548  14.1254
Weighted 2.5 % Perc.:     7.1763   -0.8534      2.8214   0.1654  16.4642
Weighted Median:         18.4333    1.3158     10.4421   0.2183  23.5013
Weighted Mean:           18.0423    1.5084     11.1828   0.2185  23.6488
Weighted Mode:           21.8762    0.1017      8.0346   0.1859  23.7741
Weighted 97.5 % Perc.:   28.7128    4.4678     20.1690   0.2731  32.2961
Max.:                    30.3855    5.2116     23.7994   0.2944  37.8058
                             t5 w_lowl_t5 admix_w_e_t6 admix_e_w_t6       t7
Min.:                    0.2320    7.1666     -22.6781     -22.2627   5.6300
Weighted 2.5 % Perc.:    0.2572   24.5255     -14.2133     -14.1896   6.5796
Weighted Median:         0.4762   43.1594      16.2478      39.6570   8.3345
Weighted Mean:           0.4602   42.5981      23.8872      42.0231   8.2991
Weighted Mode:           0.5347   48.1653       3.1615      40.1239   8.6265
Weighted 97.5 % Perc.:   0.6153   57.3682     116.1886      97.7170   9.9498
Max.:                    0.6402   61.5332     178.1282     145.2287  10.5785
                       w_anc_t7       t8  gor_anc
Min.:                   51.8909   7.5971   9.7929
Weighted 2.5 % Perc.:   62.8131   8.9328  21.0023
Weighted Median:        84.5198  13.3088  35.4170
Weighted Mean:          82.7070  13.1915  35.9544
Weighted Mode:          86.6025  13.3765  34.1359
Weighted 97.5 % Perc.:  98.2452  16.9102  48.2289
Max.:                  105.2163  18.6323  79.3707

#-----------------------------------------------------------------------------------------------------------------------
#scp -r hpawar@172.16.10.20:/scratch/devel/hpawar/admix/abc/simul/test/modelcomp/out/pls/abc_17plscomp_18feb22  /Users/harvi/Downloads/gorilla_abc/modelchoice/pls
#scp -r hpawar@172.16.10.20:/scratch/devel/hpawar/admix/abc/simul/test/modelcomp/out/pls/trget3_17plscomp_18feb22  /Users/harvi/Downloads/gorilla_abc/modelchoice/pls
# need to output the param file **
#scp -r hpawar@172.16.10.20:/scratch/devel/hpawar/admix/abc/simul/test/modelcomp/out/pls/param_17plscomp_18feb22 /Users/harvi/Downloads/gorilla_abc/modelchoice/pls



# on local
setwd("/Users/harvi/Downloads/gorilla_abc/modelchoice/pls")
require(data.table)

# A) tol 0.005
load("/Users/harvi/Downloads/gorilla_abc/modelchoice/pls/abc_17plscomp_18feb22",verbose=T)
load("/Users/harvi/Downloads/gorilla_abc/modelchoice/pls/param_17plscomp_18feb22",verbose=T)
library(abc)
# histogram of the weighted posterior sample  : https://cran.r-project.org/web/packages/abc/vignettes/abcvignette.pdf
pdf("/Users/harvi/Downloads/gorilla_abc/modelchoice/pls/abc_17plscom_posterior.18feb22.pdf") 
hist(myabc2_17)
dev.off() 
# generate the diagnostic plots of the estimation of the posterior distribution 
# a density plot of the prior distribution, a density plot of the posterior distribution estimated with and without regression- based correction, a scatter plot of the Euclidean distances as a function of the parameter values, and a normal Q-Q plot of the residuals from the regression
pdf("/Users/harvi/Downloads/gorilla_abc/modelchoice/pls/abc_17plscom_diagnostic.18feb22.pdf") 
plot(myabc2_17, param)
dev.off() 
