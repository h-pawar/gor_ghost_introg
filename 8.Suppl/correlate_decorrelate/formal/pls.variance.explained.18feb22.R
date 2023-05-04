#Fri 18 Feb 2022 14:15:58 CET

# 1) Calculate how many PLS components explain 99% of the variance
# 2) Perform ABC inference using these 17 components (rather than the optimal 10 components - as performed in post.pls.abc.18feb22.R)

#-----------------------------------------------------------------------------------------------------------------------
#module load R/4.0.1

library("pls"); library("MASS");

options(scipen=100)
require(data.table)
options(stringsAsFactors=F)
load("/scratch/devel/hpawar/admix/abc/simul/test/modelcomp/out/pls/plsr_woutvalid_10comp.17feb22", verbose=T)
cumsum(explvar(myPlsr))
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

numComp=12

# perform pls # without validation
myPlsr<-plsr(as.matrix(param)~as.matrix(test), scale=F,ncomp=numComp)
 cumsum(explvar(myPlsr))
 
numComp=15

# perform pls # without validation
myPlsr<-plsr(as.matrix(param)~as.matrix(test), scale=F,ncomp=numComp)
 cumsum(explvar(myPlsr))


numComp=20
# perform pls # without validation
myPlsr<-plsr(as.matrix(param)~as.matrix(test), scale=F,ncomp=numComp)
 cumsum(explvar(myPlsr))

#-----------------------------------------------------------------------------------------------------------------------

#to get to 99% of variance explained I have to increase the PLS components to 17.
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
# then do abc on the first 17 dimensions
myabc2_17<-abc(target=trget3[1:10],param=param,tol=0.01,sumstat=scrs[,1:10],method="neuralnet",numnet=100)


myabc2_17
summary(myabc2_17)
save(myabc2_17,file=paste("/scratch/devel/hpawar/admix/abc/simul/test/modelcomp/out/pls/abc_17plscomp_18feb22"))
save(trget3,file=paste("/scratch/devel/hpawar/admix/abc/simul/test/modelcomp/out/pls/trget3_17plscomp_18feb22"))

save(param,file=paste("/scratch/devel/hpawar/admix/abc/simul/test/modelcomp/out/pls/param_17plscomp_18feb22"))

#-----------------------------------------------------------------------------------------------------------------------

# on local
setwd("/Users/harvi/Downloads/gorilla_abc/modelchoice/pls")
require(data.table)

# A) tol 0.005
load("/Users/harvi/Downloads/gorilla_abc/modelchoice/pls/abc_17plscomp_18feb22",verbose=T)
load("/Users/harvi/Downloads/gorilla_abc/modelchoice/pls/param_17plscomp_18feb22",verbose=T)
library(abc)
pdf("/Users/harvi/Downloads/gorilla_abc/modelchoice/pls/abc_17plscom_posterior.18feb22.pdf") 
hist(myabc2_17)
dev.off() 
pdf("/Users/harvi/Downloads/gorilla_abc/modelchoice/pls/abc_17plscom_diagnostic.18feb22.pdf") 
plot(myabc2_17, param)
dev.off() 
