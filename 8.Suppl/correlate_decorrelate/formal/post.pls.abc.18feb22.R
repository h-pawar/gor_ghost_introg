# Thu 17 Feb 2022 18:19:30 CET
# Perform ABC inference using the decorrelated statistics (ie directly on the PLS components)

#-----------------------------------------------------------------------------------------------------------------------

#module load R/4.0.1

library("pls"); library("MASS");
library(abc)

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

#  read in the pls results
load("/scratch/devel/hpawar/admix/abc/simul/test/modelcomp/out/pls/plsr_woutvalid_10comp.17feb22", verbose=T)

cumsum(explvar(myPlsr)) # how much of variance is explained by 10 components used

# 2) get the scores for each simulation and each dimension, and predict the scores for the real statistics
scrs<-matrix(myPlsr$scores,ncol=ncol(myPlsr$scores))  
trget3<-predict(myPlsr,newdata=matrix(trget2,nrow=1),type="scores")

# then do abc on the first 10 dimensions
myabc2<-abc(target=trget3[1:10],param=param,tol=0.01,sumstat=scrs[,1:10],method="neuralnet",numnet=100)

myabc2[[17]][[1]]<-colnames(param)
myabc2[[17]][[2]]<-names(trget3)

myabc2
summary(myabc2)
save(myabc2,file=paste("/scratch/devel/hpawar/admix/abc/simul/test/modelcomp/out/pls/abc_pls_18feb22"))
save(trget3,file=paste("/scratch/devel/hpawar/admix/abc/simul/test/modelcomp/out/pls/trget3_18feb22"))

save(param,file=paste("/scratch/devel/hpawar/admix/abc/simul/test/modelcomp/out/pls/param_pls_18feb22"))


#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

# on local
setwd("/Users/harvi/Downloads/gorilla_abc/modelchoice/pls")
require(data.table)

# A) tol 0.005
load("/Users/harvi/Downloads/gorilla_abc/modelchoice/pls/abc_pls_18feb22",verbose=T)
load("/Users/harvi/Downloads/gorilla_abc/modelchoice/pls/param_pls_18feb22",verbose=T)
library(abc)
pdf("/Users/harvi/Downloads/gorilla_abc/modelchoice/pls/abc_pls_posterior.18feb22.pdf") 
hist(myabc2)
dev.off() 
pdf("/Users/harvi/Downloads/gorilla_abc/modelchoice/pls/abc_pls_diagnostic.18feb22.pdf") 
plot(myabc2, param)
dev.off() 

