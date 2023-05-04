## Thu 17 Feb 2022 11:07:48 CET
# Re-run find_pls.r with this optimum number of PLS components. (following https://github.com/hirzi/LSD)
# follows from Volumes/"Ultra USB 3.0"/IBE/further.analysis.feb.2020/gorillas/abc/simul_postabc/abc+ghost/modelcomp/uncorrelate.stats.26jan22.R)
#-----------------------------------------------------------------------------------------------------------------------

#module load R/4.0.1

library("pls"); library("MASS");

options(scipen=100)
require(data.table)
options(stringsAsFactors=F)

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
numComp=10

# perform pls # without validation
myPlsr<-plsr(as.matrix(param)~as.matrix(test), scale=F,ncomp=numComp)

save(myPlsr,file=paste("/scratch/devel/hpawar/admix/abc/simul/test/modelcomp/out/pls/plsr_woutvalid_10comp.17feb22"))

lambda1<-lambda[c(1:5,7:21,23:38)]
myGM1<-myGM[c(1:5,7:21,23:38)]


#write pls to a file
myPlsrDataFrame<-data.frame(comp1=myPlsr$loadings[,1]);
for(i in 2:numComp){myPlsrDataFrame<-cbind(myPlsrDataFrame, myPlsr$loadings[,i]);}
write.table(cbind(colnames(test), myMax1, myMin1, lambda1, myGM1, myBCMeans1, myBCSDs1, myPlsrDataFrame), file="/scratch/devel/hpawar/admix/abc/simul/test/modelcomp/out/pls/PLSfile_sstat.nullabc_woutvalid_10comp.17feb22.txt", col.names=F,row.names=F, sep='', quote=F); 

# write out in more readable manner
write.table(cbind(colnames(test), myMax1, myMin1, lambda1, myGM1, myBCMeans1, myBCSDs1, myPlsrDataFrame), file="/scratch/devel/hpawar/admix/abc/simul/test/modelcomp/out/pls/PLSfile_sstat.nullabc_woutvalid_10comp.17feb22.1.txt", col.names=F,row.names=F, sep='\t', quote=F); 

save(myPlsrDataFrame,file=paste("/scratch/devel/hpawar/admix/abc/simul/test/modelcomp/out/pls/plsrdf_woutvalid_10comp.17feb22"))

# pdf plots = Number of PLS components versus root mean square error of prediction (RMSEP) for the normal distribution model.
 pdf("/scratch/devel/hpawar/admix/abc/simul/test/modelcomp/out/pls/RMSEP_sstat.nullabc_woutvalid_10comp.17feb22.pdf");
plot(RMSEP(myPlsr));
dev.off();

q()

#-----------------------------------------------------------------------------------------------------------------------

