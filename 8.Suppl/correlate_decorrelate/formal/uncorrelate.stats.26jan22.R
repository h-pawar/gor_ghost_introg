#module load R/4.0.1

library("pls"); library("MASS");


# ~35,700 ms simulations for null demog

options(scipen=100)
require(data.table)
options(stringsAsFactors=F)

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

# simulated data 
simufiles <- paste("/scratch/devel/hpawar/admix/abc/simul/test/11sep21/", list.files(path = "/scratch/devel/hpawar/admix/abc/simul/test/11sep21/", pattern="abc_sim"), sep = "")

simns=list()
simns_param=list()
for (i in 1:length(simufiles)){
load(file=simufiles[[i]],verbose=T)
simns[[i]]<-simuresults
simns_param[[i]]<-inputsets
}

# simulations which failed 
probl<-grep("Error", simns)

problsimns=list()
for (i in 1:length(probl)){
problsimns[[i]]<-grep("Error", simns[[probl[i]]])
}


#-----------------------------------------------------------------------------------------------------------------------
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

out_1_10<-afun(1,10)
#11
out_12_45<-afun(12,45)
#46
out_47_102<-afun(47,102)
#103
out_104_241<-afun(104,241)
#242
out_243_246<-afun(243,246)
#247
out_248_293<-afun(248,293)
#294,295
out_296_426<-afun(296,426)
#427
out_428_465<-afun(428,465)
#466
out_467_475<-afun(467,475)
#476
out_477_482<-afun(477,482)
#483
out_484_509<-afun(484,509)
#510
out_511_563<-afun(511,563)
#564
out_565_672<-afun(565,672)
#673
out_674_700<-afun(674,700)



# apply for all 14 simn reps which failed
out_rest=list()
for (i in 1:length(probl)){
out_rest[[i]]<-afun1(i)
}


out_12_45<-out_12_45[lengths(out_12_45) != 0]
out_47_102<-out_47_102[lengths(out_47_102) != 0]
out_104_241<-out_104_241[lengths(out_104_241) != 0]
out_243_246<-out_243_246[lengths(out_243_246) != 0]
out_248_293<-out_248_293[lengths(out_248_293) != 0]
out_296_426<-out_296_426[lengths(out_296_426) != 0]
out_428_465<-out_428_465[lengths(out_428_465) != 0]
out_467_475<-out_467_475[lengths(out_467_475) != 0]
out_477_482<-out_477_482[lengths(out_477_482) != 0]
out_484_509<-out_484_509[lengths(out_484_509) != 0]
out_511_563<-out_511_563[lengths(out_511_563) != 0]
out_565_672<-out_565_672[lengths(out_565_672) != 0]
out_674_700<-out_674_700[lengths(out_674_700) != 0]


test2<-do.call(c, list(out_1_10,list(out_rest[[1]]),
out_12_45,
list(out_rest[[2]]),
out_47_102,
list(out_rest[[3]]),
out_104_241,
list(out_rest[[4]]),
out_243_246,
list(out_rest[[5]]),
out_248_293,
list(out_rest[[6]]),
list(out_rest[[7]]),
out_296_426,
list(out_rest[[8]]),
out_428_465,
list(out_rest[[9]]),
out_467_475,
list(out_rest[[10]]),
out_477_482,
list(out_rest[[11]]),
out_484_509,
list(out_rest[[12]]),
out_511_563,
list(out_rest[[13]]),
out_565_672,
list(out_rest[[14]]),
out_674_700))

out_1_700<-as.data.frame(do.call(rbind,test2)) # out df for the 700 simns

return(out_1_700)
}
#-----------------------------------------------------------------------------------------------------------------------
comb_process_fun<-function(afun,afun1,a,b) {

out_1_10<-afun(1,10,a,b)
#11
out_12_45<-afun(12,45,a,b)
#46
out_47_102<-afun(47,102,a,b)
#103
out_104_241<-afun(104,241,a,b)
#242
out_243_246<-afun(243,246,a,b)
#247
out_248_293<-afun(248,293,a,b)
#294,295
out_296_426<-afun(296,426,a,b)
#427
out_428_465<-afun(428,465,a,b)
#466
out_467_475<-afun(467,475,a,b)
#476
out_477_482<-afun(477,482,a,b)
#483
out_484_509<-afun(484,509,a,b)
#510
out_511_563<-afun(511,563,a,b)
#564
out_565_672<-afun(565,672,a,b)
#673
out_674_700<-afun(674,700,a,b)



# apply for all 14 simn reps which failed
out_rest=list()
for (i in 1:length(probl)){
out_rest[[i]]<-afun1(i,a,b)
}


out_12_45<-out_12_45[lengths(out_12_45) != 0]
out_47_102<-out_47_102[lengths(out_47_102) != 0]
out_104_241<-out_104_241[lengths(out_104_241) != 0]
out_243_246<-out_243_246[lengths(out_243_246) != 0]
out_248_293<-out_248_293[lengths(out_248_293) != 0]
out_296_426<-out_296_426[lengths(out_296_426) != 0]
out_428_465<-out_428_465[lengths(out_428_465) != 0]
out_467_475<-out_467_475[lengths(out_467_475) != 0]
out_477_482<-out_477_482[lengths(out_477_482) != 0]
out_484_509<-out_484_509[lengths(out_484_509) != 0]
out_511_563<-out_511_563[lengths(out_511_563) != 0]
out_565_672<-out_565_672[lengths(out_565_672) != 0]
out_674_700<-out_674_700[lengths(out_674_700) != 0]


test2<-do.call(c, list(out_1_10,list(out_rest[[1]]),
out_12_45,
list(out_rest[[2]]),
out_47_102,
list(out_rest[[3]]),
out_104_241,
list(out_rest[[4]]),
out_243_246,
list(out_rest[[5]]),
out_248_293,
list(out_rest[[6]]),
list(out_rest[[7]]),
out_296_426,
list(out_rest[[8]]),
out_428_465,
list(out_rest[[9]]),
out_467_475,
list(out_rest[[10]]),
out_477_482,
list(out_rest[[11]]),
out_484_509,
list(out_rest[[12]]),
out_511_563,
list(out_rest[[13]]),
out_565_672,
list(out_rest[[14]]),
out_674_700))

out_1_700<-as.data.frame(do.call(rbind,test2)) # out df for the 700 simns

return(out_1_700)
}

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
#-----------------------------------------------------------------------------------------------------------------------

#  popn-wise fixed sites
s_popfix<-comb_process_fun(out_het_fun1,het_format_problsimns_fun,2,1)
s_popfix_perkb<-((s_popfix/(2500*40000))*1000)

#  popn-wise seg sites
s_popseg<-comb_process_fun(out_het_fun1,het_format_problsimns_fun,2,2)
s_popseg_perkb<-((s_popseg/(2500*40000))*1000)

#-----------------------------------------------------------------------------------------------------------------------

sumstat<-cbind(s_het_mu,s_het_sd,
s_popfix_perkb,s_popseg_perkb,
s_ifix_mu,s_ifix_sd,
s_fst,
s_pi_mu,
s_pi_sd,
s_tajima_mu,
s_tajima_sd)

sumstat<-data.matrix(sumstat)

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

# format input parameters used to generate the simulations generating these summary stats

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

for (i in 1:length(probl)){
pn_out[[probl[[i]]]]<-pn_out[[probl[[i]]]][-c(problsimns[[i]])]
}

pn_format=list()
for (i in 1:length(pn_out)){
pn_format[[i]]<-data.frame(matrix(unlist(pn_out[[i]]), nrow=length(pn_out[[i]]), byrow=TRUE),stringsAsFactors=FALSE)
}

param_df<-as.data.frame(do.call(rbind, pn_format))

nrow(param_df) 
ncol(param_df) 

param<-data.matrix(param_df)

paraname=c("w_lowl_t0","w_cros_t0","e_lowl_t0","e_moun_t0","e_lowl_t1","e_lowl_t2","e_moun_t3","e_moun_t3.1","t4","e_anc_t4","t5","w_lowl_t5","admix_w_e_t6","admix_e_w_t6","t7","w_anc_t7","t8","gor_anc","id")
  
colnames(param)<-paraname
head(param)

param<-param[,c(1:18)]

#-----------------------------------------------------------------------------------------------------------------------
# remove same cols from the simulated vals
sumstat1<-sumstat[,-c(6,22,32,36)]

#-----------------------------------------------------------------------------------------------------------------------


##force stat in [1,2]
myMax<-c(); myMin<-c(); lambda<-c(); myGM<-c();
for(i in 1:ncol(sumstat1)){
  myMax<-c(myMax, max(sumstat1[,i])); myMin<-c(myMin, min(sumstat1[,i]));
  sumstat1[,i]<-1+(sumstat1[,i] -myMin[i])/(myMax[i]-myMin[ i]);
  print(c(range(sumstat1[,i]),length(which(sumstat1[,i]==NA)),mean(sumstat1[,i])))
}

#transform statistics via boxcox
for(i in 1:ncol(sumstat1)){
print(i)
d<-cbind(data.frame(sumstat1[,i]), data.frame(param));
mylm<-lm(as.formula(d), data=d);
myboxcox<-boxcox(mylm, lambda=seq(-20,100,1/10), interp=T, eps=1/50);
lambda<-c(lambda, myboxcox$x[myboxcox$y==max(myboxcox$y)]);
myGM<-c(myGM, mean(exp(log(sumstat[,i]))));
}


sstat<-sumstat1
#standardize the BC-stat
myBCMeans<-c(); myBCSDs<-c();
for(i in 1:ncol(sstat)){
    sstat[,i]<-(sstat[,i] ^lambda[i] - 1)/(lambda[i]*myGM[i] ^(lambda[i]-1));
    myBCSDs<-c(myBCSDs, sd(sstat[,i]));
    myBCMeans<-c(myBCMeans, mean(sstat[,i]));
    sstat[,i]<-(sstat[,i] -myBCMeans[i])/myBCSDs[i];
  }

# remove the cols with NaN (so the pls function runs)
test<-sstat[,c(1:5,7:21,23:38)]

# perform pls
myPlsr<-plsr(as.matrix(param)~as.matrix(test), scale=F, validation='LOO');

#write pls to a file
myPlsrDataFrame<-data.frame(comp1=myPlsr$loadings[,1]);
for(i in 2:numComp){myPlsrDataFrame<-cbind(myPlsrDataFrame, myPlsr$loadings[,i]);}
write.table(cbind(names(stats), myMax, myMin, lambda, myGM, myBCMeans, myBCSDs, myPlsrDataFrame), file="/scratch/devel/hpawar/admix/abc/simul/test/modelcomp/out/PLSfile_sstat.nullabc.15feb22.txt", col.names=F,row.names=F, sep=’’, quote=F); 

#make RMSEP plot
pdf("/scratch/devel/hpawar/admix/abc/simul/test/modelcomp/out/RMSEP_sstat.nullabc.15feb22.pdf");
plot(RMSEP(myPlsr));
dev.off();

q()
#-----------------------------------------------------------------------------------------------------------------------
