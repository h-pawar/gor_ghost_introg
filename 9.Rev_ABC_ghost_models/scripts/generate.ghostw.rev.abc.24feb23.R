# Fri 24 Feb 2023 10:13:26 CET
# ABC ghost parameter inference - sample all parameters from priors

# /scratch/devel/hpawar/admix/abc/simul/scripts/revisions_feb23/ghostw.rev.17feb23.arr
# /scratch/devel/hpawar/admix/abc/simul/scripts/revisions_feb23/ghostw.rev.17feb23.R

#module load R/4.0.1
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
library(abc)
options(scipen=100)
require(data.table)
options(stringsAsFactors=F)


# load in the empirical data

# directly read in target1
load("/scratch/devel/hpawar/admix/abc/results/test/11sep21simns/segrecalc/target1_6oct", verbose=T)

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


#howmanyproblsimns=list()
#for (i in 1:length(probl)){
#howmanyproblsimns[[i]]<-length(problsimns[[i]])
#}

#-----------------------------------------------------------------------------------------------------------------------



#Thu 30 Sep 2021 11:59:08 CEST
#STREAMLINED VERSION OF FUNCTIONS
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

all_fun<-function(afun){

out_test=list()
for (i in 1:length(simufiles)) {
  
  skip_to_next <- FALSE
  
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

#-----------------------------------------------------------------------------------------------------------------------

call_fun<-function(afun,a,b){

out_test=list()
for (i in 1:length(simufiles)) {
  
  skip_to_next <- FALSE
  
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


sumstat<-cbind(s_het_mu,s_het_sd,
s_popfix_perkb,s_popseg_perkb,
s_ifix_mu,s_ifix_sd,
s_fst,
s_pi_mu,
s_pi_sd,
s_taj)


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


# remove iterations with any failed rep - probl
pn_out1<-pn_out[-probl]

pn_format=list()
for (i in 1:length(pn_out1)){
pn_format[[i]]<-data.frame(matrix(unlist(pn_out1[[i]]), nrow=length(pn_out1[[i]]), byrow=TRUE),stringsAsFactors=FALSE)
}



param_df<-as.data.frame(do.call(rbind, pn_format))


param<-data.matrix(param_df)

# all parameters now being sampled from priors **
paraname=c("w_lowl_t0","w_cros_t0","e_lowl_t0","e_moun_t0","ghost_t0","e_lowl_t1","e_lowl_t2","e_moun_t3","e_moun_t3.1","t4","e_anc_t4","t_archintrog","archaicintrog","t5","w_lowl_t5","admix_w_e_t6","admix_e_w_t6","t7","w_anc_t7","t8","gor_anc","t9","gor_ghost_anc","id")

colnames(param)<-paraname
head(param)


param<-param[,c(1:23)]



#-----------------------------------------------------------------------------------------------------------------------

# remove same cols from the simulated vals
sumstat1<-sumstat[,-c(6,22,32,36)]


#-----------------------------------------------------------------------------------------------------------------------

# remove ghost_t0 from the parameter search
param1<-param[,-c(5)]


#-----------------------------------------------------------------------------------------------------------------------


param1<-data.matrix(param1)  
]
#-----------------------------------------------------------------------------------------------------------------------

# should also force these posteriors to be within the priors, same logit transf

# W priors
lowerp<-c(3,0.1,0.1,0.1,0.1,0.01,1,0.01,0.1,0.14,0.1,6.0025,0,0.19,5,0,0,0.1,10,1.5,10,15.0025,10)
upperp<-c(100,20,30,20,100,20,25,5,20,0.25,30,14.9975,100,0.6,50,100,100,6,100,15,100,50,100)

# without ghost_t0 priors as follows
lowerp1<-lowerp[-5]
upperp1<-upperp[-5]


priordf<-cbind(lowerp1,upperp1)
prior_ranges<-data.matrix(priordf)

# Which you feed into the ABC like this:
ghostwabcrev_24feb23<-abc(target=target1,param=param1,tol=0.005,sumstat=sumstat1,method="neuralnet",numnet=100,transf="logit",logit.bounds=prior_ranges)

#-----------------------------------------------------------------------------------------------------------------------

ghostwabcrev_24feb23[[17]][[1]]<-colnames(param1)
ghostwabcrev_24feb23[[17]][[2]]<-names(target1)

ghostwabcrev_24feb23
summary(ghostwabcrev_24feb23)

save(ghostwabcrev_24feb23,file=paste("/scratch/devel/hpawar/admix/abc/simul/test/ghost/revisions_8feb23/ghostw/segrecalc/ghostwabcrev_24feb23"))
save(param1,file=paste("/scratch/devel/hpawar/admix/abc/simul/test/ghost/revisions_8feb23/ghostw/segrecalc/ghostwabc_param_24feb23"))
save(sumstat1,file=paste("/scratch/devel/hpawar/admix/abc/simul/test/ghost/revisions_8feb23/ghostw/segrecalc/ghostwabc_sumstat1_24feb23"))


#-----------------------------------------------------------------------------------------------------------------------

pdf("/scratch/devel/hpawar/admix/abc/simul/test/ghost/revisions_8feb23/ghostw/segrecalc/ghostwabcrev_posterior_24feb23.pdf") 
hist(ghostwabcrev_24feb23)
dev.off()
pdf("/scratch/devel/hpawar/admix/abc/simul/test/ghost/revisions_8feb23/ghostw/segrecalc/ghostwabcrev_diagnostic_24feb23.pdf") 
plot(ghostwabcrev_24feb23, param1)
dev.off() 

