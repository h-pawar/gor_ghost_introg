# Mon 20 Feb 2023 09:18:56 CET

# ABC ghost parameter inference - sample all parameters from priors

# 1) process newly generated simulations
# add condition to retain only simulations where archaic introgression time < species split time → generate more simulations to get to ~35k (to make revised parameter inference complete)

# generated with scripts :  /scratch/devel/hpawar/admix/abc/simul/scripts/revisions_feb23/ghoste.rev.17feb23.arr
				#  /scratch/devel/hpawar/admix/abc/simul/scripts/revisions_feb23/ghoste.rev.17feb23.R


# 2) read in subset of simulations (from 15feb23) where archaic introgression time < species split time:
	# following /Volumes/"Ultra USB 3.0"/IBE/further.analysis.feb.2020/gorillas/abc/revisions_feb23/generateghoste.rev.abc.15feb23_subs.R


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
# load in subset of simulations generated 15feb23 **
	# objects: param_subs, sumstat1_subs
load(file=paste("/scratch/devel/hpawar/admix/abc/simul/test/ghost/revisions_8feb23/ghoste/segrecalc/ghosteabc_param_15feb23_subs"),verbose=T)
load(file=paste("/scratch/devel/hpawar/admix/abc/simul/test/ghost/revisions_8feb23/ghoste/segrecalc/ghosteabc_sumstat1_15feb23_subs"),verbose=T)

#-----------------------------------------------------------------------------------------------------------------------


# simulated data 
simufiles <- paste("/scratch/devel/hpawar/admix/abc/simul/test/ghost/revisions_8feb23/ghoste/rev_17feb23/", list.files(path = "/scratch/devel/hpawar/admix/abc/simul/test/ghost/revisions_8feb23/ghoste/rev_17feb23/", pattern="ghoste.rev.abc_sim"), sep = "")



simns=list()
simns_param=list()
for (i in 1:length(simufiles)){
load(file=simufiles[[i]],verbose=T)
simns[[i]]<-simuresults
simns_param[[i]]<-inputsets
}


# check if any simulations failed 
probl<-grep("Error", simns)

problsimns=list()
for (i in 1:length(probl)){
problsimns[[i]]<-grep("Error", simns[[probl[i]]])
}


#Thu 30 Sep 2021 11:59:08 CEST
#STREAMLINED VERSION OF FUNCTIONS

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
for (i in 1:length(t)){
probltest[[i]]<-cbind(
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

out_test=list()
for (i in 1:length(simufiles)) {
  
  skip_to_next <- FALSE
   
  tryCatch(

out_test[[i]]<-simu_fst_stats_fun(i)

, error = function(e) { skip_to_next <<- TRUE})
  
  if(skip_to_next) { next }     
}

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

# combine this with the other meaningful simulations from 15feb23 where tarchintrog < t8


sumstat<-data.matrix(sumstat)

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

# remove the id parameter
param<-param[,c(1:23)]

# remove same cols from the simulated vals
sumstat1<-sumstat[,-c(6,22,32,36)]

#-----------------------------------------------------------------------------------------------------------------------

# add in meaningful simn reps 

sumstat1<-data.matrix(sumstat1)  
sumstat1_subs<-data.matrix(sumstat1_subs)

csum<-rbind(sumstat1_subs, sumstat1)
csum<-data.matrix(csum)

#-----------------------------------------------------------------------------------------------------------------------
# remove ghost_t0 from the parameter search
param1<-param[,-c(5)]
#-----------------------------------------------------------------------------------------------------------------------

param1<-data.matrix(param1)  
param_subs<-data.matrix(param_subs)  


cparam<-rbind(param_subs, param1)


cparam<-data.matrix(cparam)

#-----------------------------------------------------------------------------------------------------------------------


lowerp<-c(3,0.1,0.1,0.1,0.1,0.01,1,0.01,0.1,0.14,0.1,0.2525,0,0.19,5,0,0,0.1,10,1.5,10,15.0025,10)
upperp<-c(100,20,30,20,100,20,25,5,20,0.25,30,14.9975,100,0.6,50,100,100,6,100,15,100,50,100)

# without ghost_t0 priors as follows
lowerp1<-lowerp[-5]
upperp1<-upperp[-5]


priordf<-cbind(lowerp1,upperp1)
prior_ranges<-data.matrix(priordf)

# Which you feed into the ABC like this:
ghosteabcrev_20feb23<-abc(target=target1,param=cparam,tol=0.005,sumstat=csum,method="neuralnet",numnet=100,transf="logit",logit.bounds=prior_ranges)

#-----------------------------------------------------------------------------------------------------------------------

ghosteabcrev_20feb23[[17]][[1]]<-colnames(cparam1)
ghosteabcrev_20feb23[[17]][[2]]<-names(target1)

ghosteabcrev_20feb23
summary(ghosteabcrev_20feb23)


save(ghosteabcrev_20feb23,file=paste("/scratch/devel/hpawar/admix/abc/simul/test/ghost/revisions_8feb23/ghoste/segrecalc/ghosteabcrev_20feb23"))
save(cparam,file=paste("/scratch/devel/hpawar/admix/abc/simul/test/ghost/revisions_8feb23/ghoste/segrecalc/ghosteabc_param_20feb23"))
save(csum,file=paste("/scratch/devel/hpawar/admix/abc/simul/test/ghost/revisions_8feb23/ghoste/segrecalc/ghosteabc_sumstat1_20feb23"))


#-----------------------------------------------------------------------------------------------------------------------

pdf("/scratch/devel/hpawar/admix/abc/simul/test/ghost/revisions_8feb23/ghoste/segrecalc/ghosteabcrev_posterior_20feb23.pdf") 
hist(ghosteabcrev_20feb23)

pdf("/scratch/devel/hpawar/admix/abc/simul/test/ghost/revisions_8feb23/ghoste/segrecalc/ghosteabcrev_diagnostic_20feb23.pdf") 
plot(ghosteabcrev_20feb23, cparam)
dev.off() 


