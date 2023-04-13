# Mon 20 Feb 2023 09:18:56 CET

# ABC ghost parameter inference - sample all parameters from priors

# 1) process newly generated simulations
# add condition to retain only simulations where archaic introgression time < species split time → generate more simulations to get to ~35k (to make revised parameter inference complete)

# generated with scripts :  /scratch/devel/hpawar/admix/abc/simul/scripts/revisions_feb23/ghoste.rev.17feb23.arr
				#  /scratch/devel/hpawar/admix/abc/simul/scripts/revisions_feb23/ghoste.rev.17feb23.R


# paths to simulated data: /scratch/devel/hpawar/admix/abc/simul/test/ghost/revisions_8feb23/ghoste/rev_17feb23/


# 2) read in subset of simulations (from 15feb23) where archaic introgression time < species split time:
	# following /Volumes/"Ultra USB 3.0"/IBE/further.analysis.feb.2020/gorillas/abc/revisions_feb23/generateghoste.rev.abc.15feb23_subs.R


# paths to simulated data: /scratch/devel/hpawar/admix/abc/simul/test/ghost/revisions_8feb23/ghoste


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
# load in subset of simulations generated 15feb23 **
	# objects: param_subs, sumstat1_subs
load(file=paste("/scratch/devel/hpawar/admix/abc/simul/test/ghost/revisions_8feb23/ghoste/segrecalc/ghosteabc_param_15feb23_subs"),verbose=T)
load(file=paste("/scratch/devel/hpawar/admix/abc/simul/test/ghost/revisions_8feb23/ghoste/segrecalc/ghosteabc_sumstat1_15feb23_subs"),verbose=T)

# str(sumstat1_subs)
#'data.frame':	21623 obs. of  40 variables:

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
# simulations which failed 
probl<-grep("Error", simns)

problsimns=list()
for (i in 1:length(probl)){
problsimns[[i]]<-grep("Error", simns[[probl[i]]])
}


# probl
# [1]   6  74  75  76  82  84  89 103 107 109 122 154 169 178 181 192 194 203 220
#[20] 235 239 282

# more failed reps -> will need to amend the functions below **
	# am making a mess of this...


#howmanyproblsimns=list()
#for (i in 1:length(probl)){
#howmanyproblsimns[[i]]<-length(problsimns[[i]])
#}

# sum(unlist(howmanyproblsimns))
#[1] 327

#(290*51)- sum(unlist(howmanyproblsimns))
#[1] 14463

# will need to check how dealt with this in the null scenario *
# refer back to generateabc.amendtajima.R


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

out_test=list()
for (i in 1:length(simufiles)) {
  
  skip_to_next <- FALSE
  
  # Note that print(b) fails since b doesn't exist
  
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
#[1] 18564

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

# remove same cols from the simulated vals
sumstat1<-sumstat[,-c(6,22,32,36)]

# ncol(sumstat1) #[1] 40


#-----------------------------------------------------------------------------------------------------------------------

# add in meaningful simn reps 

#nrow(sumstat1)+nrow(sumstat1_subs)
#[1] 36086

# target is 35,700 reps -> ie combine then subset to this number

sumstat1<-data.matrix(sumstat1)  
sumstat1_subs<-data.matrix(sumstat1_subs)


# explicitly check how mnay 
#nrow(sumstat1[which(param1$t_archintrog < param1$t8),])

csum<-rbind(sumstat1_subs, sumstat1)

#str(csum)
# num [1:32505, 1:40] 1.552 1.336 1.501 0.899 1.14 ...
# - attr(*, "dimnames")=List of 2
#  ..$ : chr [1:32505] "2" "4" "7" "10" ...
#  ..$ : chr [1:40] "X1" "X2" "X3" "X4" ...


#csum1<-csum[c(1:35700),]

#csum1<-data.matrix(csum1)
# ie this is the object of summary stats to query -> combine also the param objects *

csum<-data.matrix(csum)

#-----------------------------------------------------------------------------------------------------------------------


# remove ghost_t0 from the parameter search
param1<-param[,-c(5)]
#head(param1)

#-----------------------------------------------------------------------------------------------------------------------



# nrow(param1)
#[1] 18564

param1<-data.matrix(param1)  
param_subs<-data.matrix(param_subs)  


cparam<-rbind(param_subs, param1)


cparam<-data.matrix(cparam)

# str(cparam)
# num [1:35291, 1:22] 92.3 99.9 13.8 85.1 37.6 ...
# - attr(*, "dimnames")=List of 2
#  ..$ : chr [1:35291] "3" "4" "7" "9" ...
#  ..$ : chr [1:22] "w_lowl_t0" "w_cros_t0" "e_lowl_t0" "e_moun_t0" ...

#-----------------------------------------------------------------------------------------------------------------------


lowerp<-c(3,0.1,0.1,0.1,0.1,0.01,1,0.01,0.1,0.14,0.1,0.2525,0,0.19,5,0,0,0.1,10,1.5,10,15.0025,10)
upperp<-c(100,20,30,20,100,20,25,5,20,0.25,30,14.9975,100,0.6,50,100,100,6,100,15,100,50,100)

# remove ghostt0 prior from this * - remove entry 5 
# without ghost_t0 priors as follows
lowerp1<-lowerp[-5]
upperp1<-upperp[-5]


priordf<-cbind(lowerp1,upperp1)
prior_ranges<-data.matrix(priordf)

# Which you feed into the ABC like this:
ghosteabcrev_20feb23<-abc(target=target1,param=cparam,tol=0.005,sumstat=csum,method="neuralnet",numnet=100,transf="logit",logit.bounds=prior_ranges)

# The "tol" parameter is unknown and might be tested a bit.



#-----------------------------------------------------------------------------------------------------------------------

ghosteabcrev_20feb23[[17]][[1]]<-colnames(cparam1)
ghosteabcrev_20feb23[[17]][[2]]<-names(target1)

ghosteabcrev_20feb23
summary(ghosteabcrev_20feb23)


# mkdir -p /scratch/devel/hpawar/admix/abc/simul/test/ghost/revisions_8feb23/ghoste/segrecalc


save(ghosteabcrev_20feb23,file=paste("/scratch/devel/hpawar/admix/abc/simul/test/ghost/revisions_8feb23/ghoste/segrecalc/ghosteabcrev_20feb23"))
save(cparam,file=paste("/scratch/devel/hpawar/admix/abc/simul/test/ghost/revisions_8feb23/ghoste/segrecalc/ghosteabc_param_20feb23"))
save(csum,file=paste("/scratch/devel/hpawar/admix/abc/simul/test/ghost/revisions_8feb23/ghoste/segrecalc/ghosteabc_sumstat1_20feb23"))


#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

> ghosteabcrev_20feb23
Call:
abc(target = target1, param = cparam, sumstat = csum, tol = 0.005, 
    method = "neuralnet", transf = "logit", logit.bounds = prior_ranges, 
    numnet = 100)
Method:
Non-linear regression via neural networks
with correction for heteroscedasticity

Parameters:
w_lowl_t0, w_cros_t0, e_lowl_t0, e_moun_t0, e_lowl_t1, e_lowl_t2, e_moun_t3, e_moun_t3.1, t4, e_anc_t4, t_archintrog, archaicintrog, t5, w_lowl_t5, admix_w_e_t6, admix_e_w_t6, t7, w_anc_t7, t8, gor_anc, t9, gor_ghost_anc

Statistics:
het_mu_WL, het_mu_WC, het_mu_EL, het_mu_EM, het_sd_WL, het_sd_EL, het_sd_EM, fixedsites_WL, fixedsites_WC, fixedsites_EL, fixedsites_EM, segsites_WL, segsites_WC, segsites_EL, segsites_EM, fixedperid_mu_WL, fixedperid_mu_WC, fixedperid_mu_EL, fixedperid_mu_EM, fixedperid_sd_WL, fixedperid_sd_EL, fixedperid_sd_EM, fst_WL.WC, fst_WL.EL, fst_WL.EM, fst_WC.EL, fst_WC.EM, fst_EL.EM, pi_mu_WL, pi_mu_EL, pi_mu_EM, pi_sd_WL, pi_sd_EL, pi_sd_EM, tajima_mu_WL, tajima_mu_EL, tajima_mu_EM, tajima_sd_WL, tajima_sd_EL, tajima_sd_EM

Total number of simulations 35291 

Number of accepted simulations:  177 

> summary(ghosteabcrev_20feb23)
Call: 
abc(target = target1, param = cparam, sumstat = csum, tol = 0.005, 
    method = "neuralnet", transf = "logit", logit.bounds = prior_ranges, 
    numnet = 100)
Data:
 abc.out$adj.values (177 posterior samples)
Weights:
 abc.out$weights

                       w_lowl_t0 w_cros_t0 e_lowl_t0 e_moun_t0 e_lowl_t1
Min.:                    32.9331    1.2511    0.6345    0.1004    1.1174
Weighted 2.5 % Perc.:    58.1893    6.9713    2.4681    0.1069    1.1906
Weighted Median:         90.2319   15.8214   12.1790    0.1601    3.1168
Weighted Mean:           87.2338   14.8643   12.8065    0.2641    3.1870
Weighted Mode:           93.0488   17.7933    6.8437    0.1451    3.1381
Weighted 97.5 % Perc.:   99.4077   19.7522   25.7445    1.4019    5.8492
Max.:                    99.8949   19.9596   27.5437    4.8083   18.9757
                       e_lowl_t2 e_moun_t3 e_moun_t3.1      t4 e_anc_t4
Min.:                     1.0808    0.0166      3.0245  0.1462   4.2525
Weighted 2.5 % Perc.:     1.2541    0.0213      4.3852  0.2159   7.7414
Weighted Median:          4.7071    0.8084     13.8682  0.2457  16.5929
Weighted Mean:            6.1806    1.0546     13.3829  0.2419  17.2678
Weighted Mode:            2.6042    0.4323     17.5620  0.2477  15.0057
Weighted 97.5 % Perc.:   21.7987    4.0057     19.3047  0.2499  27.0713
Max.:                    24.3148    4.5542     19.7309  0.2499  29.9839
                       t_archintrog archaicintrog      t5 w_lowl_t5
Min.:                        0.7645        1.0595  0.1919    7.0381
Weighted 2.5 % Perc.:        1.9242       13.3059  0.2097   25.1404
Weighted Median:            10.7847       77.3443  0.4044   43.1105
Weighted Mean:               9.7616       70.2802  0.3976   41.1250
Weighted Mode:              12.6058       87.0123  0.4497   44.5144
Weighted 97.5 % Perc.:      14.0861       97.9366  0.5876   49.0612
Max.:                       14.6253       99.9060  0.5944   49.8853
                       admix_w_e_t6 admix_e_w_t6      t7 w_anc_t7      t8
Min.:                        0.0737       0.5827  1.2932  92.3449  5.5779
Weighted 2.5 % Perc.:        0.3273       5.8783  2.2953  94.0782  8.9344
Weighted Median:             3.1324      45.3267  4.2411  98.3656 13.1881
Weighted Mean:               6.2213      45.7968  4.2295  97.9403 12.8826
Weighted Mode:               1.6226      23.4525  4.7387  99.2182 14.0101
Weighted 97.5 % Perc.:      44.1508      94.9658  5.8868  99.8745 14.8760
Max.:                       56.7357      99.4308  5.9299  99.9888 14.9998
                       gor_anc      t9 gor_ghost_anc
Min.:                  10.0499 15.1547       14.5909
Weighted 2.5 % Perc.:  10.7771 18.8038       59.2967
Weighted Median:       17.3283 41.0872       93.0008
Weighted Mean:         19.8359 37.5536       89.7248
Weighted Mode:         14.8940 44.2662       95.3068
Weighted 97.5 % Perc.: 45.4792 49.6404       99.0966
Max.:                  96.7386 49.9038       99.4782




#-----------------------------------------------------------------------------------------------------------------------

# histogram of the weighted posterior sample  : https://cran.r-project.org/web/packages/abc/vignettes/abcvignette.pdf
pdf("/scratch/devel/hpawar/admix/abc/simul/test/ghost/revisions_8feb23/ghoste/segrecalc/ghosteabcrev_posterior_20feb23.pdf") 

hist(ghosteabcrev_20feb23)


# generate the diagnostic plots of the estimation of the posterior distribution 
# a density plot of the prior distribution, a density plot of the posterior distribution estimated with and without regression- based correction, a scatter plot of the Euclidean distances as a function of the parameter values, and a normal Q-Q plot of the residuals from the regression
pdf("/scratch/devel/hpawar/admix/abc/simul/test/ghost/revisions_8feb23/ghoste/segrecalc/ghosteabcrev_diagnostic_20feb23.pdf") 
plot(ghosteabcrev_20feb23, cparam)
dev.off() 

#~/Downloads/gor_ghost_revisions_feb23

q()

scp -r -oHostKeyAlgorithms=+ssh-rsa  hpawar@172.16.10.21:/scratch/devel/hpawar/admix/abc/simul/test/ghost/revisions_8feb23/ghoste/segrecalc/ghosteabcrev_diagnostic_20feb23.pdf ~/Downloads/gor_ghost_revisions_feb23
scp -r -oHostKeyAlgorithms=+ssh-rsa  hpawar@172.16.10.21:/scratch/devel/hpawar/admix/abc/simul/test/ghost/revisions_8feb23/ghoste/segrecalc/ghosteabcrev_posterior_20feb23.pdf ~/Downloads/gor_ghost_revisions_feb23



#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------




