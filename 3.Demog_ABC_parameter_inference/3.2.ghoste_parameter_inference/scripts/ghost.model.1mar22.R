#!/usr/bin/r

# Tue  5 Apr 2022 15:34:56 CEST
# amended with posteriors from final null abc parameter inference on 5/4/22

# Mon 28 Feb 2022 15:13:58 CET
# reinfer ABC parameter inference of ghoste model
# having reperformed ABC Parameter Inference for the null window-based model - introducing logit transform in the ABC - to force posteriors to be within the priors

# Tue  2 Nov 2021 09:37:14 CET
# generate demog model + ghost introgression into E


#for model comparison I should compare A) the null model B) archaic to E C) archaic to W.
# MK - Take the same priors for B) and C) to start from things being equal.

# following test.abc.model.v5.R (final null demog model)
#-----------------------------------------------------------------------------------------------------------------------
# new batch of simulated data (700 simn reps) & summary stats were generated with 
    # /scratch/devel/hpawar/admix/abc/simul/scripts/test.abc.model.arr - which calls
    # /scratch/devel/hpawar/admix/abc/simul/scripts/test.abc.model.v5.R # new version supercedes test.abc.model.v4.R (used to generate prev 1-2000 simns)
#-----------------------------------------------------------------------------------------------------------------------

# Fri 10 Sep 2021 15:28:07 CEST
# Need to regenerate simulations (initially 1/3 of the 2000 reps - to see if this new approach is working)
# 1) change way of calculating segregating sites 
  # to match empirical data need to remove sites which are fixed across all gorilla populations (ie 1/1 across all pops)
  # b/c empirical data being mapped to human ref introduces a large number of such sites
# 2) do not simulate where t8<t7 (& t5 or t6 > t7) - ie need (t5/t6<t7, t7<t8)
# 3) fix the time of admixture b/n the extant lineages (atm high levels of confusion b/c timing of admixture is overlapping the divergence times b/n the sp)


# need to incorporate MK suggestions *** done - send for 10 test reps & if fine then scale up - Sat 11 Sep 2021 11:29:56 CEST

#------------------------------------------------------------------------------------------------------------------------

# adapting /scratch/devel/hpawar/admix/abc/simul/scripts/test.abc.model.v4.R

# Wed 19 May 2021 10:29:24 CEST
# adapt key steps from MK script simul_abc_model.R

# adapt gorilla ms steps from simul.testdemog.R
  # generate new simulations w/out forcing a given number of segregating sites - # segregating sites depend on the demographic model
  # except here won't be applying the S* statistic & additional demog events proposed - detailed in gor.initial.priors.18mar.sh

# gorilla wide uniform priors for the parameter vals from gor.initial.priors.18mar.sh
# calc summary statistics for this ms simn - code from het.simld.15apr.R & allfunction.18may.R (heterozygosity, number of segregating sites, pairwise Fst, pi, Tajima's D)
#module load gcc/6.3.0 R/3.4.2
# for troubleshooting - Mon 31 May 2021 16:09:42 CEST
#mkdir /dev/shm/mydata
#cd /dev/shm/mydata
#------------------------------------------------------------------------------------------------------------------------
require(data.table)
options(stringsAsFactors=F)
library(mgcv)
library(gap)
library(parallel)
library('pegas')
## see below what ptype stands for
  # ptype = array number
ptype=(commandArgs(TRUE))
#ptype=1 # for troubleshooting purposes
ptype=as.numeric(as.character(ptype))
print(ptype)
options(scipen=100)

# number of iterations per run of script # Q - how to choose this? * # MK: recommend that you may use iter=102 instead to target ~100,000 simulations.  
#iter=1 # run as iter=1 to get sense of time, then for actual jobs return to iter=102

# MK: Then, try to run the script with iter=6 and measure how long that takes. 
#iter=6

# Thu  3 Jun 2021 12:49:00 CEST MK - If half the data takes 40min, you can launch jobs with 51 iterations on 4 CPUs, which should take ~8.5h, so you can set 9h for the time. 
iter=51

# random numbers
randtyp<-sample(100000:500000,10000,replace=F)
#------------------------------------------------------------------------------------------------------------------------

mura<-1.235e-08 # mu, mutation rate
rbm<-9.40e-09 # R recombination rate
gent<-19 # generation time
slen<-40000 #length of region to simulate 40kb
# simulate 40kb b/c calculated S* on a genome-wide scale with a window size of 40kb

mr<-mura*slen*4000 ## scaled mutation rate
## scaled mutation rate. 4 * N * u * length 
  # reference N = 1000
  # => mu * length of fragment * 4 * 1000
sdv=mr*0.233
# standard deviation of the mutation rate - whether this should indeed be *0.233 for gor? ask ** use as a starting point
recv<-rbm*4000*slen
# scaled recombination rate. R * 4 * N * length of fragment

# sample mutation rate from normal dist (rnorm(n, mean = 0, sd = 1))
# sample recombination rate from -ve binomial dist (sampling parameter vals for input into ms from prob dist using tbs flag)

tetaro<-cbind(rnorm(20000,mr,sdv),rnbinom(20000,mu=recv,size=0.5))
  # rnorm(n, mean = 0, sd = 1) n=sample size, ie here number of reps to simulate - 20,000 reps
  # rnbinom(n, size, prob, mu) # where prob here is probability of success in each trial - https://stat.ethz.ch/R-manual/R-devel/library/stats/html/NegBinomial.html
      # size = dispersion parameter
      # mu = mean value for the -ve binomial

tetaro[,1]<-ifelse(tetaro[,1]<0.001,0.001,tetaro[,1])
tetaro[,2]<-ifelse(tetaro[,2]<0.001,0.001,tetaro[,2])

# Q whether to add ghost pop into types? or this is only for the s* setting?
types <- data.table(
#  out_pop = c(1,1,3,3),
#  in_pop = c(3,4,1,2),
  out_pop = c(1,1,3,3)-1,
  in_pop = c(3,4,1,2)-1,
  outgroup = c("WL","WL","EL","EL"),
  ingroup = c("EL","EM","WL","WC"),
  out_chr = c(54,54,18,18),
  in_chr = c(18,24,54,2)
)
# should be 27wl instead of 22 - Thu 10 Jun 2021 14:54:59 CEST
# 27*2 = 54, changing 44 to 54 & 22 to 27

types$sum <- types$out_chr + types$in_chr

comps<-list(c(54,0,18,0),c(54,0,0,24),c(54,0,18,0),c(0,2,18,0))

# this number of haploid chr to simulate: w_lowl, w_cros, e_lowl, e_moun
#ivec<-c(54,2,18,24);isu=sum(ivec);ivec=paste(ivec,collapse=" ")

# Q - check previous ms simulations - whether can sample 0 archaics in the present
  # don't think i got to the stage of adding an archaic there - b/c issues with the extant gorilla demography

# this number of haploid chr to simulate: w_lowl (1), w_cros (2), e_lowl (3), e_moun (4), ghost (5)
  # how many are sampled in the present - check if runs if specify 0 here **
ivec<-c(54,2,18,24,0);isu=sum(ivec);ivec=paste(ivec,collapse=" ")

#------------------------------------------------------------------------------------------------------------------------
# gorilla wide uniform priors for the parameter vals from gor.initial.priors.18mar.sh:
# range of parameters, randomly sampling from uniform distribution 1000 times

# MK: 2) Yes, fix the nice ones, and then use narrowed ranges for the rest (95%)

# MK
#3. Improvements in the priors.
#I suggest to remove t1 and use a fixed value because its so narrow. 
#For t2 and t3, you may take the midpoint and also take it out of the parameter search. 
#We are not that much interested in when the Ne changes happened. 
#Same for t6. I think it might even be irrelevant if the gene flow happens before or after the split of the western populations, 
#so maybe just take the midpoint time to population 1.
#So, lets try it that way and see how the simulations look like.
# For t6, I meant that the gene flow happens between the first population (WLG or WAnc), regardless of whether there is a second population (CRG).
# me
# for t6 then would you suggest setting this at 0.44 (34kya: 34,000/(4 * 1000 * 19) - the same as mcmanus)?
# MK
#yes, setting it to the literature value is fine.

#------------------------------------------------------------------------------------------------------------------------

#Mon 28 Feb 2022 15:16:25 CET *******
# need to go through & adjust these priors again?? **** YES - NEED TO BE UPDATED ***********
# adding weighted median posteriors from /Volumes/"Ultra USB 3.0"/IBE/further.analysis.feb.2020/gorillas/abc/generateabc.amendtajima.R - Tue  5 Apr 2022 14:40:10 CEST


#------------------------------------------------------------------------------------------------------------------------
# taking priors from: possible.priors.ghoste.28feb22.R # have checked w MK
#--------------------------------------------------------------------------------------------------------------------
# annotate below with 1) weighted median posteriors 2) visually what makes sense from the posterior curves 

# 1) priors for current Ne 
w_lowl_t0=64.7488 # w_lowl_t0 # could fix?
w_cros_t0=14.5582 # w_cros_t0 # fix
e_lowl_t0=20.8943 # e_lowl_t0 # fix
e_moun_t0=2.1578 # e_moun_t0 # fix
ghost_t0=runif(1000,min=0.1,max=100)
# 2) rapid & severe population decline of E lowland over last 20y (0.2kya) of ~80% - van der valk 2018, plumptr 2016
t1=0.0003947368
e_lowl_t1=0.2431 # e_lowl_t1 # fix
# 3)  Grauer’s gorillas went through a period of population growth and range expansion 5,000–10,000 years ago [22, 23]' van der valk 2019 citing van der valk 2018 & Tocheri 2016
t2=0.095
e_lowl_t2=21.9824 # e_lowl_t2 # fix
## 4) severe population bottleneck experienced by mountain gorillas - xue et al - this was likely after their split from eastern lowland
t3=0.12
e_moun_t3=0.1150 # e_moun_t3 # fix
# MK: 148 You may model this bottleneck as a temporal change, where at t3 it goes down to this value, and t3 minus 5 or 10 generations it goes back to the previous/terminal Ne.
t3.1=t3+0.0025 # 10 generations
e_moun_t3.1=runif(1000,min=2.9191,max=19.1638)
#e_moun_t3.1=9.9114 # e_moun_t3.1 - priors - take credible interval # 2.9191, 19.1638
## 5) E subspecies split # merge eas_moun into eas_lowl ** (adding parameters from xue et al. 2015 psmc into the mcmanus model)
t4=runif(1000,min=0.1446,max=0.2420)
#t4=0.1839 # t4  - priors t4 # 0.1446, 0.2420 
e_anc_t4=5.3253 # e_anc_t4 # fix 
t_archintrog=runif(1000,min=0.2445,max=14.3616) # this is taking the upper bound for t4 + 0.0025, but could take the lower bound & then filter out such reps?
# if take the lower bound woudl be min=0.1471
archaicintrog=runif(1000,min=0,max=100)
t_archintrog.1=t_archintrog+0.00025 # 1 generation
## 6) recent contraction in population size 23 ka in w lowland - mcmanus et al. 2015
t5=0.5381 # t5 # fix
w_lowl_t5=48.2880 # w_lowl_t5 # fix
# 7) migration pulse W > E 3% - for admixture need to give as a proportion not percent** 
t6=0.4473684
admix_w_e_t6=runif(1000,min=6.7955,max=99.2159)
#admix_w_e_t6=61.0343 # admix_w_e_t6 - priors #6.7955, 99.2159
admix_e_w_t6=33.0968 # admix_e_w_t6 # fix
t6.1=t6+0.00025 # 1 generation
# 8) W subspecies split ie merge wes_cros into wes_lowl
t7=5.9778 # t7 # fix
w_anc_t7=98.1352 # w_anc_t7 # fix
# 9) gor split (W - E) # merge eastern into western at time 261,000/(4 * 1000 * 19) = 3.434211
t8=12.7037 # t8 # fix
gor_anc=14.3641 # gor_anc # fix
# 11) split between common anc of all gorillas + ghost pop
t9=runif(1000,min=15,max=50)
gor_ghost_anc=runif(1000,min=10,max=100)


# whether to take upper bound of t4 for timing of ghost introg?
# or take within the range & => 
# GO BACK TO TAKING ONLY SIMNS WHERE X < Y, b/c now t_archintrog dependent on a non-fixed parameter of t4
# but on the other hand the event is unlikely to be so recent, so i think its fine if taking the upper bound for t4


#------------------------------------------------------------------------------------------------------------------------


inputsets=
  as.data.frame(t(cbind(
  sample(ghost_t0,iter,replace=T),
  sample(e_moun_t3.1,iter,replace=T),
  sample(t4,iter,replace=T),
  sample(t_archintrog,iter,replace=T),
  sample(archaicintrog,iter,replace=T),
  sample(admix_w_e_t6,iter,replace=T),
  sample(t9,iter,replace=T),
  sample(gor_ghost_anc,iter,replace=T),randtyp[1:iter])))



# fix parameters - t1, t2,t3,t6 => no longer subsample from the parameter ranges of priors for these


# subsample from the parameter ranges, to get iterations; add random number as identifier
#inputsets=as.list(as.data.frame(t(cbind(
#  sample(w_lowl_t0,iter,replace=T),
#  sample(w_cros_t0,iter,replace=T),
#  sample(e_lowl_t0,iter,replace=T),
#  sample(e_moun_t0,iter,replace=T),
#  #sample(t1,iter,replace=T),
#  sample(e_lowl_t1,iter,replace=T),
#  #sample(t2,iter,replace=T),
#  sample(e_lowl_t2,iter,replace=T),
#  #sample(t3,iter,replace=T),
#  sample(e_moun_t3,iter,replace=T),
#  #sample(t3.1,iter,replace=T),
#  sample(e_moun_t3.1,iter,replace=T),
#  sample(t4,iter,replace=T),
#  sample(e_anc_t4,iter,replace=T),
#  sample(t5,iter,replace=T),
#  sample(w_lowl_t5,iter,replace=T),
#  #sample(t6,iter,replace=T),
#  sample(admix_w_e_t6,iter,replace=T),
#  sample(admix_e_w_t6,iter,replace=T),
#  #sample(t6.1,iter,replace=T),
#  sample(t7,iter,replace=T),
#  sample(w_anc_t7,iter,replace=T),
#  sample(t8,iter,replace=T),
#  sample(gor_anc,iter,replace=T),randtyp[1:iter]))))


#------------------------------------------------------------------------------------------------------------------------
# MK:
#2. Improvements in the accepted simulations
#This could be done like that. You first create a generous set of values, 
#and only take iter numbers of iterations out of these that fulfill the criteria:

# priors shoudl be in the right directions -> not generating extra sets of values first
##inputsets=
 ## as.data.frame(t(cbind(
  ##sample(w_lowl_t0,iter*100,replace=T),
   ##sample(w_cros_t0,iter*100,replace=T),
   ##sample(e_lowl_t0,iter*100,replace=T),
  #sample(e_moun_t0,iter*100,replace=T),
   ##sample(ghost_t0,iter*100,replace=T),
  #sample(t1,iter*100,replace=T),
  #sample(e_lowl_t1,iter*100,replace=T),
  #sample(t2,iter*100,replace=T),
   ##sample(e_lowl_t2,iter*100,replace=T),
  #sample(t3,iter*100,replace=T),
  #sample(e_moun_t3,iter*100,replace=T),
  #sample(t3.1,iter,replace=T),
  #sample(e_moun_t3.1,iter*100,replace=T),
  #sample(t4,iter*100,replace=T),
  #sample(e_anc_t4,iter*100,replace=T),
   ##sample(t_archintrog,iter*100,replace=T),
   ##sample(archintrog,iter*100,replace=T),
  #sample(t_archintrog.1,iter,replace=T),
   ##sample(t5,iter*100,replace=T),
  #sample(w_lowl_t5,iter*100,replace=T),
  #sample(t6,iter*100,replace=T),
   ##sample(admix_w_e_t6,iter*100,replace=T),
   ##sample(admix_e_w_t6,iter*100,replace=T),
  #sample(t6.1,iter,replace=T),
  #sample(t7,iter*100,replace=T),
  #sample(w_anc_t7,iter*100,replace=T),
  #sample(t8,iter*100,replace=T),
  #sample(gor_anc,iter*100,replace=T),
   ##sample(t9,iter*100,replace=T),
   ##sample(gor_ghost_anc,iter*100,replace=T),randtyp[1:iter*100])))

## here, you would consider divergence times, but also directionality of population size changes, 
#all with adding a 5% margin.
# This is a suggestion, you can think about which parameters are the important ones to rule out.
#selec<-which(inputsets[21,]>inputsets[19,]*1.05 & inputsets[14,]*1.05<inputsets[19,] & inputsets[10,]*1.05<inputsets[11,] & inputsets[6,]*1.05<inputsets[8,] & inputsets[15,]*1.05<inputsets[20,] & inputsets[15,]>inputsets[1,]*1.05)

# amend to new numbering b/c have now fixed some parameters
#selec<-which(inputsets[21,]>inputsets[19,]*1.05 &  # time8=input[21], time7=input[19] time8 > time 7
 # inputsets[14,]*1.05<inputsets[19,] & # time5=input[14], time7=input[19] # time5 < time 7
 # inputsets[10,]*1.05<inputsets[11,] & # time3_em=input[10], time3.1_em=input[11] # size_em_afterbottleneck < size em_beforebottleneck
 # inputsets[6,]*1.05<inputsets[8,] & # time1_el=input[6], time2_el=input[8] # size_el_afterdecline < size_el_beforedecline
 # inputsets[15,]*1.05<inputsets[20,] & # time5_wl=input[15], time7_ancw=input[20] # size_wl_aftercontrac < size_wl_before
 # inputsets[15,]<inputsets[1,]*1.05) # time5_wl=input[15], initial_wl=input[1] # size_wl_ # should be < not > **


#selec<-which(inputsets[21,]>inputsets[19,]*1.05 &  
# time8=input[21], time7=input[19] time8 > time 7
# 17>15
  #inputsets[14,]*1.05<inputsets[19,] & 
 # time5=input[14], time7=input[19] # time5 < time 7
  # 11 < 15
 # inputsets[10,]*1.05<inputsets[11,] & 
  # time3_em=input[10], time3.1_em=input[11] # size_em_afterbottleneck < size em_beforebottleneck
  #  7 < 8
  #inputsets[6,]*1.05<inputsets[8,] & 
  # time1_el=input[6], time2_el=input[8] # size_el_afterdecline < size_el_beforedecline
  # 5 < 6
  #inputsets[15,]*1.05<inputsets[20,] & 
  # time5_wl=input[15], time7_ancw=input[20] # size_wl_aftercontrac < size_wl_before
  # 12 < 16
#  inputsets[15,]<inputsets[1,]*1.05) 
  # time5_wl=input[15], initial_wl=input[1] # size_wl_ # should be < not > **
 # 12 < 1

## here, you would consider divergence times, but also directionality of population size changes, 
#all with adding a 5% margin.
# This is a suggestion, you can think about which parameters are the important ones to rule out.
#selec<-which(inputsets[17,]>inputsets[15,]*1.05 & inputsets[11,]*1.05<inputsets[15,] & inputsets[7,]*1.05<inputsets[8,] & inputsets[5,]*1.05<inputsets[6,] & inputsets[12,]*1.05<inputsets[16,] &  inputsets[12,]<inputsets[1,]*1.05)              
#inputsets<-as.list(inputsets[,selec][,1:iter])

# update selec **
# to check orders of the parameters **
  # may need to uncomment & adjust lines 428 & 429

  
#------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------

#selec<-which(inputsets[17,]>inputsets[15,]*1.05 & inputsets[11,]*1.05<inputsets[15,] & inputsets[7,]*1.05<inputsets[8,] & inputsets[5,]*1.05<inputsets[6,] & inputsets[12,]*1.05<inputsets[16,] &  inputsets[12,]<inputsets[1,]*1.05) 
# i think now the selec step is not needed - b/c i'm adjusting the priors - Wed  3 Nov 2021 11:29:03 CET
    # Q - whether this is the way to go?? **

#------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------
# code from het.simld.15apr.R - MK amended - to be more efficient
# nested function to calculate heterozygosity & seg sites & fst: 
  # bhatiafst function from - readme_fstcomponents (script written with JS)

config<-c(54,2,18,24)
cofig<-cbind(c(0,cumsum(config/2)[-length(config)]),cumsum(config/2))
#sum(config)
len=40000


all_function<-function(nput) {
  msout <-read.ms.output(paste("/dev/shm/mydata/ghost.abc.sim_",nput,sep=""))
  ## create empty vector
    # bind pairs of cols, ie 1_2, 3_4, 5_6 etc (the haploids -> diploid ids)
    # MK: You would divide the number of hets by all sites for which you have data, in this case 40kbp (times the number of simulated windows)
    # this is for the 1 simn rep - will need to add another layer for going through all the simn reps
    # loop through every i with step size of 2 again
      # bind pairs of cols, ie 1_2, 3_4, 5_6 etc (the haploids -> diploid ids)
      # MK: To make it more efficient, you could even paste without the separator, then you get 00 01 10 and 11.
    # remove null elements of the list - https://stackoverflow.com/questions/33004238/r-removing-null-elements-from-a-list
    
  # generalise to any number of reps : # msout[grep("nreps",names(msout))] # https://stackoverflow.com/questions/4220631/how-do-i-grep-in-r
  ## apply the function for all simn reps from 1:msout[grep("nreps",names(msout))][[1]]
#    return(count_het)
#  }

  # MK: get genotype tables  
  vec<-seq(1,(ncol(msout$gametes[[1]])-1),2)
  ids.loop<-function(y){
    get_gt<-function(i) {     cbind(as.data.frame(paste0(msout$gametes[[y]][,c(i)], msout$gametes[[y]][,c(i+1)]))) }
    lapply(vec,get_gt)
  }
  ids=lapply(1:length(msout$gametes),ids.loop)

  ## MK: output Heterozygosity
  count_fun<-function(y,ids) {
    count_fun2<-function(i,y,ids) { length(which(ids[[y]][[i]]=="01")) + length(which(ids[[y]][[i]]=="10")) }
    unlist(lapply(1:length(vec),count_fun2,ids=ids,y=y))
    }  

  full=lapply(1:msout[grep("nreps",names(msout))][[1]],count_fun,ids=ids)
  full_het=rowSums(do.call(cbind,full))/(msout[grep("nreps",names(msout))][[1]]*len)
  mean_fun<-function(typ,fhets) { c(mean(fhets[(typ[1]+1):typ[2]]),sd(fhets[(typ[1]+1):typ[2]])) }
  het_op<-apply(cofig,1,mean_fun,fhets=full_het)*1000
  
  # estimating heterozygosity per simn rep - is this the way to go? # add "*1000" to this, to get heterozygosity per kbp, which is a more accesible value
  
  ## MK: output segsites # new way to calculate below
#  segfun<-function(y,ids) { 
#    idx<-do.call(cbind,ids[[y]])
#    gt_fun<-function(nput) { length(which(nput!="00")) }
#    pofun<-function(typ,idx) { pop=c((typ[1]+1):typ[2]); sum(apply(idx[,pop,drop=F],1,gt_fun)) }
#    apply(cofig,1,pofun,idx=idx)
#    }
#  full_seg=colSums(do.call(rbind,lapply(1:msout[grep("nreps",names(msout))][[1]],segfun,ids=ids)))

#------------------------------------------------------------------------------------------------------------------------
# Fri 10 Sep 2021 16:54:46 CEST
# MK:
#lets try to improve the simulations.
#I have the following suggestions.

#1. Improvements in the statistics. I found that there seem to be no fully fixed sites in the simulated data. 
#I suggest to calculate the following statistics: 
#population-wise fixed sites & segregating sites, individual fixed sites (like heterozygosity). 
#That could be done this way, building on your code:

## you need to define beforehands that cofig2 is:
cofig2<-cbind(c(0,cumsum(config)[-length(config)]),cumsum(config))

 #it makes sense to take the population-wise segregating sites as non-fixed sites. This seems like a more informative measure. 
 #This would exclude the CRG population, because this value is the same as heterozygosity with just one individua

#good that you are checking the things carefully! It seems that you need to say
#        op2=sum(op[which(as.numeric(names(op))>0 & as.numeric(names(op))<length(pop))]);op2=ifelse(length(op2)==0,0,op2)
#Only add "as.numeric()", the typical R problem... For safety also do that in op1:         
#op1=op[which(as.numeric(names(op))==length(pop))];op1=ifelse(length(op1)==0,0,op1)

 ## population fixed sites & segregating sites
segfun<-function(y,ids) {
      idx<-msout$gametes[[y]]
      pofun<-function(typ,idx) { pop=c((typ[1]+1):typ[2]);op<-table(rowSums(idx[,pop,drop=F]))
      op1=op[which(as.numeric(names(op))==length(pop))];op1=ifelse(length(op1)==0,0,op1)
      op2=sum(op[which(as.numeric(names(op))>0 & as.numeric(names(op))<length(pop))]);op2=ifelse(length(op2)==0,0,op2)
      c(op1,op2)}
      apply(cofig2,1,pofun,idx=idx)
  }
  full_segs=do.call(rbind,lapply(1:msout[grep("nreps",names(msout))][[1]],segfun,ids=ids))
  full_segsites=rbind(colSums(full_segs[seq(1,4999,2),]),colSums(full_segs[seq(2,5000,2),])) 

# MK
#1) The output of each iteration in full_segs is a data frame with two rows, each representing a different statistic. 
#This line takes the sum of each other row to give the final value for each statistic.

  ## individual fixed sites
  count_fun_f<-function(y,ids) {
    count_fun2_f<-function(i,y,ids) { length(which(ids[[y]][[i]]=="11")) }
    unlist(lapply(1:length(vec),count_fun2_f,ids=ids,y=y))
  } 
 
  full_f=lapply(1:msout[grep("nreps",names(msout))][[1]],count_fun_f,ids=ids)
  full_fix=rowSums(do.call(cbind,full_f))/(msout[grep("nreps",names(msout))][[1]]*len)
  mean_fun_f<-function(typ,ffix) { c(mean(ffix[(typ[1]+1):typ[2]]),sd(ffix[(typ[1]+1):typ[2]])) }
  fix_op<-apply(cofig,1,mean_fun_f,ffix=full_fix)*1000

# to output 
#population-wise fixed sites & segregating sites, individual fixed sites (like heterozygosity). 
# full_segsites, fix_op (yes not full_segs)
 
# MK: full_segsites is the # of population-wise fixed sites and the # of population-wise segregating sites; 
 # fix_op is the fixed sites per individual (same format as heterozygotes in your script).

#------------------------------------------------------------------------------------------------------------------------
  
  ## you may add fst here
  alinds<-c(rep(1,config[1]),rep(2,config[2]),rep(3,config[3]),rep(4,config[4]))
algam<-do.call(rbind,msout$gametes)
alfrfun<-function(nput) {
   inds<-which(alinds==nput)
   op<-rowSums(algam[,inds])/(length(inds))
   return(op)
   }
alfre<-lapply(1:4,alfrfun)

# frequencies needed for input into the fst function: 
bhatiaFST<-function(pop1f,pop2f,nchr1,nchr2) {
  altAllelePop1 <- round(pop1f*nchr1) # count of alt (1) alleles
  altAllelePop2 <- round(pop2f*nchr2)
  
  h1 <- (altAllelePop1*(nchr1-altAllelePop1))/(nchr1*(nchr1-1))
  h2 <- (altAllelePop2*(nchr2-altAllelePop2))/(nchr2*(nchr2-1))
  nHat <- (pop1f - pop2f)^2  - (h1/nchr1) - (h2/nchr2)
  dHat <- nHat+h1+h2
  fst <- nHat/dHat
  return(fst)
}

fst.all<-function(j,i){
fst.hold<-bhatiaFST(alfre[[j]],alfre[[i]], config[j], config[i])
# add this line into the function for calculating pairwise fst, set negative fst vals to 0 (considered loci with no popn differentiation)
fst.hold<-replace(fst.hold, fst.hold<0, 0) 
mean.fst<-mean(fst.hold, na.rm=TRUE) # remove nan vals
return(mean.fst)
}

fst.out<-rbind(
fst.all(1,2),
fst.all(1,3),
fst.all(1,4),
fst.all(2,3),
fst.all(2,4),
fst.all(3,4)
)
 
#----------------------------------------------------------------------------------------
# pi & tajima's d
refun<-function(input) { op=gsub(0,"C",input);gsub(1,"T",op) }
vec1<-seq(1,ncol(msout$gametes[[1]]))
  haploid.loop<-function(y){
    get_gt2<-function(i) {     as.data.frame((msout$gametes[[y]][,c(i)])) }
    lapply(vec1,get_gt2)
  }
  haploids=lapply(1:length(msout$gametes),haploid.loop)
format.function2<-function(y, i){
lapply(haploids[[y]][[i]],refun)
}

loop.pi<-function(x,y){
# add another layer  - loop through simn reps
reps_pop1=list()
for (j in 1:length(msout$gametes)){
pop1=list()
for (i in x: y) {
pop1[[i]]<-unlist(format.function2(j,i))
}
reps_pop1[[j]]<-pop1
}

hold.pi=list()
for (j in 1:length(msout$gametes)){
reps_pop1[[j]]<-reps_pop1[[j]][lengths(reps_pop1[[j]]) != 0]
hold.pi[[j]]<-nuc.div(as.DNAbin(reps_pop1[[j]]))
}
# output mean & sd
out.pi<-cbind(mean(unlist(hold.pi)), sd(unlist(hold.pi)))
return(out.pi)
}


# Wed 23 Mar 2022 14:24:36 CET
# updated tajimas d stat below
#loop.taj<-function(x,y){
# add another layer  - loop through simn reps
#reps_pop1=list()
#for (j in 1:length(msout$gametes)){
#pop1=list()
#for (i in x: y) {
#pop1[[i]]<-unlist(format.function2(j,i))
#}
#reps_pop1[[j]]<-pop1
#}
# add tajima's d to the function
#hold.tajima=list()
#for (j in 1:length(msout$gametes)){
#reps_pop1[[j]]<-reps_pop1[[j]][lengths(reps_pop1[[j]]) != 0]
#hold.tajima[[j]]<-tajima.test(as.DNAbin(reps_pop1[[j]]))
#}

#out.tajima<-cbind(mean(unlist(hold.tajima)), sd(unlist(hold.tajima)))
#return(out.tajima)
#}

fout.pi<-rbind(
loop.pi(1,54),
loop.pi(55,56),
loop.pi(57,74),
loop.pi(75,98)
)

#fout.taj<-rbind(
#loop.taj(1,54),
#loop.pi(45,46)[[2]],
#loop.taj(57,74),
#loop.taj(75,98)
#)

#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
# amended Thu 24 Mar 2022 15:28:41 CET


loop.taj<-function(x,y){
reps_pop1=list()
for (j in 1:length(msout$gametes)){
pop1=list()
for (i in x: y) {
pop1[[i]]<-unlist(format.function2(j,i))
}
reps_pop1[[j]]<-pop1
}

# add tajima's d to the function
hold.tajima=list()
for (j in 1:length(msout$gametes)){
reps_pop1[[j]]<-reps_pop1[[j]][lengths(reps_pop1[[j]]) != 0]
hold.tajima[[j]]<-tajima.test(as.DNAbin(reps_pop1[[j]]))
}
# need to add this step, b/c output of tajimas differs from pi
onlytajima=list()
for (j in 1:length(hold.tajima) ) {
onlytajima[[j]]<-hold.tajima[[j]][1]
}
#out.tajima<-cbind(mean(unlist(onlytajima)), sd(unlist(onlytajima)))
perwind<-unlist(onlytajima)
# tajimas d per window, set any nas to 0 (windows where no seg sites), then take mean & sd & output this
#perwind[is.na(perwind)] = 0 
#out.tajima<-cbind(mean(perwind), sd(perwind))
#return(out.tajima)
#return(list(out.tajima,perwind))
# MK: Tue 29 Mar 2022 10:26:15 CEST:  I guess its good if you also remove NAs instead of turning them to 0 (i.e. calculate mean and sd with na.rm=T),  as 0 by 0 might just be meaningless instead of "truly zero"
# Tue  5 Apr 2022 15:25:44 CEST - null model parameter inference, taking posteriors from calcn of taj d where set nas to 0 (ie method 2 -> so only output this here)
#out.tajima<-cbind(mean(perwind,na.rm=T), sd(perwind,na.rm=T))
# output both ways
perwind[is.na(perwind)] = 0 
out.tajima1<-cbind(mean(perwind), sd(perwind))
# first is removing na in calculation, 2nd list sets nas to 0
#return(list(out.tajima,out.tajima1))
return(out.tajima1)
}

# per pop (3: WL, EL, EM), output: mean tajimas d, sd tajimas d, per window tajimas d

allout.taj<-list(
loop.taj(1,54),
loop.taj(57,74),
loop.taj(75,98)
)


#------------------------------------------------------------------------------------------------------------------------

  ## output
  #return(list(het_op,full_seg,fst.out,fout.pi,fout.taj))
 # return(list(het_op,full_segsites,fix_op,fst.out,fout.pi,fout.taj))
 return(list(het_op,full_segsites,fix_op,fst.out,fout.pi,allout.taj))

}

#------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

wholefunction=function(input) {
# 24 parameter vals - now 22 as t3.1, t6.1 are not independent variables
# 22 - 4 = 18 parameters now (t1, t2,t3,t6 now fixed)
# 12 parameters now -  have now fixed - e_moun_t0,t1,e_lowl_t1,t2,t3,e_moun_t3,e_moun_t3.1,t4,e_anc_t4,w_lowl_t5,t6,t7,w_anc_t7,t8,gor_anc

initial_wl=w_lowl_t0
initial_wc=w_cros_t0
initial_el=e_lowl_t0
initial_em=e_moun_t0
initial_ghost=input[1]
time1=t1
time1_el=e_lowl_t1
time2=t2
time2_el=e_lowl_t2
time3=t3
time3_em=e_moun_t3
time3.1=time3+0.0025
time3.1_em=input[2]
time4=input[3]
time4_ance=e_anc_t4
timeI=input[4]
timeI_archgeneflow=input[5]
timeI.1=timeI+0.00025
time5=t5
time5_wl=w_lowl_t5
time6=t6
time6_geneflow_w_e=input[6]
time6_geneflow_e_w=admix_e_w_t6
time6.1=time6+0.00025
time7=t7
time7_ancw=w_anc_t7
time8=t8
time8_ancew=gor_anc
time9=input[7]
time9_allanc=input[8]
#------------------------------------------------------------------------------------------------------------------------
# generate simulation values, i.e. a list of complete parameters to input to the simulation
simval<-list(
c("-en",time1,paste(3,time1_el,sep=" ")),
c("-en",time2,paste(3,time2_el,sep=" ")),
c("-en",time3,paste(4,time3_em,sep=" ")),
c("-en",time3.1,paste(4,time3.1_em,sep=" ")),
c("-ej",time4,paste(4,3,sep=" ")),
c("-en",time4,paste(3,time4_ance,sep=" ")),
c("-em",timeI,paste(3,5,timeI_archgeneflow,sep=" ")),
c("-em",timeI.1,paste(3,5,0,sep=" ")),
c("-en",time5,paste(1,time5_wl,sep=" ")),
c("-em",time6,paste(3,1,time6_geneflow_w_e,sep=" ")),
c("-em",time6,paste(1,3,time6_geneflow_e_w,sep=" ")),
c("-em",time6.1,paste(3,1,0,sep=" ")),
c("-em",time6.1,paste(1,3,0,sep=" ")),
c("-ej",time7,paste(2,1,sep=" ")),
c("-en",time7,paste(1,time7_ancw,sep=" ")),
c("-ej",time8,paste(3,1,sep=" ")),
c("-en",time8,paste(1,time8_ancew,sep=" ")),
c("-ej",time9,paste(5,1,sep=" ")),
c("-en",time9,paste(1,time9_allanc,sep=" "))
)

# get an ordered list and merged part of the simulation parameters - ordering the events is necessary
simval<-do.call(rbind,simval);simval<-simval[order(as.numeric(simval[,2])),];simval=paste(t(simval),collapse=" ")

 
# generate the ms simulation 
# store data temporally - otherwise you flood the cluster - need to create tmpdir
# Q - how many simn reps to generate here? **
# MK: How many windows do you create? 20k? 
  #For this kind of simulation, that is clearly too much. I suggest 5000, which is 200Mbp of data: .. 5000 seeds ..  
  #That is already a lot of data, and for the measures we use, representative enough. Please try this and see how long it takes.
# 3/6/21 - generate 2500 windows - 100Mb of data with iter=51 instead of iter=6
#system(paste("cat /scratch/devel/hpawar/admix/sstar/simul/tetaroAr_corr.txt | /home/devel/hpawar/ms/msdir/ms ",isu," 2500 -seeds ",paste(sample(2000:6000,3),collapse=" ")," -t tbs -r tbs 40000.0 -I 4 54 2 18 24 -n 1 ",initial_wl," -n 2 ",initial_wc," -n 3 ",initial_el," -n 4 ",initial_em,"  ",simval," > /dev/shm/mydata/abc.sim_",input[23],sep=""),intern=F)

# ABC + ghost
system(paste("cat /scratch/devel/hpawar/admix/sstar/simul/tetaroAr_corr.txt | /home/devel/hpawar/ms/msdir/ms ",isu," 2500 -seeds ",paste(sample(2000:6000,3),collapse=" ")," -t tbs -r tbs 40000.0 -I 5 54 2 18 24 0 -n 1 ",initial_wl," -n 2 ",initial_wc," -n 3 ",initial_el," -n 4 ",initial_em," -n 5 ",initial_ghost,"  ",simval," > /dev/shm/mydata/ghost.abc.sim_",input[9],sep=""),intern=F)

#------------------------------------------------------------------------------------------------------------------------
# now calc summary statistics for this ms simn - code from het.simld.15apr.R
# Q - whether this shoudl be part of the 'wholefunction' also? * - or the all_function shoudl be outside the wholefunction? but called within it once?
  # that would make more sense 

#system(paste("rm /dev/shm/mydata/abc.sim_",input[23],sep=""),intern=F) # this was throwing errors here

return(all_function(input[9]))

}

# whether this wholefunction will work..

#remove temporary data (important when doing many!!)
# MK: 4) Also, VERY important, delete the things you create on /dev/shm/ !
  # Put this after all_function(input[25]) is done: system(paste("rm /dev/shm/mydata/abc.sim_",input[25],sep=""),intern=F)
  # If you dont do this, the temporary simulation files will just accumulate and accumulate, flooding the file system! 
#system(paste("rm /dev/shm/mydata/abc.sim_",input[23],sep=""),intern=F)

# everything above is done as one function instead of a continuous process in R
# this allow parallelization within the same launching of the script
# the following command generates simulated data for the "iter" number of times and interprets the output; it divides this task to 6 cores on the cluster, and collects the data of these "iter" iterations into a single object
simuresults<-mclapply(inputsets,wholefunction,mc.cores=6,mc.silent=F)
# this object, (importantly) together with the input parameters, is then stored as one object
# here, x simulations are done and being processed
# save new simulations to new dir: /scratch/devel/hpawar/admix/abc/simul/test/11sep21
# instead of /scratch/devel/hpawar/admix/abc/simul/test/out/
  # change dir - done
 # mkdir -p /scratch/devel/hpawar/admix/abc/simul/test/ghost/3nov21
#save(simuresults,inputsets,file=paste("/scratch/devel/hpawar/admix/abc/simul/test/11sep21/abc_sim",ptype,sep=""))

# need to change paths for output **
# mkdir -p /scratch/devel/hpawar/admix/abc/simul/test/ghost/final_1mar22

save(simuresults,inputsets,file=paste("/scratch/devel/hpawar/admix/abc/simul/test/ghost/final_1mar22/ghost.abc_sim",ptype,sep=""))

# this works - remove input ms files for this run here 
for (i in 1:iter){
#system(paste("rm /dev/shm/mydata/abc.sim_",inputsets[[i]][[23]],sep=""),intern=F) # have now fixed 4 more parameters
#system(paste("rm /dev/shm/mydata/abc.sim_",inputsets[[i]][[19]],sep=""),intern=F)
# ABC + ghost
system(paste("rm /dev/shm/mydata/ghost.abc.sim_",inputsets[[i]][[9]],sep=""),intern=F)
}
