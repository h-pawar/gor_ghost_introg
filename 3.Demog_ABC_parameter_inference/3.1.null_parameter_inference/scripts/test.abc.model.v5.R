#!/usr/bin/r

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

# Q - is the types no longer needed here?
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
ivec<-c(54,2,18,24);isu=sum(ivec);ivec=paste(ivec,collapse=" ")

#------------------------------------------------------------------------------------------------------------------------
# gorilla wide uniform priors for the parameter vals from gor.initial.priors.18mar.sh:
# range of parameters, randomly sampling from uniform distribution 1000 times


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


# 1) priors for current Ne 
w_lowl_t0=runif(1000,min=3,max=100) #-n 1 25.161 \  # wes_lowl pop size = 25.161 *10^3 - (parameters mcmanus 2015, using 12my divergence for humans & gorillas)
w_cros_t0=runif(1000,min=0.1,max=20) #-n 2 3.080 \  # wes_cros pop size = 3.080 *10^3 - mcmanus
e_lowl_t0=runif(1000,min=0.1,max=30) #-n 3 4.280 \  # eas_lowl pop size = 4.280 *10^3 - mcmanus
e_moun_t0=runif(1000,min=0.1,max=20) #-n 4 0.800 \  # eas_moun pop size = 800 - Xue et al. 2015 

# 2) rapid & severe population decline of E lowland over last 20y (0.2kya) of ~80% - van der valk 2018, plumptr 2016
# time of 20y: 20/(4 * 1000 * 19) [1] 0.0002631579 # take range of times around this (10-100y)
#t1=runif(1000,min=0.0001315789,max=0.0006578947) # 10 - 50ya
#(0.0001315789+0.0006578947)/2 # set t1 as midpoint (fixed val)
t1=0.0003947368
e_lowl_t1=runif(1000,min=0.01,max=20) 
  # here can i set the min val as the e_lowl_t0 Ne value?

# For t2 and t3, you may take the midpoint and also take it out of the parameter search. 

# 3)  Grauer’s gorillas went through a period of population growth and range expansion 5,000–10,000 years ago [22, 23]' van der valk 2019 citing van der valk 2018 & Tocheri 2016
#5000/(4 * 1000 * 19) [1] 0.06578947
#10000/(4 * 1000 * 19) [1] 0.1315789
#t2=runif(1000,min=0.05,max=0.14) 
# whether to instead have the pop size change of e lowland here?
  # but if i'm already overestimating the e diversity w/out adding these expansions - not sure how mcuh sense it makes 
# (0.05+0.14)/2 # midpoint
t2=0.095
e_lowl_t2=runif(1000,min=1,max=25) # 4280/0.2 =[1] 21400

## 4) severe population bottleneck experienced by mountain gorillas - xue et al - this was likely after their split from eastern lowland
#t3=runif(1000,min=0.1,max=0.14) 
# (0.1+0.14)/2 # midpoint
t3=0.12
e_moun_t3=runif(1000,min=0.01,max=5)

# MK: 148 You may model this bottleneck as a temporal change, where at t3 it goes down to this value, 
#and t3 minus 5 or 10 generations it goes back to the previous/terminal Ne.
# 5 generations = 5*19 [1] 95y # 10 gen = 190y (95/(4 * 1000 * 19) = [1] 0.00125)  190/(4 * 1000 * 19) [1] 0.0025
t3.1=t3+0.0025 # 10 generations
e_moun_t3.1=runif(1000,min=0.1,max=20)


## 5) E subspecies split # merge eas_moun into eas_lowl ** (adding parameters from xue et al. 2015 psmc into the mcmanus model)
# 15kya: 15,000/(4 * 1000 * 19) =  0.1973684
t4=runif(1000,min=0.14,max=0.25) 
e_anc_t4=runif(1000,min=0.1,max=30) #-n 3 4.280 # set size of e_anc

## 20kya = gene flow stops b/n e & w clades - xue et al
# Q - to still model as migration pulses? or as continuous?

## 6) recent contraction in population size 23 ka in w lowland - mcmanus et al. 2015
# 23kya: 23,000/(4 * 1000 * 19) =  [1] 0.3026316
# The second size change event occurred 22,800 years ago (95% CI: 16,457– 30,178) 
# and decreased the effective population size to 7,900 (95% CI: 6,433–9,240) individuals 
t5=runif(1000,min=0.19,max=0.6) 
w_lowl_t5=runif(1000,min=5,max=50)


# 7) migration pulse W > E 3% - for admixture need to give as a proportion not percent** 
# 34kya: 34,000/(4 * 1000 * 19) = 0.4473684
#t6=runif(1000,min=0.35,max=0.5)  # set as the literature val
t6=0.4473684
admix_w_e_t6=runif(1000,min=0,max=100) # admix w_lowl -> e_anc
admix_e_w_t6=runif(1000,min=0,max=100) # admix e_anc -> w_lowl

# 1 generation = 19 [1] 19y # (19/(4 * 1000 * 19) = [1]0.00025 
t6.1=t6+0.00025 # 1 generation
#\  #The migration pulse should happen for only one generation and then stop (that is, adding 19/4000 years as one generation):
#\  # ie set migrations back to 0
#-em time6.1 3 1 0 \
#-em time6.1 3 2 0 \

# 8) W subspecies split ie merge wes_cros into wes_lowl
# estimates of the separation of eastern gorillas from the western lowland/Cross River gorillas range from about 
#100 to 450 ka - cited by mcmanus et al 2015
  #450000/(4 * 1000 * 19) [1] 5.921053
# 68kya: 68,000/(4 * 1000 * 19) = 0.8947368 # was prev using
t7=runif(1000,min=0.1,max=6) 
w_anc_t7=runif(1000,min=10,max=100) # set size of w_anc

## long-term population decline over the last 100kya - decline was particularly pronounced in the eastern species ** (add in as pop size change from xue et al. 2015)
## 'backdrop of slow long-term reduction in population size' - van der valk et al. 2018
# Q - how to implement this?

#MK:182: The long-term decline could be also be interpreted as a long-term lower Ne. 
#In that sense, the difference between Ne at the node of each lineage (e_anc and w_anc) 
#and the ancestral Ne.
# q - how to specify this


# 9) gor split (W - E) # merge eastern into western at time 261,000/(4 * 1000 * 19) = 3.434211
# v large range for this 150kya - 1mya
#150000/(4 * 1000 * 19) [1] 1.973684
#1000000/(4 * 1000 * 19)[1] 13.15789
t8=runif(1000,min=1.5,max=15)  
gor_anc=runif(1000,min=10,max=100)

# MK: 190: I would even go 1.5-15. You run into a problem here which is that this split 
#could be younger (1.5) than the subsequent split in the west (6). 
#One solution would be to generate 10,000 random combinations, 
#create a dataframe with these, and remove the ones where t7 > t8; 
#and only use 1,000 (or another number) random combinations in the next step.


# Q - mcmanus also identifies this for w lowland finescale pop history - whether to include?
## history of western lowland gorillas that includes an ancestral population expansion of 1.4-fold around 970 ka - mcmanus 
## an ancestral effective population size of 31,800 (95% CI: 30,690–32,582) (table 2). 
## the first size change event occurred 969,000 years ago (95% CI: 764,074– 1,221,403) 
## and increased the effective population size to 44,200 (95% CI: 42,424–46,403) individuals. 

#------------------------------------------------------------------------------------------------------------------------

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

inputsets=
  as.data.frame(t(cbind(
  sample(w_lowl_t0,iter*100,replace=T),
  sample(w_cros_t0,iter*100,replace=T),
  sample(e_lowl_t0,iter*100,replace=T),
  sample(e_moun_t0,iter*100,replace=T),
  #sample(t1,iter*100,replace=T),
  sample(e_lowl_t1,iter*100,replace=T),
  #sample(t2,iter*100,replace=T),
  sample(e_lowl_t2,iter*100,replace=T),
  #sample(t3,iter*100,replace=T),
  sample(e_moun_t3,iter*100,replace=T),
  #sample(t3.1,iter,replace=T),
  sample(e_moun_t3.1,iter*100,replace=T),
  sample(t4,iter*100,replace=T),
  sample(e_anc_t4,iter*100,replace=T),
  sample(t5,iter*100,replace=T),
  sample(w_lowl_t5,iter*100,replace=T),
  #sample(t6,iter*100,replace=T),
  sample(admix_w_e_t6,iter*100,replace=T),
  sample(admix_e_w_t6,iter*100,replace=T),
  #sample(t6.1,iter,replace=T),
  sample(t7,iter*100,replace=T),
  sample(w_anc_t7,iter*100,replace=T),
  sample(t8,iter*100,replace=T),
  sample(gor_anc,iter*100,replace=T),randtyp[1:iter*100])))

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
selec<-which(inputsets[17,]>inputsets[15,]*1.05 & inputsets[11,]*1.05<inputsets[15,] & inputsets[7,]*1.05<inputsets[8,] & inputsets[5,]*1.05<inputsets[6,] & inputsets[12,]*1.05<inputsets[16,] &  inputsets[12,]<inputsets[1,]*1.05) 
             
inputsets<-as.list(inputsets[,selec][,1:iter])
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
  msout <-read.ms.output(paste("/dev/shm/mydata/abc.sim_",nput,sep=""))
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



loop.taj<-function(x,y){
# add another layer  - loop through simn reps
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

out.tajima<-cbind(mean(unlist(hold.tajima)), sd(unlist(hold.tajima)))
return(out.tajima)
}

fout.pi<-rbind(
loop.pi(1,54),
loop.pi(55,56),
loop.pi(57,74),
loop.pi(75,98)
)

fout.taj<-rbind(
loop.taj(1,54),
#loop.pi(45,46)[[2]],
loop.taj(57,74),
loop.taj(75,98)
)

#----------------------------------------------------------------------------------------

  ## output
  #return(list(het_op,full_seg,fst.out,fout.pi,fout.taj))
  return(list(het_op,full_segsites,fix_op,fst.out,fout.pi,fout.taj))

}

#------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

wholefunction=function(input) {
# 24 parameter vals - now 22 as t3.1, t6.1 are not independent variables
# 22 - 4 = 18 parameters now (t1, t2,t3,t6 now fixed)
initial_wl=input[1]
initial_wc=input[2]
initial_el=input[3]
initial_em=input[4]
time1=t1
time1_el=input[5]
time2=t2
time2_el=input[6]
time3=t3
time3_em=input[7]
time3.1=time3+0.0025
time3.1_em=input[8]
time4=input[9]
time4_ance=input[10]
time5=input[11]
time5_wl=input[12]
time6=t6
time6_geneflow_w_e=input[13]
time6_geneflow_e_w=input[14]
time6.1=time6+0.00025
time7=input[15]
time7_ancw=input[16]
time8=input[17]
time8_ancew=input[18]

#------------------------------------------------------------------------------------------------------------------------
# generate simulation values, i.e. a list of complete parameters to input to the simulation
simval<-list(
c("-en",time1,paste(3,time1_el,sep=" ")),
c("-en",time2,paste(3,time2_el,sep=" ")),
c("-en",time3,paste(4,time3_em,sep=" ")),
c("-en",time3.1,paste(4,time3.1_em,sep=" ")),
c("-ej",time4,paste(4,3,sep=" ")),
c("-en",time4,paste(3,time4_ance,sep=" ")),
c("-en",time5,paste(1,time5_wl,sep=" ")),
c("-em",time6,paste(3,1,time6_geneflow_w_e,sep=" ")),
c("-em",time6,paste(1,3,time6_geneflow_e_w,sep=" ")),
c("-em",time6.1,paste(3,1,0,sep=" ")),
c("-em",time6.1,paste(1,3,0,sep=" ")),
c("-ej",time7,paste(2,1,sep=" ")),
c("-en",time7,paste(1,time7_ancw,sep=" ")),
c("-ej",time8,paste(3,1,sep=" ")),
c("-en",time8,paste(1,time8_ancew,sep=" "))
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
system(paste("cat /scratch/devel/hpawar/admix/sstar/simul/tetaroAr_corr.txt | /home/devel/hpawar/ms/msdir/ms ",isu," 2500 -seeds ",paste(sample(2000:6000,3),collapse=" ")," -t tbs -r tbs 40000.0 -I 4 54 2 18 24 -n 1 ",initial_wl," -n 2 ",initial_wc," -n 3 ",initial_el," -n 4 ",initial_em,"  ",simval," > /dev/shm/mydata/abc.sim_",input[19],sep=""),intern=F)

#------------------------------------------------------------------------------------------------------------------------
# now calc summary statistics for this ms simn - code from het.simld.15apr.R
# Q - whether this shoudl be part of the 'wholefunction' also? * - or the all_function shoudl be outside the wholefunction? but called within it once?
  # that would make more sense 

all_function(input[19])

#system(paste("rm /dev/shm/mydata/abc.sim_",input[23],sep=""),intern=F) # this was throwing errors here

return(all_function(input[19]))

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
save(simuresults,inputsets,file=paste("/scratch/devel/hpawar/admix/abc/simul/test/11sep21/abc_sim",ptype,sep=""))

# this works - remove input ms files for this run here 
for (i in 1:iter){
#system(paste("rm /dev/shm/mydata/abc.sim_",inputsets[[i]][[23]],sep=""),intern=F) # have now fixed 4 more parameters
system(paste("rm /dev/shm/mydata/abc.sim_",inputsets[[i]][[19]],sep=""),intern=F)
}

#-----------------------------------------------------------------------------------------------------------------------
# MK:
#following lines need correction:
#183-206: you sample iter times from each value : "sample(w_lowl_t0,iter,replace=T)". 
#Otherwise, you just shuffle the 1000 random values, but then cbind with only 92 iterations, 
#which does not make sense and may cause problems. 
#That number 92 is kind of random, as I was targeting to get ~90,000 simulations, and some of the jobs occasionally fail. 
#I recommend that you may use  iter=102 instead to  target ~100,000 simulations. 
#As a result, the length of your inputsets object is 102. 
#Running the function will iterate through these 102 random sets and save the results together in one object.

# done

#199,202: This way times don't correspond! 
#t6 is the start of gene flow, randomly sampled, t6.1 is another random set of times. 
#If you have a closer look at the "simval" later on, you see that it doesn't match: 
#the gene flow will be 1000s of generations long instead of 1 generation. 
#In that sense, t6.1 is not an independent variable here to be used in the input vector. 
#Same is true for t3 and t3.1. Better define t6.1 and t3.1 at lines 334 and 325!

#t3.1=t3+0.0025 # 10 generations
#t6.1=t6+0.00025 # 1 generation

#time3.1=time3+0.0025
#timet6.1=time6+0.00025

# done

#366: You should use this file name: /dev/shm/mydata/abc.sim_",input[25],sep=""; "

#otherwise you will create files with the name 1, 2,..92 over and over again. 
#The random number ensures you create unique files for each simulation, which is very important, 
#while "ptype" is the name for a subset of (in this case 92) simulations.
#You will run the script 1000 independent times, each creating 102 random combinations of values, 
#each with a random number assigned. 
#Before doing so, I recommend running 1 job (or start with 1 full iteration) to see how long it takes, and adjusting the time in the submission script. 
#It should be approximately the time of one iteration*102/6.

# input[23] now (as have removed 2 parameters from sampling)
# instead of: /dev/shm/mydata/abc.sim_",ptype,sep="") "
# done - replaced ptype with input[23]

#373: all_function(input[25]) to analyze the proper simulation and nothing else
# input[23]

#220: for clarity, use different identifiers in nested functions, for example "nput" instead of "input" here
# done

#215: this should be named "cofig", as this object is used later under this name
# done

#253: you may add "*1000" to this, to get heterozygosity per kbp, which is a more accesible value

#Apart from that, in principle, the code seems to work, so it's on a good track to go for the modeling!'


#het_op
#            [,1]     [,2]        [,3]         [,4]
#[1,] 0.002152727 0.002075 0.001954444 0.0009300000
#[2,] 0.000756544       NA 0.000636948 0.0005741872
#> het_op*1000
#         [,1]  [,2]     [,3]      [,4]
#[1,] 2.152727 2.075 1.954444 0.9300000
#[2,] 0.756544    NA 0.636948 0.5741872

#> apply(cofig,1,mean_fun,fhets=full_het)
#            [,1]     [,2]        [,3]         [,4]
#[1,] 0.002152727 0.002075 0.001954444 0.0009300000
#[2,] 0.000756544       NA 0.000636948 0.0005741872

#> apply(cofig,1,mean_fun,fhets=full_het)*1000
#         [,1]  [,2]     [,3]      [,4]
#[1,] 2.152727 2.075 1.954444 0.9300000
#[2,] 0.756544    NA 0.636948 0.5741872

# done
#-----------------------------------------------------------------------------------------------------------------------

